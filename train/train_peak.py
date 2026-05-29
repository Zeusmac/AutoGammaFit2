"""
train_peak.py — Train the peak detector MLP and export to ONNX.

Architecture:
    Linear(51→32) → ReLU → Linear(32→16) → ReLU → Linear(16→1) → Sigmoid

Loss: BCELoss with positive class weight = 20
      (class imbalance: ~1–5 % of bins are peaks)
Optimiser: Adam, lr=1e-3, 20 epochs

Input data: data/combined.npz  (or data/synthetic.npz if combined not present)
Output:     ml/peak_model.onnx

Usage:
    python train/train_peak.py [--data data/combined.npz] [--epochs 20]
"""

import argparse
import os
import sys
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split


# ── model ─────────────────────────────────────────────────────────────────────
def build_model():
    return nn.Sequential(
        nn.Linear(51, 32), nn.ReLU(),
        nn.Linear(32, 16), nn.ReLU(),
        nn.Linear(16,  1), nn.Sigmoid(),
    )


# ── metrics ───────────────────────────────────────────────────────────────────
def compute_metrics(model, loader, device, threshold=0.5):
    """Return (loss, precision, recall, F1)."""
    model.eval()
    total_loss = 0.0
    tp = fp = fn = 0
    criterion = nn.BCELoss()
    with torch.no_grad():
        for X, y in loader:
            X, y  = X.to(device), y.to(device)
            prob  = model(X).squeeze(-1)
            total_loss += criterion(prob, y).item() * len(X)
            pred   = (prob >= threshold).float()
            tp += ((pred == 1) & (y == 1)).sum().item()
            fp += ((pred == 1) & (y == 0)).sum().item()
            fn += ((pred == 0) & (y == 1)).sum().item()
    loss = total_loss / len(loader.dataset)
    prec = tp / (tp + fp + 1e-9)
    rec  = tp / (tp + fn + 1e-9)
    f1   = 2 * prec * rec / (prec + rec + 1e-9)
    return loss, prec, rec, f1


# ── training loop ─────────────────────────────────────────────────────────────
def train_epoch(model, loader, criterion, optimiser, device):
    model.train()
    total_loss = 0.0
    for X, y in loader:
        X, y = X.to(device), y.to(device)
        optimiser.zero_grad()
        prob = model(X).squeeze(-1)
        loss = criterion(prob, y)
        loss.backward()
        optimiser.step()
        total_loss += loss.item() * len(X)
    return total_loss / len(loader.dataset)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Train peak detector MLP")
    parser.add_argument("--data",       default=None)
    parser.add_argument("--epochs",     type=int,   default=20)
    parser.add_argument("--lr",         type=float, default=1e-3)
    parser.add_argument("--batch",      type=int,   default=4096)
    parser.add_argument("--pos-weight", type=float, default=20.0,
                        help="BCELoss positive class weight (default 20)")
    parser.add_argument("--val",        type=float, default=0.1)
    parser.add_argument("--out",        default="ml/peak_model.onnx")
    args = parser.parse_args()

    # Locate dataset
    if args.data:
        data_path = args.data
    elif os.path.exists("data/combined.npz"):
        data_path = "data/combined.npz"
    elif os.path.exists("data/synthetic.npz"):
        data_path = "data/synthetic.npz"
    else:
        print("ERROR: no dataset found. Run generate_spectra.py first.")
        sys.exit(1)

    print(f"Loading {data_path} …")
    d     = np.load(data_path)
    X_all = torch.tensor(d["X_peak"], dtype=torch.float32)
    y_all = torch.tensor(d["y_peak"], dtype=torch.float32)
    peak_frac = y_all.mean().item()
    print(f"  {len(X_all):,} windows  peak fraction = {peak_frac:.4f}")

    n_val    = int(len(X_all) * args.val)
    n_train  = len(X_all) - n_val
    ds       = TensorDataset(X_all, y_all)
    train_ds, val_ds = random_split(ds, [n_train, n_val],
                                    generator=torch.Generator().manual_seed(42))
    train_loader = DataLoader(train_ds, batch_size=args.batch, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch * 2, shuffle=False,
                              num_workers=0)

    device    = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    model     = build_model().to(device)
    pos_w     = torch.tensor([args.pos_weight], device=device)
    criterion = nn.BCELoss(weight=None)   # weight applied manually below
    # BCEWithLogitsLoss is more numerically stable; since Sigmoid is baked in we
    # use BCELoss but scale positive-class losses via sample weights.
    criterion = nn.BCELoss(reduction="none")

    def weighted_bce(prob, y):
        raw = criterion(prob, y)
        w   = torch.where(y == 1, pos_w, torch.ones_like(y))
        return (raw * w).mean()

    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, patience=3, factor=0.5, verbose=True)

    best_f1    = 0.0
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        for X, y in train_loader:
            X, y = X.to(device), y.to(device)
            optimiser.zero_grad()
            prob = model(X).squeeze(-1)
            loss = weighted_bce(prob, y)
            loss.backward()
            optimiser.step()
            total_loss += loss.item() * len(X)
        tr_loss = total_loss / n_train

        vl_loss, prec, rec, f1 = compute_metrics(model, val_loader, device)
        scheduler.step(vl_loss)
        print(f"Epoch {epoch:3d}/{args.epochs}  "
              f"train={tr_loss:.5f}  val={vl_loss:.5f}  "
              f"P={prec:.3f}  R={rec:.3f}  F1={f1:.3f}")

        if f1 > best_f1:
            best_f1    = f1
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    print(f"\nBest validation F1: {best_f1:.4f}")
    model.load_state_dict(best_state)

    # ── Export to ONNX ────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    dummy = torch.zeros(1, 51, dtype=torch.float32, device=device)
    torch.onnx.export(
        model.cpu(), dummy.cpu(), args.out,
        input_names=["input"],
        output_names=["output"],
        dynamic_axes={"input": {0: "batch"}, "output": {0: "batch"}},
        opset_version=17,
    )
    print(f"Exported → {args.out}")


if __name__ == "__main__":
    main()
