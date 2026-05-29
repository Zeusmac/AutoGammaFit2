"""
train_bg.py — Train the background estimation MLP and export to ONNX.

Architecture:
    Linear(101→64) → ReLU → Linear(64→32) → ReLU → Linear(32→1)

Loss: MSELoss (regression — predicts normalised background at centre bin)
Optimiser: Adam, lr=1e-3, 20 epochs

Input data: data/combined.npz  (or data/synthetic.npz if combined not present)
Output:     ml/bg_model.onnx

Usage:
    python train/train_bg.py [--data data/combined.npz] [--epochs 20]
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
        nn.Linear(101, 64), nn.ReLU(),
        nn.Linear(64,  32), nn.ReLU(),
        nn.Linear(32,   1),
    )


# ── training loop ─────────────────────────────────────────────────────────────
def train(model, loader, criterion, optimiser, device):
    model.train()
    total_loss = 0.0
    for X, y in loader:
        X, y = X.to(device), y.to(device)
        optimiser.zero_grad()
        pred = model(X).squeeze(-1)
        loss = criterion(pred, y)
        loss.backward()
        optimiser.step()
        total_loss += loss.item() * len(X)
    return total_loss / len(loader.dataset)


def evaluate(model, loader, criterion, device):
    model.eval()
    total_loss = 0.0
    with torch.no_grad():
        for X, y in loader:
            X, y = X.to(device), y.to(device)
            pred = model(X).squeeze(-1)
            total_loss += criterion(pred, y).item() * len(X)
    return total_loss / len(loader.dataset)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Train background estimation MLP")
    parser.add_argument("--data",    default=None,
                        help="Path to .npz dataset (default: data/combined.npz or data/synthetic.npz)")
    parser.add_argument("--epochs",  type=int,   default=20)
    parser.add_argument("--lr",      type=float, default=1e-3)
    parser.add_argument("--batch",   type=int,   default=4096)
    parser.add_argument("--val",     type=float, default=0.1,
                        help="Validation fraction (default 0.1)")
    parser.add_argument("--out",     default="ml/bg_model.onnx")
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
    X_all = torch.tensor(d["X_bg"], dtype=torch.float32)
    y_all = torch.tensor(d["y_bg"], dtype=torch.float32)
    print(f"  {len(X_all):,} windows  (shape {tuple(X_all.shape)})")

    # Train / validation split
    n_val   = int(len(X_all) * args.val)
    n_train = len(X_all) - n_val
    ds      = TensorDataset(X_all, y_all)
    train_ds, val_ds = random_split(ds, [n_train, n_val],
                                    generator=torch.Generator().manual_seed(42))
    train_loader = DataLoader(train_ds, batch_size=args.batch, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch * 2, shuffle=False,
                              num_workers=0)

    device    = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    model     = build_model().to(device)
    criterion = nn.MSELoss()
    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, patience=3, factor=0.5, verbose=True)

    best_val  = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        tr_loss = train(model, train_loader, criterion, optimiser, device)
        vl_loss = evaluate(model, val_loader, criterion, device)
        scheduler.step(vl_loss)
        print(f"Epoch {epoch:3d}/{args.epochs}  "
              f"train MSE={tr_loss:.5f}  val MSE={vl_loss:.5f}")
        if vl_loss < best_val:
            best_val   = vl_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    print(f"\nBest validation MSE: {best_val:.5f}")
    model.load_state_dict(best_state)

    # ── Export to ONNX ────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    dummy = torch.zeros(1, 101, dtype=torch.float32, device=device)
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
