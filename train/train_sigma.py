"""
train_sigma.py — Train the peak-sigma regression MLP and export to ONNX.

Architecture:
    Linear(51→32) → ReLU → Linear(32→16) → ReLU → Linear(16→1) → Softplus

The output is sigma in bins (always positive).  Multiply by bin width (keV)
in C++ to get sigma in keV.  Constraining MIGRAD's sigma to ±8% of this
prediction breaks the A–σ correlation and reduces area uncertainties.

Loss: Huber (SmoothL1Loss, delta=0.5 bins) — robust to outlier peaks
Optimiser: Adam, lr=1e-3, 40 epochs

Input data: data/combined.npz or data/synthetic.npz  (X_sigma / y_sigma keys)
Output:     ml/sigma_model.onnx

Usage:
    python train/train_sigma.py [--data data/synthetic.npz] [--epochs 40]
"""

import argparse
import os
import sys
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split


# ── model ─────────────────────────────────────────────────────────────────────
class SigmaNet(nn.Module):
    """51-bin window → peak sigma in bins (positive scalar)."""

    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(51, 32), nn.ReLU(),
            nn.Linear(32, 16), nn.ReLU(),
            nn.Linear(16,  1), nn.Softplus(),   # guarantees positive output
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


# ── training loop ─────────────────────────────────────────────────────────────
def run_epoch(model, loader, criterion, optimiser, device, train=True):
    model.train(train)
    total_loss = 0.0
    with torch.set_grad_enabled(train):
        for X, y in loader:
            X, y = X.to(device), y.to(device)
            if train:
                optimiser.zero_grad()
            pred = model(X)
            loss = criterion(pred, y)
            if train:
                loss.backward()
                optimiser.step()
            total_loss += loss.item() * len(X)
    return total_loss / len(loader.dataset)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Train sigma regression MLP")
    parser.add_argument("--data",   default=None,
                        help="Path to .npz dataset (default: data/combined.npz or data/synthetic.npz)")
    parser.add_argument("--epochs", type=int,   default=40)
    parser.add_argument("--lr",     type=float, default=1e-3)
    parser.add_argument("--batch",  type=int,   default=2048)
    parser.add_argument("--val",    type=float, default=0.1,
                        help="Validation fraction (default 0.1)")
    parser.add_argument("--out",    default="ml/sigma_model.onnx")
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
    d = np.load(data_path)
    if "X_sigma" not in d:
        print("ERROR: dataset missing X_sigma / y_sigma keys.")
        print("       Re-run generate_spectra.py to regenerate with sigma labels.")
        sys.exit(1)

    X_all = torch.tensor(d["X_sigma"], dtype=torch.float32)
    y_all = torch.tensor(d["y_sigma"], dtype=torch.float32)
    print(f"  {len(X_all):,} sigma windows  "
          f"(sigma mean={y_all.mean().item():.4f}  std={y_all.std().item():.4f} bins)")

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

    model     = SigmaNet().to(device)
    criterion = nn.SmoothL1Loss(beta=0.5)   # Huber loss, robust to outlier peaks
    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, patience=5, factor=0.5, verbose=True)

    best_val   = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        tr_loss = run_epoch(model, train_loader, criterion, optimiser, device, train=True)
        vl_loss = run_epoch(model, val_loader,   criterion, optimiser, device, train=False)
        scheduler.step(vl_loss)
        print(f"Epoch {epoch:3d}/{args.epochs}  "
              f"train Huber={tr_loss:.6f}  val Huber={vl_loss:.6f}", flush=True)
        if vl_loss < best_val:
            best_val   = vl_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    print(f"\nBest validation Huber loss: {best_val:.6f}")
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
