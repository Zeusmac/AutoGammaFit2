"""
train_centroid.py — Train the sub-bin centroid regression MLP and export to ONNX.

Architecture:
    Linear(51→32) → ReLU → Linear(32→16) → ReLU → Linear(16→1) → Tanh × 0.5

The output is in bins ∈ [-0.5, 0.5].  Multiply by bin width (keV) in C++ to get keV offset.

Loss: MSELoss
Optimiser: Adam, lr=1e-3, 30 epochs

Input data: data/combined.npz or data/synthetic.npz  (X_centroid / y_centroid keys)
Output:     ml/centroid_model.onnx

Usage:
    python train/train_centroid.py [--data data/synthetic.npz] [--epochs 30]
"""

import argparse
import os
import sys
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split


# ── model ─────────────────────────────────────────────────────────────────────
class CentroidNet(nn.Module):
    """51-bin window → sub-bin centroid offset in [-0.5, 0.5] bins."""

    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(51, 32), nn.ReLU(),
            nn.Linear(32, 16), nn.ReLU(),
            nn.Linear(16,  1),
        )

    def forward(self, x):
        return torch.tanh(self.net(x)) * 0.5


# ── training loop ─────────────────────────────────────────────────────────────
def run_epoch(model, loader, criterion, optimiser, device, train=True):
    model.train(train)
    total_loss = 0.0
    with torch.set_grad_enabled(train):
        for X, y in loader:
            X, y = X.to(device), y.to(device)
            if train:
                optimiser.zero_grad()
            pred = model(X).squeeze(-1)
            loss = criterion(pred, y)
            if train:
                loss.backward()
                optimiser.step()
            total_loss += loss.item() * len(X)
    return total_loss / len(loader.dataset)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Train centroid regression MLP")
    parser.add_argument("--data",   default=None,
                        help="Path to .npz dataset (default: data/combined.npz or data/synthetic.npz)")
    parser.add_argument("--epochs", type=int,   default=30)
    parser.add_argument("--lr",     type=float, default=1e-3)
    parser.add_argument("--batch",  type=int,   default=2048)
    parser.add_argument("--val",    type=float, default=0.1,
                        help="Validation fraction (default 0.1)")
    parser.add_argument("--out",    default="ml/centroid_model.onnx")
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
    if "X_centroid" not in d:
        print("ERROR: dataset missing X_centroid / y_centroid keys.")
        print("       Re-run generate_spectra.py to regenerate with centroid labels.")
        sys.exit(1)

    X_all = torch.tensor(d["X_centroid"], dtype=torch.float32)
    y_all = torch.tensor(d["y_centroid"], dtype=torch.float32)
    print(f"  {len(X_all):,} centroid windows  "
          f"(offset rms: {y_all.std().item():.4f} bins)")

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

    model     = CentroidNet().to(device)
    criterion = nn.MSELoss()
    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, patience=4, factor=0.5, verbose=True)

    best_val   = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        tr_loss = run_epoch(model, train_loader, criterion, optimiser, device, train=True)
        vl_loss = run_epoch(model, val_loader,   criterion, optimiser, device, train=False)
        scheduler.step(vl_loss)
        print(f"Epoch {epoch:3d}/{args.epochs}  "
              f"train MSE={tr_loss:.6f}  val MSE={vl_loss:.6f}  "
              f"(RMSE={vl_loss**0.5:.4f} bins)")
        if vl_loss < best_val:
            best_val   = vl_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    print(f"\nBest validation MSE: {best_val:.6f}  "
          f"(RMSE={best_val**0.5:.4f} bins)")
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
