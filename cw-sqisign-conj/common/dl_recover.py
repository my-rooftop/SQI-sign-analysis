#!/usr/bin/env python3
"""
Deep-learning side-channel attack on Algorithm 3.18 line 2.

Why try this after templates failed
-----------------------------------
The template and LDA attacks assumed a leakage model -- Hamming weight of a limb,
or of a byte of it -- and under that assumption the leakage is symmetric in the
bits, so no individual bit is recoverable however many traces are taken. A neural
network assumes no model: it learns whatever function of the samples separates
the classes, including non-linear combinations across many points. If anything
beyond the Hamming weight is present in the trace, this is the method that finds
it.

What is trained
---------------
Per target, a profiling split trains an MLP on a window cropped around where that
limb leaks, and a disjoint held-out split scores it.

  hw     Hamming weight of a limb, 0..32. A sanity target: it is known to leak,
         so if the network cannot learn this the pipeline is broken and the
         negative results below mean nothing.
  bit    one bit of a limb, binary. The decisive target -- this is what the
         Hamming-weight model says is impossible.
  byte   one byte of a limb, 256 classes. Reported for the bits that matter if
         bit recovery works at all.

Every number is compared against the majority-class baseline on the same split,
because a skewed label makes a useless classifier look accurate.

Usage:
  ~/.cwvenv/bin/python3 dl_recover.py --npz /tmp/sqisign_line2.npz
  ~/.cwvenv/bin/python3 dl_recover.py --npz ... --targets hw,bit --coord 1 --limb 0
"""

import argparse
import pathlib
import sys
import time

import numpy as np

try:
    import torch
    import torch.nn as nn
except ImportError:
    print("torch is required", file=sys.stderr)
    raise


def find_poi(T_prof, labels, n=1):
    """Sample(s) most correlated with the label, on the profiling split only."""
    x = np.asarray(labels, dtype=np.float64)
    x = x - x.mean()
    xn = np.sqrt((x * x).sum())
    Tc = T_prof - T_prof.mean(axis=0, keepdims=True)
    Tn = np.sqrt((Tc * Tc).sum(axis=0))
    Tn = np.where(Tn == 0, np.finfo(np.float64).tiny, Tn)
    r = np.abs((Tc * x[:, None]).sum(axis=0) / (Tn * xn))
    return int(np.argmax(r)), r


class MLP(nn.Module):
    def __init__(self, n_in, n_out, width=512, depth=3, p_drop=0.2):
        super().__init__()
        layers, d = [], n_in
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.BatchNorm1d(width), nn.ReLU(), nn.Dropout(p_drop)]
            d = width
        layers.append(nn.Linear(d, n_out))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def train_eval(Xp, yp, Xv, yv, Xa, ya, n_cls, epochs, lr, batch, device, seed=0,
               shuffle_labels=False):
    """Train on (Xp, yp), pick the epoch on (Xv, yv), report on (Xa, ya).

    Taking the best epoch by held-out accuracy would be selecting the model on
    the test set: over 60 epochs the maximum of a 50%-baseline binary accuracy
    drifts to about 52% on 4000 traces by chance alone, which is the same size as
    the effect being looked for. The epoch is chosen on a separate validation
    split instead.

    shuffle_labels permutes the training labels as a control: anything above
    baseline there is a flaw in the pipeline, not leakage.
    """

    torch.manual_seed(seed)
    if shuffle_labels:
        yp = np.random.default_rng(seed).permutation(yp)
    mu, sd = Xp.mean(0, keepdims=True), Xp.std(0, keepdims=True) + 1e-8
    Xp = (Xp - mu) / sd
    Xv = (Xv - mu) / sd
    Xa = (Xa - mu) / sd

    Xp_t = torch.tensor(Xp, dtype=torch.float32, device=device)
    yp_t = torch.tensor(yp, dtype=torch.long, device=device)
    Xv_t = torch.tensor(Xv, dtype=torch.float32, device=device)
    yv_t = torch.tensor(yv, dtype=torch.long, device=device)
    Xa_t = torch.tensor(Xa, dtype=torch.float32, device=device)
    ya_t = torch.tensor(ya, dtype=torch.long, device=device)

    model = MLP(Xp.shape[1], n_cls).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=lr, total_steps=epochs * max(1, len(Xp_t) // batch))
    lossf = nn.CrossEntropyLoss()

    best_val, test_at_best = -1.0, 0.0
    n = len(Xp_t)
    for ep in range(epochs):
        model.train()
        perm = torch.randperm(n, device=device)
        for i in range(0, n - batch + 1, batch):
            idx = perm[i:i + batch]
            opt.zero_grad()
            loss = lossf(model(Xp_t[idx]), yp_t[idx])
            loss.backward()
            opt.step()
            sched.step()
        model.eval()
        with torch.no_grad():
            va = (model(Xv_t).argmax(1) == yv_t).float().mean().item()
            ta = (model(Xa_t).argmax(1) == ya_t).float().mean().item()
        if va > best_val:
            best_val, test_at_best = va, ta
    return test_at_best


def main():
    p = argparse.ArgumentParser(description="Deep-learning SCA on Alg 3.18 line 2")
    p.add_argument("--npz", required=True)
    p.add_argument("--coord", type=int, default=1)
    p.add_argument("--limbs", default="0,1,2")
    p.add_argument("--targets", default="hw,bit,byte")
    p.add_argument("--bits", default="0,1,7,15,31", help="bit positions to try")
    p.add_argument("--window", type=int, default=1000, help="half-width of the crop")
    p.add_argument("--epochs", type=int, default=60)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--batch", type=int, default=256)
    p.add_argument("--traces", type=int, default=None)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    device = "cuda" if torch.cuda.is_available() else "cpu"
    d = np.load(args.npz, allow_pickle=False)
    T = d["traces"]
    if args.traces:
        T = T[:args.traces]
    n = T.shape[0]
    T = T.astype(np.float32)

    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(n)
    n_prof, n_val = int(0.7 * n), int(0.1 * n)
    prof = perm[:n_prof]
    val = perm[n_prof:n_prof + n_val]
    att = perm[n_prof + n_val:]

    print(f"\n{'=' * 78}")
    print(f"{args.npz}")
    print(f"  {n} traces, {T.shape[1]} samples -> train {len(prof)}, val {len(val)}, test {len(att)}")
    print(f"  device {device}, crop +/-{args.window} samples, {args.epochs} epochs")
    print(f"  every figure is against the majority-class baseline of the same split\n")

    c = args.coord
    targets = args.targets.split(",")
    rows = []
    for k in [int(x) for x in args.limbs.split(",")]:
        key = f"c{c}_limb{k}"
        if key not in d:
            continue
        lv = d[key][:n].astype(np.uint32)
        hw = d[f"{key}_hw"][:n].astype(np.int64)

        poi, r = find_poi(T[prof].astype(np.float64), hw[prof])
        lo = max(0, poi - args.window)
        hi = min(T.shape[1], poi + args.window)
        Xp, Xv, Xa = T[prof][:, lo:hi], T[val][:, lo:hi], T[att][:, lo:hi]
        print(f"  coord[{c}] limb {k}: POI {poi} (|corr(HW)|={r[poi]:.3f}), "
              f"crop [{lo},{hi}) = {hi - lo} samples")

        if "hw" in targets:
            t0 = time.time()
            acc = train_eval(Xp, hw[prof], Xv, hw[val], Xa, hw[att], 33, args.epochs,
                             args.lr, args.batch, device, args.seed)
            v, cnt = np.unique(hw[att], return_counts=True)
            base = cnt.max() / cnt.sum()
            se = np.sqrt(base * (1 - base) / len(att))
            z = (acc - base) / se
            print(f"    HW (33 classes)     acc={acc * 100:5.2f}%  base={base * 100:5.2f}%  "
                  f"z={z:+6.1f}  {'learned' if z > 3 else 'baseline'}  ({time.time() - t0:.0f}s)")
            rows.append(("hw", k, -1, acc, base, z))

        if "bit" in targets:
            for b in [int(x) for x in args.bits.split(",")]:
                bit = ((lv >> b) & 1).astype(np.int64)
                if bit.min() == bit.max():
                    continue
                acc = train_eval(Xp, bit[prof], Xv, bit[val], Xa, bit[att], 2, args.epochs,
                                 args.lr, args.batch, device, args.seed)
                ctl = train_eval(Xp, bit[prof], Xv, bit[val], Xa, bit[att], 2, args.epochs,
                                 args.lr, args.batch, device, args.seed, shuffle_labels=True)
                base = max(bit[att].mean(), 1 - bit[att].mean())
                se = np.sqrt(base * (1 - base) / len(att))
                z = (acc - base) / se
                zc = (ctl - base) / se
                print(f"    bit {b:2d}   acc={acc * 100:5.2f}%  base={base * 100:5.2f}%  z={z:+6.1f}"
                      f"   [shuffled {ctl * 100:5.2f}% z={zc:+5.1f}]  "
                      f"{'recovered' if z > 3 and zc < 3 else 'baseline'}")
                rows.append(("bit", k, b, acc, base, z))

        if "byte" in targets:
            byte = (lv & 0xFF).astype(np.int64)
            acc = train_eval(Xp, byte[prof], Xv, byte[val], Xa, byte[att], 256, args.epochs,
                             args.lr, args.batch, device, args.seed)
            v, cnt = np.unique(byte[att], return_counts=True)
            base = cnt.max() / cnt.sum()
            se = np.sqrt(base * (1 - base) / len(att))
            z = (acc - base) / se
            print(f"    byte 0 (256 classes) acc={acc * 100:5.2f}%  base={base * 100:5.2f}%  "
                  f"z={z:+6.1f}  {'recovered' if z > 3 else 'baseline'}")
            rows.append(("byte", k, 0, acc, base, z))
        print()

    print("  summary")
    for kind in ("hw", "bit", "byte"):
        sel = [r for r in rows if r[0] == kind]
        if not sel:
            continue
        ok = [r for r in sel if r[5] > 3]
        print(f"    {kind:<5} {len(ok)}/{len(sel)} significantly above baseline "
              f"(max z = {max((r[5] for r in sel), default=0):+.1f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
