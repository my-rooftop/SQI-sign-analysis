#!/usr/bin/env python3
"""
Can the single-limb HW recovery extend to all 10 limbs (the whole 320-bit coord0)?

Reuses the cached whole-coord0 windows (from cnn_limb_hw.py) and runs a fast
SNR-POI + nearest-centroid single-trace template per limb. NC is a lower bound
vs the CNN (CNN gave corr 0.94 on limb0), but it is seconds/limb instead of hours,
so it sweeps all 10 limbs cheaply. Reports per-limb top1 / within±1 / corr and the
size-confound check corr(HW, bit-length).

  ~/.cwvenv/bin/python3 scripts/sweep_all_limbs.py \
      --cache /home/ugrp/sqisign_traces/coord0_whole_c0_l0 \
      --labels /home/ugrp/sqisign_traces/prof02/labels.jsonl --coord 0
"""
import argparse, json
import numpy as np


def snr(tr, lab):
    lab = np.asarray(lab); cl = np.unique(lab)
    m = np.empty((len(cl), tr.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = lab == c; w[i] = s.sum(); sub = tr[s]; m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


def nc(X, y, ntr, npoi=40, seed=0):
    rng = np.random.default_rng(seed); idx = rng.permutation(len(y))
    tr, te = idx[:ntr], idx[ntr:]
    poi = np.argsort(snr(X[tr], y[tr]))[-npoi:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    yte = y[te]
    top1 = (pred == yte).mean(); w1 = (np.abs(pred - yte) <= 1).mean()
    corr = np.corrcoef(pred, yte)[0, 1] if len(np.unique(pred)) > 1 else 0.0
    maj = np.max(np.bincount(yte, minlength=33)) / len(yte)
    return top1, w1, corr, maj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, help="prefix; loads <cache>_X.npy")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--coord", type=int, default=0)
    ap.add_argument("--n", type=int, default=120000)
    ap.add_argument("--limbs", type=int, default=10)
    a = ap.parse_args()

    X = np.load(a.cache + "_X.npy", mmap_mode="r")
    recs = [json.loads(l) for l in open(a.labels)]
    N = min(a.n, X.shape[0], len(recs))
    X = np.asarray(X[:N], dtype=np.float32)
    vals = [abs(int(recs[i]["vals"][f"coord{a.coord}"], 16)) for i in range(N)]
    bl = np.array([v.bit_length() for v in vals])
    print(f"coord{a.coord}: {N:,} traces, {X.shape[1]} features (whole-coord window)\n")
    print(f"{'limb':5} {'bits':9} {'top1%':7} {'within1%':9} {'corr':7} {'maj%':6} "
          f"{'corr(HW,size)':13}")
    ntr = int(N * 0.7)
    for lb in range(a.limbs):
        y = np.array([bin((v >> (32 * lb)) & 0xFFFFFFFF).count("1") for v in vals], dtype=np.int64)
        if len(np.unique(y)) < 5:
            print(f"limb{lb:<2} (limb mostly zero / too few HW classes, skip)"); continue
        top1, w1, corr, maj = nc(X, y, ntr)
        csize = np.corrcoef(y, bl)[0, 1]
        lo, hi = 32 * lb, 32 * lb + 31
        print(f"limb{lb:<2} {lo}-{hi:<4} {top1*100:6.1f} {w1*100:8.1f} {corr:+6.3f} "
              f"{maj*100:5.1f} {csize:+12.3f}")
    print("\ncorr = corr(predicted HW, true HW), single-trace, held-out test.")
    print("corr(HW,size) near 0 => genuine data leakage; large => contaminated by the size channel.")
    print("NC is a lower bound; the CNN gave corr 0.94 on limb0.")


if __name__ == "__main__":
    main()
