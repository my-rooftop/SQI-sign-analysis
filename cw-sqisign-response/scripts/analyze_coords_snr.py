#!/usr/bin/env python3
"""
Per-coordinate SNR on a collect_coords corpus (traces.npy + labels.jsonl).

Two leakage channels, kept separate:
  SIZE  -- pulse width vs bit-length. Pure timing; the strong, already-known leak.
  DATA  -- LS 32-bit limb Hamming weight, size-independent. This is the number
           that decides whether coordinate VALUES (not just their size) leak.

Alignment: coordinates are variable-length and ~50% negative (leading '-'), so a
rise-aligned window smears. The LS hex digits sit just before the FALL edge
regardless of sign or length, so DATA analysis aligns on the fall edge.

Safe to run while collection is in progress: opens traces.npy read-only and uses
only --n rows (default: well under the current label count).

  ~/.cwvenv/bin/python3 scripts/analyze_coords_snr.py /home/ugrp/sqisign_traces/prof02 --n 140000
"""
import argparse, json, sys
import numpy as np


def snr(traces, labels):
    labels = np.asarray(labels)
    cl = np.unique(labels)
    m = np.empty((len(cl), traces.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = labels == c; w[i] = s.sum(); sub = traces[s]
        m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


def nc_cv(X, y, frac=0.7, npoi=40, seed=0):
    rng = np.random.default_rng(seed); idx = rng.permutation(len(y))
    ntr = int(len(y) * frac); tr, te = idx[:ntr], idx[ntr:]
    poi = np.argsort(snr(X[tr], y[tr]))[-npoi:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    return (pred == y[te]).mean(), len(cls)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--n", type=int, default=0, help="rows to consider (0 = all safe)")
    ap.add_argument("--subset", type=int, default=30000, help="traces used for SNR est.")
    ap.add_argument("--ls-window", type=int, default=12000, help="samples before fall edge")
    a = ap.parse_args()

    meta = json.load(open(f"{a.base}/meta.json"))
    L = meta["crop_len"]; iv = meta["interval_ns"]
    recs = [json.loads(l) for l in open(f"{a.base}/labels.jsonl")]
    navail = len(recs)
    N = a.n if a.n else max(0, navail - 500)      # stay off the freshest rows
    N = min(N, navail)
    T = np.load(f"{a.base}/traces.npy", mmap_mode="r")
    N = min(N, T.shape[0])
    print(f"{a.base}: {navail:,} labels available, using {N:,} rows, "
          f"{iv:.1f} ns/sample, source={meta.get('source')}")

    rng = np.random.default_rng(0)
    use = rng.choice(N, size=min(a.subset, N), replace=False)
    use.sort()
    W = a.ls_window

    print(f"\n{'coord':6} {'bitlen':10} | {'SIZE: w~bl r':12} {'1tr bl%':8} | "
          f"{'DATA: LShw SNR':14} {'null':6} {'1tr LShw%':10} {'base%':6}")
    for c in range(4):
        pi = 2 + c
        rise = np.array([recs[i]["rise"][pi] for i in use])
        fall = np.array([recs[i]["fall"][pi] for i in use])
        vals = np.array([int(recs[i]["vals"][f"coord{c}"], 16) for i in use])  # signed
        av = np.abs(vals)
        bl = np.array([int(x).bit_length() for x in av])
        lshw = np.array([bin(int(x) & 0xFFFFFFFF).count("1") for x in av], dtype=np.int16)
        width = fall - rise

        good = (rise >= 0) & (fall > rise) & (fall <= L) & (fall - W >= 0)
        gi = np.where(good)[0]

        # SIZE channel: width vs bit-length
        wc = np.corrcoef(width[gi], bl[gi])[0, 1]
        acc_bl, _ = nc_cv(width[gi].reshape(-1, 1).astype(np.float32), bl[gi], npoi=1)

        # DATA channel: fall-aligned LS window
        seg = np.empty((len(gi), W), dtype=np.float32)
        for r, i in enumerate(gi):
            f = int(fall[i]); seg[r] = T[use[i], f - W:f]
        s_hw = snr(seg, lshw[gi]).max()
        s_null = snr(seg, np.random.default_rng(1).permutation(lshw[gi])).max()
        acc_hw, _ = nc_cv(seg, lshw[gi])
        base_hw = np.max(np.bincount(lshw[gi] - lshw[gi].min())) / len(gi)

        print(f"coord{c:<1} {bl.min()}-{bl.max():<5} | {wc:+11.3f} {acc_bl*100:7.0f}% | "
              f"{s_hw:14.3f} {s_null:6.3f} {acc_hw*100:9.1f}% {base_hw*100:5.1f}%")

    print(f"\nSIZE 1tr bl% = single-trace bit-length from pulse width.")
    print(f"DATA LShw SNR vs null: >null by a clear margin => LS-limb VALUE leaks; "
          f"~null => data at noise floor (only size leaks).")


if __name__ == "__main__":
    main()
