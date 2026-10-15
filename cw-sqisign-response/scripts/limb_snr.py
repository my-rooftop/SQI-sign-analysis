#!/usr/bin/env python3
"""
Limb-level SNR over the whole coord0 region (runs alongside DL training).

At the true 32-bit granularity (not hex nibbles). For each 32-bit limb of coord0,
computes partition SNR keyed on that limb's Hamming weight, over the fall-aligned,
lightly-pooled coord0 window. A clear SNR peak above the permutation null means
that limb's value leaks somewhere in the region.

Note on TVLA: a fixed-vs-random Welch t-test needs a FIXED input group, which
prof02 does not contain (every value is unique/random by design). So the
first-order detector here is SNR keyed on the limb HW, which is the appropriate
analog for an all-random corpus. (A proper TVLA would need a separate fixed-input
capture.)

  ~/.cwvenv/bin/python3 scripts/limb_snr.py /home/ugrp/sqisign_traces/prof02 \
      --coord 0 --limbs 10 --n 40000 --pool 4
"""
import argparse, json, time
import numpy as np


def snr(tr, lab):
    lab = np.asarray(lab); cl = np.unique(lab)
    m = np.empty((len(cl), tr.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = lab == c; w[i] = s.sum(); sub = tr[s]; m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--coord", type=int, default=0)
    ap.add_argument("--limbs", type=int, default=10)
    ap.add_argument("--window", type=int, default=66000)
    ap.add_argument("--pool", type=int, default=4)
    ap.add_argument("--n", type=int, default=40000)
    a = ap.parse_args()

    meta = json.load(open(f"{a.base}/meta.json")); L = meta["crop_len"]; pi = 2 + a.coord
    recs = [json.loads(l) for l in open(f"{a.base}/labels.jsonl")]
    T = np.load(f"{a.base}/traces.npy", mmap_mode="r")
    N = min(a.n, T.shape[0], len(recs))
    Lp = a.window // a.pool
    t0 = time.time()

    X = np.empty((N, Lp), dtype=np.float32)
    vals = np.empty(N, dtype=object)
    k = 0
    for i in range(N):
        f = recs[i]["fall"][pi]; s = f - a.window
        if s < 0 or f > L:
            continue
        seg = np.asarray(T[i, s:f], dtype=np.float32)
        X[k] = seg[:Lp * a.pool].reshape(Lp, a.pool).mean(1)
        vals[k] = abs(int(recs[i]["vals"][f"coord{a.coord}"], 16)); k += 1
    X = X[:k]; vals = vals[:k]
    print(f"coord{a.coord}: {k:,} traces, window {a.window} pooled {a.pool}x -> {Lp}, "
          f"{time.time()-t0:.0f}s to load", flush=True)

    rngp = np.random.default_rng(1)
    print(f"\n{'limb':5} {'HW SNR peak':12} {'null':8} {'ratio':6}  {'@sample(pooled)':16}")
    for lb in range(a.limbs):
        hw = np.array([bin((int(v) >> (32 * lb)) & 0xFFFFFFFF).count("1") for v in vals],
                      dtype=np.int16)
        if len(np.unique(hw)) < 3:
            print(f"limb{lb}: <3 HW classes (limb mostly zero), skip"); continue
        s = snr(X, hw); sn = snr(X, rngp.permutation(hw))
        peak = s.max(); npeak = sn.max(); j = int(np.argmax(s))
        flag = "  <-- above null" if peak > 3 * npeak and peak > 0.02 else ""
        print(f"limb{lb:<2} {peak:12.4f} {npeak:8.4f} {peak/(npeak+1e-9):6.1f}  "
              f"{j:>10}{flag}", flush=True)
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
