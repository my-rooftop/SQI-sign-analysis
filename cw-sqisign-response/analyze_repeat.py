#!/usr/bin/env python3
"""
Averaging gain: with m repeats of the SAME coord0 value averaged, how does
  (a) DATA-leakage SNR (LS-limb Hamming weight), and
  (b) full-value identification accuracy (among the K known candidates)
improve?  m=1 == the real single-trace limit (SQIsign coords are never repeatable).
"""
import sys
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sqisign_repeat_avg"
d = np.load(BASE + ".npz", allow_pickle=True)
vidx = np.asarray(d["vidx"]); edges = np.asarray(d["edges"])
values = [int(v) for v in d["values"]]
K = int(d["K"]); R = int(d["R"]); interval = float(d["interval_ns"])
T = np.lib.format.open_memmap(BASE + "_traces.npy", mode="r")
Ntot = min(T.shape[0], len(vidx))
print(f"K={K} values, R={R} repeats, N={Ntot}, interval={interval}ns")


def snr(tr, lb):
    lb = np.asarray(lb); cl = np.unique(lb)
    m = np.empty((len(cl), tr.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = lb == c; w[i] = s.sum(); sub = tr[s]; m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


# coord0 = pulse index 2
p = 2
rise = edges[:Ntot, p]; fall = edges[:Ntot, 6 + p]
good = (rise > 0) & (fall > rise)
a0 = int(np.median(rise[good])); Lwin = int(np.median((fall - rise)[good]) * 0.95)
seg = np.empty((Ntot, Lwin), dtype=np.float32)
for i in range(Ntot):
    r = int(rise[i]) if rise[i] > 0 else a0
    seg[i] = T[i, r:r + Lwin]

# group trace-rows by value index
rows_by_val = [np.where(vidx == k)[0] for k in range(K)]
ls_hw = np.array([bin(values[k] & 0xFFFFFFFF).count("1") for k in range(K)])

rng = np.random.default_rng(0)
levels = [1, 2, 4, 8, 16, 32]
print(f"\n{'m':>3} {'dataSNR(LS-HW)':>15} {'null':>7} {'value-ID acc':>13} {'(chance=1/K)':>12}")

# templates for value-ID: use first half of repeats per value
half = R // 2
train_tmpl = np.stack([seg[rows_by_val[k][:half]].mean(0) for k in range(K)])  # (K,Lwin)
# POIs = most value-discriminating samples (variance of templates across values)
poi = np.argsort(train_tmpl.var(0))[-60:]
tm = train_tmpl[:, poi]; mu = tm.mean(0); sd = tm.std(0) + 1e-6
tmz = (tm - mu) / sd

for m in levels:
    # (a) SNR of LS-limb HW using m-averaged traces (floor(R/m) per value)
    avg_traces = []; avg_hw = []
    for k in range(K):
        rws = rows_by_val[k]
        nchunk = len(rws) // m
        for c in range(nchunk):
            avg_traces.append(seg[rws[c * m:(c + 1) * m]].mean(0))
            avg_hw.append(ls_hw[k])
    A = np.stack(avg_traces); Ahw = np.array(avg_hw)
    s = snr(A, Ahw); snull = snr(A, rng.permutation(Ahw)).max()

    # (b) value identification: average m of the held-out (second half) repeats, match template
    correct = 0; total = 0
    for k in range(K):
        rws = rows_by_val[k][half:]
        for _ in range(8):
            pick = rng.choice(rws, size=min(m, len(rws)), replace=False)
            q = seg[pick].mean(0)[poi]; qz = (q - mu) / sd
            pred = int(np.argmin(((tmz - qz) ** 2).sum(1)))
            correct += (pred == k); total += 1
    acc = correct / total
    print(f"{m:>3} {s.max():>15.3f} {snull:>7.3f} {acc*100:>11.1f}% {100.0/K:>11.1f}%")

print("\nm=1 row = real single-trace limit. Higher m needs replaying the same value.")
