#!/usr/bin/env python3
"""
Train per-coordinate templates on the large known-value profiling set and measure
single-trace recovery (train/test split) at scale. Targets per coord:
  - bit-length (overall magnitude)
  - most-significant 32-bit limb: Hamming weight (33 cls) and top byte value (256 cls)
mpz_get_str emits hex MS-first, so the MS limb sits at a fixed offset from each
coord's pulse rise -> well aligned across variable-length values.
"""
import json, sys
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sqisign_profile_known"
meta = np.load(BASE + ".npz", allow_pickle=True)
recs = json.loads(str(meta["records_json"]))
edges = np.asarray(meta["edges"]); N = int(meta["n"]); S = int(meta["samples"])
interval_ns = float(meta["interval_ns"])
T = np.lib.format.open_memmap(BASE + "_traces.npy", mode="r")
N = min(N, T.shape[0], len(recs), len(edges))
print(f"N={N} S={S} interval={interval_ns}ns")


def snr(tr, lb):
    lb = np.asarray(lb); cl = np.unique(lb)
    m = np.empty((len(cl), tr.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = lb == c; w[i] = s.sum(); sub = tr[s]; m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


def recover(X, y, npoi=40, frac=0.7, seed=0):
    rng = np.random.default_rng(seed); idx = rng.permutation(len(y))
    ntr = int(len(y) * frac); tr, te = idx[:ntr], idx[ntr:]
    s = snr(X[tr], y[tr]); poi = np.argsort(s)[-npoi:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    acc = (pred == y[te]).mean()
    # majority-class baseline on test
    base = np.max(np.bincount(y[te] - y[te].min())) / len(y[te])
    return acc, base, s.max(), len(cls)


def mslimb(v, bitlen):
    v = abs(int(v)); topbit = bitlen  # number of bits
    sh = max(0, topbit - 32)
    return (v >> sh) & 0xFFFFFFFF


print(f"\n{'coord':7} {'bitlen':9} {'SNRbl':6} {'1tr_bl':7}| "
      f"{'SNRmsHW':7} {'msHW%':6} {'base':5}| {'msTopByte%':10} {'base':5}")
for ci, ck in enumerate(["coord0", "coord1", "coord2", "coord3"]):
    p = 2 + ci
    rise = edges[:N, p]; fall = edges[:N, 6 + p]
    good = (rise > 0) & (fall > rise)
    gi = np.where(good)[0]
    vals = np.array([int(recs[i][ck], 16) for i in gi])
    bl = np.array([abs(int(x)).bit_length() for x in vals])
    ms = np.array([mslimb(vals[k], bl[k]) for k in range(len(gi))], dtype=np.uint32)
    ms_hw = np.array([bin(int(x)).count("1") for x in ms], dtype=np.int16)
    ms_topbyte = ((ms >> 24) & 0xFF).astype(np.int16)

    width = (fall - rise)[gi]
    Lwin = max(2000, min(int(np.median(width) * 0.95), 200000))
    seg = np.empty((len(gi), Lwin), dtype=np.float32)
    for r, i in enumerate(gi):
        seg[r] = T[i, int(rise[i]):int(rise[i]) + Lwin]

    acc_bl, _, snr_bl, ncl_bl = recover(seg, bl)
    acc_hw, base_hw, snr_hw, _ = recover(seg, ms_hw)
    acc_tb, base_tb, _, _ = recover(seg, ms_topbyte)
    print(f"{ck:7} {bl.min()}-{bl.max():<4} {snr_bl:6.2f} {acc_bl*100:5.0f}% ({100/ncl_bl:2.0f})| "
          f"{snr_hw:7.2f} {acc_hw*100:5.0f}% {base_hw*100:4.0f}%| "
          f"{acc_tb*100:8.0f}%  {base_tb*100:4.0f}%")

print("\n(1tr_bl parenthesis = random baseline; msHW base / msTopByte base = majority-class)")
print("done.")
