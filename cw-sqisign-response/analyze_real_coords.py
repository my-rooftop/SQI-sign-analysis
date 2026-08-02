#!/usr/bin/env python3
"""
Analyse real-response serialization traces (firmware mode 6, per-value pulses).
Answers: how much of each real coord0..3 is recoverable, single-trace?

Pulses per trace (Ch B): [content, denom, coord0, coord1, coord2, coord3].
For coord i we use pulse index 2+i (rise/fall edges stored in `edges`).

Metrics per coord:
  - pulse WIDTH vs bit-length  (direct timing leak of the coordinate's size)
  - power SNR keyed on bit-length (size) and on LS-limb Hamming weight (data, size-independent)
  - single-trace recovery accuracy (template / nearest-centroid, train/test split) for
    bit-length and LS-limb HW, with a permutation-null baseline.
"""
import json, sys
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sqisign_real_coords"
meta = np.load(BASE + ".npz", allow_pickle=True)
recs = json.loads(str(meta["records_json"]))
edges = np.asarray(meta["edges"])            # (N,12): [rise0..5, fall0..5]
interval_ns = float(meta["interval_ns"])
N = int(meta["n"]); S = int(meta["samples"])
T = np.lib.format.open_memmap(BASE + "_traces.npy", mode="r")
N = min(N, T.shape[0], len(recs), len(edges))
print(f"N={N} S={S} interval={interval_ns}ns")


def snr(traces, labels):
    labels = np.asarray(labels); cl = np.unique(labels)
    m = np.empty((len(cl), traces.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        sel = labels == c; w[i] = sel.sum(); sub = traces[sel]
        m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)


def nc_cv(X, y, ntr, npoi=15, seed=0):
    """nearest-centroid single-trace recovery, train/test split; returns test accuracy."""
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y)); tr, te = idx[:ntr], idx[ntr:]
    s = snr(X[tr], y[tr]); poi = np.argsort(s)[-npoi:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    return (pred == y[te]).mean(), len(te), len(cls)


def limbs(v, k):
    v = abs(int(v)); return [(v >> (32 * j)) & 0xFFFFFFFF for j in range(k)]


coord_keys = ["coord0", "coord1", "coord2", "coord3"]
print("\n=== per-coordinate recovery (real responses) ===")
for ci, ck in enumerate(coord_keys):
    p = 2 + ci                                   # pulse index for this coord
    rise = edges[:N, p]; fall = edges[:N, 6 + p]
    ok = (rise > 0) & (fall > rise)
    vals = np.array([int(recs[i][ck], 16) for i in range(N)])
    bitlen = np.array([abs(int(x)).bit_length() for x in vals])
    width = (fall - rise)                         # samples the pulse is high

    # (1) timing: pulse width vs bit-length
    good = ok & (bitlen > 0)
    cw = np.corrcoef(width[good], bitlen[good])[0, 1]
    # single-trace bit-length from width alone (round nearest observed mean per class)
    # nearest-centroid on 1-D width
    accw, nte, ncl = nc_cv(width[good].reshape(-1, 1).astype(np.float32), bitlen[good],
                           ntr=int(good.sum() * 0.7), npoi=1)

    # (2) power leakage: align segment at pulse rise, fixed window
    Lwin = int(np.median(width[good]) * 0.9)
    Lwin = max(2000, min(Lwin, 120000))
    seg = np.empty((good.sum(), Lwin), dtype=np.float32)
    gi = np.where(good)[0]
    for r, i in enumerate(gi):
        a = int(rise[i]); seg[r] = T[i, a:a + Lwin]
    bl = bitlen[gi]
    lslimb_hw = np.array([bin(limbs(vals[i], 1)[0]).count("1") for i in gi], dtype=np.int16)

    snr_bl = snr(seg, bl); snr_hw = snr(seg, lslimb_hw)
    # null for LS-limb HW
    rng = np.random.default_rng(3)
    snr_hw_null = snr(seg, rng.permutation(lslimb_hw)).max()

    ntr = int(len(gi) * 0.7)
    acc_bl, nte2, ncl2 = nc_cv(seg, bl, ntr)
    acc_hw, _, nclhw = nc_cv(seg, lslimb_hw, ntr)
    # random baselines
    rb_bl = 1.0 / ncl2; rb_hw = np.max(np.bincount(lslimb_hw - lslimb_hw.min())) / len(lslimb_hw)

    print(f"\n{ck}: n={good.sum()} bitlen {bl.min()}-{bl.max()} ({ncl2} distinct)")
    print(f"  [size]  width~bitlen corr={cw:+.3f}   single-trace bitlen from WIDTH: {accw*100:.0f}% ({ncl} cls)")
    print(f"  [power] SNR bitlen peak={snr_bl.max():.3f}   SNR LSlimb-HW peak={snr_hw.max():.3f} (null {snr_hw_null:.3f})")
    print(f"  [1-trace] bitlen(power)={acc_bl*100:.0f}% (rand {rb_bl*100:.0f}%)   "
          f"LSlimb-HW={acc_hw*100:.0f}% (baseline {rb_hw*100:.0f}%)")

print("\ndone.")
