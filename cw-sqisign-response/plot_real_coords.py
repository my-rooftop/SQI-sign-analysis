#!/usr/bin/env python3
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/tmp/sqisign_real_coords"
meta = np.load(BASE + ".npz", allow_pickle=True)
recs = json.loads(str(meta["records_json"]))
edges = np.asarray(meta["edges"]); N = int(meta["n"])
T = np.lib.format.open_memmap(BASE + "_traces.npy", mode="r")
N = min(N, T.shape[0], len(recs), len(edges))

def snr(tr, lb):
    lb = np.asarray(lb); cl = np.unique(lb)
    m = np.empty((len(cl), tr.shape[1])); v = np.empty_like(m); w = np.empty(len(cl))
    for i, c in enumerate(cl):
        s = lb == c; w[i] = s.sum(); sub = tr[s]; m[i] = sub.mean(0); v[i] = sub.var(0)
    w /= w.sum(); g = (w[:, None] * m).sum(0)
    return (w[:, None] * (m - g) ** 2).sum(0) / ((w[:, None] * v).sum(0) + 1e-12)

def nc_cv(X, y, npoi=25, frac=0.7, seed=0):
    rng = np.random.default_rng(seed); idx = rng.permutation(len(y))
    ntr = int(len(y) * frac); tr, te = idx[:ntr], idx[ntr:]
    s = snr(X[tr], y[tr]); poi = np.argsort(s)[-npoi:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    return (pred == y[te]).mean(), 1.0 / len(cls)

fig, ax = plt.subplots(4, 2, figsize=(15, 13))
summary = []
for ci, ck in enumerate(["coord0", "coord1", "coord2", "coord3"]):
    p = 2 + ci; rise = edges[:N, p]; fall = edges[:N, 6 + p]
    good = (rise > 0) & (fall > rise)
    gi = np.where(good)[0]
    vals = np.array([int(recs[i][ck], 16) for i in gi])
    bl = np.array([abs(int(x)).bit_length() for x in vals])
    lshw = np.array([bin(abs(int(x)) & 0xFFFFFFFF).count("1") for x in vals], dtype=np.int16)
    width = (fall - rise)[gi]
    Lwin = max(2000, min(int(np.median(width) * 0.9), 120000))
    seg = np.empty((len(gi), Lwin), dtype=np.float32)
    for r, i in enumerate(gi):
        seg[r] = T[i, int(rise[i]):int(rise[i]) + Lwin]

    s_bl = snr(seg, bl); s_hw = snr(seg, lshw)
    # multi-perm null
    rng = np.random.default_rng(7)
    nulls = [snr(seg, rng.permutation(lshw)).max() for _ in range(8)]
    null_hw = float(np.mean(nulls))
    acc_bl, rb_bl = nc_cv(seg, bl); acc_hw, rb_hw = nc_cv(seg, lshw)
    summary.append((ck, len(gi), bl.min(), bl.max(), s_bl.max(), s_hw.max(), null_hw,
                    acc_bl, rb_bl, acc_hw, rb_hw))

    us = 4.0 / 1e3  # us per sample
    x = np.arange(Lwin) * us
    ax[ci, 0].plot(x, s_bl, lw=0.5, color="C3", label=f"SNR(bit-length)  peak={s_bl.max():.2f}")
    ax[ci, 0].plot(x, s_hw, lw=0.5, color="C0", alpha=0.8, label=f"SNR(LS-limb HW)  peak={s_hw.max():.2f}")
    ax[ci, 0].axhline(null_hw, color="gray", ls="--", lw=0.8, label=f"HW null={null_hw:.2f}")
    ax[ci, 0].set_title(f"{ck}: power SNR over its serialization pulse"); ax[ci, 0].set_ylabel("SNR")
    ax[ci, 0].legend(fontsize=7); ax[ci, 0].set_xlabel("us from pulse start")

    ax[ci, 1].scatter(width * us, bl, s=6, alpha=0.4)
    ax[ci, 1].set_title(f"{ck}: pulse width vs bit-length  (corr={np.corrcoef(width,bl)[0,1]:+.2f})")
    ax[ci, 1].set_xlabel("pulse width (us)"); ax[ci, 1].set_ylabel("bit-length")

fig.tight_layout(); fig.savefig("/tmp/sqisign_real_coords_analysis.png", dpi=120); plt.close(fig)

print(f"{'coord':7} {'n':4} {'bitlen':10} {'SNRbl':6} {'SNRhw':6} {'null':6} "
      f"{'1tr_bl':8} {'rand':5} {'1tr_hw':7} {'base':5}")
for (ck, n, blo, bhi, sbl, shw, nul, abl, rbl, ahw, rhw) in summary:
    print(f"{ck:7} {n:4} {blo}-{bhi:<6} {sbl:6.3f} {shw:6.3f} {nul:6.3f} "
          f"{abl*100:6.0f}%  {rbl*100:4.0f}% {ahw*100:5.0f}%  {rhw*100:4.0f}%")
print("\nsaved /tmp/sqisign_real_coords_analysis.png")
