#!/usr/bin/env python3
"""
Analyse the fixed-key REAL-CORE dataset (collect_fixedkey.py, mode 2, command c).

This is the attack-target measurement (secret consumed in the core), unlike the
serialization proxy (prof02). Only ~2,100 traces, so a deep CNN would overfit; the
primary readout here is training-free SNR (limb Hamming weight), plus a small
nearest-centroid single-trace check. A CNN is optionally run with train/test both
reported so overfit is visible.

Window: the mode-2 pulse wraps materialization (basis*x, where coords are first
created -- the primary target) then normalization. We align at the trigger rise
(constant) and take a fixed window; the materialization sits at the start.

  ~/.cwvenv/bin/python3 scripts/analyze_fixedkey.py /home/ugrp/sqisign_traces/fixedkey01 \
      --coord 0 --limb 0 --pool 8
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


def load(base, coord, W, pool, nmax):
    meta = json.load(open(f"{base}/meta.json"))
    recs = [json.loads(l) for l in open(f"{base}/labels.jsonl")]
    T = np.load(f"{base}/traces.npy", mmap_mode="r")
    N = min(nmax or len(recs), T.shape[0], len(recs))
    Lp = W // pool
    X = np.empty((N, Lp), dtype=np.float32)
    vals = []
    k = 0
    for i in range(N):
        r = recs[i]["trig_rise"]
        if r + W > T.shape[1]:
            continue
        seg = np.asarray(T[i, r:r + W], dtype=np.float32)
        X[k] = seg[:Lp * pool].reshape(Lp, pool).mean(1)
        vals.append(abs(int(recs[i]["vals"][f"coord{coord}"], 16)))
        k += 1
    return meta, X[:k], vals, k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--coord", type=int, default=0)
    ap.add_argument("--limb", type=int, default=0)
    ap.add_argument("--window", type=int, default=2_180_000, help="samples from trig rise")
    ap.add_argument("--pool", type=int, default=8)
    ap.add_argument("--nmax", type=int, default=0)
    ap.add_argument("--cnn", action="store_true", help="also train a small CNN")
    ap.add_argument("--epochs", type=int, default=40)
    a = ap.parse_args()

    t0 = time.time()
    meta, X, vals, N = load(a.base, a.coord, a.window, a.pool, a.nmax)
    iv = meta["interval_ns"]
    print(f"{a.base}: {N:,} traces, window {a.window:,} pooled {a.pool}x -> {X.shape[1]:,}, "
          f"{iv:.1f} ns/sample, load {time.time()-t0:.0f}s", flush=True)

    def hw(v, lb): return bin((int(v) >> (32 * lb)) & 0xFFFFFFFF).count("1")
    y = np.array([hw(v, a.limb) for v in vals], dtype=np.int64)
    bl = np.array([int(v).bit_length() for v in vals])
    print(f"coord{a.coord} limb{a.limb} HW: mean {y.mean():.1f} std {y.std():.1f}; "
          f"corr(HW,bitlen)={np.corrcoef(y,bl)[0,1]:+.3f} (near 0 => data, not size)", flush=True)

    # --- training-free SNR (limb HW) ---
    s = snr(X, y)
    sn = snr(X, np.random.default_rng(1).permutation(y))
    j = int(np.argmax(s))
    print(f"\n[SNR] limb{a.limb} HW peak={s.max():.4f}  null={sn.max():.4f}  "
          f"ratio={s.max()/(sn.max()+1e-9):.1f}  @pooled {j} "
          f"(~{j*a.pool*iv/1e6:.2f} ms into finalize)", flush=True)

    # --- small single-trace nearest-centroid on top SNR points ---
    rng = np.random.default_rng(0); idx = rng.permutation(N)
    ntr = int(N * 0.85); tr, te = idx[:ntr], idx[ntr:]
    poi = np.argsort(snr(X[tr], y[tr]))[-40:]
    mu = X[tr][:, poi].mean(0); sd = X[tr][:, poi].std(0) + 1e-6
    Ztr = (X[tr][:, poi] - mu) / sd; Zte = (X[te][:, poi] - mu) / sd
    cls = np.unique(y[tr]); cent = np.stack([Ztr[y[tr] == c].mean(0) for c in cls])
    pred = cls[np.argmin(((Zte[:, None] - cent[None]) ** 2).sum(2), 1)]
    top1 = (pred == y[te]).mean(); within1 = (np.abs(pred - y[te]) <= 1).mean()
    corr = np.corrcoef(pred, y[te])[0, 1] if len(np.unique(pred)) > 1 else 0.0
    maj = np.max(np.bincount(y[te], minlength=33)) / len(te)
    print(f"[NC]  single-trace: top1={top1*100:.1f}%  within±1={within1*100:.1f}%  "
          f"corr={corr:+.3f}  majority={maj*100:.1f}%  chance=3.0%  (test n={len(te)})", flush=True)

    if a.cnn:
        cnn_run(X, y, a.epochs)
    print(f"\ntotal {time.time()-t0:.0f}s", flush=True)


def cnn_run(X, y, epochs):
    import torch, torch.nn as nn
    torch.manual_seed(0); torch.set_num_threads(12)
    N = len(X); rng = np.random.default_rng(0); perm = rng.permutation(N)
    ntr = int(N * 0.85); tr, te = perm[:ntr], perm[ntr:]
    Xtr, ytr, Xte, yte = X[tr], y[tr], X[te], y[te]
    mu, sd = Xtr.mean(), Xtr.std() + 1e-6; Xtr = (Xtr - mu) / sd; Xte = (Xte - mu) / sd

    class Net(nn.Module):
        def __init__(s, win):
            super().__init__()
            s.f = nn.Sequential(
                nn.Conv1d(1, 16, 11, padding=5), nn.BatchNorm1d(16), nn.ReLU(), nn.AvgPool1d(8),
                nn.Conv1d(16, 32, 11, padding=5), nn.BatchNorm1d(32), nn.ReLU(), nn.AvgPool1d(8),
                nn.Conv1d(32, 64, 7, padding=3), nn.BatchNorm1d(64), nn.ReLU(), nn.AdaptiveAvgPool1d(8))
            s.c = nn.Sequential(nn.Flatten(), nn.Dropout(0.4), nn.Linear(64 * 8, 33))
        def forward(s, x): return s.c(s.f(x))

    net = Net(X.shape[1]); opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-3)
    crit = nn.CrossEntropyLoss()
    Xt = torch.tensor(Xtr).unsqueeze(1); yt = torch.tensor(ytr); n = len(Xt)
    for ep in range(epochs):
        net.train(); pm = torch.randperm(n)
        for b in range(0, n, 128):
            i = pm[b:b+128]
            if len(i) < 8: continue
            loss = crit(net(Xt[i]), yt[i]); opt.zero_grad(); loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        ptr = net(Xt).argmax(1).numpy()
        pte = net(torch.tensor(Xte).unsqueeze(1)).argmax(1).numpy()
    tr_acc = (ptr == ytr).mean(); te_acc = (pte == yte).mean()
    w1 = (np.abs(pte - yte) <= 1).mean()
    corr = np.corrcoef(pte, yte)[0, 1] if len(np.unique(pte)) > 1 else 0.0
    print(f"[CNN] train_top1={tr_acc*100:.1f}%  test_top1={te_acc*100:.1f}%  "
          f"within±1={w1*100:.1f}%  corr={corr:+.3f}  "
          f"(train>>test => overfit from small n)", flush=True)


if __name__ == "__main__":
    main()
