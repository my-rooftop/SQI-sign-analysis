#!/usr/bin/env python3
"""
Stable CNN classifier for a coord limb's Hamming weight (replaces the diverged
triplet run). One trace = one sample; whole coord region -> one 32-bit limb HW.

Why a plain softmax CNN instead of triplet here: batch-hard triplet on 180k
CPU samples diverged (loss 0.40 -> 0.49, top1 ~ chance). Direct cross-entropy
classification is far more stable and gives the accuracy number directly. The
leakage question is unchanged; only the estimator is more robust.

  - Input : whole coord region, fall-aligned, average-pooled. One window/trace.
  - Label : Hamming weight (0..32, 33 classes) of one 32-bit limb.
  - Split : train / test by trace.
  - Extras: class-weighted loss (HW is binomial, very unbalanced), cosine LR,
            reports top1, within +/-1, macro-F1-ish balanced acc, vs majority.

  ~/.cwvenv/bin/python3 scripts/cnn_limb_hw.py /home/ugrp/sqisign_traces/prof02 \
      --coord 0 --limb 0 --n-train 180000 --n-test 20000 --epochs 30
"""
import argparse, json, time
import numpy as np


def build(base, coord, limb, n_total, W, pool, cache):
    import os
    tag = f"{cache}_c{coord}_l{limb}"
    if os.path.exists(tag + "_X.npy"):
        return np.load(tag + "_X.npy"), np.load(tag + "_y.npy")
    meta = json.load(open(f"{base}/meta.json")); L = meta["crop_len"]; pi = 2 + coord
    recs = [json.loads(l) for l in open(f"{base}/labels.jsonl")]
    T = np.load(f"{base}/traces.npy", mmap_mode="r")
    n_total = min(n_total, T.shape[0], len(recs)); Lp = W // pool
    X = np.empty((n_total, Lp), dtype=np.float32); y = np.empty(n_total, dtype=np.int64); k = 0
    for i in range(n_total):
        f = recs[i]["fall"][pi]; a = f - W
        if a < 0 or f > L:
            continue
        seg = np.asarray(T[i, a:f], dtype=np.float32)
        X[k] = seg[:Lp * pool].reshape(Lp, pool).mean(1)
        v = abs(int(recs[i]["vals"][f"coord{coord}"], 16))
        y[k] = bin((v >> (32 * limb)) & 0xFFFFFFFF).count("1"); k += 1
        if (i + 1) % 40000 == 0:
            print(f"  built {i+1}/{n_total}", flush=True)
    X, y = X[:k], y[:k]
    np.save(tag + "_X.npy", X); np.save(tag + "_y.npy", y)
    print(f"built X={X.shape} pooled {pool}x", flush=True)
    return X, y


def run(X, y, n_train, n_test, epochs, batch, lr, seed=0):
    import torch, torch.nn as nn
    torch.manual_seed(seed); torch.set_num_threads(12)
    N = len(X)
    if n_train + n_test > N:
        n_test = max(1, N // 10); n_train = N - n_test
    rng = np.random.default_rng(seed); perm = rng.permutation(N)
    tr, te = perm[:n_train], perm[n_train:n_train + n_test]
    Xtr, ytr, Xte, yte = X[tr], y[tr], X[te], y[te]
    mu, sd = Xtr.mean(), Xtr.std() + 1e-6
    Xtr = (Xtr - mu) / sd; Xte = (Xte - mu) / sd
    ncls = 33

    class Net(nn.Module):
        def __init__(s, win):
            super().__init__()
            s.f = nn.Sequential(
                nn.Conv1d(1, 32, 11, padding=5), nn.BatchNorm1d(32), nn.ReLU(), nn.AvgPool1d(4),
                nn.Conv1d(32, 64, 11, padding=5), nn.BatchNorm1d(64), nn.ReLU(), nn.AvgPool1d(4),
                nn.Conv1d(64, 128, 7, padding=3), nn.BatchNorm1d(128), nn.ReLU(), nn.AdaptiveAvgPool1d(8))
            s.c = nn.Sequential(nn.Flatten(), nn.Linear(128 * 8, 128), nn.ReLU(),
                                nn.Dropout(0.2), nn.Linear(128, ncls))
        def forward(s, x): return s.c(s.f(x))

    net = Net(X.shape[1])
    # class weights (HW ~ binomial, extremes rare)
    cnt = np.bincount(ytr, minlength=ncls).astype(np.float64)
    with np.errstate(divide="ignore"):
        w = np.where(cnt > 0, 1.0 / np.where(cnt > 0, cnt, 1), 0.0)
    w = w / w[cnt > 0].mean()
    crit = nn.CrossEntropyLoss(weight=torch.tensor(w, dtype=torch.float32))
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    Xtr_t = torch.tensor(Xtr).unsqueeze(1); ytr_t = torch.tensor(ytr); n = len(Xtr_t)
    for ep in range(epochs):
        net.train(); pm = torch.randperm(n); tot = 0.0; nb = 0
        for b in range(0, n, batch):
            idx = pm[b:b + batch]
            if len(idx) < 8: continue
            out = net(Xtr_t[idx]); loss = crit(out, ytr_t[idx])
            opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item(); nb += 1
        sched.step()
        if (ep + 1) % 5 == 0 or ep == 0:
            print(f"  epoch {ep+1}/{epochs} loss={tot/max(1,nb):.4f} lr={sched.get_last_lr()[0]:.1e}",
                  flush=True)

    net.eval()
    preds = []
    with torch.no_grad():
        for b in range(0, len(Xte), 8192):
            preds.append(net(torch.tensor(Xte[b:b+8192]).unsqueeze(1)).argmax(1).numpy())
    pred = np.concatenate(preds)
    top1 = (pred == yte).mean()
    within1 = (np.abs(pred - yte) <= 1).mean()
    majority = np.max(np.bincount(yte, minlength=ncls)) / len(yte)
    # balanced acc = mean per-class recall (robust to the binomial imbalance)
    bacc = np.mean([ (pred[yte==c]==c).mean() for c in np.unique(yte) if (yte==c).sum()>0 ])
    corr = np.corrcoef(pred, yte)[0, 1]
    return dict(top1=top1, within1=within1, majority=majority, bacc=bacc, corr=corr,
                n_train=n_train, n_test=len(yte))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--coord", type=int, default=0)
    ap.add_argument("--limb", type=int, default=0)
    ap.add_argument("--window", type=int, default=66000)
    ap.add_argument("--pool", type=int, default=16)
    ap.add_argument("--n-train", type=int, default=180000)
    ap.add_argument("--n-test", type=int, default=20000)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--cache", default="/home/ugrp/sqisign_traces/cnn_whole")
    a = ap.parse_args()
    t0 = time.time()
    X, y = build(a.base, a.coord, a.limb, a.n_train + a.n_test, a.window, a.pool, a.cache)
    print(f"class balance HW: min/max/mode = {np.bincount(y,minlength=33).min()}/"
          f"{np.bincount(y,minlength=33).max()}/{np.bincount(y).argmax()}", flush=True)
    r = run(X, y, a.n_train, a.n_test, a.epochs, a.batch, a.lr)
    print(f"\n== coord{a.coord} limb{a.limb} (32-bit) HW, single-trace CNN ==")
    print(f"   train {r['n_train']:,} / test {r['n_test']:,}")
    print(f"   top1={r['top1']*100:.2f}%  within±1={r['within1']*100:.1f}%  "
          f"balanced_acc={r['bacc']*100:.2f}%  majority={r['majority']*100:.2f}%  "
          f"chance=3.03%  corr(pred,true)={r['corr']:+.3f}")
    print(f"total {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
