#!/usr/bin/env python3
"""Figures for RESPONSE_EXPERIMENT.md (Alg 3.18 line 5).

Run from the package root after both captures exist:
  ~/.cwvenv/bin/python3 make_resp_figures.py
"""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False
OUT = "figures/"
d  = np.load('/tmp/sqisign_resp.npz'); T = d['traces']
dt = np.load('/tmp/sqisign_resp_tvla.npz'); Tt = dt['traces'].astype(np.float64); g = dt['group']
spc = float(d['samples_per_cycle'])

def welch(a, b):
    na, nb = len(a), len(b); ma, mb = a.mean(0), b.mean(0)
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    s = np.sqrt(va/na + vb/nb); s[s == 0] = np.finfo(float).tiny
    return (ma - mb) / s

# ---- figure 1: signal quality -------------------------------------------
mean, std = T.mean(0), T.std(0)
fig, ax = plt.subplots(3, 1, figsize=(13, 8))
w = slice(0, int(spc*20))
for i in range(4): ax[0].plot(T[i][w], lw=0.8, alpha=0.75)
ax[0].set_title("4 single traces, first 20 clock cycles -- they nearly "
                "superimpose, so the traces are well aligned")
ax[1].plot(mean[w], lw=1.1, color="#2a78d6")
for k in range(21): ax[1].axvline(k*spc, color="#e34948", lw=0.4, alpha=0.4)
ax[1].set_title(f"Mean of 20,000, same window -- the clock structure "
                f"resolves ({spc:.1f} samples/cycle)")
ax[2].plot(std, lw=0.4, color="#1baf7a")
ax[2].set_title("Per-sample standard deviation across traces -- structured, "
                "not flat, and the peaks are where the data shows")
ax[2].set_xlabel("sample")
for a in ax: a.grid(alpha=0.25)
fig.suptitle("Alg 3.18 line 5 formation -- signal quality "
             "(deterministic component is 197x the noise floor)",
             fontsize=12, y=1.002)
fig.tight_layout(); fig.savefig(OUT+"resp_signal_quality.png", dpi=130,
                                bbox_inches="tight"); plt.close(fig)

# ---- figure 2: TVLA, two panels -----------------------------------------
t_fr = welch(Tt[g == 1], Tt[g == 0])
f = np.flatnonzero(g == 1); h = len(f)//2   # same split the analyzer uses
t_ff = welch(Tt[f[:h]], Tt[f[h:2*h]])
fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
ax[0].plot(t_fr, lw=0.5, color="#e34948")
ax[0].set_title(f"Fixed vs random  (537 / 463)   max|t| = {np.abs(t_fr).max():.2f}"
                f"   --  {(np.abs(t_fr)>4.5).mean()*100:.1f}% of samples over threshold")
ax[1].plot(t_ff, lw=0.5, color="#eda100")
ax[1].set_title(f"Control: fixed vs fixed   max|t| = {np.abs(t_ff).max():.2f}  "
                f"-- must stay inside the dashed lines")
for a in ax:
    for s_ in (4.5, -4.5): a.axhline(s_, color="k", ls="--", lw=0.8)
    a.set_ylabel("t"); a.grid(alpha=0.25)
ax[1].set_xlabel("sample")
ax[0].set_ylim(-28, 28); ax[1].set_ylim(-28, 28)
fig.suptitle("Alg 3.18 line 5 (alpha_rsp) -- TVLA", fontsize=12)
fig.tight_layout(); fig.savefig(OUT+"response_tvla.png", dpi=130); plt.close(fig)
print("wrote resp_signal_quality.png, response_tvla.png")
