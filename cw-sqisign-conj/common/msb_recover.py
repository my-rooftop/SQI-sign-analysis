#!/usr/bin/env python3
"""
How many most significant bits of alpha_sk,chl can be recovered from one trace of
Algorithm 3.18 line 2?

Context
-------
Full recovery of alpha_sk,chl from this operation is out of reach: the leakage is
the Hamming weight of each 32-bit limb and nothing finer, which caps the
extractable information at H(HW) = 3.55 bits per limb, about 64 of the 610 bits.
Partial recovery of the leading bits is a different and much more useful
question, because an algebraic reconstruction can take MSBs as input.

alpha is redrawn every signature, so the attack phase is single-trace; many
traces only buy a better profile. This is a profiled attack with templates built
on a disjoint split, and profiling happens on the same board, so read the numbers
as an upper bound.

What is measured
----------------
For each coordinate and each width T, the top T bits of the magnitude are taken
as a class label,

    MSB_T(v) = floor(|v| / 2^(L - T)),   L = the coordinate's magnitude bound,

so leading zeros count, which is the form a lattice reconstruction wants. A
Gaussian template per class is fitted on the profiling split and scored on
held-out traces. Reported per T: exact-match rate, the rank of the truth, and the
rate within the top few candidates, since an algebraic step can usually enumerate
a small list.

A guess is only worth calling recovered if it beats always answering the most
common class, so the majority-class baseline is reported alongside.

Usage:
  ~/.cwvenv/bin/python3 msb_recover.py --npz /tmp/sqisign_line2.npz
  ~/.cwvenv/bin/python3 msb_recover.py --npz ... --coords 1,2,3 --widths 1,2,3,4,6,8
"""

import argparse
import pathlib
import sys

import numpy as np


def corr_per_sample(T_c, T_n, model):
    m = np.asarray(model, dtype=np.float64)
    m = m - m.mean()
    mn = np.sqrt((m * m).sum())
    if mn == 0:
        return np.zeros(T_c.shape[1])
    return np.abs((T_c * m[:, None]).sum(axis=0) / (T_n * mn))


def pick_pois(T_c, T_n, labels, n_poi, min_gap):
    """Samples most correlated with the label, kept a clock cycle apart.

    Neighbouring samples inside one cycle carry nearly the same measurement, so
    an unguarded top-k would stack copies of one point and make the pooled
    covariance singular.
    """
    r = corr_per_sample(T_c, T_n, labels.astype(np.float64))
    chosen = []
    for idx in np.argsort(r)[::-1]:
        if all(abs(int(idx) - c) >= min_gap for c in chosen):
            chosen.append(int(idx))
            if len(chosen) == n_poi:
                break
    return np.array(sorted(chosen)), r


def fit_templates(X, y, classes, reg=1e-6):
    means, counts, resid = [], [], []
    for c in classes:
        sel = X[y == c]
        m = sel.mean(axis=0)
        means.append(m)
        counts.append(sel.shape[0])
        resid.append(sel - m)
    means = np.asarray(means)
    cov = np.atleast_2d(np.cov(np.vstack(resid), rowvar=False))
    cov += reg * np.trace(cov) / cov.shape[0] * np.eye(cov.shape[0])
    return means, np.linalg.inv(cov), np.asarray(counts, dtype=np.float64)


def log_post(X, means, prec, priors):
    out = np.empty((X.shape[0], means.shape[0]))
    for i, m in enumerate(means):
        dx = X - m
        out[:, i] = -0.5 * np.einsum("ij,jk,ik->i", dx, prec, dx)
    out += np.log(priors / priors.sum())[None, :]
    out -= out.max(axis=1, keepdims=True)
    return out


def msb_labels(mags, values, L, T):
    """Top T bits of |v| inside a fixed L-bit window, leading zeros included."""
    shift = max(L - T, 0)
    return np.array([int(abs(int(v)) >> shift) for v in values], dtype=np.int64)


def run_coord(T_all, values, mags, c, widths, n_poi, min_gap, prof_frac, min_class, seed):
    n = T_all.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_prof = int(prof_frac * n)
    prof, att = perm[:n_prof], perm[n_prof:]

    L = int(mags.max())
    print(f"\n  coord[{c}]: magnitude {int(mags.min())}-{L} bits, "
          f"window L={L}, profiling {len(prof)}, attack {len(att)}")
    print(f"    {'T':>3}{'classes':>9}{'exact':>9}{'baseline':>10}{'top-2':>8}"
          f"{'top-8':>8}{'rank':>8}{'verdict':>12}")

    rows = []
    for T in widths:
        y = msb_labels(mags, values, L, T)
        yp, ya = y[prof], y[att]

        Xp_full = T_all[prof]
        mu = Xp_full.mean(axis=0, keepdims=True)
        Xc = Xp_full - mu
        Xn = np.sqrt((Xc * Xc).sum(axis=0))
        Xn = np.where(Xn == 0, np.finfo(np.float64).tiny, Xn)
        pois, _ = pick_pois(Xc, Xn, yp, n_poi, min_gap)

        Xp = Xp_full[:, pois]
        Xa = T_all[att][:, pois]

        classes = np.array([c_ for c_ in np.unique(yp) if (yp == c_).sum() >= min_class])
        if classes.size < 2:
            print(f"    {T:>3}{'-':>9}  too few classes to measure")
            continue
        keep = np.isin(yp, classes)
        means, prec, counts = fit_templates(Xp[keep], yp[keep], classes)

        lp = log_post(Xa, means, prec, counts)
        order = np.argsort(lp, axis=1)[:, ::-1]
        pred = classes[order[:, 0]]

        # traces whose true class was never profiled count as failures, not as
        # samples to be quietly dropped
        in_set = np.isin(ya, classes)
        exact = float((pred == ya).mean())
        ranks = np.full(len(ya), classes.size + 1, dtype=int)
        for i in np.flatnonzero(in_set):
            ranks[i] = int(np.flatnonzero(classes[order[i]] == ya[i])[0]) + 1
        top2 = float((ranks <= 2).mean())
        top8 = float((ranks <= 8).mean())

        vals, cnts = np.unique(ya, return_counts=True)
        base = float(cnts.max() / cnts.sum())
        verdict = "recovered" if exact > base + 3 * np.sqrt(base * (1 - base) / len(ya)) else "at baseline"

        print(f"    {T:>3}{classes.size:>9}{exact * 100:>8.1f}%{base * 100:>9.1f}%"
              f"{top2 * 100:>7.1f}%{top8 * 100:>7.1f}%{np.mean(ranks):>8.1f}   {verdict}")
        rows.append({"coord": c, "T": T, "exact": exact, "base": base,
                     "top2": top2, "top8": top8, "classes": int(classes.size)})
    return rows


def main():
    p = argparse.ArgumentParser(description="MSB recovery for Alg 3.18 line 2")
    p.add_argument("--npz", required=True)
    p.add_argument("--coords", default="0,1,2,3")
    p.add_argument("--widths", default="1,2,3,4,5,6,8")
    p.add_argument("--n-poi", type=int, default=16)
    p.add_argument("--min-gap", type=int, default=34)
    p.add_argument("--prof-frac", type=float, default=0.8)
    p.add_argument("--min-class", type=int, default=20)
    p.add_argument("--traces", type=int, default=None)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    d = np.load(args.npz, allow_pickle=False)
    T_all = d["traces"].astype(np.float64)
    if args.traces:
        T_all = T_all[:args.traces]
    n = T_all.shape[0]

    print(f"\n{'=' * 78}")
    print(f"{args.npz}")
    print(f"  {n} traces, {T_all.shape[1]} samples, "
          f"POIs {args.n_poi} at >= {args.min_gap} samples apart")
    print(f"  profiled, single-trace attack phase; same-board profiling, so an upper bound")

    widths = [int(x) for x in args.widths.split(",")]
    rows = []
    for c in [int(x) for x in args.coords.split(",")]:
        values = [int(h, 16) for h in d[f"c{c}_value_hex"][:n]]
        # stored as a 416-bit two's-complement pattern; recover the signed value
        values = [v - (1 << 416) if (v >> 415) & 1 else v for v in values]
        mags = d[f"c{c}_magnitude_bits"][:n]
        rows += run_coord(T_all, values, mags, c, widths, args.n_poi, args.min_gap,
                          args.prof_frac, args.min_class, args.seed)

    if rows:
        print(f"\n  summary: largest T judged recovered, per coordinate")
        for c in sorted({r["coord"] for r in rows}):
            ok = [r["T"] for r in rows if r["coord"] == c
                  and r["exact"] > r["base"] + 3 * np.sqrt(r["base"] * (1 - r["base"]) / max(n // 5, 1))]
            print(f"    coord[{c}]: {max(ok) if ok else 0} bits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
