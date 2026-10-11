#!/usr/bin/env python3
"""
Step 5: profiled single-trace recovery of alpha_sk,chl from Algorithm 3.18 line 2.

Why a template attack and not CPA
---------------------------------
alpha_sk,chl is redrawn for every signature (Alg 3.18 line 1 samples it), so there
is no fixed secret to average over. Many traces only buy a better profile; the
attack itself has to work on the one trace belonging to the signature under
attack. That is a profiled template attack with a single-trace attack phase, not
CPA.

This profiles and attacks on the same board with known ground truth, which is the
easiest case an attacker could face -- treat the numbers as an upper bound, not
as what a cross-device attack would achieve.

What it measures
----------------
  HW templates   Per limb, Gaussian templates over the Hamming weight of the
                 limb, built on the profiling split and scored on held-out
                 traces. Reports exact-match accuracy, the rank of the true
                 class, and the information actually gained, H(HW) - H(HW|trace),
                 estimated from the posteriors.

  bit LDA        Per limb and bit position, a linear discriminant on the same
                 points of interest, scored on held-out traces. This is the
                 decisive test: HW alone caps what is recoverable, so full
                 recovery needs the leakage to separate values inside one HW
                 class. Accuracy is compared against the majority-class baseline,
                 not against 50%, because individual bits of a limb are balanced
                 but finite samples are not.

  budget         Bits of the secret the measured leakage actually removes,
                 against the ceiling perfect HW knowledge would give. The secret
                 size is summed from the capture's own magnitudes: about 610 bits
                 for alpha_sk,chl at line 2, about 746 for alpha_rsp at line 5.

Usage:
  ~/.cwvenv/bin/python3 recover_conj_template.py --npz /tmp/sqisign_conj_prof.npz
  ~/.cwvenv/bin/python3 recover_conj_template.py --npz ... --curve 100,200,400,800
"""

import argparse
import json
import pathlib
import sys

import numpy as np

def clean_limbs(d, c, n):
    """Limbs of coordinate c that lie entirely below the magnitude.

    Derived from the data rather than hardcoded: alpha_sk,chl at line 2 is about
    234/234/71/71 bits wide, alpha_rsp at line 5 about 266/267/104/103, so the
    counts differ per target.
    """
    key = f"c{c}_magnitude_bits"
    if key not in d:
        return 0
    m = int(np.median(d[key][:n]))
    return max(0, m // 32)


# ----------------------------------------------------------------- utilities

def corr_per_sample(traces, values):
    v = np.asarray(values, dtype=np.float64)
    v = v - v.mean()
    vn = np.sqrt((v * v).sum())
    if vn == 0:
        return np.zeros(traces.shape[1])
    t = traces - traces.mean(axis=0, keepdims=True)
    tn = np.sqrt((t * t).sum(axis=0))
    tn = np.where(tn == 0, np.finfo(np.float64).tiny, tn)
    return (t * v[:, None]).sum(axis=0) / (tn * vn)


def pick_pois(traces, labels, n_poi, min_gap):
    """Top |corr| samples against `labels`, keeping them min_gap apart.

    Neighbouring samples inside one clock cycle are nearly the same measurement,
    so taking the raw top-k would stack a dozen copies of one point and make the
    template covariance singular.
    """
    r = np.abs(corr_per_sample(traces, labels))
    order = np.argsort(r)[::-1]
    chosen = []
    for idx in order:
        if all(abs(idx - c) >= min_gap for c in chosen):
            chosen.append(int(idx))
            if len(chosen) == n_poi:
                break
    return np.array(sorted(chosen)), r


def gaussian_templates(X, y, classes, reg=1e-6):
    """Class means plus one pooled covariance (linear discriminant)."""
    means, counts = [], []
    resid = []
    for c in classes:
        sel = X[y == c]
        m = sel.mean(axis=0)
        means.append(m)
        counts.append(sel.shape[0])
        resid.append(sel - m)
    means = np.asarray(means)
    R = np.vstack(resid)
    cov = np.cov(R, rowvar=False)
    cov = np.atleast_2d(cov)
    cov += reg * np.trace(cov) / cov.shape[0] * np.eye(cov.shape[0])
    return means, np.linalg.inv(cov), np.asarray(counts, dtype=np.float64)


def log_posteriors(X, means, prec, priors):
    """Log posterior over classes for each row of X, up to a constant."""
    out = np.empty((X.shape[0], means.shape[0]))
    for i, m in enumerate(means):
        dx = X - m
        out[:, i] = -0.5 * np.einsum("ij,jk,ik->i", dx, prec, dx)
    out += np.log(priors / priors.sum())[None, :]
    out -= out.max(axis=1, keepdims=True)
    return out


def entropy_bits(p):
    p = np.clip(p, 1e-300, 1.0)
    return float(-(p * np.log2(p)).sum())


# ------------------------------------------------------------- HW templates

def attack_limb_hw(traces, hw, prof_idx, att_idx, n_poi, min_gap):
    Xp_all = traces[prof_idx]
    yp = hw[prof_idx]
    pois, _ = pick_pois(Xp_all, yp, n_poi, min_gap)
    Xp = Xp_all[:, pois]
    Xa = traces[att_idx][:, pois]
    ya = hw[att_idx]

    classes = np.unique(yp)
    keep = np.array([c for c in classes if (yp == c).sum() >= 2])
    if keep.size < 2:
        return None
    mask = np.isin(yp, keep)
    means, prec, counts = gaussian_templates(Xp[mask], yp[mask], keep)

    lp = log_posteriors(Xa, means, prec, counts)
    post = np.exp(lp)
    post /= post.sum(axis=1, keepdims=True)

    pred = keep[np.argmax(post, axis=1)]
    # rank of the true class; traces whose true class was never profiled are
    # scored as a miss at the worst rank rather than silently dropped
    ranks, hcond = [], []
    for i, t in enumerate(ya):
        order = keep[np.argsort(post[i])[::-1]]
        where = np.flatnonzero(order == t)
        ranks.append(int(where[0]) + 1 if where.size else keep.size)
        hcond.append(entropy_bits(post[i]))

    prior = counts / counts.sum()
    return {
        "n_poi": len(pois),
        "pois": pois.tolist(),
        "acc": float((pred == ya).mean()),
        "mae": float(np.abs(pred.astype(int) - ya.astype(int)).mean()),
        "rank_mean": float(np.mean(ranks)),
        "n_classes": int(keep.size),
        "h_prior": entropy_bits(prior),
        "h_cond": float(np.mean(hcond)),
    }


# ----------------------------------------------------------------- bit LDA

def attack_limb_bits(traces, limb_vals, prof_idx, att_idx, n_poi, min_gap, hw):
    """Binary LDA per bit position, on POIs chosen for that bit."""
    res = []
    for b in range(32):
        bp = ((limb_vals[prof_idx] >> b) & 1).astype(np.int8)
        ba = ((limb_vals[att_idx] >> b) & 1).astype(np.int8)
        if bp.min() == bp.max() or ba.min() == ba.max():
            continue
        pois, _ = pick_pois(traces[prof_idx], bp, n_poi, min_gap)
        Xp = traces[prof_idx][:, pois]
        Xa = traces[att_idx][:, pois]
        means, prec, counts = gaussian_templates(Xp, bp, np.array([0, 1]))
        lp = log_posteriors(Xa, means, prec, counts)
        pred = np.argmax(lp, axis=1)
        # baseline: always answer the majority class of the attack set
        base = max(ba.mean(), 1 - ba.mean())
        res.append({"bit": b, "acc": float((pred == ba).mean()), "base": float(base)})
    return res


# -------------------------------------------------------------------- main

def run(npz, n_poi, min_gap, prof_n, seed, do_bits, quiet=False, shuffle=False):
    d = np.load(npz, allow_pickle=False)
    T = d["traces"].astype(np.float64)
    N = T.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(N)
    prof_n = min(prof_n, N - 20)
    prof_idx, att_idx = perm[:prof_n], perm[prof_n:]

    if not quiet:
        print(f"\n{'=' * 78}")
        print(f"{npz}")
        print(f"  {N} traces, {T.shape[1]} samples -> profiling {len(prof_idx)}, "
              f"attack {len(att_idx)} (held out)")
        print(f"  POIs per target: {n_poi}, minimum separation {min_gap} samples")

    hw_rows, bit_rows = [], []
    for c in (0, 1, 2, 3):
        for k in range(clean_limbs(d, c, N)):
            key = f"c{c}_limb{k}"
            if key not in d:
                continue
            lv = d[key].astype(np.uint32)
            hw = d[f"{key}_hw"].astype(np.int16)
            if shuffle:
                hw = hw[rng.permutation(N)]
                lv = lv[rng.permutation(N)]
            r = attack_limb_hw(T, hw, prof_idx, att_idx, n_poi, min_gap)
            if r is None:
                continue
            r.update(coord=c, limb=k)
            hw_rows.append(r)
            if do_bits:
                for br in attack_limb_bits(T, lv, prof_idx, att_idx, n_poi, min_gap, hw):
                    br.update(coord=c, limb=k)
                    bit_rows.append(br)

    return hw_rows, bit_rows, len(prof_idx), len(att_idx)


def report(hw_rows, bit_rows, n_prof, n_att):
    print(f"\n-- Hamming-weight templates, scored on {n_att} held-out traces --")
    print(f"  {'limb':<12}{'classes':>8}{'exact':>8}{'MAE':>7}{'rank':>7}"
          f"{'H(HW)':>8}{'H(HW|t)':>9}{'gain':>7}")
    tot_gain = 0.0
    for r in hw_rows:
        gain = r["h_prior"] - r["h_cond"]
        tot_gain += gain
        print(f"  c{r['coord']} limb {r['limb']:<5}{r['n_classes']:>8}"
              f"{r['acc'] * 100:>7.1f}%{r['mae']:>7.2f}{r['rank_mean']:>7.2f}"
              f"{r['h_prior']:>8.2f}{r['h_cond']:>9.2f}{gain:>7.2f}")
    if hw_rows:
        print(f"  {'sum over ' + str(len(hw_rows)) + ' limbs':<48}{tot_gain:>31.2f} bits")

    if bit_rows:
        acc = np.array([b["acc"] for b in bit_rows])
        base = np.array([b["base"] for b in bit_rows])
        # one-sided binomial z against each bit's own majority baseline
        z = (acc - base) / np.sqrt(np.clip(base * (1 - base), 1e-9, None) / n_att)
        print(f"\n-- per-bit LDA, {len(bit_rows)} bit positions --")
        print(f"  accuracy  mean={acc.mean() * 100:.1f}%  max={acc.max() * 100:.1f}%")
        print(f"  majority baseline mean={base.mean() * 100:.1f}%")
        print(f"  bits beating their baseline by z > 3: {int((z > 3).sum())} of {len(bit_rows)}")
        if (z > 3).sum() == 0:
            print("  -> no bit is individually recoverable: the leakage carries the "
                  "Hamming weight, not which bits are set")
    return tot_gain


def main():
    p = argparse.ArgumentParser(description="Profiled single-trace recovery experiment")
    p.add_argument("--npz", required=True)
    p.add_argument("--n-poi", type=int, default=12)
    p.add_argument("--min-gap", type=int, default=34,
                   help="minimum POI separation in samples (one target clock cycle)")
    p.add_argument("--prof", type=int, default=None, help="profiling traces (default: 80%%)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-bits", action="store_true")
    p.add_argument("--no-control", action="store_true",
                   help="skip the shuffled-label control (halves the runtime)")
    p.add_argument("--curve", default=None,
                   help="comma-separated profiling sizes for a learning curve")
    p.add_argument("--secret-bits", type=int, default=None,
                   help="size of the secret; default: summed from the capture's "
                        "own magnitudes, so line 2 and line 5 each report their own")
    p.add_argument("--name", default=None,
                   help="what the secret is called in the report; default: inferred "
                        "from the trigger mode")
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    d0 = np.load(args.npz, allow_pickle=False)
    n_total = d0["traces"].shape[0]
    prof = args.prof if args.prof else int(0.8 * n_total)

    # The two targets differ in size, so read both the name and the bit count off
    # the capture instead of carrying line 2's numbers into a line 5 report.
    mode = int(d0["trig_mode"]) if "trig_mode" in d0 else 0
    name = args.name or ("alpha_rsp" if mode >= 15 else "alpha_sk,chl")
    secret_bits = args.secret_bits
    if secret_bits is None:
        secret_bits = 4 + sum(int(np.median(d0[f"c{c}_magnitude_bits"]))
                              for c in range(4) if f"c{c}_magnitude_bits" in d0)

    hw_rows, bit_rows, np_, na_ = run(args.npz, args.n_poi, args.min_gap, prof,
                                      args.seed, not args.no_bits)
    gain = report(hw_rows, bit_rows, np_, na_)

    # Picking points of interest by correlation and fitting templates to them
    # finds structure in pure noise, so the estimate carries an upward bias.
    # Measure it on this very capture with the labels shuffled.
    floor = None
    if not args.no_control:
        ctrl, _, _, _ = run(args.npz, args.n_poi, args.min_gap, prof,
                            args.seed + 1, False, quiet=True, shuffle=True)
        floor = sum(max(0.0, r["h_prior"] - r["h_cond"]) for r in ctrl)

    # What perfect Hamming-weight knowledge would be worth, as a ceiling.
    from math import comb, log2
    pk = np.array([comb(32, h) for h in range(33)], dtype=float)
    pk /= pk.sum()
    h_hw = entropy_bits(pk)
    print(f"\n-- recovery budget --")
    print(f"  secret: {secret_bits} bits of {name} across "
          f"{len(hw_rows)} fully-informative limbs (+ boundary limbs and 4 signs)")
    print(f"  ceiling if every limb's Hamming weight were known exactly: "
          f"{len(hw_rows)} x {h_hw:.2f} = {len(hw_rows) * h_hw:.0f} bits")
    print(f"  measured: {gain:.1f} bits")
    if floor is not None:
        print(f"  shuffled-label control on the same traces: {floor:.1f} bits "
              f"-- this is the estimator's own bias, not leakage")
        print(f"  leakage above the control: {gain - floor:+.1f} bits")
    print(f"  remaining search space: about {secret_bits - gain:.0f} bits")

    if args.curve:
        print(f"\n-- learning curve --")
        print(f"  {'profiling':>10}{'HW exact':>10}{'bits gained':>13}")
        for n in [int(x) for x in args.curve.split(",")]:
            if n >= n_total - 20:
                continue
            h, b, np2, na2 = run(args.npz, args.n_poi, args.min_gap, n, args.seed,
                                 False, quiet=True)
            if not h:
                continue
            acc = np.mean([r["acc"] for r in h])
            g = sum(r["h_prior"] - r["h_cond"] for r in h)
            print(f"  {np2:>10}{acc * 100:>9.1f}%{g:>13.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
