#!/usr/bin/env python3
"""
Step 4: quantify how much Algorithm 3.18 line 2 leaks. No key recovery here --
this only answers "is there exploitable leakage, where, and how strong".

Three measurements:

  TVLA   Welch t-test between the fixed and random groups of a --mode tvla
         capture. |t| > 4.5 at any sample is the usual threshold for declaring
         leakage. Also reported as a half-split sanity check: fixed vs fixed must
         NOT cross 4.5, otherwise the rig itself is drifting and the fixed-vs-random
         result means nothing.

  SNR    Per limb, over a --mode random capture, using the Hamming weight of the
         limb as the class label:
             SNR(t) = Var_over_classes(mean_c signal) / mean_c(Var_within_class)
         An SNR of ~0.1 is already workable for CPA; below ~0.01 expect to need
         many traces. Reported per (coordinate, limb) so the 25 sign-extension
         limbs -- which carry no secret -- can be told apart from the 14
         informative ones.

  CORR   Pearson correlation of each sample against HW(limb) and HW(~limb), which
         is the CPA distinguisher evaluated at the correct hypothesis. Its peak
         says how well a real attack would separate, and where in the window each
         limb lands.

Usage:
  ~/.cwvenv/bin/python3 analyze_conj_leakage.py --tvla /tmp/sqisign_conj_tvla.npz
  ~/.cwvenv/bin/python3 analyze_conj_leakage.py --snr  /tmp/sqisign_conj_random.npz
  ~/.cwvenv/bin/python3 analyze_conj_leakage.py --tvla /tmp/sqisign_conj_tvla.npz \
                                                --snr  /tmp/sqisign_conj_random.npz
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TVLA_THRESHOLD = 4.5


def welch_t(a, b):
    """Welch t statistic per sample between two trace sets."""
    na, nb = a.shape[0], b.shape[0]
    if na < 2 or nb < 2:
        raise ValueError(f"need >=2 traces per group, got {na} and {nb}")
    ma, mb = a.mean(axis=0), b.mean(axis=0)
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    denom = np.sqrt(va / na + vb / nb)
    denom[denom == 0] = np.finfo(np.float64).tiny
    return (ma - mb) / denom


def bin_hw(hw, nbins):
    """Group Hamming weights into equal-width bins.

    HW of a uniform 32-bit word is binomial(32, 1/2): HW=16 is common, HW<=6 and
    HW>=26 are rare. Using the raw HW as the class label therefore leaves the
    extreme classes with a handful of traces each, and their noisy between-class
    means push the SNR argmax onto noise. Binning trades a little of the signal
    for classes that are actually populated.
    """
    if nbins <= 0:
        return np.asarray(hw)
    h = np.asarray(hw, dtype=np.float64)
    edges = np.linspace(0, 32, nbins + 1)[1:-1]
    return np.digitize(h, edges).astype(np.int16)


def snr_by_class(traces, labels, min_class=5):
    """Signal-to-noise ratio per sample, classes given by `labels`.

    Classes with fewer than `min_class` traces are dropped: their sample means are
    too noisy to contribute signal and instead inflate the numerator.
    """
    classes = np.unique(labels)
    if classes.size < 2:
        return None
    means, variances, weights = [], [], []
    for c in classes:
        sel = traces[labels == c]
        if sel.shape[0] < max(2, min_class):
            continue
        means.append(sel.mean(axis=0))
        variances.append(sel.var(axis=0, ddof=1))
        weights.append(sel.shape[0])
    if len(means) < 2:
        return None
    means = np.asarray(means)
    variances = np.asarray(variances)
    weights = np.asarray(weights, dtype=np.float64)
    signal = means.var(axis=0, ddof=1)
    noise = np.average(variances, axis=0, weights=weights)
    noise = np.where(noise == 0, np.finfo(np.float64).tiny, noise)
    return signal / noise


def classify_limb(values, sign=None):
    """What kind of label a limb's observed values actually form.

    A limb above the value's magnitude is sign extension: it is 0x00000000 or
    0xffffffff depending only on the sign of the coordinate. Across traces it is
    not constant, but it is a one-bit sign label wearing a 32-bit limb's clothes,
    and every sign-extension limb of the same coordinate carries the identical
    label. Correlating against it finds the sign leak -- strongly -- and reports
    the same peak sample for all of them, which would otherwise be mistaken for
    per-limb Hamming-weight recovery.

    Returns one of "constant", "sign", "informative".
    """
    v = np.asarray(values, dtype=np.uint64)
    distinct = np.unique(v)
    if distinct.size <= 1:
        return "constant"
    if np.all((v == 0) | (v == 0xFFFFFFFF)):
        return "sign"
    # a handful of distinct values is a near-sign label too, not a HW label
    if distinct.size <= 3:
        return "sign"
    # The boundary limb -- the one straddling the top of the magnitude -- is the
    # subtle case. It holds a few real value bits under a run of sign-extension
    # bits, so it has thousands of distinct values and slips past every test
    # above, yet its Hamming weight is still almost entirely decided by the sign:
    # for alpha_rsp, HW averages 5.5 when the coordinate is positive and 26.5
    # when it is negative. Correlating against it re-finds the sign leak at
    # |corr| ~ 0.95 and would be read as the strongest per-limb value leak in the
    # capture, which it is not -- that one bit is already counted as the sign.
    if sign is not None:
        hw = np.array([int(x).bit_count() for x in v], dtype=np.float64)
        sg = np.asarray(sign, dtype=np.float64)
        if hw.std() > 0 and sg.std() > 0 and abs(np.corrcoef(hw, sg)[0, 1]) >= 0.9:
            return "sign"
    return "informative"


def corr_per_sample(traces, values):
    """Pearson correlation of every sample against a per-trace scalar."""
    v = np.asarray(values, dtype=np.float64)
    v = v - v.mean()
    vn = np.sqrt((v * v).sum())
    if vn == 0:
        return np.zeros(traces.shape[1])
    t = traces.astype(np.float64)
    t = t - t.mean(axis=0, keepdims=True)
    tn = np.sqrt((t * t).sum(axis=0))
    tn = np.where(tn == 0, np.finfo(np.float64).tiny, tn)
    return (t * v[:, None]).sum(axis=0) / (tn * vn)


# ============================================================
# TVLA
# ============================================================

def do_tvla(path, plot_path):
    d = np.load(path, allow_pickle=False)
    traces = d["traces"].astype(np.float64)
    group = d["group"]
    spc = float(d["samples_per_cycle"]) if "samples_per_cycle" in d else 0.0

    fixed = traces[group == 1]
    rand = traces[group == 0]
    print(f"\n{'=' * 72}\nTVLA  {path}")
    print(f"  traces: fixed={fixed.shape[0]} random={rand.shape[0]} "
          f"samples={traces.shape[1]}"
          + (f"  ({spc:.1f} samples per target cycle)" if spc else ""))
    if fixed.shape[0] < 2 or rand.shape[0] < 2:
        print("  not enough traces in both groups; capture with --mode tvla")
        return None

    t = welch_t(fixed, rand)
    peak = int(np.argmax(np.abs(t)))
    n_over = int((np.abs(t) > TVLA_THRESHOLD).sum())
    print(f"  fixed vs random : max|t| = {np.abs(t).max():.2f} at sample {peak}"
          + (f" (cycle ~{peak / spc:.0f})" if spc else ""))
    print(f"                    {n_over} of {t.size} samples exceed "
          f"{TVLA_THRESHOLD} ({100 * n_over / t.size:.2f}%)")

    # sanity: splitting the fixed group in half must NOT show leakage
    half = fixed.shape[0] // 2
    t_null = None
    if half >= 2:
        t_null = welch_t(fixed[:half], fixed[half:2 * half])
        n_null = int((np.abs(t_null) > TVLA_THRESHOLD).sum())
        print(f"  fixed vs fixed  : max|t| = {np.abs(t_null).max():.2f}, "
              f"{n_null} samples over threshold  "
              f"{'OK' if np.abs(t_null).max() < TVLA_THRESHOLD else 'RIG DRIFT -- result not trustworthy'}")

    verdict = ("LEAKAGE DETECTED" if np.abs(t).max() > TVLA_THRESHOLD
               else "no leakage above threshold at this trace count")
    print(f"  verdict: {verdict}")

    fig, ax = plt.subplots(2 if t_null is None else 3, 1,
                           figsize=(14, 4 + 3 * (2 if t_null is None else 3)))
    ax[0].plot(fixed.mean(axis=0), lw=0.5, label="fixed mean")
    ax[0].plot(rand.mean(axis=0), lw=0.5, label="random mean")
    ax[0].set_title("Group means")
    ax[0].set_ylabel("mV")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    ax[1].plot(t, color="C3", lw=0.5)
    for s in (TVLA_THRESHOLD, -TVLA_THRESHOLD):
        ax[1].axhline(s, color="k", ls="--", lw=0.7)
    ax[1].set_title(f"TVLA fixed vs random, max|t| = {np.abs(t).max():.2f}")
    ax[1].set_ylabel("t")
    ax[1].grid(alpha=0.3)

    if t_null is not None:
        ax[2].plot(t_null, color="gray", lw=0.5)
        for s in (TVLA_THRESHOLD, -TVLA_THRESHOLD):
            ax[2].axhline(s, color="k", ls="--", lw=0.7)
        ax[2].set_title(f"Sanity: fixed vs fixed, max|t| = {np.abs(t_null).max():.2f} "
                        f"(should stay inside the dashed lines)")
        ax[2].set_ylabel("t")
        ax[2].grid(alpha=0.3)

    ax[-1].set_xlabel("sample")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=130)
    plt.close(fig)
    print(f"  saved {plot_path}")
    return t


# ============================================================
# SNR / correlation per limb
# ============================================================

def do_snr(path, plot_path, coords, top_n, hw_bins, min_class):
    d = np.load(path, allow_pickle=False)
    traces = d["traces"].astype(np.float64)
    nlimbs = int(d["nlimbs"])
    spc = float(d["samples_per_cycle"]) if "samples_per_cycle" in d else 0.0
    recs = json.loads(str(d["records_json"])) if "records_json" in d else []

    ntr = traces.shape[0]
    nsamp = traces.shape[1]
    # Threshold for a |corr| MAXIMISED OVER SAMPLES, which is what is reported
    # below. Under the null each sample's r has sd ~1/sqrt(N), and the maximum of
    # nsamp of them sits near sqrt(2 ln nsamp) sd -- for 35k samples that is 4.6
    # sd, far above the 4/sqrt(N) figure appropriate to a single pre-chosen
    # sample. Comparing a max against the per-sample threshold declares leakage
    # on pure noise.
    corr_floor_pointwise = 4.0 / np.sqrt(max(ntr, 2))
    corr_floor = np.sqrt(2.0 * np.log(max(nsamp, 2))) / np.sqrt(max(ntr, 2))

    print(f"\n{'=' * 72}\nSNR / correlation  {path}")
    print(f"  traces={ntr} samples={traces.shape[1]} nlimbs={nlimbs}"
          + (f"  ({spc:.1f} samples per target cycle)" if spc else ""))
    print(f"  HW class binning: {hw_bins if hw_bins > 0 else 'none (raw HW)'}, "
          f"min traces per class: {min_class}")
    print(f"  |corr| thresholds: {corr_floor_pointwise:.3f} per-sample, "
          f"{corr_floor:.3f} for a max over all {nsamp} samples (the figure used here)")
    if ntr < 20:
        print("  WARNING: fewer than 20 traces; SNR estimates will be very noisy")

    rows = []
    for c in coords:
        nl_key = f"c{c}_nlimbs"
        se_key = f"c{c}_signext_limbs"
        nl = int(np.median(d[nl_key])) if nl_key in d else nlimbs
        se = int(np.median(d[se_key])) if se_key in d else 0
        # Count the limbs the way the ranking below does, so the header cannot
        # disagree with the summary: the stored sign-extension count only covers
        # the all-0x00/all-0xff limbs, and misses the boundary limb whose
        # Hamming weight the sign still decides.
        sgn = d[f"c{c}_sign"] if f"c{c}_sign" in d else None
        kinds = [classify_limb(d[f"c{c}_limb{k}"].astype(np.uint32), sgn)
                 for k in range(nl) if f"c{c}_limb{k}" in d]
        info_limbs = sum(1 for k in kinds if k == "informative")
        sign_limbs = sum(1 for k in kinds if k == "sign")
        print(f"\n  coord[{c}]: {nl} limbs, {sign_limbs} decided by the sign "
              f"({se} of them pure sign extension), {info_limbs} informative"
              + ("   (copied by ibz_copy, not negated)" if c == 0 else ""))

        for k in range(nl):
            hw_key = f"c{c}_limb{k}_hw"
            lv_key = f"c{c}_limb{k}"
            if hw_key not in d:
                continue
            hw = d[hw_key]
            lv = d[lv_key]
            kind = classify_limb(lv, d[f"c{c}_sign"] if f"c{c}_sign" in d else None)
            s = snr_by_class(traces, bin_hw(hw, hw_bins), min_class=min_class)
            if kind == "constant" or s is None:
                print(f"    limb {k:2d}: constant across traces "
                      f"(0x{int(lv[0]):08x}), no leakage to measure")
                rows.append({"coord": c, "limb": k, "snr": np.nan, "corr": np.nan,
                             "peak": -1, "kind": "constant"})
                continue
            r_hw = corr_per_sample(traces, hw)
            # HW(~limb) is the complementary hypothesis a CPA attacker also has
            # to test. Older captures stored it; derive it when they did not.
            neg_key = f"c{c}_limb{k}_neg_hw"
            neg = d[neg_key] if neg_key in d else np.array(
                [32 - int(v).bit_count() for v in lv], dtype=np.uint8)
            r_neg = corr_per_sample(traces, neg)

            # The correlation argmax is the more stable locator, so the reported
            # point of interest is where |corr(HW)| peaks; SNR is then read at the
            # same sample. Two independent statistics peaking together is the
            # evidence we want -- an SNR spike on its own is often just a noisy
            # sparse class.
            poi = int(np.argmax(np.abs(r_hw)))
            snr_at_poi = float(s[poi])
            snr_max = float(s.max())
            snr_argmax = int(np.argmax(s))
            n_over = int((np.abs(r_hw) > corr_floor).sum())

            rows.append({
                "coord": c, "limb": k,
                "snr": snr_at_poi,
                "snr_max": snr_max,
                "snr_argmax": snr_argmax,
                "corr": float(np.abs(r_hw)[poi]),
                "corr_neg": float(np.abs(r_neg).max()),
                "peak": poi,
                "n_over_floor": n_over,
                "kind": kind,
                "distinct": int(np.unique(lv).size),
            })
            agree = abs(snr_argmax - poi) <= max(4, int(spc)) if spc else abs(snr_argmax - poi) <= 4
            tag = "  <- SIGN label, not a HW label" if kind == "sign" else ""
            print(f"    limb {k:2d}: |corr(HW)|={np.abs(r_hw)[poi]:.3f} at sample {poi}"
                  + (f" (cycle ~{poi / spc:5.0f})" if spc else "")
                  + f"  SNR@poi={snr_at_poi:7.4f}"
                  f"  distinct={np.unique(lv).size:3d}"
                  f"  over_floor={n_over}{tag}"
                  + ("" if agree else "   [SNR peak elsewhere: %d]" % snr_argmax))

    real = [r for r in rows if r.get("kind") == "informative" and np.isfinite(r["snr"])]
    signs = [r for r in rows if r.get("kind") == "sign"]

    if signs:
        sc = np.array([r["corr"] for r in signs])
        print(f"\n  {len(signs)} sign-extension limbs (excluded from the ranking): "
              f"|corr| median={np.median(sc):.3f}. These all carry the same one-bit "
              f"label and peak at the same samples, so they measure the sign leak, "
              f"not per-limb Hamming weight.")

    if real:
        real_sorted = sorted(real, key=lambda r: -r["corr"])
        print(f"\n  strongest {min(top_n, len(real_sorted))} informative limbs by |corr(HW)|:")
        for r in real_sorted[:top_n]:
            print(f"    coord[{r['coord']}] limb {r['limb']:2d}  |corr|={r['corr']:.3f}  "
                  f"SNR@poi={r['snr']:.4f}  sample {r['peak']}  "
                  f"distinct={r['distinct']}  over floor={r['n_over_floor']}")
        snrs = np.array([r["snr"] for r in real])
        corrs = np.array([r["corr"] for r in real])
        n_sig = int((corrs > corr_floor).sum())
        print(f"\n  across {len(real)} informative limbs: |corr| median={np.median(corrs):.3f} "
              f"max={corrs.max():.3f}; {n_sig}/{len(real)} peak above the "
              f"{corr_floor:.3f} noise floor")
        print(f"  SNR@poi: median={np.median(snrs):.4f} max={snrs.max():.4f}")
        base = np.median(corrs)
        # traces needed for CPA to separate the right hypothesis, the usual
        # Mangard estimate n ~ 28 / ln^2((1+r)/(1-r)) at 99.99% confidence
        rr = min(base, 0.999)
        n_needed = 28.0 / (np.log((1 + rr) / (1 - rr)) ** 2)
        print(f"  at the median informative-limb correlation r={base:.3f}, a CPA "
              f"against a FIXED value would need about {n_needed:.0f} traces per limb")
        # That number does not describe an attack on this target. CPA averages
        # many traces against one fixed secret; both alpha_sk,chl and alpha_rsp
        # are freshly sampled for every signature, so there is nothing to average
        # over and more traces only buy a better profile. The attack that counts
        # is profiled, with a single-trace attack phase -- see
        # recover_conj_template.py and msb_recover.py. What the correlation
        # really says here is how well a template can be built, not how fast CPA
        # would converge.
        print("  note: these operands are redrawn per signature, so CPA does not "
              "apply; the figure above is a leakage-strength yardstick only.")
        if snrs.max() >= 0.1:
            print("  -> strong per-limb leakage; a profiled single-trace attack is "
                  "worth attempting. Whether it recovers bits is a separate "
                  "question, answered by the template and MSB experiments.")
        elif snrs.max() >= 0.01:
            print("  -> weak but present; expect a profiled attack to struggle")
        else:
            print("  -> below the useful range at this trace count; improve the "
                  "measurement (probe, range, averaging) before attempting recovery")

    # ---- plot ----
    fig, ax = plt.subplots(len(coords) + 1, 1, figsize=(14, 3.2 * (len(coords) + 1)))
    ax[0].plot(traces.mean(axis=0), color="C0", lw=0.6)
    ax[0].set_title(f"Mean trace ({traces.shape[0]} traces)")
    ax[0].set_ylabel("mV")
    ax[0].grid(alpha=0.3)

    for idx, c in enumerate(coords, start=1):
        nl = int(np.median(d[f"c{c}_nlimbs"])) if f"c{c}_nlimbs" in d else nlimbs
        drawn = 0
        for k in range(nl):
            hw_key = f"c{c}_limb{k}_hw"
            if hw_key not in d:
                continue
            if np.unique(d[f"c{c}_limb{k}"]).size <= 1:
                continue
            s = snr_by_class(traces, bin_hw(d[hw_key], hw_bins), min_class=min_class)
            if s is None:
                continue
            ax[idx].plot(s, lw=0.7, label=f"limb {k}")
            drawn += 1
        ax[idx].set_title(f"SNR vs HW(limb), coord[{c}] ({drawn} varying limbs)")
        ax[idx].set_ylabel("SNR")
        ax[idx].grid(alpha=0.3)
        if drawn:
            ax[idx].legend(ncol=7, fontsize=7)

    ax[-1].set_xlabel("sample")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=130)
    plt.close(fig)
    print(f"\n  saved {plot_path}")
    return rows


def main():
    p = argparse.ArgumentParser(description="TVLA and SNR for Alg 3.18 line 2 Conjugate")
    p.add_argument("--tvla", help="npz from --mode tvla")
    p.add_argument("--snr", help="npz from --mode random")
    p.add_argument("--coords", default="1,2,3",
                   help="coordinates to analyse (0 is copied rather than negated)")
    p.add_argument("--top", type=int, default=10)
    p.add_argument("--hw-bins", type=int, default=8,
                   help="bin Hamming weights into this many classes (0 = raw HW). "
                        "Raw HW leaves the extreme classes nearly empty.")
    p.add_argument("--min-class", type=int, default=5,
                   help="drop SNR classes with fewer traces than this")
    p.add_argument("--tvla-plot", default="/tmp/sqisign_conj_tvla_plot.png")
    p.add_argument("--snr-plot", default="/tmp/sqisign_conj_snr_plot.png")
    args = p.parse_args()

    if not args.tvla and not args.snr:
        p.error("give --tvla and/or --snr")

    coords = [int(x) for x in args.coords.split(",") if x.strip() != ""]

    if args.tvla:
        if not pathlib.Path(args.tvla).exists():
            print(f"missing: {args.tvla}", file=sys.stderr)
            return 1
        do_tvla(args.tvla, args.tvla_plot)

    if args.snr:
        if not pathlib.Path(args.snr).exists():
            print(f"missing: {args.snr}", file=sys.stderr)
            return 1
        do_snr(args.snr, args.snr_plot, coords, args.top, args.hw_bins, args.min_class)

    return 0


if __name__ == "__main__":
    sys.exit(main())
