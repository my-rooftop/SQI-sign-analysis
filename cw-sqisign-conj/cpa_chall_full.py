#!/usr/bin/env python3
"""
Recover a whole 128-bit entry of Msk from Algorithm 3.23 line 12.

cpa_chall.py recovers limb 0. This continues through the remaining limbs by
following what the schoolbook multiply actually does (mp.c:452):

    for i in range(n1):
        carry = 0
        for j in range(n3):
            prod = a[i] * b[j]
            acc  = res[i+j] + low32(prod) + carry
            res[i+j] = low32(acc)
            carry = high32(prod) + high32(acc)
        res[i+n3] = carry

The outer loop runs over the secret's limbs. When iteration i reaches j = 0 it
writes

    res[i] = low32( res_i_before + low32(a_i * b_0) )

and res_i_before depends only on a_0..a_{i-1}, which are already recovered, and
on b, which is public. So each limb in turn has an intermediate that is a known
offset plus an unknown limb times a known limb -- the same shape limb 0 had, and
attackable the same way.

Each limb is recovered in two rounds, low 16 bits then high 16 bits, keeping a
candidate list from the first round: a 16-bit model predicts only half the
product, which is enough to put the truth in the top few but not to rank it
first, and that does not improve with more traces because it is a limit of the
model rather than of the statistics.

Usage:
  ~/.cwvenv/bin/python3 cpa_chall_full.py --npz /tmp/sqisign_chall_cpa.npz
  ~/.cwvenv/bin/python3 cpa_chall_full.py --npz ... --traces 2000 --keep 16
"""

import argparse
import pathlib
import sys

import numpy as np

MASK32 = np.uint64(0xFFFFFFFF)
POPCNT = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def hw_rows(v, nbytes):
    out = np.zeros(v.shape, dtype=np.int16)
    for k in range(nbytes):
        out += POPCNT[((v >> np.uint64(8 * k)) & np.uint64(0xFF)).astype(np.uint8)]
    return out


def res_before(limbs_known, b, idx, nlimb):
    """res[idx] just before outer iteration `idx` touches it.

    Simulates the multiply with the already-recovered low limbs in place and the
    unknown ones zero, which is exactly the state the real loop is in when it
    arrives at (i = idx, j = 0).
    """
    n = b.shape[0]
    res = np.zeros((n, 2 * nlimb + 2), dtype=np.uint64)
    for i in range(idx):
        ai = np.uint64(limbs_known[i])
        carry = np.zeros(n, dtype=np.uint64)
        for j in range(nlimb):
            prod = ai * b[:, j]
            acc = res[:, i + j] + (prod & MASK32) + carry
            res[:, i + j] = acc & MASK32
            carry = (prod >> np.uint64(32)) + (acc >> np.uint64(32))
        res[:, i + nlimb] = res[:, i + nlimb] + carry
    return res[:, idx].copy()


def score(hyps, base, b0, x, xn, nbits, batch=1 << 13):
    """|corr| at one sample for each hypothesis of the limb's low `nbits`."""
    mask = np.uint64((1 << nbits) - 1)
    nby = (nbits + 7) // 8
    out = np.empty(len(hyps))
    for s in range(0, len(hyps), batch):
        h = np.asarray(hyps[s:s + batch], dtype=np.uint64)
        pred = (base[None, :] + ((h[:, None] * b0[None, :]) & MASK32)) & mask
        M = hw_rows(pred, nby).astype(np.float64)
        M -= M.mean(axis=1, keepdims=True)
        mn = np.sqrt((M * M).sum(axis=1))
        mn = np.where(mn == 0, np.finfo(np.float64).tiny, mn)
        out[s:s + len(h)] = np.abs(M @ x) / (mn * xn)
    return out


def score_all_samples(hyps, base, b0, Tc, Tn, nbits, batch=2048):
    """Same, but maximised over every sample -- no POI and no key needed."""
    mask = np.uint64((1 << nbits) - 1)
    nby = (nbits + 7) // 8
    best = np.empty(len(hyps), np.float32)
    arg = np.empty(len(hyps), np.int64)
    for s in range(0, len(hyps), batch):
        h = np.asarray(hyps[s:s + batch], dtype=np.uint64)
        pred = (base[None, :] + ((h[:, None] * b0[None, :]) & MASK32)) & mask
        M = hw_rows(pred, nby).astype(np.float32)
        M -= M.mean(1, keepdims=True)
        mn = np.sqrt((M * M).sum(1))
        mn = np.where(mn == 0, np.float32(1e-30), mn)
        R = np.abs((M @ Tc) / (mn[:, None] * Tn[None, :]))
        j = R.argmax(1)
        best[s:s + len(h)] = R[np.arange(len(h)), j]
        arg[s:s + len(h)] = j
    return best, arg


def recover_one_limb(base, b0, Tc32, Tn32, Xc, keep, truth=None, idx=0):
    # round 1: low 16 bits, searched over every sample so the attack picks its
    # own point of interest
    s1, a1 = score_all_samples(np.arange(1 << 16, dtype=np.uint64), base, b0,
                               Tc32, Tn32, 16)
    order = np.argsort(s1)[::-1]
    poi = int(a1[order[0]])
    r1 = None
    if truth is not None:
        r1 = int(np.flatnonzero(order == (truth & 0xFFFF))[0]) + 1

    x = Xc[:, poi].astype(np.float64).copy()
    x -= x.mean()
    xn = np.sqrt((x * x).sum())

    hi = np.arange(1 << 16, dtype=np.uint64)
    best_v, best_s = None, -1.0
    for c in order[:keep].astype(np.uint64):
        full = c | (hi << np.uint64(16))
        sc = score(full, base, b0, x, xn, 32)
        j = int(np.argmax(sc))
        if sc[j] > best_s:
            best_s, best_v = float(sc[j]), int(full[j])

    tag = ""
    if truth is not None:
        tag = "  RECOVERED" if best_v == truth else f"  (true 0x{truth:08x})"
    print(f"  limb {idx}: round 1 rank {r1 if r1 else '-'}/65536, POI {poi}; "
          f"round 2 -> 0x{best_v:08x} |corr|={best_s:.4f}{tag}", flush=True)
    return best_v, best_s


def main():
    p = argparse.ArgumentParser(description="Full 128-bit Msk entry recovery")
    p.add_argument("--npz", required=True)
    p.add_argument("--traces", type=int, default=2000)
    p.add_argument("--keep", type=int, default=16)
    p.add_argument("--entry", type=int, default=1,
                   help="which Msk entry the capture triggered on (1 = Msk[0][1])")
    p.add_argument("--limbs", type=int, default=4)
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    d = np.load(args.npz, allow_pickle=False)
    T = d["traces"]
    n = min(args.traces, T.shape[0])
    T = T[:n].astype(np.float64)
    b = d["coeff_limbs"][:n].astype(np.uint64)
    secret = [int(x, 16) for x in d["secret_hex"]]
    target = secret[args.entry]

    Xc = T - T.mean(axis=0, keepdims=True)
    Tc32 = Xc.astype(np.float32)
    Tn32 = np.sqrt((Tc32 * Tc32).sum(0))
    Tn32 = np.where(Tn32 == 0, np.float32(1e-30), Tn32)

    nlimb = b.shape[1]
    print(f"\n{'=' * 76}")
    print(f"{args.npz}")
    print(f"  {n} traces, {T.shape[1]} samples, {nlimb} known limbs per coefficient")
    print(f"  target Msk[{args.entry // 2}][{args.entry % 2}] = 0x{target:032x}")
    print(f"  non-profiled: round 1 searches all samples and nominates its own POI\n")

    got = []
    for i in range(args.limbs):
        truth = (target >> (32 * i)) & 0xFFFFFFFF
        base = res_before(got, b, i, nlimb)
        v, _ = recover_one_limb(base, b[:, 0], Tc32, Tn32, Xc, args.keep,
                                truth=truth, idx=i)
        got.append(v)
        if v != truth:
            print(f"  -> limb {i} wrong; later limbs build on it, stopping", flush=True)
            break

    rec = 0
    for i, v in enumerate(got):
        rec |= v << (32 * i)
    bits_ok = sum(1 for k in range(32 * len(got))
                  if ((rec >> k) & 1) == ((target >> k) & 1))
    print(f"\n  recovered 0x{rec:032x}")
    print(f"  truth     0x{target & ((1 << (32 * len(got))) - 1):032x}")
    print(f"  limbs correct: {sum(1 for i, v in enumerate(got) if v == ((target >> (32 * i)) & 0xFFFFFFFF))}"
          f"/{args.limbs}   bits correct: {bits_ok}/{32 * args.limbs} of the entry's 128")
    return 0


if __name__ == "__main__":
    sys.exit(main())
