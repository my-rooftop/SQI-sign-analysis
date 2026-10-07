#!/usr/bin/env python3
"""
Recover the whole 512-bit secret matrix Msk from Algorithm 3.23 line 12.

ibz_mat_2x2_eval (dim2.c:42) with vec = [1, achl] does

    prod = m00 * 1    ; matvec0 = prod
    prod = m01 * achl ; matvec0 = matvec0 + prod
    prod = m10 * 1    ; matvec1 = prod
    prod = m11 * achl ; matvec1 = matvec1 + prod

Two of the entries meet the known coefficient in a multiply and two do not, so
they need different intermediates:

  m01, m11   attacked through ibz_mul (mp.c:452). Outer iteration i writes
             res[i] = low32(res_i_before + low32(a_i * b_0)), where res_i_before
             follows from the already-recovered lower limbs and the public
             coefficient. Four limbs, low 16 bits then high 16, candidate list
             carried between rounds.

  m00, m10   multiplied by the constant 1, so their multiply has no varying
             input and leaks only a Hamming weight. The addition that follows is
             the way in: once m01 is known, prod = m01 * achl is known and varies
             per trace, and ibz_add (mp.c:383) writes
                 sum[k] = low32(matvec[k] + prod[k] + carry)
             with matvec[k] being the unknown limb of m00. Same shape, same
             attack. m10 follows from m11 the same way.

So the order is forced: m01 before m00, m11 before m10.

Non-profiled throughout -- round 1 of each limb scans every sample and the
winning hypothesis nominates the point of interest.

Usage:
  ~/.cwvenv/bin/python3 cpa_chall_matrix.py --npz /tmp/sqisign_chall_eval.npz
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


def mul_state(limbs_known, b, idx, nlimb, width):
    """res[idx] as ibz_mul leaves it just before outer iteration idx."""
    n = b.shape[0]
    res = np.zeros((n, width), dtype=np.uint64)
    for i in range(idx):
        ai = np.uint64(limbs_known[i])
        carry = np.zeros(n, dtype=np.uint64)
        for j in range(nlimb):
            prod = ai * b[:, j]
            acc = res[:, i + j] + (prod & MASK32) + carry
            res[:, i + j] = acc & MASK32
            carry = (prod >> np.uint64(32)) + (acc >> np.uint64(32))
        res[:, i + nlimb] += carry
    return res[:, idx].copy()


def full_product(a_limbs, b, nlimb, width):
    """All limbs of a * b, for a already recovered."""
    n = b.shape[0]
    res = np.zeros((n, width), dtype=np.uint64)
    for i, ai_ in enumerate(a_limbs):
        ai = np.uint64(ai_)
        carry = np.zeros(n, dtype=np.uint64)
        for j in range(nlimb):
            prod = ai * b[:, j]
            acc = res[:, i + j] + (prod & MASK32) + carry
            res[:, i + j] = acc & MASK32
            carry = (prod >> np.uint64(32)) + (acc >> np.uint64(32))
        res[:, i + nlimb] += carry
    return res


def add_carry_in(known_limbs, P, idx):
    """Carry entering limb idx of ibz_add, given the recovered lower limbs."""
    n = P.shape[0]
    carry = np.zeros(n, dtype=np.uint64)
    for k in range(idx):
        s = np.uint64(known_limbs[k]) + P[:, k] + carry
        carry = s >> np.uint64(32)
    return carry


def score_all_samples(models, Tc, Tn, batch=2048):
    best = np.empty(models.shape[0], np.float32)
    arg = np.empty(models.shape[0], np.int64)
    for s in range(0, models.shape[0], batch):
        M = models[s:s + batch].astype(np.float32)
        M -= M.mean(1, keepdims=True)
        mn = np.sqrt((M * M).sum(1))
        mn = np.where(mn == 0, np.float32(1e-30), mn)
        R = np.abs((M @ Tc) / (mn[:, None] * Tn[None, :]))
        j = R.argmax(1)
        best[s:s + M.shape[0]] = R[np.arange(M.shape[0]), j]
        arg[s:s + M.shape[0]] = j
    return best, arg


def score_at(models, x, xn, batch=1 << 13):
    out = np.empty(models.shape[0])
    for s in range(0, models.shape[0], batch):
        M = models[s:s + batch].astype(np.float64)
        M -= M.mean(1, keepdims=True)
        mn = np.sqrt((M * M).sum(1))
        mn = np.where(mn == 0, np.finfo(np.float64).tiny, mn)
        out[s:s + M.shape[0]] = np.abs(M @ x) / (mn * xn)
    return out


def recover_limb(predict, Tc32, Tn32, Xc, keep, truth=None, label="", batch=1 << 13,
                 stride=1, window=None):
    """Two rounds: low 16 bits over all samples, then full 32 at the chosen POI.

    `predict(hyps, nbits)` returns the modelled intermediate, masked to nbits,
    for each hypothesis of the limb's low bits.
    """
    hyps = np.arange(1 << 16, dtype=np.uint64)
    models = np.empty((len(hyps), Tc32.shape[0]), dtype=np.int16)
    for s in range(0, len(hyps), batch):
        models[s:s + batch] = hw_rows(predict(hyps[s:s + batch], 16), 2)
    # Round 1 on a strided subset: the correlation peak spans a clock cycle, about
    # 17 samples here, so sampling every few loses nothing and cuts the search
    # over a 99k-sample window from minutes to seconds per limb.
    # Each entry has to be searched inside its own stretch of the trace. The
    # model a_i * b_0 has the same shape for Msk[0][1] and Msk[1][1], so a sweep
    # over the whole window simply finds whichever of the two multiplies leaks
    # harder and returns that entry's value for both.
    w0, w1 = (0, Tc32.shape[1]) if window is None else window
    sub = Tc32[:, w0:w1:stride]
    subn = Tn32[w0:w1:stride]
    s1, a1 = score_all_samples(models, sub, subn)
    order = np.argsort(s1)[::-1]
    poi_coarse = w0 + int(a1[order[0]]) * stride
    r1 = int(np.flatnonzero(order == (truth & 0xFFFF))[0]) + 1 if truth is not None else None

    lo = max(w0, poi_coarse - 2 * stride)
    hi_s = min(w1, poi_coarse + 2 * stride + 1)
    W = Xc[:, lo:hi_s].astype(np.float64)
    W = W - W.mean(axis=0, keepdims=True)
    Wn = np.sqrt((W * W).sum(axis=0))
    Wn = np.where(Wn == 0, np.finfo(np.float64).tiny, Wn)

    hi = np.arange(1 << 16, dtype=np.uint64)
    best_v, best_s, poi = None, -1.0, poi_coarse
    for c in order[:keep].astype(np.uint64):
        full = c | (hi << np.uint64(16))
        for s in range(0, len(full), batch):
            M = hw_rows(predict(full[s:s + batch], 32), 4).astype(np.float64)
            M -= M.mean(axis=1, keepdims=True)
            mn = np.sqrt((M * M).sum(axis=1))
            mn = np.where(mn == 0, np.finfo(np.float64).tiny, mn)
            R = np.abs((M @ W) / (mn[:, None] * Wn[None, :]))
            j = int(np.argmax(R))
            jr, jc = divmod(j, R.shape[1])
            if R[jr, jc] > best_s:
                best_s = float(R[jr, jc])
                best_v = int(full[s + jr])
                poi = lo + jc

    tag = ""
    if truth is not None:
        tag = "  OK" if best_v == truth else f"  (true 0x{truth:08x})"
    print(f"    {label}: r1 rank {r1 if r1 else '-'}/65536, POI {poi}; "
          f"-> 0x{best_v:08x} |corr|={best_s:.4f}{tag}", flush=True)
    return best_v


def main():
    p = argparse.ArgumentParser(description="Full Msk recovery from Alg 3.23 line 12")
    p.add_argument("--npz", required=True)
    p.add_argument("--traces", type=int, default=3000)
    p.add_argument("--keep", type=int, default=16)
    p.add_argument("--limbs", type=int, default=4)
    p.add_argument("--stride", type=int, default=4,
                   help="sample stride for the round 1 sweep; the peak is a clock "
                        "cycle wide so a few samples apart loses nothing")
    args = p.parse_args()

    if not pathlib.Path(args.npz).exists():
        print(f"missing: {args.npz}", file=sys.stderr)
        return 1

    d = np.load(args.npz, allow_pickle=False)
    T = d["traces"]
    n = min(args.traces, T.shape[0])
    T = T[:n].astype(np.float64)
    b = d["coeff_limbs"][:n].astype(np.uint64)
    secret = [int(x, 16) for x in d["secret_hex"]]  # m00, m01, m10, m11

    Xc = T - T.mean(axis=0, keepdims=True)
    Tc32 = Xc.astype(np.float32)
    Tn32 = np.sqrt((Tc32 * Tc32).sum(0))
    Tn32 = np.where(Tn32 == 0, np.float32(1e-30), Tn32)

    nlimb = b.shape[1]
    width = 2 * nlimb + 4
    print(f"\n{'=' * 76}")
    print(f"{args.npz}")
    print(f"  {n} traces, {T.shape[1]} samples; non-profiled\n")

    got = {}
    # Operation boundaries from the firmware's own cycle profile (command 'f'):
    # mul0 905, copy, mul1 1560, add, mul2 905, copy, mul3 1533, add, eval 5894.
    # The copies and adds share the remaining 991 cycles. Nothing here depends on
    # the key -- these are properties of the implementation.
    S = T.shape[1]
    cyc = S / 5894.0
    def win(c0, c1, pad=0.04):
        a = int(max(0, (c0 - pad * 5894) * cyc))
        b = int(min(S, (c1 + pad * 5894) * cyc))
        return (a, b)
    windows = {"m01": win(1153, 2713), "m11": win(4114, 5647),
               "m00": win(2713, 2961), "m10": win(5647, 5894)}
    print(f"  연산 구간 (샘플): " + ", ".join(f"{k}={v}" for k, v in windows.items()) + "\n")

    for mul_idx, add_idx, mul_name, add_name in ((1, 0, "m01", "m00"), (3, 2, "m11", "m10")):
        # --- the entry that multiplies achl ---
        print(f"  {mul_name} = Msk[{mul_idx // 2}][{mul_idx % 2}]  (multiply by achl)")
        limbs = []
        for i in range(args.limbs):
            truth = (secret[mul_idx] >> (32 * i)) & 0xFFFFFFFF
            base = mul_state(limbs, b, i, nlimb, width)
            b0 = b[:, 0]

            def predict(h, nbits, base=base, b0=b0):
                mask = np.uint64((1 << nbits) - 1)
                return (base[None, :] + ((h[:, None] * b0[None, :]) & MASK32)) & mask

            v = recover_limb(predict, Tc32, Tn32, Xc, args.keep, truth, f"limb {i}",
                                 stride=args.stride, window=windows[mul_name])
            limbs.append(v)
            if v != truth:
                break
        got[mul_name] = limbs

        # --- the entry that multiplies 1, reached through the addition ---
        if len(limbs) == args.limbs:
            P = full_product(limbs, b, nlimb, width)
            print(f"  {add_name} = Msk[{add_idx // 2}][{add_idx % 2}]  "
                  f"(constant multiplier; attacked through the add)")
            alimbs = []
            for k in range(args.limbs):
                truth = (secret[add_idx] >> (32 * k)) & 0xFFFFFFFF
                cin = add_carry_in(alimbs, P, k)
                Pk = P[:, k]

                def predict(h, nbits, Pk=Pk, cin=cin):
                    mask = np.uint64((1 << nbits) - 1)
                    return (h[:, None] + Pk[None, :] + cin[None, :]) & mask

                v = recover_limb(predict, Tc32, Tn32, Xc, args.keep, truth, f"limb {k}",
                                     stride=args.stride, window=windows[add_name])
                alimbs.append(v)
                if v != truth:
                    break
            got[add_name] = alimbs
        print()

    names = {"m00": 0, "m01": 1, "m10": 2, "m11": 3}
    total_bits = 0
    print(f"  {'entry':<6}{'limbs':>8}{'bits':>8}   value")
    print("  " + "-" * 70)
    for nm, idx in names.items():
        ls = got.get(nm, [])
        ok = sum(1 for i, v in enumerate(ls) if v == ((secret[idx] >> (32 * i)) & 0xFFFFFFFF))
        bits = 32 * ok
        total_bits += bits
        val = "".join(f"{v:08x}" for v in reversed(ls)) if ls else "-"
        print(f"  {nm:<6}{ok:>5}/{args.limbs}{bits:>8}   0x{val}")
    print("  " + "-" * 70)
    print(f"  전체: {total_bits}/512 bit 복구")
    return 0


if __name__ == "__main__":
    sys.exit(main())
