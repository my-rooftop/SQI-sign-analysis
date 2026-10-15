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


def public_projector(b):
    """Projector onto the span of the public coefficient's Hamming weights.

    The hypothesis a = 1 predicts low32(1 * b0) = b0, so its model is just the
    Hamming weight of the public coefficient -- and HW(b0) leaks hard on its own
    (|corr| 0.43 at one sample), carrying no secret. Left in, that degenerate
    hypothesis beats the real one on roughly a third of trace ranges. Projecting
    the public Hamming weights out of both the traces and the models removes it:
    anything proportional to them scores zero.
    """
    n = b.shape[0]
    cols = [np.ones(n)]
    for j in range(b.shape[1]):
        cols.append(hw_rows(b[:, j][None, :], 4)[0].astype(np.float64))
    C = np.vstack(cols).T
    return C @ np.linalg.pinv(C)


def residualise(M, P):
    return M - (M @ P.T)


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
                 stride=1, window=None, P=None):
    """Two rounds: low 16 bits over all samples, then full 32 at the chosen POI.

    `predict(hyps, nbits)` returns the modelled intermediate, masked to nbits,
    for each hypothesis of the limb's low bits.
    """
    hyps = np.arange(1 << 16, dtype=np.uint64)
    models = np.empty((len(hyps), Tc32.shape[0]), dtype=np.float32)
    for s in range(0, len(hyps), batch):
        M = hw_rows(predict(hyps[s:s + batch], 16), 2).astype(np.float32)
        models[s:s + batch] = residualise(M, P) if P is not None else M
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
            if P is not None:
                M = residualise(M, P)
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
    p.add_argument("--keep", type=int, default=32)
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
    # Operation boundaries measured directly, from the per-multiply trigger
    # pulses of CONJ_TRIG_MODE=11, rather than estimated from the cycle profile.
    # The estimate had the m10 window overlapping m11's multiply, so the m10
    # attack picked up m11's leak instead of its own.
    S = T.shape[1]
    if "pulse_bounds" in d and len(d["pulse_bounds"]) == 4:
        pb = d["pulse_bounds"]
        mul = [(int(a), int(b_)) for a, b_ in pb]
        windows = {
            "m01": mul[1],                       # ibz_mul(m01, achl)
            "m11": mul[3],                       # ibz_mul(m11, achl)
            "m00": (mul[1][1], mul[2][0]),       # add matvec0 += m01*achl
            "m10": (mul[3][1], S),               # add matvec1 += m11*achl
        }
        src = "measured from the trigger pulses"
    else:
        cyc = S / 5894.0
        windows = {"m01": (int(1153 * cyc), int(2713 * cyc)),
                   "m11": (int(4114 * cyc), int(5647 * cyc)),
                   "m00": (int(2713 * cyc), int(2961 * cyc)),
                   "m10": (int(5647 * cyc), S)}
        src = "estimated from cycle counts"
    print(f"  operation windows ({src}): " + ", ".join(f"{k}={v}" for k, v in windows.items()))
    bad = [k for k, (a, b_) in windows.items() if b_ - a < 100]
    if bad:
        print(f"  warning: windows too narrow {bad}")
    print(f"  public-value HW projection removed: on\n")
    PROJ = public_projector(b)  # projector for the public value (named apart from the product P)

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
                                 stride=args.stride, window=windows[mul_name], P=PROJ)
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
                                     stride=args.stride, window=windows[add_name], P=PROJ)
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
    print(f"  total: {total_bits}/512 bits recovered")
    return 0


if __name__ == "__main__":
    sys.exit(main())
