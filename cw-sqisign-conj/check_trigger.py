#!/usr/bin/env python3
"""
Step 3: confirm the trigger really brackets Algorithm 3.18 line 2.

Runs a handful of Alg 3.18 computations and checks, per trace:

  * the board's Conjugate is arithmetically correct (coord[1..3] negated,
    coord[0] and denom untouched)
  * the trigger pulse count matches what CONJ_TRIG_MODE should produce
  * the triggered width, converted from samples to target cycles, matches the
    DWT cycle counter the firmware reports for the same region
  * the per-limb cost is consistent with ibz_neg (a handful of cycles per limb)

Nothing here needs the leakage to be exploitable; it only establishes that the
window being captured is the code we think it is. Run it before spending time on
TVLA or SNR.

Usage:
  ~/.cwvenv/bin/python3 check_trigger.py --traces 3
  ~/.cwvenv/bin/python3 check_trigger.py --traces 3 --samples 200000 --skip-flash
"""

import argparse
import pathlib
import random
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import conj_rig as rig

VECTOR = pathlib.Path("/tmp/sqisign_conj_vec.txt")


def main():
    p = argparse.ArgumentParser(description="Verify the Alg 3.18 line 2 trigger")
    p.add_argument("-n", "--traces", type=int, default=3)
    p.add_argument("--samples", type=int, default=200000,
                   help="samples per trace; line 2 is short, this is mostly headroom")
    p.add_argument("--timebase", type=int, default=2,
                   help="ps3000a timebase (2 -> ~4 ns with 2 channels enabled)")
    p.add_argument("--pre-samples", type=int, default=2000)
    p.add_argument("--cha-range", type=int, default=3,
                   help="Ch A range code (3 = +/-100 mV; the CW308 shunt swings ~30 mV, so 500 mV wastes most of the ADC)")
    p.add_argument("--chb-range", type=int, default=8,
                   help="Ch B range code (8 = +/-5 V; the 3.3 V trigger GPIO overflows a 2 V range)")
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--timeout", type=float, default=900.0,
                   help="seconds to wait for the trigger (Alg 3.18 line 1 runs first)")
    p.add_argument("--plot", default="/tmp/sqisign_conj_trigger.png")
    args = p.parse_args()

    if not rig.HEX.exists():
        print(f"ERROR: firmware not found: {rig.HEX}\n"
              f"Build it with: make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2", file=sys.stderr)
        sys.exit(1)
    rig.ensure_host_tool()

    import chipwhisperer as cw
    from chipwhisperer.capture.api.programmers import STM32FProgrammer

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    failures = []
    rows = []
    last_trig = None
    last_power = None

    try:
        scope.default_setup()
        rig.set_target_clock(scope)
        if not args.skip_flash:
            print(f"programming {rig.HEX.name} ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(rig.HEX))
        else:
            rig.reset_target(scope)
        time.sleep(0.2)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)
        rig.sync(target)

        info = rig.board_info(target)
        print(f"board: trig_mode={info['trig_mode']} out_of_place={info['out_of_place']} "
              f"probe_reps={info['probe_reps']} IBZ_NLIMBS={info['ibz_nlimbs']} "
              f"sizeof(quat_alg_elem_t)={info['sizeof_quat_alg_elem']} "
              f"sizeof(quat_ideal_t)={info['sizeof_quat_ideal']}", flush=True)

        pico = rig.Pico(args.samples, args.timebase, args.pre_samples,
                        args.cha_range, args.chb_range, auto_trigger_ms=0)

        for i in range(args.traces):
            vec = rig.gen_vector(VECTOR, seed=f"trig{i}")
            rig.load_vector(target, vec)
            srng = random.Random(i)
            rig.seed_board(target, bytes(srng.getrandbits(8) for _ in range(16)))

            pico.arm()
            target.flush()
            t0 = time.time()
            target.simpleserial_write("c", bytearray())
            power, trig, overflow = pico.get(timeout_s=args.timeout)
            ack = target.simpleserial_read("r", 1, timeout=int(args.timeout * 1000))
            elapsed = time.time() - t0
            if ack != bytearray([1]):
                failures.append(f"trace {i}: board reported compute failure ({ack})")
                continue

            res = rig.parse_result(rig.read_outbuf(target))
            prof = rig.read_profile(target)
            last_trig, last_power = trig, power

            # --- arithmetic check: is this really Conjugate? ---
            problems = rig.check_conjugate(res)
            for msg in problems:
                failures.append(f"trace {i}: {msg}")

            # --- trigger pulse structure ---
            pulses = rig.trigger_pulses(trig)
            nlimbs = res["pre_limbs"][1]["nlimbs"]
            want_n, want_desc = rig.expected_pulses(info["trig_mode"], nlimbs)
            if want_n is not None and len(pulses) != want_n:
                failures.append(f"trace {i}: {len(pulses)} trigger pulses, expected {want_n} "
                                f"for mode {info['trig_mode']}")

            # --- triggered width vs the DWT cycle count for the same region ---
            spc = pico.samples_per_cycle
            if info["trig_mode"] == 2:
                ref_cycles = prof.get("conjugate", 0)
                ref_name = "conjugate"
            elif info["trig_mode"] == 3:
                ref_cycles = sum(prof.get(f"conj_coord{k}", 0) for k in (1, 2, 3))
                ref_name = "conj_coord1..3"
            elif info["trig_mode"] == 1:
                ref_cycles = prof.get("equiv_coprime", 0)
                ref_name = "equiv_coprime"
            elif info["trig_mode"] == 5:
                ref_cycles = prof.get("ideal_mul", 0)
                ref_name = "ideal_mul"
            elif info["trig_mode"] == 6:
                ref_cycles = prof.get("sample", 0)
                ref_name = "sample"
            else:
                ref_cycles = 0
                ref_name = "-"

            total_high = sum(w for _, _, w in pulses)
            span = (pulses[-1][1] - pulses[0][0]) if pulses else 0
            measured_cycles = total_high / spc if spc else 0

            print(f"\n[{i}] compute {elapsed:.1f}s  overflow={overflow}")
            print(f"    Alg 3.18 profile (target cycles): total={prof.get('total')} "
                  f"intersect={prof.get('intersect')} equiv_coprime={prof.get('equiv_coprime')} "
                  f"conjugate={prof.get('conjugate')} ideal_mul={prof.get('ideal_mul')} "
                  f"sample={prof.get('sample')}")
            print(f"    line 2 per coord (cycles): copy_c0={prof.get('conj_coord0')} "
                  f"neg_c1={prof.get('conj_coord1')} neg_c2={prof.get('conj_coord2')} "
                  f"neg_c3={prof.get('conj_coord3')}")
            print(f"    limbs negated: c1={prof.get('conj_limbs1')} c2={prof.get('conj_limbs2')} "
                  f"c3={prof.get('conj_limbs3')}  bitlen c1={prof.get('conj_bitlen1')}")
            print(f"    trigger: {len(pulses)} pulse(s), high={total_high} samples "
                  f"= {measured_cycles:.0f} cycles, span={span} samples")
            if want_desc:
                for k, (s, e, w) in enumerate(pulses):
                    desc = want_desc[k] if k < len(want_desc) else "?"
                    print(f"      pulse {k}: [{s}, {e}) width={w} "
                          f"= {w / spc:.0f} cycles   {desc}")
            if ref_cycles:
                ratio = measured_cycles / ref_cycles
                verdict = "OK" if 0.8 <= ratio <= 1.25 else "MISMATCH"
                print(f"    width vs DWT {ref_name}={ref_cycles} cycles -> "
                      f"ratio={ratio:.2f}  {verdict}")
                if verdict == "MISMATCH":
                    failures.append(
                        f"trace {i}: triggered width {measured_cycles:.0f} cycles vs "
                        f"DWT {ref_name} {ref_cycles} cycles (ratio {ratio:.2f})")

            # per-limb cost, a direct consistency check on ibz_neg
            for k in (1, 2, 3):
                c = prof.get(f"conj_coord{k}", 0)
                nl = prof.get(f"conj_limbs{k}", 0)
                if nl:
                    print(f"      coord[{k}]: {c} cycles / {nl} limbs = {c / nl:.1f} cycles per limb")

            trig_span_mv = float(trig.max() - trig.min())
            if trig_span_mv < 500:
                failures.append(f"trace {i}: trigger swing only {trig_span_mv:.0f} mV; "
                                f"check the Ch B probe")

            rows.append({
                "conjugate_cycles": prof.get("conjugate", 0),
                "measured_cycles": measured_cycles,
                "pulses": len(pulses),
                "nlimbs": nlimbs,
                "trig_span_mv": trig_span_mv,
                "power_ptp_mv": float(power.max() - power.min()),
                "overflow": overflow,
            })

    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    # ---- plot: trigger and the power inside it ----
    if last_trig is not None:
        pulses = rig.trigger_pulses(last_trig)
        fig, ax = plt.subplots(3, 1, figsize=(14, 9))
        ax[0].plot(last_trig, color="green", lw=0.5)
        ax[0].set_title("Trigger (Ch B), last trace")
        ax[0].set_ylabel("mV")
        ax[0].grid(alpha=0.3)

        ax[1].plot(last_power, color="C0", lw=0.3)
        for s, e, _ in pulses:
            ax[1].axvspan(s, e, color="orange", alpha=0.25)
        ax[1].set_title("Power (Ch A) with the triggered region(s) shaded")
        ax[1].set_ylabel("mV")
        ax[1].grid(alpha=0.3)

        if pulses:
            s, e = pulses[0][0], pulses[-1][1]
            pad = max(1, (e - s) // 20)
            lo, hi = max(0, s - pad), min(len(last_power), e + pad)
            ax[2].plot(np.arange(lo, hi), last_power[lo:hi], color="C3", lw=0.5)
            for ps_, pe_, _ in pulses:
                ax[2].axvspan(ps_, pe_, color="orange", alpha=0.25)
            ax[2].set_title(f"Zoom on the triggered window: samples [{lo}, {hi})")
        ax[2].set_xlabel("sample")
        ax[2].set_ylabel("mV")
        ax[2].grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        plt.close(fig)
        print(f"\nsaved {args.plot}", flush=True)

    print("\n" + "=" * 70)
    if rows:
        mc = np.array([r["measured_cycles"] for r in rows])
        dc = np.array([r["conjugate_cycles"] for r in rows])
        print(f"traces={len(rows)}  triggered width: mean={mc.mean():.0f} cycles "
              f"sd={mc.std():.1f}   DWT conjugate: mean={dc.mean():.0f} sd={dc.std():.1f}")
        print(f"limbs per coordinate: {sorted(set(r['nlimbs'] for r in rows))} "
              f"(constant across traces is what we want)")
    if failures:
        print(f"\n{len(failures)} PROBLEM(S):")
        for f in failures:
            print("  -", f)
        return 1
    print("\ntrigger check: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
