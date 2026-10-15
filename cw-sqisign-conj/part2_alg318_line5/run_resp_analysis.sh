#!/bin/bash
# Run every section 5.3.2 analysis on one line 5 capture, in the order the
# report presents them. Each stage writes its own log so a later stage failing
# does not cost the earlier results.
#
#   ./run_resp_analysis.sh [npz] [outdir]
set -u
NPZ=${1:-/tmp/sqisign_resp.npz}
OUT=${2:-/tmp/resp_analysis}
PY=~/.cwvenv/bin/python3
HERE=$(cd "$(dirname "$0")" && pwd)
COMMON=$HERE/../common
TVLA=${TVLA_NPZ:-/tmp/sqisign_resp_tvla.npz}

mkdir -p "$OUT"
[ -f "$NPZ" ] || { echo "no capture at $NPZ" >&2; exit 1; }

$PY - "$NPZ" <<'PY'
import sys, numpy as np
d = np.load(sys.argv[1], allow_pickle=False)
T = d["traces"]
print(f"capture: {T.shape[0]} traces x {T.shape[1]} samples, "
      f"trig_mode={int(d['trig_mode'])}, {float(d['interval_ns']):.3f} ns/sample")
for c in range(4):
    m = d[f"c{c}_magnitude_bits"]
    print(f"  coord[{c}]: magnitude {m.min()}-{m.max()} bits "
          f"(median {int(np.median(m))}), {int(np.median(d[f'c{c}_nlimbs']))} limbs stored")
tot = sum(int(np.median(d[f"c{c}_magnitude_bits"])) for c in range(4)) + 4
print(f"  secret to infer: about {tot} bits (magnitudes + 4 signs)")
PY

run () {  # run <name> <command...>
  local name=$1 log="$OUT/$1.log"; shift
  echo "=== $name ==="
  if "$@" > "$log" 2>&1; then
    echo "  ok -> $log"
  else
    echo "  FAILED (exit $?) -> $log"
    tail -5 "$log" | sed 's/^/    /'
  fi
}

# 1. leakage detection. TVLA needs its own fixed-vs-random capture; SNR and the
#    per-limb correlation come off the random campaign.
if [ -f "$TVLA" ]; then
  run tvla $PY "$COMMON/analyze_conj_leakage.py" --tvla "$TVLA" \
      --tvla-plot "$OUT/tvla.png"
else
  echo "=== tvla ==="
  echo "  skipped: no fixed-vs-random capture at $TVLA"
  echo "  make one with: capture_resp_pico.py --records ... --mode tvla -n 1000 --out $TVLA"
fi
run snr $PY "$COMMON/analyze_conj_leakage.py" --snr "$NPZ" --snr-plot "$OUT/snr.png"

# 2. what the leakage is worth: HW templates, per-bit LDA, and the bit budget.
run template $PY "$COMMON/recover_conj_template.py" --npz "$NPZ"

# 3. the question the paper actually asks: how many leading bits come out.
run msb $PY "$COMMON/msb_recover.py" --npz "$NPZ" --widths 1,2,3,4,6,8

# 4. does a learned model beat the Gaussian templates.
run dl $PY "$COMMON/dl_recover.py" --npz "$NPZ"

echo
echo "logs and figures in $OUT"
