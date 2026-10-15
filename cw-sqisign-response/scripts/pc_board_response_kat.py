#!/usr/bin/env python3
import argparse
import pathlib
import subprocess
import time

import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer


APP = pathlib.Path(__file__).resolve().parent
PKG = APP.parent
FW_DIR = PKG / "firmware"
HOST = PKG / "host" / "response_split_host"
HOST_BUILD = PKG / "host" / "build_host_tool.sh"
HEX = FW_DIR / "simpleserial-sqisign-response-CW308_STM32F4.hex"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
BOARD_RESPONSE = pathlib.Path("/tmp/sqisign_response_board.txt")
CHUNK = 48


def run(cmd, **kwargs):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def ss_cmd(target, cmd, data=b"", pay_len=1, timeout_ms=250000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", pay_len, timeout=timeout_ms)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-flash", action="store_true")
    parser.add_argument("--platform", default="CW308_STM32F4")
    parser.add_argument("--firmware", type=pathlib.Path, default=HEX)
    args = parser.parse_args()
    firmware = args.firmware

    if not HOST.exists():
        run([HOST_BUILD])

    with VECTOR.open("wb") as f:
        run([HOST, "gen"], stdout=f)
    vector = VECTOR.read_bytes()
    print(f"vector bytes={len(vector)} file={VECTOR}", flush=True)

    if not firmware.exists():
        run(["make", f"PLATFORM={args.platform}"], cwd=FW_DIR)

    scope = cw.scope()
    target = cw.target(scope)
    try:
        scope.default_setup()
        if not args.skip_flash:
            cw.program_target(scope, STM32FProgrammer, str(firmware))

        time.sleep(0.2)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)

        r = ss_cmd(target, "x", b"", 1, timeout_ms=5000)
        print("reset:", r, flush=True)
        if r != bytearray([1]):
            raise RuntimeError("board reset command failed")

        padded = vector
        if len(padded) % CHUNK:
            padded += b"\n" * (CHUNK - (len(padded) % CHUNK))
        for off in range(0, len(padded), CHUNK):
            r = ss_cmd(target, "a", padded[off : off + CHUNK], 1, timeout_ms=5000)
            if r != bytearray([1]):
                raise RuntimeError(f"write chunk failed at offset {off}: {r}")
        print(f"wrote chunks={len(padded) // CHUNK}", flush=True)

        t0 = time.time()
        r = ss_cmd(target, "c", b"", 1, timeout_ms=900000)
        elapsed = time.time() - t0
        print(f"compute: {r} elapsed={elapsed:.3f}s", flush=True)
        if r != bytearray([1]):
            raise RuntimeError("compute command failed")

        r = ss_cmd(target, "p", b"", 1, timeout_ms=5000)
        if r == bytearray([1]):
            r = ss_cmd(target, "l", b"", 2, timeout_ms=5000)
            prof_len = (r[0] << 8) | r[1]
            prof = bytearray()
            for idx in range((prof_len + CHUNK - 1) // CHUNK):
                r = ss_cmd(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, timeout_ms=5000)
                prof.extend(r)
            del prof[prof_len:]
            print("profile:", flush=True)
            print(prof.decode("ascii", errors="replace"), flush=True)
            r = ss_cmd(target, "s", b"", 1, timeout_ms=5000)
            if r != bytearray([1]):
                raise RuntimeError("serialize-last command failed after profile read")

        r = ss_cmd(target, "l", b"", 2, timeout_ms=5000)
        outlen = (r[0] << 8) | r[1]
        print(f"outlen={outlen}", flush=True)

        out = bytearray()
        for idx in range((outlen + CHUNK - 1) // CHUNK):
            r = ss_cmd(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, timeout_ms=5000)
            out.extend(r)
        del out[outlen:]
        BOARD_RESPONSE.write_bytes(out)
        print(f"response bytes={len(out)} file={BOARD_RESPONSE}", flush=True)
        print(out.decode("ascii", errors="replace"), flush=True)
    finally:
        try:
            target.dis()
        finally:
            scope.dis()

    run([HOST, "verify", VECTOR, BOARD_RESPONSE])


if __name__ == "__main__":
    main()
