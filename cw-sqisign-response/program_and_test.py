#!/usr/bin/env python3
import time
from pathlib import Path

import chipwhisperer as cw


FW = Path(__file__).with_name("simpleserial-sqisign-response-CW308_STM32F4.hex")


def main() -> int:
    scope = cw.scope()
    try:
        scope.default_setup()
        cw.program_target(scope, cw.programmers.STM32FProgrammer, str(FW))

        target = cw.target(scope)
        scope.io.nrst = "low"
        time.sleep(0.05)
        scope.io.nrst = "high_z"
        time.sleep(0.5)
        target.flush()

        target.simpleserial_write("i", bytearray())
        resp = target.simpleserial_read("r", 8, timeout=5000)
        print(f"i: {None if resp is None else list(resp)}")
    finally:
        try:
            scope.dis()
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
