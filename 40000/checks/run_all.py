#!/usr/bin/env python3
"""Run every check in this folder, one after another, and print one summary.

    python run_all.py            # fast suite (about two minutes)
    python run_all.py --full     # also the exhaustive constant-prime scan (~3 min)
    python run_all.py --arb      # also the Arb audit (needs python-flint)
    python run_all.py --all      # everything
    python run_all.py --list     # only list the steps

Each step runs in a separate process, so a missing optional dependency
(``python-flint``) only fails that one step.  The exit status is 0 when every
step that was run succeeded.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

FAST = [
    ("artifacts", ["verify_40000.py"],
     "certificate, digits file, b-file, OEIS prefix"),
    ("checkpoints", ["verify_checkpoints.py"],
     "1209 chunks: structure, coverage, exact sums"),
    ("formulas", ["verify_nodes_decimal.py", "--methods"],
     "three independent formulas for E1"),
    ("constants", ["verify_nodes_decimal.py", "--plan"],
     "B_disc, tau_R, tau_L and the certificate corrections"),
    ("chunks", ["verify_nodes_decimal.py", "--chunks", "12"],
     "12 spread-out chunks recomputed at 100 decimals"),
    ("mathworld", ["verify_mathworld.py", "--prime-limit", "60"],
     "every statement of the MathWorld page (fast prime scan)"),
]

FULL = [
    ("primes-to-3000", ["verify_mathworld.py", "--prime-limit", "3000"],
     "exhaustive constant-prime search to 3000 digits"),
]

ARB = [
    ("arb-audit", ["verify_checkpoints_arb.py", "--audit"],
     "Arb: constants and certificate (needs python-flint)"),
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true", help="add the exhaustive prime scan")
    parser.add_argument("--arb", action="store_true", help="add the Arb audit")
    parser.add_argument("--all", action="store_true", help="add both")
    parser.add_argument("--list", action="store_true", help="list the steps and exit")
    args = parser.parse_args(argv)

    steps = list(FAST)
    if args.full or args.all:
        steps += FULL
    if args.arb or args.all:
        steps += ARB

    if args.list:
        for name, cmd, note in steps:
            print(f"  {name:<16} python {' '.join(cmd):<42} {note}")
        return 0

    print(f"bundle root: {HERE.parent}", flush=True)
    print(f"python:      {sys.version.split()[0]}", flush=True)
    print()

    results = []
    for name, cmd, note in steps:
        print(f"=== {name} :: {note}", flush=True)
        started = time.perf_counter()
        code = subprocess.run([sys.executable, *cmd], cwd=str(HERE)).returncode
        elapsed = time.perf_counter() - started
        results.append((name, code, elapsed))
        print(f"--- {name}: {'OK' if code == 0 else 'FAILED'} in {elapsed:.1f} s", flush=True)
        print()

    print("summary")
    for name, code, elapsed in results:
        print(f"  {name:<16} {'OK' if code == 0 else 'FAILED':<7} {elapsed:6.1f} s")

    failed = [name for name, code, _ in results if code != 0]
    if failed:
        print()
        print("FAILED steps: " + ", ".join(failed), file=sys.stderr)
        return 1
    print()
    print(f"OK: all {len(results)} steps passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
