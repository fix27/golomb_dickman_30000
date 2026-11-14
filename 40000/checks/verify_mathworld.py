#!/usr/bin/env python3
"""Check every numerical statement of the MathWorld page
"Golomb-Dickman Constant Digits" against the certified 40,000-digit expansion.

Source of the statements
------------------------
https://mathworld.wolfram.com/Golomb-DickmanConstantDigits.html  (accessed 2026-11-15)

Cross-references used only to label the claims: OEIS A084945 (digits),
A174974 / A174975 (constant primes), A225242 (Earls sequence), A229195
(first occurrence of each digit).

Only the Python standard library is used, so the check runs anywhere:

    python verify_mathworld.py
    python verify_mathworld.py --root .. --prime-limit 3000 --json report.json
    python verify_mathworld.py --prime-limit 60          # fast run

Exit status is 0 when no statement is contradicted (NOTEs and SKIPs are
allowed) and 1 as soon as one statement FAILs.

Status codes
------------
PASS  the certified digits reproduce the statement exactly
FAIL  the certified digits contradict the statement
NOTE  historical, contextual or open statement; recorded, not machine-checkable
SKIP  machine-checkable but not run at the configured limit (reason given)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Callable

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(400_000)

CERT_NAME = "golomb_dickman_40000_certified.txt"
DIGITS_NAME = "golomb_dickman_40000_digits.txt"

# ---------------------------------------------------------------- MathWorld data

MW_PREFIX = "6243299885"                      # "lambda = 0.6243299885..."
MW_COMPUTED_2013 = 15000                      # "computed to 15000 decimal digits"

# "The starting positions of the first occurrence of 0, 1, 2, ... are ..."
MW_FIRST_OCCURRENCE = {0: 15, 1: 28, 2: 2, 3: 4, 4: 3, 5: 10, 6: 1, 7: 17, 8: 8, 9: 6}

# "the last strings for n = 1, 2, ... are 1, 33, 821, ..., ending at 28, 587, 6322"
MW_LAST_STRING = {
    1: ("1", 28),
    2: ("33", 587),
    3: ("821", 6322),
}

# Earls sequence, A225242: first occurrence of n copies of the digit n
MW_EARLS = {1: 28, 2: 256, 3: 1967, 4: 387}

# Lengths (in decimal digits) for which the truncation of lambda is prime,
# A174974.  The two smallest primes themselves are quoted by MathWorld.
MW_PRIME_LENGTHS = [6, 27, 57, 60, 1659, 2508]
MW_PRIME_VALUES = {
    6: "624329",
    27: "624329988543550870992936383",
}

# "reported counts of each decimal digit d in the first 10^n decimal digits"
MW_DIGIT_COUNTS = {
    10: [0, 0, 2, 1, 1, 1, 1, 0, 2, 2],
    100: [9, 7, 12, 10, 9, 12, 10, 5, 14, 12],
    1000: [89, 108, 93, 94, 100, 98, 104, 96, 109, 109],
    10000: [987, 999, 996, 989, 1021, 983, 1042, 995, 993, 995],
}

# "0123456789 and 9876543210 do not occur in the first 15000 digits"
MW_ABSENT = ("0123456789", "9876543210")

# 102 digits of OEIS A084945 (the page attributes the expansion to it).
OEIS_PREFIX = (
    "6243299885435508709929363831008372441796426201805292869735519024"
    "95638088855113254462460276195539868869"
)

MR_DETERMINISTIC = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
MR_DETERMINISTIC_LIMIT = 3_317_044_064_679_887_385_961_981  # 3.3 * 10^24
MR_MANY = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
)


# ---------------------------------------------------------------- small helpers

def bundle_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / CERT_NAME).is_file():
            return candidate
    return start


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def load_digits(root: Path) -> str:
    """The certified digits, or the standalone digits file if that is absent."""
    digits_path = root / DIGITS_NAME
    if digits_path.is_file():
        return digits_path.read_text(encoding="ascii").strip()
    cert = root / CERT_NAME
    if not cert.is_file():
        fail(f"missing {CERT_NAME} in {root}")
    lines = cert.read_text(encoding="ascii").splitlines()
    if len(lines) < 4 or not lines[3].strip().startswith("0."):
        fail("cannot find the digit line in the certificate")
    return lines[3].strip()[2:]


def sieve(limit: int) -> list[int]:
    flags = bytearray([1]) * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for i in range(2, int(limit**0.5) + 1):
        if flags[i]:
            flags[i * i:: i] = bytearray(len(flags[i * i:: i]))
    return [i for i in range(limit + 1) if flags[i]]


def is_probable_prime(n: int, bases: tuple[int, ...], small: list[int]) -> bool:
    """Miller-Rabin.  Exact when n < MR_DETERMINISTIC_LIMIT and bases cover it."""
    if n < 2:
        return False
    for p in small:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in bases:
        a %= n
        if a < 2:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


class Report:
    def __init__(self) -> None:
        self.rows: list[dict[str, Any]] = []

    def add(self, code: str, claim: str, status: str, detail: str,
            extra: list[str] | None = None) -> None:
        self.rows.append(
            {"code": code, "claim": claim, "status": status,
             "detail": detail, "extra": extra or []}
        )

    def failures(self) -> int:
        return sum(1 for row in self.rows if row["status"] == "FAIL")


# ---------------------------------------------------------------- the checks

def check_c1(report: Report, digits: str) -> None:
    """lambda = 0.6243299885..., OEIS A084945."""
    claim = "decimal expansion begins 0.6243299885... (OEIS A084945)"
    ok = digits.startswith(MW_PREFIX)
    detail = f"first 20 digits: 0.{digits[:20]}"
    report.add("C1", claim, "PASS" if ok else "FAIL", detail)

    ok53 = digits[:53] == OEIS_PREFIX[:53]
    ok102 = digits[:102] == OEIS_PREFIX
    report.add(
        "C1b", "the expansion is the one of OEIS A084945",
        "PASS" if ok102 else "FAIL",
        f"102-digit prefix {'matches' if ok102 else 'differs'}; "
        f"first 53 digits {'match' if ok53 else 'differ'}",
    )


def check_c2(report: Report, digits: str) -> None:
    """Positions count from 1 at the first decimal digit."""
    claim = "positions count from 1 at the first digit after the point"
    bad = []
    for digit, position in MW_FIRST_OCCURRENCE.items():
        text = str(digit)
        if digits[position - 1] != text or text in digits[: position - 1]:
            bad.append(digit)
    report.add(
        "C2", claim, "PASS" if not bad else "FAIL",
        "1-based indexing reproduces all ten first-occurrence positions"
        if not bad else f"inconsistent for digits {bad}",
    )


def check_c3(report: Report, digits: str) -> None:
    """Mitchell (1968): 53 decimal places."""
    claim = "Mitchell (1968) computed lambda to 53 decimal places"
    report.add(
        "C3", claim, "NOTE",
        "historical attribution; not machine-checkable. "
        f"Our digits 1..53 are 0.{digits[:53]} and agree with OEIS A084945",
    )


def check_c4(report: Report, digits: str) -> None:
    """15000 digits computed by Weisstein (2013)."""
    claim = f"lambda computed to {MW_COMPUTED_2013} decimal digits (2013)"
    report.add(
        "C4", claim, "NOTE",
        f"historical; this bundle certifies {len(digits)} digits, so"
        f" {len(digits) - MW_COMPUTED_2013} digits beyond that range are new"
        if len(digits) > MW_COMPUTED_2013 else "historical; not machine-checkable",
    )


def check_c5(report: Report, digits: str) -> None:
    """Earls sequence A225242."""
    claim = "Earls sequence begins 28, 256, 1967, 387 (A225242)"
    bad = []
    found: dict[int, int] = {}
    for n, expected in MW_EARLS.items():
        pattern = str(n) * n
        index = digits.find(pattern)
        found[n] = index + 1 if index >= 0 else 0
        if found[n] != expected:
            bad.append(n)
    extra = []
    for n in (5, 6):
        pattern = str(n) * n
        index = digits.find(pattern)
        if index < 0:
            extra.append(f"n={n}: {pattern} not found in {len(digits)} digits")
        else:
            extra.append(f"n={n}: {pattern} first at position {index + 1}"
                         f" (beyond the 15000 digits known in 2013)")
    report.add(
        "C5", claim, "PASS" if not bad else "FAIL",
        "all four positions confirmed: "
        + ", ".join(f"n={n}: {found[n]}" for n in sorted(found))
        if not bad else f"wrong position for n={bad}",
        extra,
    )


def check_c6(report: Report, digits: str, prime_limit: int, small: list[int],
             blocks: dict[str, Any]) -> None:
    """Lambda-constant primes, A174974 / A174975."""
    claim = "constant primes at 6, 27, 57, 60, 1659, 2508 digits (A174974)"
    detail_bits = []
    ok = True
    for k in MW_PRIME_LENGTHS:
        value = int(digits[:k])
        bases = MR_MANY if value >= MR_DETERMINISTIC_LIMIT else MR_DETERMINISTIC
        prime = is_probable_prime(value, bases, small)
        label = "prime" if value < MR_DETERMINISTIC_LIMIT else "probable prime (64 bases)"
        detail_bits.append(f"{k}: {label}" if prime else f"{k}: COMPOSITE")
        if not prime:
            ok = False
        if k in MW_PRIME_VALUES:
            if digits[:k] != MW_PRIME_VALUES[k]:
                ok = False
                detail_bits.append(f"{k}:value differs from MathWorld")
    report.add("C6", claim, "PASS" if ok else "FAIL",
               "; ".join(detail_bits))

    # Exhaustive search below the configured limit.
    if prime_limit <= 0:
        report.add("C6b", "no other prime truncation below the limit", "SKIP",
                   "exhaustive search disabled")
        return
    found, seconds = blocks["scan"](prime_limit)
    expected = [k for k in MW_PRIME_LENGTHS if k <= prime_limit]
    status = "PASS" if found == expected else "FAIL"
    note = "" if found == expected else f" (MathWorld: {expected})"
    report.add(
        "C6b", f"exhaustively, the prime truncations up to {prime_limit} digits are {expected}",
        status,
        f"found {found} in {seconds:.1f} s{note}",
    )


def check_c7(report: Report, digits: str) -> None:
    """A229195."""
    claim = "first occurrence of 0,1,...,9 at 15, 28, 2, 4, 3, 10, 1, 17, 8, 6"
    got = {}
    for digit in range(10):
        position = digits.find(str(digit))
        got[digit] = position + 1 if position >= 0 else 0
    ok = all(got[d] == MW_FIRST_OCCURRENCE[d] for d in range(10))
    report.add(
        "C7", claim, "PASS" if ok else "FAIL",
        "recomputed: " + ", ".join(str(got[d]) for d in range(10)),
    )


def check_c8(report: Report, digits: str) -> None:
    """Coverage of all strings of length n."""
    claim = "last strings 1, 33, 821 end at 28, 587, 6322"
    lines = []
    ok = True
    for n, (last_expected, end_expected) in MW_LAST_STRING.items():
        need = 10**n
        seen: set[str] = set()
        last = ""
        end = -1
        for i in range(len(digits) - n + 1):
            piece = digits[i:i + n]
            if piece not in seen:
                seen.add(piece)
                if len(seen) == need:
                    last = piece
                    end = i + n
                    break
        lines.append(f"n={n}: last {last!r} at {end}")
        if (last, end) != (last_expected, end_expected):
            ok = False
    # n = 4 is reported for context only: 10^4 strings need ~9.2e4 digits.
    n = 4
    seen4: set[str] = set()
    for i in range(len(digits) - n + 1):
        seen4.add(digits[i:i + n])
    lines.append(f"n=4: {len(seen4)} of 10000 strings seen (not complete)")
    report.add("C8", claim, "PASS" if ok else "FAIL", "; ".join(lines))


def check_c9(report: Report, digits: str) -> None:
    """Absence of 0123456789 and 9876543210."""
    claim = f"{MW_ABSENT[0]} and {MW_ABSENT[1]} do not occur in {MW_COMPUTED_2013} digits"
    window = digits[:MW_COMPUTED_2013]
    lines = []
    ok = True
    for pattern in MW_ABSENT:
        pos = window.find(pattern)
        pos_all = digits.find(pattern)
        if pos >= 0:
            lines.append(f"{pattern}: found at {pos + 1} in {len(window)} digits")
            ok = False
        elif pos_all >= 0:
            lines.append(f"{pattern}: absent in {len(window)} digits, "
                         f"first found at {pos_all + 1} in all {len(digits)}")
        elif len(digits) > len(window):
            lines.append(f"{pattern}: absent in {len(window)} and in all {len(digits)} digits")
        else:
            lines.append(f"{pattern}: absent in {len(window)} digits")
    report.add("C9", claim, "PASS" if ok else "FAIL", "; ".join(lines))


def check_c10(report: Report, digits: str) -> None:
    """Normality is open; give the empirical picture."""
    claim = "it is not known whether lambda is normal in base 10"
    counts = [digits.count(str(d)) for d in range(10)]
    expected = len(digits) / 10
    chi2 = sum((c - expected) ** 2 / expected for c in counts)
    report.add(
        "C10", claim, "NOTE",
        f"open problem; digit counts in {len(digits)} digits: {counts}, "
        f"chi-squared = {chi2:.2f} with 9 degrees of freedom "
        f"(mean 4000, 5% critical value 16.92)",
    )


def check_c11(report: Report, digits: str) -> None:
    """The count table for 10, 100, 1000, 10000 digits."""
    claim = "digit counts in the first 10, 100, 1000, 10000 digits"
    mismatches = []
    total = 0
    for size, expected in MW_DIGIT_COUNTS.items():
        window = digits[:size]
        for digit in range(10):
            total += 1
            if window.count(str(digit)) != expected[digit]:
                mismatches.append(f"d={digit}, n={size}: "
                                  f"{window.count(str(digit))} != {expected[digit]}")
    ok = not mismatches
    report.add(
        "C11", claim, "PASS" if ok else "FAIL",
        f"all {total} table entries reproduced exactly" if ok
        else "; ".join(mismatches[:8]),
    )


def check_c12(report: Report, digits: str) -> None:
    """A174975: the primes themselves."""
    claim = ("first constant primes are 624329 and "
             "624329988543550870992936383")
    ok = all(digits[:k] == value for k, value in MW_PRIME_VALUES.items())
    report.add(
        "C12", claim, "PASS" if ok else "FAIL",
        f"6-digit and 27-digit truncations "
        f"{'match' if ok else 'differ'}: {digits[:6]}, {digits[:27]}",
    )


def check_c13(report: Report, digits: str, prime_limit: int) -> None:
    """Search for primes through 15000 digits."""
    claim = f"completed search for primes through {MW_COMPUTED_2013} digits"
    if prime_limit >= MW_COMPUTED_2013:
        report.add("C13", claim, "PASS", "exhaustive search reached the stated bound")
        return
    report.add(
        "C13", claim, "SKIP",
        f"exhaustive search run to {prime_limit} digits "
        f"(--prime-limit {MW_COMPUTED_2013} to reproduce the full claim; "
        f"the 15000-digit range needs hours)",
    )


def check_c14(report: Report, digits: str) -> None:
    """The discoverer table."""
    claim = "largest reported: 1659 (Broadhurst 2010), 2508 (Weisstein 2010)"
    report.add(
        "C14", claim, "NOTE",
        "attribution is historical; the two lengths are prime in this "
        "expansion (see C6)",
    )


# ---------------------------------------------------------------- driver

def build_scan(digits: str, small: list[int]) -> Callable[[int], tuple[list[int], float]]:
    """Return a function that finds all prime truncations up to a limit."""
    def scan(limit: int) -> tuple[list[int], float]:
        import time
        start = time.perf_counter()
        trial = sieve(5000)
        residues = [0] * len(trial)
        found: list[int] = []
        value = 0
        for k, ch in enumerate(digits[:limit], start=1):
            digit = ord(ch) - 48
            value = value * 10 + digit
            for i, p in enumerate(trial):
                residues[i] = (residues[i] * 10 + digit) % p
            if value < 5000:
                if value > 9 and is_probable_prime(value, MR_DETERMINISTIC, small):
                    found.append(k)
                continue
            if digit not in (1, 3, 7, 9):
                continue
            if 0 in residues:
                continue
            if not is_probable_prime(value, (2,), small):
                continue
            bases = MR_MANY if value >= MR_DETERMINISTIC_LIMIT else MR_DETERMINISTIC
            if is_probable_prime(value, bases, small):
                found.append(k)
        return found, time.perf_counter() - start
    return scan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify the MathWorld 'Golomb-Dickman Constant Digits' page "
                    "against the certified expansion.",
    )
    parser.add_argument("--root", type=Path, default=bundle_root(Path(__file__).resolve().parent))
    parser.add_argument("--prime-limit", type=int, default=3000,
                        help="exhaustive constant-prime search up to this many digits (default 3000)")
    parser.add_argument("--json", type=Path, default=None, help="write a JSON report")
    args = parser.parse_args(argv)

    digits = load_digits(args.root)
    small = sieve(1000)
    blocks = {"scan": build_scan(digits, small)}

    print('MathWorld, "Golomb-Dickman Constant Digits"')
    print("  https://mathworld.wolfram.com/Golomb-DickmanConstantDigits.html")
    print(f"  bundle: {args.root}")
    print(f"  digits checked after the decimal point: {len(digits)}")
    print()

    report = Report()
    check_c1(report, digits)
    check_c2(report, digits)
    check_c3(report, digits)
    check_c4(report, digits)
    check_c5(report, digits)
    check_c6(report, digits, args.prime_limit, small, blocks)
    check_c7(report, digits)
    check_c8(report, digits)
    check_c9(report, digits)
    check_c10(report, digits)
    check_c11(report, digits)
    check_c12(report, digits)
    check_c13(report, digits, args.prime_limit)
    check_c14(report, digits)

    for row in report.rows:
        print(f"  [{row['code']:>4s}] {row['status']:<4s}  {row['claim']}")
        print(f"          {row['detail']}")
        for line in row["extra"]:
            print(f"          {line}")

    counts: dict[str, int] = {}
    for row in report.rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print()
    print("  summary: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))

    if args.json is not None:
        args.json.write_text(
            json.dumps({"source": "MathWorld Golomb-Dickman Constant Digits",
                        "digits": len(digits), "claims": report.rows}, indent=2) + "\n",
            encoding="ascii",
        )
        print(f"  report written to {args.json}")

    if report.failures():
        print(f"FAIL: {report.failures()} statement(s) contradicted", file=sys.stderr)
        return 1
    print("OK: no statement of the MathWorld page is contradicted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
