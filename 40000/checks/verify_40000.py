#!/usr/bin/env python3
"""Verify the 40,000-digit Golomb-Dickman certificate and OEIS b-file.

The check uses only the Python standard library and exact integer arithmetic.
It verifies that:

1. the certified interval does not cross a decimal boundary;
2. the digits in the certificate equal floor(10^D * lambda);
3. the standalone digits file equals the certificate digits;
4. the b-file has exactly D terms, indexed 0..D-1, and each term matches;
5. the first 102 digits match OEIS A084945.

Run from any directory:

    python verify_40000.py
"""

from __future__ import annotations

import hashlib
import math
import re
import sys
from pathlib import Path

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(200_000)

EXPECTED_DIGITS = 40_000
KNOWN_PREFIX = (
    "6243299885435508709929363831008372441796426201805292869735519024"
    "95638088855113254462460276195539868869"
)


def _bundle_root(start: Path) -> Path:
    """Locate the bundle root: the directory that holds the certificate.

    This lets every script be run from anywhere, including from ``checks/``.
    """
    for candidate in (start, *start.parents):
        if (candidate / "golomb_dickman_40000_certified.txt").is_file():
            return candidate
    return start

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def read_ascii(path: Path) -> str:
    raw = path.read_bytes()
    if b"\xef\xbb\xbf" in raw[:3]:
        fail(f"{path.name}: UTF-8 BOM is not allowed")
    try:
        return raw.decode("ascii")
    except UnicodeDecodeError as exc:
        fail(f"{path.name}: not plain ASCII ({exc})")


def main() -> None:
    base = _bundle_root(Path(__file__).resolve().parent)

    cert_path = base / "golomb_dickman_40000_certified.txt"
    digits_path = base / "golomb_dickman_40000_digits.txt"
    bfile_path = base / "b084945.txt"

    for path in (cert_path, digits_path, bfile_path):
        if not path.is_file():
            fail(f"missing file: {path}")

    cert_lines = read_ascii(cert_path).splitlines()
    if len(cert_lines) < 4:
        fail("certificate must contain at least 4 lines")

    match_s = re.search(r"\bS=(\d+)\b", cert_lines[1])
    match_lo = re.search(r"\blam_lo=(\d+)\b", cert_lines[1])
    match_hi = re.search(r"\blam_hi=(\d+)\b", cert_lines[2])
    if not (match_s and match_lo and match_hi):
        fail("cannot parse S, lam_lo, or lam_hi from the certificate")

    scale_bits = int(match_s.group(1))
    lam_lo = int(match_lo.group(1))
    lam_hi = int(match_hi.group(1))
    if not (0 < lam_lo < lam_hi):
        fail("expected 0 < lam_lo < lam_hi")

    cert_digits = cert_lines[3].strip()
    if not cert_digits.startswith("0."):
        fail("fourth certificate line must start with '0.'")
    cert_digits = cert_digits[2:]
    digits_count = len(cert_digits)
    if digits_count != EXPECTED_DIGITS:
        fail(f"expected {EXPECTED_DIGITS} digits, found {digits_count}")

    scale = 10**digits_count
    lower = (lam_lo * scale) >> scale_bits
    upper = (lam_hi * scale) >> scale_bits
    if lower != upper:
        fail("the certified interval crosses a decimal boundary")

    proven_digits = str(lower).zfill(digits_count)
    if proven_digits != cert_digits:
        fail("certificate digits do not equal floor(10^D * lambda)")

    standalone_digits = read_ascii(digits_path).strip()
    if standalone_digits != cert_digits:
        fail("digits file and certificate disagree")

    bfile_text = read_ascii(bfile_path)
    if "\t" in bfile_text:
        fail("b-file contains tab characters")
    if "\r" in bfile_text:
        fail("b-file contains CR characters; OEIS expects plain ASCII with LF")
    if not bfile_text.endswith("\n\n"):
        fail("b-file should end with a blank line")

    terms: list[tuple[int, int]] = []
    for line_number, line in enumerate(bfile_text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) != 2:
            fail(f"b-file line {line_number}: expected 'n value'")
        try:
            index, value = int(parts[0]), int(parts[1])
        except ValueError:
            fail(f"b-file line {line_number}: non-integer value")
        terms.append((index, value))

    if len(terms) != digits_count:
        fail(f"b-file has {len(terms)} terms, expected {digits_count}")

    for position, (index, value) in enumerate(terms):
        if index != position:
            fail(f"b-file index {index} at position {position}; expected {position}")
        if not 0 <= value <= 9:
            fail(f"b-file term a({index}) = {value} is not a decimal digit")
        if value != int(cert_digits[position]):
            fail(f"b-file term a({index}) does not match the certificate")

    if cert_digits[: len(KNOWN_PREFIX)] != KNOWN_PREFIX:
        fail("first 102 digits do not match OEIS A084945")

    width_bits = lam_hi - lam_lo
    # log10(width) = (log2(width) - S) * log10(2); log2() on a large int is exact
    # enough for one decimal place (it does not go through float conversion).
    width_log10 = (math.log2(width_bits) - scale_bits) * math.log10(2)

    print(f"OK: {digits_count} digits certified")
    print(f"  interval: lam_lo/2^{scale_bits} .. lam_hi/2^{scale_bits}")
    print(f"  interval width ~ 10^{width_log10:.1f}")
    print(f"  b-file: {len(terms)} terms, indices 0..{len(terms) - 1}")
    print(f"  first 102 digits match OEIS A084945")
    print("  SHA-256:")
    for path in (cert_path, digits_path, bfile_path):
        print(f"    {sha256(path)}  {path.name}")


if __name__ == "__main__":
    main()
