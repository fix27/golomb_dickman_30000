#!/usr/bin/env python3
"""Independent auditor for the 40,000-digit Golomb-Dickman checkpoint set.

This script uses only the Python standard library.  It does NOT trust the
stored chunk intervals as proof of the constant: it re-reads every pickle,
validates the chunk structure and coverage, recomputes the exact integer sums
``Lsum`` and ``Hsum``, and checks that those sums, together with the certificate
corrections, reproduce the published decimal digits.

What this script proves / checks
--------------------------------
1. ``meta.pkl`` is structurally valid and consistent with the production plan
   (``D``, ``G``, ``Ebits``, ``K_L``, ``K_R``, ``chunk``, ``h``).
2. All 1,209 chunk files are present, parse cleanly, have the expected
   half-open ranges ``[lo, hi)``, and cover ``[-K_L, K_R + 1)`` exactly once.
3. Every stored interval satisfies ``L <= H`` and has a width compatible with
   the per-node error budget ``2^-Ebits`` and the fixed-point scale ``2^-S``.
4. The exact integer sums ``Lsum``/``Hsum`` reproduce the certificate endpoints
   ``lam_lo``/``lam_hi`` up to the non-negative corrections ``B`` and
   ``TR + TL + B``.
5. ``floor(10^D * lam_lo / 2^S) == floor(10^D * lam_hi / 2^S)`` and equals
   the published digits; the standalone digits file and the OEIS b-file agree.
6. The first 102 digits agree with OEIS A084945.

What this script does NOT prove
-------------------------------
It does not recompute the Arb node balls.  A corrupted checkpoint that is
internally consistent could still pass this audit.  The deep verifier
``verify_checkpoints_arb.py`` recomputes the node intervals with Arb and checks
containment; the two tools are meant to be used together.

The corrections ``B`` (discretisation) and ``TR``/``TL`` (tails) are only
checked for non-negativity and exact endpoint consistency here.  Their
numerical validity is checked in Arb by ``verify_checkpoints_arb.py``.

Usage
-----
    python verify_checkpoints.py
    python verify_checkpoints.py --root C:\\path\\to\\40000
    python verify_checkpoints.py --ckpt gd_rigorous_ckpt_40000 --json audit.json

Exit code 0 means every check in this file passed.
"""

from __future__ import annotations

import argparse
import json
import math
import pickle
import re
import sys
from pathlib import Path
from typing import Any, Iterable

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

EXPECTED_D = 40_000
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

CERT_NAME = "golomb_dickman_40000_certified.txt"
DIGITS_NAME = "golomb_dickman_40000_digits.txt"
BFILE_NAME = "b084945.txt"
CKPT_NAME = "gd_rigorous_ckpt_40000"
META_NAME = "meta.pkl"


class AuditError(Exception):
    """A failed structural or arithmetic check."""


def fail(message: str) -> None:
    raise AuditError(message)


def log10_int(value: int) -> float:
    """log10 of a positive integer without converting the whole integer to float."""
    if value <= 0:
        raise ValueError("log10_int expects a positive integer")
    text = str(value)
    head = int(text[:16])
    return math.log10(head) + (len(text) - 16)


def ceil_div(a: int, b: int) -> int:
    """Ceiling of a / b for b > 0."""
    return -((-a) // b)


def read_ascii(path: Path, *, allow_empty: bool = False) -> str:
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        fail(f"{path.name}: UTF-8 BOM is not allowed")
    if b"\x00" in raw:
        fail(f"{path.name}: NUL byte found")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        fail(f"{path.name}: not plain ASCII ({exc})")
    if not text and not allow_empty:
        fail(f"{path.name}: file is empty")
    return text


def check_lf_only(path: Path, text: str) -> None:
    if "\r" in text:
        fail(f"{path.name}: contains CR (expected LF-only)")
    if "\t" in text:
        fail(f"{path.name}: contains a tab character")


class SafeUnpickler(pickle.Unpickler):
    """Refuse to import anything.  The checkpoint data is primitive only."""

    def find_class(self, module: str, name: str) -> Any:  # noqa: D401
        raise pickle.UnpicklingError(
            f"refusing to load class {module}.{name} from checkpoint data"
        )


def load_pickle(path: Path) -> Any:
    try:
        with path.open("rb") as handle:
            return SafeUnpickler(handle).load()
    except AuditError:
        raise
    except Exception as exc:  # pragma: no cover - diagnostic path
        fail(f"{path.name}: cannot be unpickled ({type(exc).__name__}: {exc})")


def as_int(value: Any, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        fail(f"{what}: expected int, got {type(value).__name__}")
    return value


def as_pair(value: Any, what: str) -> tuple[int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        fail(f"{what}: expected a 2-tuple, got {type(value).__name__}")
    lo = as_int(value[0], f"{what}[0]")
    hi = as_int(value[1], f"{what}[1]")
    if not (lo <= hi):
        fail(f"{what}: lower bound {lo} exceeds upper bound {hi}")
    return lo, hi


def parse_meta(meta_path: Path) -> dict[str, Any]:
    meta = load_pickle(meta_path)
    if not isinstance(meta, (tuple, list)) or len(meta) != 8:
        fail("meta.pkl: expected an 8-tuple (D, G, h, Ebits, K_L, K_R, chunk, version)")
    D = as_int(meta[0], "meta.D")
    G = as_int(meta[1], "meta.G")
    h = meta[2]
    if not isinstance(h, float):
        fail(f"meta.h: expected float, got {type(h).__name__}")
    if not math.isfinite(h) or h <= 0.0:
        fail(f"meta.h: expected a positive finite float, got {h!r}")
    Ebits = as_int(meta[3], "meta.Ebits")
    K_L = as_int(meta[4], "meta.K_L")
    K_R = as_int(meta[5], "meta.K_R")
    chunk = as_int(meta[6], "meta.chunk")
    version = meta[7]
    if version != "v2":
        fail(f"meta.version: expected 'v2', got {version!r}")
    if D != EXPECTED_D:
        fail(f"meta.D: expected {EXPECTED_D}, got {D}")
    if G < 0:
        fail(f"meta.G: expected non-negative guard, got {G}")
    if Ebits <= 0 or K_L <= 0 or K_R <= 0 or chunk <= 0:
        fail("meta: Ebits, K_L, K_R and chunk must be positive")

    # h must be a dyadic rational; the exact integer certificate arithmetic
    # relies on this.  ``as_integer_ratio`` returns the exact double value.
    hn, hd = h.as_integer_ratio()
    if hd <= 0 or (hd & (hd - 1)) != 0:
        fail(f"meta.h: denominator {hd} is not a power of two")
    if hn <= 0:
        fail("meta.h: numerator must be positive")

    # Reproduce the production formulas for the plan parameters.  These are
    # double-precision selections; the rigorous inequalities are checked in
    # verify_checkpoints_arb.py with Arb.
    N = D + G
    ln10 = math.log(10.0)
    ln2 = math.log(2.0)
    euler = 0.5772156649015329
    Ebits_expected = int(math.ceil((N + 8) * (ln10 / ln2))) + 6
    if Ebits != Ebits_expected:
        fail(f"meta.Ebits: expected {Ebits_expected} from the plan formula, got {Ebits}")
    S_expected = Ebits + 64
    if S_expected <= 0:
        fail("meta: invalid S")

    UR = int(math.ceil((N + 8) * ln10)) + 2
    tR = solve_t(math.log(UR), 0.0, 40.0)
    K_R_expected = int(math.ceil(tR / h)) + 1
    lnuL = 0.5 * (-(N + 8) * ln10 + ln2 - euler) - 1.0
    tL = solve_t(lnuL, -600.0, 0.0)
    K_L_expected = int(math.ceil(-tL / h)) + 1
    if K_R != K_R_expected:
        fail(f"meta.K_R: expected {K_R_expected} from the tail formula, got {K_R}")
    if K_L != K_L_expected:
        fail(f"meta.K_L: expected {K_L_expected} from the tail formula, got {K_L}")

    return {
        "D": D,
        "G": G,
        "h": h,
        "hn": hn,
        "hd": hd,
        "Ebits": Ebits,
        "S": S_expected,
        "K_L": K_L,
        "K_R": K_R,
        "chunk": chunk,
        "version": version,
        "nodes": K_L + K_R + 1,
        "n_chunks": ceil_div(K_L + K_R + 1, chunk),
    }


def solve_t(target_lnu: float, lo: float, hi: float) -> float:
    """Double-precision solve of t - exp(-t) = target_lnu."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - math.exp(-mid) < target_lnu:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def expected_ranges(K_L: int, K_R: int, chunk: int) -> list[tuple[int, int]]:
    n0 = -K_L
    n1 = K_R + 1
    return [(lo, min(lo + chunk, n1)) for lo in range(n0, n1, chunk)]


def scan_chunks(ckpt: Path, meta: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    expected = expected_ranges(meta["K_L"], meta["K_R"], meta["chunk"])
    expected_by_lo = {lo: hi for lo, hi in expected}
    errors: list[str] = []

    files: dict[int, Path] = {}
    for path in sorted(ckpt.glob("c_*.pkl")):
        suffix = path.stem[2:]
        try:
            key = int(suffix)
        except ValueError:
            errors.append(f"unexpected chunk file name: {path.name}")
            continue
        if key in files:
            errors.append(f"duplicate chunk key {key}: {files[key].name} and {path.name}")
        files[key] = path

    missing = sorted(set(expected_by_lo) - set(files))
    extra = sorted(set(files) - set(expected_by_lo))
    if missing:
        errors.append(f"missing {len(missing)} chunk(s), first: {missing[:5]}")
    if extra:
        errors.append(f"unexpected {len(extra)} chunk(s), first: {extra[:5]}")

    chunks: list[dict[str, Any]] = []
    previous_hi: int | None = None
    for lo, hi in expected:
        path = files.get(lo)
        if path is None:
            continue
        value = load_pickle(path)
        if not isinstance(value, (tuple, list)) or len(value) != 4:
            errors.append(f"{path.name}: expected a 4-tuple (lo, hi, even, odd)")
            continue
        try:
            got_lo = as_int(value[0], f"{path.name}.lo")
            got_hi = as_int(value[1], f"{path.name}.hi")
            even = as_pair(value[2], f"{path.name}.even")
            odd = as_pair(value[3], f"{path.name}.odd")
        except AuditError as exc:
            errors.append(str(exc))
            continue
        if got_lo != lo or got_hi != hi:
            errors.append(
                f"{path.name}: range ({got_lo}, {got_hi}) does not match expected ({lo}, {hi})"
            )
            continue
        if previous_hi is not None and got_lo != previous_hi:
            errors.append(
                f"{path.name}: gap/overlap at {previous_hi} -> {got_lo}"
            )
        previous_hi = got_hi
        chunks.append(
            {
                "lo": got_lo,
                "hi": got_hi,
                "even": even,
                "odd": odd,
                "path": path,
            }
        )

    if chunks:
        if chunks[0]["lo"] != -meta["K_L"]:
            errors.append(
                f"first chunk starts at {chunks[0]['lo']}, expected {-meta['K_L']}"
            )
        if chunks[-1]["hi"] != meta["K_R"] + 1:
            errors.append(
                f"last chunk ends at {chunks[-1]['hi']}, expected {meta['K_R'] + 1}"
            )
    else:
        errors.append("no chunks could be read")

    return chunks, errors


def check_chunk_widths(chunks: list[dict[str, Any]], meta: dict[str, Any]) -> list[str]:
    """Sanity-check the saved interval widths against the per-node error budget.

    A saved chunk interval is the outward-rounded sum of ``hi - lo`` node balls,
    each of radius < 2^-Ebits.  In units of 2^-S with S = Ebits + 64, each node
    contributes at most 2 * 2^64 = 2^65 units of width.  We allow a factor two
    for the outward rounding and the two parities.
    """
    errors: list[str] = []
    max_width = 0
    for chunk in chunks:
        n = chunk["hi"] - chunk["lo"]
        width = (chunk["even"][1] - chunk["even"][0]) + (chunk["odd"][1] - chunk["odd"][0])
        if width < 0:
            errors.append(f"chunk {chunk['lo']}: negative aggregate width")
            continue
        # Generous but still meaningful: 2^66 units per node plus 16 for rounding.
        limit = n * (1 << 66) + 16
        if width > limit:
            errors.append(
                f"chunk {chunk['lo']}: width {width} exceeds sanity limit {limit}"
            )
        max_width = max(max_width, width)
    return errors


def parse_certificate(path: Path, meta: dict[str, Any]) -> dict[str, Any]:
    text = read_ascii(path)
    check_lf_only(path, text)
    lines = text.splitlines()
    if len(lines) < 4:
        fail(f"{path.name}: expected at least 4 lines")
    if "Golomb-Dickman constant" not in lines[0]:
        fail(f"{path.name}: unexpected header line")
    match_s = re.search(r"\bS=(\d+)\b", lines[1])
    match_lo = re.search(r"\blam_lo=(\d+)\b", lines[1])
    match_hi = re.search(r"\blam_hi=(\d+)\b", lines[2])
    if not (match_s and match_lo and match_hi):
        fail(f"{path.name}: cannot parse S / lam_lo / lam_hi")
    S = int(match_s.group(1))
    lam_lo = int(match_lo.group(1))
    lam_hi = int(match_hi.group(1))
    if S != meta["S"]:
        fail(f"{path.name}: S={S}, expected {meta['S']} from meta.pkl")
    if not (0 < lam_lo < lam_hi):
        fail(f"{path.name}: expected 0 < lam_lo < lam_hi")
    digit_line = lines[3].strip()
    if not digit_line.startswith("0."):
        fail(f"{path.name}: fourth line must start with '0.'")
    digits = digit_line[2:]
    if len(digits) != EXPECTED_D:
        fail(f"{path.name}: expected {EXPECTED_D} digits, got {len(digits)}")
    if not digits.isdigit():
        fail(f"{path.name}: digit line contains non-digits")
    return {
        "S": S,
        "lam_lo": lam_lo,
        "lam_hi": lam_hi,
        "digits": digits,
    }


def check_digits_file(path: Path, cert_digits: str) -> None:
    text = read_ascii(path)
    check_lf_only(path, text)
    value = text.strip()
    if len(value) != EXPECTED_D:
        fail(f"{path.name}: expected {EXPECTED_D} digits, got {len(value)}")
    if not value.isdigit():
        fail(f"{path.name}: contains non-digits")
    if value != cert_digits:
        fail(f"{path.name}: digits differ from the certificate")


def check_bfile(path: Path, cert_digits: str) -> int:
    text = read_ascii(path)
    check_lf_only(path, text)
    if not text.endswith("\n\n"):
        fail(f"{path.name}: expected a final blank line (LF LF)")
    lines = text.splitlines()
    terms: list[tuple[int, int]] = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) != 2:
            fail(f"{path.name}: line {line_number}: expected 'n value'")
        try:
            index, value = int(parts[0]), int(parts[1])
        except ValueError:
            fail(f"{path.name}: line {line_number}: non-integer value")
        terms.append((index, value))
    if len(terms) != EXPECTED_D:
        fail(f"{path.name}: expected {EXPECTED_D} terms, got {len(terms)}")
    for position, (index, value) in enumerate(terms):
        if index != position:
            fail(f"{path.name}: index {index} at position {position}, expected {position}")
        if not 0 <= value <= 9:
            fail(f"{path.name}: a({index})={value} is not a decimal digit")
        if value != int(cert_digits[position]):
            fail(f"{path.name}: a({index}) differs from the certificate")
    return len(terms)


def audit(root: Path, ckpt_name: str, cert_name: str, digits_name: str, bfile_name: str) -> dict[str, Any]:
    ckpt = root / ckpt_name
    if not ckpt.is_dir():
        fail(f"missing checkpoint directory: {ckpt}")
    meta_path = ckpt / META_NAME
    if not meta_path.is_file():
        fail(f"missing meta file: {meta_path}")

    meta = parse_meta(meta_path)
    chunks, errors = scan_chunks(ckpt, meta)
    errors.extend(check_chunk_widths(chunks, meta))

    cert = parse_certificate(root / cert_name, meta)
    check_digits_file(root / digits_name, cert["digits"])
    n_terms = check_bfile(root / bfile_name, cert["digits"])

    Lsum = sum(c["even"][0] + c["odd"][0] for c in chunks)
    Hsum = sum(c["even"][1] + c["odd"][1] for c in chunks)
    hn, hd = meta["hn"], meta["hd"]
    lo_fin = (hn * Lsum) // hd
    hi_fin = -((-hn * Hsum) // hd)
    B = lo_fin - cert["lam_lo"]
    C = cert["lam_hi"] - hi_fin
    if B < 0:
        errors.append(f"certificate lam_lo is above floor(h*Lsum) by {-B}")
    if C < 0:
        errors.append(f"certificate lam_hi is below ceil(h*Hsum) by {-C}")
    if C < B:
        errors.append(f"upper correction C={C} is smaller than lower correction B={B}")

    S = meta["S"]
    scale = 10 ** EXPECTED_D
    a = (cert["lam_lo"] * scale) >> S
    b = (cert["lam_hi"] * scale) >> S
    if a != b:
        errors.append("certified interval crosses a decimal cell")
    proven = str(a).zfill(EXPECTED_D)
    if proven != cert["digits"]:
        errors.append("certificate digits differ from floor(10^D * lam_lo / 2^S)")
    if cert["digits"][: len(KNOWN_PREFIX)] != KNOWN_PREFIX:
        errors.append("first 102 digits do not match OEIS A084945")

    width = cert["lam_hi"] - cert["lam_lo"]
    chunk_width = hi_fin - lo_fin
    width_log10 = log10_int(width) - S * math.log10(2.0)
    chunk_width_log10 = (
        log10_int(chunk_width) - S * math.log10(2.0) if chunk_width > 0 else float("-inf")
    )
    correction_log10 = (
        log10_int(C) - S * math.log10(2.0) if C > 0 else float("-inf")
    )
    guard_digits = -width_log10

    if errors:
        raise AuditError("\n".join(errors))

    return {
        "root": str(root),
        "meta": {
            "D": meta["D"],
            "G": meta["G"],
            "h": meta["h"],
            "h_ratio": [meta["hn"], meta["hd"]],
            "Ebits": meta["Ebits"],
            "S": meta["S"],
            "K_L": meta["K_L"],
            "K_R": meta["K_R"],
            "chunk": meta["chunk"],
            "nodes": meta["nodes"],
            "n_chunks": meta["n_chunks"],
        },
        "chunks_read": len(chunks),
        "bfile_terms": n_terms,
        "Lsum_digits": len(str(Lsum)),
        "Hsum_digits": len(str(Hsum)),
        "B": B,
        "C_upper_correction": C,
        "chunk_width_log10": chunk_width_log10,
        "certificate_width_log10": width_log10,
        "correction_log10": correction_log10,
        "guard_digits": guard_digits,
        "digits": cert["digits"],
    }


def print_report(report: dict[str, Any]) -> None:
    meta = report["meta"]
    print("OK: checkpoint audit passed")
    print(f"  root: {report['root']}")
    print(
        "  meta: D={D} G={G} Ebits={Ebits} S={S} K_L={K_L} K_R={K_R} "
        "chunk={chunk}".format(**meta)
    )
    print(f"  h = {meta['h_ratio'][0]} / {meta['h_ratio'][1]}")
    print(f"  chunks: {report['chunks_read']} (expected {meta['n_chunks']})")
    print(f"  nodes covered: {meta['nodes']}")
    print(f"  Lsum: {report['Lsum_digits']} digits; Hsum: {report['Hsum_digits']} digits")
    print(f"  lower correction B = {report['B']}")
    print(f"  upper correction C = TR + TL + B = {report['C_upper_correction']}")
    print(f"  chunk-sum width: 10^{report['chunk_width_log10']:.6f}")
    print(f"  certificate width: 10^{report['certificate_width_log10']:.6f}")
    print(f"  guard digits: {report['guard_digits']:.2f}")
    print(f"  b-file terms: {report['bfile_terms']}")
    print("  first 102 digits match OEIS A084945")
    print("")
    print("NOTE: this audit does not recompute the Arb node balls.")
    print("      Run verify_checkpoints_arb.py for the deep Arb verification.")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=_bundle_root(Path(__file__).resolve().parent),
        help="directory containing the certificate, digits, b-file and checkpoint folder",
    )
    parser.add_argument("--ckpt", default=CKPT_NAME, help="checkpoint directory name")
    parser.add_argument("--cert", default=CERT_NAME, help="certificate file name")
    parser.add_argument("--digits", default=DIGITS_NAME, help="digits file name")
    parser.add_argument("--bfile", default=BFILE_NAME, help="OEIS b-file name")
    parser.add_argument("--json", type=Path, default=None, help="write a JSON report")
    args = parser.parse_args(argv)

    try:
        report = audit(args.root, args.ckpt, args.cert, args.digits, args.bfile)
    except AuditError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    print_report(report)
    if args.json is not None:
        payload = dict(report)
        payload.pop("digits", None)
        args.json.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"  JSON report: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
