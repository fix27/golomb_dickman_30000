#!/usr/bin/env python3
"""Deep, independent Arb verification of the 40,000-digit checkpoint set.

This is the companion of ``verify_checkpoints.py``.  The stdlib auditor checks
the pickle structure, coverage, exact integer sums and published artifacts.
This script additionally recomputes the Arb node balls from the defining
formulas, without importing ``golomb_dickman.py`` and without using the saved
interval endpoints as proof.

What is recomputed
------------------
* the constants M, B_disc, tau_R, tau_L and the plan indices K_L, K_R;
* the node values G(kh) = F(u) u (1 + e^-t),  u = exp(t - e^-t),
  F(u) = exp(-u - E1(u)), using three independent formulas:
    - the alternating series for E1 (small u),
    - the Stieltjes continued fraction for e^u E1(u) (moderate u),
    - the asymptotic series for E1 (large u);
* the outward-rounded fixed-point sums for every selected chunk;
* the final interval and its decimal digits.

Checks
------
* every recomputed chunk interval is contained in the saved checkpoint
  interval (per parity);
* the recomputed integer corrections B, TR, TL reproduce the corrections
  implied by the certificate;
* the final interval proves exactly the published 40,000 digits.

Modes
-----
    python verify_checkpoints_arb.py --audit
        Constants, certificate corrections and exact integer sums only.
        No node recomputation.  Fast.

    python verify_checkpoints_arb.py --sample 20
        Recompute 20 chunks spread over the whole index range (including the
        first and last chunk).  This is the recommended first deep check.

    python verify_checkpoints_arb.py --all --jobs 8
        Recompute all 1,209 chunks.  This is the full independent verification
        and can take hours.  Results are cached in ``verify_arb_ckpt_40000/``;
        rerun with the same command to resume.

    python verify_checkpoints_arb.py --nodes 20
        Self-test: recompute single nodes by two different formulas and check
        that the resulting intervals overlap.

Requirements
------------
``python-flint`` (Arb).  The script deliberately does not install it and does
not fall back to floating point.  If it is missing, the script exits with code
2.

Trust base
----------
Arb's ball arithmetic, the Stieltjes continued-fraction theorem, the standard
alternating/asymptotic remainder bounds for E1, and the elementary estimates
implemented in ``strip_constants``.  This is not a proof-assistant
formalisation.

Exit codes: 0 = all checks passed, 1 = a check failed, 2 = missing dependency
or usage error.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

try:
    from flint import arb, ctx
except ImportError:  # pragma: no cover - handled in main()
    arb = None  # type: ignore[assignment]
    ctx = None  # type: ignore[assignment]

# Import the stdlib auditor from the same directory.  It is used for parsing
# and structural checks; all numerical node work below is independent.
try:
    import verify_checkpoints as vc
except ImportError:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import verify_checkpoints as vc  # type: ignore[no-redef]

LN2 = math.log(2.0)
LN10 = math.log(10.0)
LOG2_10 = LN10 / LN2
EULER_F = 0.5772156649015329  # only for float estimates, never for proof

ASYM_KMAX = 3000
RS_M = 8
def _bundle_root(start: Path) -> Path:
    """Locate the bundle root: the directory that holds the certificate.

    This lets every script be run from anywhere, including from ``checks/``.
    """
    for candidate in (start, *start.parents):
        if (candidate / "golomb_dickman_40000_certified.txt").is_file():
            return candidate
    return start

U_SERIES_DEFAULT = 2000.0
EXTRA_BITS_DEFAULT = 32


class Escalate(Exception):
    """Raised when a formula needs more precision or a different method."""


_EG_CACHE: dict[str, Any] = {"bits": 0, "val": None}


def require_flint() -> None:
    if arb is None:
        print(
            "FAIL: python-flint is not installed.\n"
            "Install it with:  python -m pip install python-flint",
            file=sys.stderr,
        )
        raise SystemExit(2)


def get_eg(bits: int):
    """e^gamma with at least ``bits`` bits of precision (cached per process)."""
    if _EG_CACHE["val"] is None or _EG_CACHE["bits"] < bits:
        old = ctx.prec
        work = int(bits * 1.2) + 256
        ctx.prec = work
        _EG_CACHE["val"] = arb.const_euler().exp()
        _EG_CACHE["bits"] = work
        ctx.prec = old
    return _EG_CACHE["val"]


def log10_arb(x) -> int:
    """Floor of log10(x) for a positive Arb ball, computed rigorously."""
    old = ctx.prec
    ctx.prec = max(old, 512)
    try:
        value = (x.upper().log() / arb(10).log()).floor().lower().unique_fmpz()
        return int(value)
    finally:
        ctx.prec = old


def arb_ceil_int(x) -> int:
    """Smallest integer >= x for a ball x."""
    return int(x.upper().floor().upper().unique_fmpz()) + 1


def to_fixed(x, S: int) -> tuple[int, int]:
    """Outward-rounded fixed-point interval (L, H) with x in [L,H] * 2^-S."""
    old = ctx.prec
    ctx.prec = max(old, S + 128)
    try:
        scale = arb(1 << S)
        lower = int((x.lower() * scale).floor().lower().unique_fmpz())
        upper = int((x.upper() * scale).floor().upper().unique_fmpz()) + 1
        return lower, upper
    finally:
        ctx.prec = old


def strip_constants(prec: int = 768) -> dict[str, Any]:
    """Lemma 4: M >= sup_{|y| <= d} integral |G(x+iy)| dx, d = 7/5.

    Independent re-implementation of the three explicit regions.
    """
    old = ctx.prec
    ctx.prec = prec
    try:
        one = arb(1)
        pi = arb.pi()
        gam = arb.const_euler()
        e = one.exp()
        d = arb(7) / 5
        c = d.cos()
        if not (c > 0):
            raise ArithmeticError("strip_constants: cos(7/5) is not positive")
        a = 2 * c

        # t <= 0: |u| <= 1, |F| <= C1 |u|
        C1 = (gam + 2 + (e - 2) / 2).exp()
        M1 = C1 * (-a).exp() * (2 / a + 1 / (a * a))

        # 0 <= t <= 2: |u| <= e^2
        e2 = arb(2).exp()
        em2 = (-arb(2)).exp()
        B2 = 2 * gam.exp() * (e2 * e2) * (2 * e2 + e2.exp() / 2).exp()
        M2 = 2 * B2

        # t >= 2: |arg u| <= theta_max < pi/2, Re u >= kappa |u|
        th = d + em2 * d.sin()
        if not (th < pi / 2):
            raise ArithmeticError("strip_constants: theta_max is not below pi/2")
        kap = th.cos()
        if not (kap > 0):
            raise ArithmeticError("strip_constants: kappa is not positive")
        kp = kap * (-em2).exp()
        facE = ((em2 - 2).exp()).exp()
        M3 = facE * (1 + em2) * (-kp * e2).exp() / kp

        return {
            "d": d,
            "M": M1 + M2 + M3,
            "M1": M1,
            "M2": M2,
            "M3": M3,
            "pi": pi,
            "gam": gam,
        }
    finally:
        ctx.prec = old


def solve_t(target_lnu: float, lo: float, hi: float) -> float:
    """Double-precision solve of t - exp(-t) = target_lnu (monotone)."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - math.exp(-mid) < target_lnu:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def estimate(k: int, h: float) -> tuple[float, float, float]:
    """Float estimates (xf, uf, lng) used only to choose a method."""
    tf = k * h
    emf = math.exp(-tf) if tf > -745.0 else math.inf
    if math.isinf(emf):
        xf = -math.inf
        uf = 0.0
        lng = -math.inf
        return xf, uf, lng
    xf = tf - emf
    if xf < -745.0:
        uf = 0.0
    elif xf > 700.0:
        uf = math.inf
    else:
        uf = math.exp(xf)
    l1p = math.log1p(emf) if emf < 1e300 else math.log(emf) + math.log1p(1.0 / emf)
    if xf < 0:
        lng = EULER_F + 2 * xf + l1p
    else:
        lng = xf - uf + l1p
    return xf, uf, lng


# ---------------------------------------------------------------------------
# Node evaluation
# ---------------------------------------------------------------------------

def _make_u(k: int, h: float):
    """Return (t, em, u) as Arb balls: t = kh, em = e^-t, u = e^(t-em)."""
    t = arb(h) * k
    em = (-t).exp()
    u = (t - em).exp()
    return t, em, u


def eval_series(k: int, h: float, thr_bits: int, P: int):
    """E1 series for small u, with a rigorous alternating-tail bound.

    E1(u) = -gamma - log(u) - sum_{j>=1} (-u)^j/(j j!), so
    G = e^gamma u^2 (1+e^-t) exp(-2u + S + delta), |delta| <= rem.
    The finite sum is evaluated by rectangular block splitting; the remainder
    is 2|a_{K+1}| once K+1 >= 2u.
    """
    old = ctx.prec
    ctx.prec = P
    try:
        _, em, u = _make_u(k, h)
        m = RS_M
        mu = -u
        xp = [arb(1), mu]
        for _j in range(2, m + 1):
            xp.append(xp[-1] * mu)
        xm = xp[m]

        S = arb(0)
        pref = xp[2] / 4  # a_2 = (-u)^2/(2*2!)
        kb = 2
        thr = arb(1) / arb(1 << thr_bits)

        blocks = 0
        while True:
            blocks += 1
            if blocks > 20_000_000:
                raise Escalate("series: too many blocks")
            inner = arb(1)
            dprod = 1
            for j in range(1, m):
                dprod *= kb + j
                inner = inner + (xp[j] * kb) / ((kb + j) * dprod)
            S = S + pref * inner
            dprod *= kb + m
            pref = (pref * xm) * kb / ((kb + m) * dprod)
            kb += m
            if kb >= 2 * u.upper() and 2 * abs(pref).upper() < thr:
                break

        rem = 2 * abs(pref).upper()
        expo = -2 * u + S + arb(0, rem)
        return get_eg(P) * u * u * (1 + em) * expo.exp()
    finally:
        ctx.prec = old


def eval_cf(k: int, h: float, gap_bits: int, P: int, xf: float, uf: float, lng: float, extra: int):
    """Stieltjes continued fraction for f(u) = e^u E1(u).

    f lies between consecutive convergents, so |f - f_n| <= |f_n - f_{n-1}|.
    """
    old = ctx.prec
    try:
        ctx.prec = P
        _, em, u = _make_u(k, h)

        X = max(extra * LN2 + lng - uf + xf, 1.0)
        nmax = int(3.0 * X * (0.22 + 0.14 * X / max(uf, 1.0))) + 4000
        nmax = min(max(nmax, 4000), 2_000_000)

        ctx.prec = gap_bits + 96 + int(math.log2(nmax + 2))
        A2, A1 = arb(1), u  # A_{-1}, A_0
        B2, B1 = arb(0), arb(1)  # B_{-1}, B_0
        n = 0
        converged = False
        thr = arb(1) / arb(1 << gap_bits)
        while n < nmax:
            for _ in range(32):
                n += 1
                a = (n + 1) // 2
                if n & 1:
                    A2, A1 = A1, A1 + a * A2
                    B2, B1 = B1, B1 + a * B2
                else:
                    A2, A1 = A1, u * A1 + a * A2
                    B2, B1 = B1, u * B1 + a * B2
            fn = B1 / A1
            D = abs(fn - B2 / A2).upper()
            if D < thr:
                converged = True
                break
        if not converged:
            raise Escalate("continued fraction: nmax exceeded")

        f = fn + arb(0, D)
        ctx.prec = P
        E1 = (-u).exp() * f
        return (-u - E1).exp() * u * (1 + em)
    finally:
        ctx.prec = old


def choose_asym_K(xf: float, uf: float, lng: float, L_eps: float) -> int | None:
    """Smallest K with e^{-2u} K!/u^K <= e^{-L_eps-3} (float selection only)."""
    if not math.isfinite(uf) or uf <= 1.0:
        return None
    top = min(ASYM_KMAX, int(uf) - 1)
    if top < 0:
        return None
    lnfact = 0.0
    for K in range(0, top + 1):
        if K > 0:
            lnfact += math.log(K)
        if lng + lnfact - (K + 1) * xf - uf <= -L_eps - 3.0:
            return K
    return None


def eval_asym(k: int, h: float, K: int, P: int):
    """Asymptotic series for E1(u) = e^-u (S/u + R), |R| <= K!/u^(K+1)."""
    old = ctx.prec
    ctx.prec = P
    try:
        _, em, u = _make_u(k, h)
        w = (-u).exp()
        if K == 0:
            Ssum = arb(0)
        else:
            Ssum = arb(1)
            for j in range(K - 1, 0, -1):
                Ssum = 1 - Ssum * j / u
        f = arb(1)
        for j in range(1, K + 1):
            f = f * j / u
        rem = f / u
        E1 = w * (Ssum / u + arb(0, rem))
        return (-u - E1).exp() * u * (1 + em)
    finally:
        ctx.prec = old


def method_ball(k: int, params: dict[str, Any], method: str, extra: int):
    """Evaluate one node with a forced method and a tighter-by-extra target."""
    h = params["h"]
    Ebits = params["Ebits"]
    S = params["S"]
    xf, uf, lng = estimate(k, h)
    L_eps = (Ebits + extra) * LN2

    if method == "series":
        need_ln = L_eps + lng + 2.0
        thr_bits = max(64, int(need_ln / LN2) + 8)
        P = max(256, thr_bits + int(max(0.0, math.log2(uf + 2.0))) + 128)
        return eval_series(k, h, thr_bits, P)

    if method == "cf":
        gap_bits = max(32, int((L_eps + lng - uf + 4.0) / LN2) + 2)
        P = max(256, Ebits + extra + 128)
        return eval_cf(k, h, gap_bits, P, xf, uf, lng, extra)

    if method == "asym":
        K = choose_asym_K(xf, uf, lng, L_eps)
        if K is None:
            raise Escalate("asymptotic: no admissible K")
        P = max(256, Ebits + extra + 128)
        return eval_asym(k, h, K, P)

    raise ValueError(f"unknown method: {method}")


def node_ball(k: int, params: dict[str, Any]):
    """Recompute a rigorous Arb ball for G(kh) with radius below 2^-(Ebits+extra)."""
    h = params["h"]
    Ebits = params["Ebits"]
    u_series = params.get("u_series", U_SERIES_DEFAULT)
    extra = params.get("extra", EXTRA_BITS_DEFAULT)
    xf, uf, lng = estimate(k, h)

    for attempt in range(12):
        L_eps = (Ebits + extra) * LN2
        if uf < u_series:
            method = "series"
        else:
            K = choose_asym_K(xf, uf, lng, L_eps)
            method = "asym" if K is not None else "cf"
        if attempt >= 8:
            method = "series"
        try:
            g = method_ball(k, params, method, extra)
            target = arb(1) / arb(1 << (Ebits + extra))
            if g.is_finite() and g.rad() < target:
                return g
        except Escalate:
            pass
        extra += 64

    raise ArithmeticError(
        f"node k={k}: could not reach radius target (last method {method})"
    )


# ---------------------------------------------------------------------------
# Plan, chunks and caching
# ---------------------------------------------------------------------------

def recompute_plan(meta: dict[str, Any]) -> dict[str, Any]:
    """Recompute M, B_disc, tails and the integer corrections with Arb."""
    D = meta["D"]
    G = meta["G"]
    h = meta["h"]
    Ebits = meta["Ebits"]
    K_L = meta["K_L"]
    K_R = meta["K_R"]
    S = meta["S"]
    N = D + G

    K = strip_constants(768)
    M, d, pi, gam = K["M"], K["d"], K["pi"], K["gam"]

    old = ctx.prec
    ctx.prec = 512
    try:
        ten = arb(10)
        hA = arb(h)
        Bdisc = 2 * M / ((2 * pi * d / hA).exp() - 1)

        tR = K_R * hA
        uR = (tR - (-tR).exp()).exp()
        tailR = (-uR).exp()

        tL = -(K_L * hA)
        uL = (tL - (-tL).exp()).exp()
        tailL = gam.exp() * uL * uL / 2

        # h selection formula (float); h itself is a valid explicit parameter
        # even if this reproduction differs in the last bit.
        Lint = arb_ceil_int((2 * M).log() + (N + 6) * ten.log())
        h_calc = math.floor(0.97 * 2 * math.pi * 1.4 / Lint * 2**40) / 2**40

        checks = {
            "Bdisc_small": bool(Bdisc * ten ** (N + 6) < 1),
            "tailR_small": bool(tailR * ten ** (N + 8) < 1),
            "tailL_small": bool(tailL * ten ** (N + 8) < 1),
            "uR_ge_1": bool(uR >= 1),
            "uR_large_enough": bool(uR >= math.ceil((N + 8) * LN10) + 2),
            "uL_positive": bool(uL > 0),
        }

        B = to_fixed(Bdisc, S)[1]
        TR = to_fixed(tailR, S)[1]
        TL = to_fixed(tailL, S)[1]
    finally:
        ctx.prec = old

    return {
        "M": M,
        "Bdisc": Bdisc,
        "tailR": tailR,
        "tailL": tailL,
        "B": B,
        "TR": TR,
        "TL": TL,
        "h_calc": h_calc,
        "h_reproduced": (h_calc == h),
        "checks": checks,
        "log10_M": log10_arb(M),
        "log10_Bdisc": log10_arb(Bdisc),
        "log10_tailR": log10_arb(tailR),
        "log10_tailL": log10_arb(tailL),
    }


def recompute_chunk(task: tuple[int, int, dict[str, Any]]) -> tuple[int, int, int, int, int, int]:
    """Worker: recompute one half-open chunk and return even/odd fixed intervals."""
    lo, hi, params = task
    S = params["S"]
    Ebits = params["Ebits"]
    extra = params["extra"]
    ctx.prec = max(256, Ebits + extra + 256)
    even = arb(0)
    odd = arb(0)
    for k in range(lo, hi):
        g = node_ball(k, params)
        ctx.prec = max(256, Ebits + extra + 256)
        if k & 1:
            odd = odd + g
        else:
            even = even + g
    Le, He = to_fixed(even, S)
    Lo, Ho = to_fixed(odd, S)
    return lo, hi, Le, He, Lo, Ho


def cache_signature(meta: dict[str, Any], extra: int, u_series: float) -> str:
    return (
        f"v1|{meta['h']!r}|{meta['Ebits']}|{meta['S']}|"
        f"{extra}|{u_series!r}"
    )


def load_cached(cache_dir: Path, lo: int, signature: str):
    path = cache_dir / f"c_{lo}.pkl"
    if not path.is_file():
        return None
    try:
        with path.open("rb") as handle:
            payload = pickle.load(handle)
        if isinstance(payload, tuple) and len(payload) == 2 and payload[0] == signature:
            result = payload[1]
            if isinstance(result, tuple) and len(result) == 6 and result[0] == lo:
                return result
    except Exception:
        return None
    return None


def save_cached(cache_dir: Path, signature: str, result: tuple[int, ...]) -> None:
    path = cache_dir / f"c_{result[0]}.pkl"
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as handle:
        pickle.dump((signature, result), handle, protocol=4)
    os.replace(tmp, path)


def sample_chunks(chunks: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    """Evenly spaced chunks, always including the first and the last."""
    if count <= 0:
        return []
    if count >= len(chunks):
        return list(chunks)
    if count == 1:
        return [chunks[0]]
    indices = sorted({round(i * (len(chunks) - 1) / (count - 1)) for i in range(count)})
    return [chunks[i] for i in indices]


def compare_chunk(saved: dict[str, Any], result: tuple[int, ...]) -> dict[str, Any]:
    """Check that the recomputed interval is inside the saved interval."""
    lo, hi, Le, He, Lo, Ho = result
    se = saved["even"]
    so = saved["odd"]
    even_ok = se[0] <= Le and He <= se[1]
    odd_ok = so[0] <= Lo and Ho <= so[1]
    # Positive values mean the recomputed interval sticks out of the saved one.
    left_even = se[0] - Le
    right_even = He - se[1]
    left_odd = so[0] - Lo
    right_odd = Ho - so[1]
    return {
        "lo": lo,
        "hi": hi,
        "ok": bool(even_ok and odd_ok),
        "even_ok": bool(even_ok),
        "odd_ok": bool(odd_ok),
        "even_saved": se,
        "even_new": (Le, He),
        "odd_saved": so,
        "odd_new": (Lo, Ho),
        "left_even": left_even,
        "right_even": right_even,
        "left_odd": left_odd,
        "right_odd": right_odd,
    }


def run_chunks(
    selected: list[dict[str, Any]],
    params: dict[str, Any],
    signature: str,
    cache_dir: Path | None,
    jobs: int,
    use_cache: bool,
) -> dict[int, tuple[int, ...]]:
    """Recompute selected chunks in parallel, with optional on-disk cache."""
    results: dict[int, tuple[int, ...]] = {}
    todo: list[dict[str, Any]] = []

    for chunk in selected:
        lo = chunk["lo"]
        cached = None
        if use_cache and cache_dir is not None:
            cached = load_cached(cache_dir, lo, signature)
        if cached is not None:
            results[lo] = cached
        else:
            todo.append(chunk)

    if not todo:
        print(f"  all {len(selected)} selected chunks loaded from cache")
        return results

    total = len(todo)
    done = 0
    started = time.time()
    print(f"  recomputing {total} chunk(s) with {jobs} worker(s) ...", flush=True)

    if jobs <= 1:
        for chunk in todo:
            task = (chunk["lo"], chunk["hi"], params)
            result = recompute_chunk(task)
            results[result[0]] = result
            if use_cache and cache_dir is not None:
                save_cached(cache_dir, signature, result)
            done += 1
            if done % 5 == 0 or done == total:
                elapsed = time.time() - started
                print(f"    [{done}/{total}] {elapsed:.1f}s", flush=True)
        return results

    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(recompute_chunk, (chunk["lo"], chunk["hi"], params)): chunk
            for chunk in todo
        }
        for future in as_completed(futures):
            chunk = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                raise ArithmeticError(
                    f"chunk {chunk['lo']}: recomputation failed ({type(exc).__name__}: {exc})"
                ) from exc
            results[result[0]] = result
            if use_cache and cache_dir is not None:
                save_cached(cache_dir, signature, result)
            done += 1
            if done % 5 == 0 or done == total:
                elapsed = time.time() - started
                rate = done / elapsed if elapsed > 0 else 0.0
                print(
                    f"    [{done}/{total}] {elapsed:.1f}s ({rate:.2f} chunks/s)",
                    flush=True,
                )
    return results


# ---------------------------------------------------------------------------
# Self-test and final certificate checks
# ---------------------------------------------------------------------------

def build_params(meta: dict[str, Any], extra: int, u_series: float) -> dict[str, Any]:
    return {
        "h": meta["h"],
        "Ebits": meta["Ebits"],
        "S": meta["S"],
        "extra": extra,
        "u_series": u_series,
    }


def overlaps(g1, g2) -> bool:
    """True if the two Arb balls intersect (both are rigorous enclosures)."""
    return not (g1.upper() < g2.lower() or g2.upper() < g1.lower())


def node_selftest(meta: dict[str, Any], extra: int, u_series: float, count: int) -> tuple[int, list[str]]:
    """Recompute a few nodes by two formulas and check interval overlap."""
    params = build_params(meta, extra, u_series)
    k0 = -meta["K_L"]
    k1 = meta["K_R"]
    if count <= 0:
        return 0, []
    if count == 1:
        ks = [0]
    else:
        ks = sorted({round(k0 + i * (k1 - k0) / (count - 1)) for i in range(count)})

    errors: list[str] = []
    tested = 0
    print(f"  node self-test on {len(ks)} node(s):")
    for k in ks:
        xf, uf, lng = estimate(k, params["h"])
        if uf < u_series:
            m1, m2 = "series", "cf"
        else:
            K = choose_asym_K(xf, uf, lng, (meta["Ebits"] + extra) * LN2)
            m1 = "cf"
            m2 = "asym" if K is not None else "cf"
        try:
            g1 = method_ball(k, params, m1, extra)
            if m1 == m2:
                # Precision self-check: the tighter ball must lie inside the wider one.
                g2 = method_ball(k, params, m2, extra + 64)
                ok = g2.lower() >= g1.lower() and g2.upper() <= g1.upper()
                relation = "tighter-inside"
            else:
                g2 = method_ball(k, params, m2, extra)
                ok = overlaps(g1, g2)
                relation = "overlap"
            tested += 1
            if not ok:
                errors.append(
                    f"node k={k}: {m1} vs {m2} {relation} check failed"
                )
                print(f"    k={k:7d}  u~{uf:.6g}  {m1} vs {m2}: FAIL")
            else:
                print(f"    k={k:7d}  u~{uf:.6g}  {m1} vs {m2}: OK ({relation})")
        except Escalate as exc:
            errors.append(f"node k={k}: {m1}/{m2} escalated: {exc}")
            print(f"    k={k:7d}  u~{uf:.6g}  {m1} vs {m2}: ESCALATE ({exc})")
    print(f"  node self-test: {tested}/{len(ks)} pairs completed")
    return tested, errors


def final_interval(meta: dict[str, Any], Lsum: int, Hsum: int, B: int, TR: int, TL: int):
    hn, hd = meta["hn"], meta["hd"]
    lo_fin = (hn * Lsum) // hd
    hi_fin = -((-hn * Hsum) // hd)
    return {
        "lo_fin": lo_fin,
        "hi_fin": hi_fin,
        "lam_lo": lo_fin - B,
        "lam_hi": hi_fin + TR + TL + B,
    }


def check_final(meta: dict[str, Any], cert: dict[str, Any], Lsum: int, Hsum: int, plan: dict[str, Any]) -> dict[str, Any]:
    interval = final_interval(meta, Lsum, Hsum, plan["B"], plan["TR"], plan["TL"])
    lam_lo = interval["lam_lo"]
    lam_hi = interval["lam_hi"]
    S = meta["S"]
    D = meta["D"]
    scale = 10 ** D
    a = (lam_lo * scale) >> S
    b = (lam_hi * scale) >> S
    digits_proven = str(a).zfill(D) if a == b else None
    digits_ok = bool(a == b and digits_proven == cert["digits"])
    overlap = bool(max(lam_lo, cert["lam_lo"]) <= min(lam_hi, cert["lam_hi"]))
    subset_saved = bool(cert["lam_lo"] <= lam_lo and lam_hi <= cert["lam_hi"])
    contains_saved = bool(lam_lo <= cert["lam_lo"] and cert["lam_hi"] <= lam_hi)
    width = lam_hi - lam_lo
    width_log10 = vc.log10_int(width) - S * math.log10(2.0) if width > 0 else float("-inf")
    return {
        "lam_lo": lam_lo,
        "lam_hi": lam_hi,
        "width": width,
        "width_log10": width_log10,
        "digits_ok": digits_ok,
        "digits_proven_head": (digits_proven or "")[:120],
        "overlap": overlap,
        "subset_saved": subset_saved,
        "contains_saved": contains_saved,
        "a": a,
        "b": b,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def serialize_cmp(cmp: dict[str, Any]) -> dict[str, Any]:
    return {
        "lo": cmp["lo"],
        "hi": cmp["hi"],
        "ok": cmp["ok"],
        "even_ok": cmp["even_ok"],
        "odd_ok": cmp["odd_ok"],
        "even_saved": [str(v) for v in cmp["even_saved"]],
        "even_new": [str(v) for v in cmp["even_new"]],
        "odd_saved": [str(v) for v in cmp["odd_saved"]],
        "odd_new": [str(v) for v in cmp["odd_new"]],
        "left_even": str(cmp["left_even"]),
        "right_even": str(cmp["right_even"]),
        "left_odd": str(cmp["left_odd"]),
        "right_odd": str(cmp["right_odd"]),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Deep Arb verification of the 40,000-digit Golomb-Dickman checkpoint set.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python verify_checkpoints_arb.py --audit\n"
            "  python verify_checkpoints_arb.py --sample 20\n"
            "  python verify_checkpoints_arb.py --all --jobs 8\n"
            "  python verify_checkpoints_arb.py --nodes 20\n"
        ),
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=_bundle_root(Path(__file__).resolve().parent),
        help="directory containing the certificate, digits, b-file and checkpoint folder",
    )
    parser.add_argument("--ckpt", default=vc.CKPT_NAME, help="checkpoint directory name")
    parser.add_argument("--cert", default=vc.CERT_NAME, help="certificate file name")
    parser.add_argument("--digits", default=vc.DIGITS_NAME, help="digits file name")
    parser.add_argument("--bfile", default=vc.BFILE_NAME, help="OEIS b-file name")

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--audit", action="store_true", help="constants and certificate only (fast)")
    mode.add_argument("--all", action="store_true", help="recompute every chunk (slow)")
    mode.add_argument(
        "--sample",
        nargs="?",
        const=20,
        type=int,
        default=None,
        metavar="N",
        help="recompute N evenly spaced chunks (default 20)",
    )
    mode.add_argument(
        "--nodes",
        nargs="?",
        const=20,
        type=int,
        default=None,
        metavar="N",
        help="self-test N individual nodes by two formulas (default 20)",
    )

    parser.add_argument(
        "--extra",
        type=int,
        default=EXTRA_BITS_DEFAULT,
        help=f"extra target bits beyond Ebits (default {EXTRA_BITS_DEFAULT})",
    )
    parser.add_argument(
        "--u-series",
        type=float,
        default=U_SERIES_DEFAULT,
        help=f"series/continued-fraction threshold (default {U_SERIES_DEFAULT})",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="worker processes; 0 means all CPUs",
    )
    parser.add_argument("--cache-dir", type=Path, default=None, help="verifier cache directory")
    parser.add_argument("--no-cache", action="store_true", help="do not read or write the verifier cache")
    parser.add_argument("--json", type=Path, default=None, help="write a JSON report")

    args = parser.parse_args(argv)
    require_flint()

    root = args.root
    ckpt = root / args.ckpt
    failures: list[str] = []

    try:
        meta = vc.parse_meta(ckpt / vc.META_NAME)
        chunks, structural_errors = vc.scan_chunks(ckpt, meta)
        structural_errors.extend(vc.check_chunk_widths(chunks, meta))
        if structural_errors:
            raise vc.AuditError("\n".join(structural_errors))
        cert = vc.parse_certificate(root / args.cert, meta)
        vc.check_digits_file(root / args.digits, cert["digits"])
        vc.check_bfile(root / args.bfile, cert["digits"])
    except vc.AuditError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    print("Structural audit: OK")
    print(
        "  {n} chunks, {nodes} nodes, K_L={kl}, K_R={kr}, chunk={ch}, Ebits={eb}, S={s}".format(
            n=len(chunks),
            nodes=meta["nodes"],
            kl=meta["K_L"],
            kr=meta["K_R"],
            ch=meta["chunk"],
            eb=meta["Ebits"],
            s=meta["S"],
        )
    )

    Lsum_saved = sum(c["even"][0] + c["odd"][0] for c in chunks)
    Hsum_saved = sum(c["even"][1] + c["odd"][1] for c in chunks)
    hn, hd = meta["hn"], meta["hd"]
    lo_fin_saved = (hn * Lsum_saved) // hd
    hi_fin_saved = -((-hn * Hsum_saved) // hd)
    B_cert = lo_fin_saved - cert["lam_lo"]
    C_cert = cert["lam_hi"] - hi_fin_saved
    if B_cert < 0:
        failures.append(f"certificate lam_lo exceeds floor(h*Lsum) by {-B_cert}")
    if C_cert < 0:
        failures.append(f"certificate lam_hi is below ceil(h*Hsum) by {-C_cert}")

    print("Recomputing constants with Arb ...")
    plan = recompute_plan(meta)
    print(
        "  M ~ 10^{m}, B_disc ~ 10^{b}, tau_R ~ 10^{r}, tau_L ~ 10^{l}".format(
            m=plan["log10_M"],
            b=plan["log10_Bdisc"],
            r=plan["log10_tailR"],
            l=plan["log10_tailL"],
        )
    )
    for name, ok in plan["checks"].items():
        if not ok:
            failures.append(f"plan check failed: {name}")
    if not plan["h_reproduced"]:
        print(
            "  warning: h does not reproduce the float selection formula exactly; "
            "the explicit h in meta.pkl is still used and the inequalities are checked"
        )
    print(
        "  corrections from Arb: B={b}, TR={tr}, TL={tl}, TR+TL+B={c}".format(
            b=plan["B"], tr=plan["TR"], tl=plan["TL"], c=plan["TR"] + plan["TL"] + plan["B"]
        )
    )
    print(f"  corrections implied by certificate: B={B_cert}, C={C_cert}")
    if plan["B"] != B_cert:
        failures.append(f"recomputed B={plan['B']} differs from certificate B={B_cert}")
    if plan["TR"] + plan["TL"] + plan["B"] != C_cert:
        failures.append(
            "recomputed TR+TL+B={} differs from certificate C={}".format(
                plan["TR"] + plan["TL"] + plan["B"], C_cert
            )
        )

    if not plan["checks"].get("Bdisc_small", False):
        print("  FAIL: B_disc is not below 10^-(D+G+6)")
    if not plan["checks"].get("tailR_small", False):
        print("  FAIL: right tail is not below 10^-(D+G+8)")
    if not plan["checks"].get("tailL_small", False):
        print("  FAIL: left tail is not below 10^-(D+G+8)")

    selected: list[dict[str, Any]] = []
    mode_name = "audit"
    if args.nodes is not None:
        mode_name = "nodes"
    elif args.all:
        mode_name = "all"
    elif args.sample is not None:
        mode_name = "sample"

    if args.nodes is not None:
        tested, node_errors = node_selftest(meta, args.extra, args.u_series, args.nodes)
        failures.extend(node_errors)
        if tested == 0:
            failures.append("node self-test: no nodes were tested")
    elif args.all:
        selected = list(chunks)
    elif args.sample is not None:
        selected = sample_chunks(chunks, args.sample)
        print(f"  sample: {len(selected)} of {len(chunks)} chunks selected")

    params = build_params(meta, args.extra, args.u_series)
    signature = cache_signature(meta, args.extra, args.u_series)

    cache_dir: Path | None = None
    if not args.no_cache:
        if args.cache_dir is not None:
            cache_dir = args.cache_dir
        else:
            base = Path(
                os.environ.get("LOCALAPPDATA")
                or os.environ.get("TEMP")
                or os.environ.get("TMP")
                or "."
            )
            cache_dir = base / "gd_verify_cache" / f"D{meta['D']}"
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"  warning: cannot create cache directory {cache_dir}: {exc}")
            cache_dir = None

    results: dict[int, tuple[int, ...]] = {}
    if selected:
        jobs = args.jobs
        if jobs <= 0:
            jobs = os.cpu_count() or 1
        jobs = max(1, min(jobs, len(selected)))
        try:
            results = run_chunks(
                selected,
                params,
                signature,
                cache_dir,
                jobs,
                not args.no_cache,
            )
        except ArithmeticError as exc:
            print(f"FAIL: {exc}", file=sys.stderr)
            return 1

    chunk_failures: list[dict[str, Any]] = []
    for chunk in selected:
        result = results.get(chunk["lo"])
        if result is None:
            failures.append(f"chunk {chunk['lo']}: no result")
            continue
        cmp = compare_chunk(chunk, result)
        if not cmp["ok"]:
            chunk_failures.append(cmp)
            failures.append(
                "chunk {lo}: recomputed interval is not contained in the saved interval".format(
                    lo=chunk["lo"]
                )
            )

    if chunk_failures:
        print(f"  chunk containment failures: {len(chunk_failures)}")
        for cmp in chunk_failures[:10]:
            print(
                "    chunk {lo}: even_ok={eo} odd_ok={oo} "
                "left_even={le} right_even={re} left_odd={lo_} right_odd={ro}".format(
                    lo=cmp["lo"],
                    eo=cmp["even_ok"],
                    oo=cmp["odd_ok"],
                    le=cmp["left_even"],
                    re=cmp["right_even"],
                    lo_=cmp["left_odd"],
                    ro=cmp["right_odd"],
                )
            )
    elif selected:
        print(f"  chunk containment: OK for {len(selected)}/{len(chunks)} chunk(s)")

    new_intervals: dict[int, tuple[tuple[int, int], tuple[int, int]]] = {
        c["lo"]: (c["even"], c["odd"]) for c in chunks
    }
    for lo, result in results.items():
        new_intervals[lo] = ((result[2], result[3]), (result[4], result[5]))
    Lsum_new = sum(pair[0][0] + pair[1][0] for pair in new_intervals.values())
    Hsum_new = sum(pair[0][1] + pair[1][1] for pair in new_intervals.values())

    final = check_final(meta, cert, Lsum_new, Hsum_new, plan)
    if not final["digits_ok"]:
        failures.append("recomputed interval does not prove the published digits")
    if not final["overlap"]:
        failures.append("recomputed interval and certificate interval do not overlap")
    if not final["subset_saved"]:
        failures.append("recomputed interval is not contained in the certificate interval")
    if mode_name in ("audit", "nodes") and not selected:
        if final["lam_lo"] != cert["lam_lo"] or final["lam_hi"] != cert["lam_hi"]:
            failures.append("reconstructed audit interval differs from the certificate endpoints")

    print(
        "  final interval width: 10^{:.6f}; digits: {}".format(
            final["width_log10"], "OK" if final["digits_ok"] else "FAIL"
        )
    )
    print(
        "  recomputed interval inside certificate: {}; overlap: {}".format(
            final["subset_saved"], final["overlap"]
        )
    )

    report = {
        "mode": mode_name,
        "root": str(root),
        "meta": {
            "D": meta["D"],
            "G": meta["G"],
            "h": meta["h"],
            "h_ratio": [str(hn), str(hd)],
            "Ebits": meta["Ebits"],
            "S": meta["S"],
            "K_L": meta["K_L"],
            "K_R": meta["K_R"],
            "chunk": meta["chunk"],
            "nodes": meta["nodes"],
            "n_chunks": meta["n_chunks"],
        },
        "plan": {
            "log10_M": plan["log10_M"],
            "log10_Bdisc": plan["log10_Bdisc"],
            "log10_tailR": plan["log10_tailR"],
            "log10_tailL": plan["log10_tailL"],
            "B": str(plan["B"]),
            "TR": str(plan["TR"]),
            "TL": str(plan["TL"]),
            "TR_TL_B": str(plan["TR"] + plan["TL"] + plan["B"]),
            "h_reproduced": plan["h_reproduced"],
            "checks": plan["checks"],
        },
        "certificate_corrections": {"B": str(B_cert), "C": str(C_cert)},
        "chunks_total": len(chunks),
        "chunks_checked": len(selected),
        "chunk_failures": [serialize_cmp(c) for c in chunk_failures[:20]],
        "final": {
            "lam_lo": str(final["lam_lo"]),
            "lam_hi": str(final["lam_hi"]),
            "width": str(final["width"]),
            "width_log10": final["width_log10"],
            "digits_ok": final["digits_ok"],
            "digits_proven_head": final["digits_proven_head"],
            "overlap": final["overlap"],
            "subset_saved": final["subset_saved"],
        },
        "pass": not failures,
    }

    if args.json is not None:
        try:
            args.json.write_text(
                json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            print(f"  JSON report: {args.json}")
        except OSError as exc:
            print(f"  warning: cannot write JSON report: {exc}")

    if failures:
        print("", file=sys.stderr)
        print(f"FAIL: {len(failures)} check(s) failed", file=sys.stderr)
        for message in failures[:20]:
            print(f"  - {message}", file=sys.stderr)
        return 1

    print("")
    print("OK: Arb checkpoint verification passed")
    if mode_name == "sample":
        print(
            f"  note: only {len(selected)}/{len(chunks)} chunks were recomputed; "
            "run --all for the complete check"
        )
    if mode_name == "audit":
        print("  note: no node balls were recomputed; run --sample or --all for the deep check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())




