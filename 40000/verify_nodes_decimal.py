#!/usr/bin/env python3
"""Arb-free consistency check of saved chunk sums (standard library only).

This is an optional companion to ``verify_checkpoints.py`` and
``verify_checkpoints_arb.py``.  It recomputes selected chunk sums with the
``decimal`` module at user-chosen precision and compares the leading decimal
digits with the interval saved in the checkpoint.

It is deliberately NOT a proof.  Its job is to give an independent numerical
opinion on the saved intervals, so that a systematic formula or indexing error
would be caught even on a machine where Arb (``python-flint``) cannot be
installed.  The rigorous check remains ``verify_checkpoints_arb.py``.

Method
------
For each node ``k`` with ``t = k h`` the code recomputes

    u  = exp(t - exp(-t))
    S  = sum_{j>=2} (-u)^j / (j j!)
    G  = e^gamma u^2 (1 + exp(-t)) exp(-2u + S)

and sums ``G`` over the chunk.  ``gamma`` is taken from OEIS A001620.  Nodes
with ``u`` above ``--u-cut`` (default 350) are provably below ``10^-140`` for
the guard used here and are skipped; at ``--decimals 100`` that cannot affect
the compared digits.

The comparison is done in integers: with ``S`` the certificate scale and
``[L, H]`` the saved chunk interval, the script checks

    floor(10^P L / 2^S) == floor(10^P H / 2^S) == floor(10^P * recomputed).

If the saved interval is wider than ``10^-P`` (it is not, for the published
checkpoints), the check is relaxed to ``[lo, hi]`` containment.

Usage
-----
    python verify_nodes_decimal.py                 # 6 spread-out chunks
    python verify_nodes_decimal.py --chunks 12
    python verify_nodes_decimal.py --lo 0 --lo 20000
    python verify_nodes_decimal.py --methods       # formula self-test only

Exit code 0 = every selected chunk agreed with its comparison tolerance.
"""

from __future__ import annotations

import argparse
import math
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Iterable

import verify_checkpoints as vc

# Euler-Mascheroni constant, first 250 decimals, OEIS A001620.
GAMMA = Decimal(
    "0.577215664901532860606512090082402431042159335939923598805767234884"
    "867726777664670936947063291746749514631447249807082480960504014486"
    "542836224173997644923536253500333742937337737673942792595258247094"
    "916008735203948165670853233151776611528621"
)

LN10 = Decimal(10).ln()


def _guard(u: Decimal) -> int:
    """Extra digits needed because the series partial sums grow like e^u."""
    return int(u / LN10) + 60


def g_series(u: Decimal, t: Decimal, digits: int) -> Decimal:
    """G(u) via the alternating series, with generous guard digits."""
    guard = _guard(u)
    with localcontext() as ctx:
        ctx.prec = digits + guard
        em = (-t).exp()
        a = -u                                   # a_1 = -u
        total = Decimal(0)                       # sum from j = 2
        j = 1
        eps = Decimal(10) ** (-(digits + guard - 20))
        while True:
            j += 1
            a = -a * u * (j - 1) / (j * j)       # a_j -> a_{j+1}
            total += a
            if j + 1 >= 2 * u and abs(a) < eps:
                break
            if j > 400_000:
                raise RuntimeError("series did not converge")
        exponent = -2 * u + total
        ctx.prec = digits + 20
        return +(GAMMA.exp() * u * u * (1 + em) * exponent.exp())


def g_cf(u: Decimal, t: Decimal, digits: int, max_n: int = 60_000) -> Decimal:
    """G(u) via the Stieltjes continued fraction (slower for small u)."""
    with localcontext() as ctx:
        ctx.prec = digits + 40
        tail = Decimal(0)
        n = max_n
        while n >= 1:
            a = Decimal((n + 1) // 2)
            b = Decimal(1) if (n & 1) else u
            tail = a / (b + tail)
            n -= 1
        f = Decimal(1) / (u + tail)
        em = (-t).exp()
        e1 = (-u).exp() * f
        ctx.prec = digits + 20
        return +((-u - e1).exp() * u * (1 + em))


def g_asym(u: Decimal, t: Decimal, k: int, digits: int) -> Decimal:
    """G(u) via the asymptotic expansion of E1 (large u only)."""
    guard = _guard(u)
    k = max(1, k)
    with localcontext() as ctx:
        ctx.prec = digits + guard
        s = Decimal(0)
        for j in range(k - 1, 0, -1):
            s = 1 - s * j / u
        if k == 1:
            s = Decimal(1)
        em = (-t).exp()
        e1 = (-u).exp() * (s / u)
        ctx.prec = digits + 20
        return +((-u - e1).exp() * u * (1 + em))


def t_from_u(u: Decimal) -> Decimal:
    """Solve t - exp(-t) = ln u for t (Newton, monotone function)."""
    target = u.ln()
    t = target
    for _ in range(200):
        e = (-t).exp()
        step = (t - e - target) / (1 + e)
        t -= step
        if abs(step) < Decimal(10) ** -40:
            break
    return t


def methods_selftest(digits: int) -> list[str]:
    """Compare series, continued fraction and asymptotics on a grid of u."""
    errors: list[str] = []
    # The continued fraction converges slowly as u -> 0, so the comparison
    # starts at u = 0.1; below that the series is the method of choice.
    table = ["0.1", "0.3679", "1", "5", "20", "50", "100", "200"]
    print(f"  formula self-test ({digits} digits, series vs continued fraction):")
    for text in table:
        u = Decimal(text)
        t = t_from_u(u)
        gs = g_series(u, t, digits)
        gc = g_cf(u, t, digits)
        tol = Decimal(10) ** (-(digits - 8))
        ok = abs(gc - gs) <= tol * max(abs(gs), Decimal(1))
        print(f"    u={text:>7s}  G={gs:.20e}  |cf-series|={abs(gc - gs):.2e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            errors.append(f"self-test: series and continued fraction disagree at u={text}")

    print("  formula self-test (continued fraction vs asymptotics):")
    for text, k in (("2000", 2000), ("5000", 5000), ("20000", 20000)):
        u = Decimal(text)
        t = t_from_u(u)
        gc = g_cf(u, t, 40, 20000)
        ga = g_asym(u, t, k, 40)
        rel = abs(ga - gc) / abs(gc)
        ok = rel < Decimal(10) ** -30
        print(f"    u={text:>7s}  K={k:<6d}  relative difference={rel:.2e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            errors.append(f"self-test: asymptotics disagree at u={text}")
    return errors


def chunk_value(chunk: dict[str, Any], h: Decimal, digits: int, u_cut: Decimal) -> tuple[Decimal, int, int]:
    """Sum G over a chunk.  Returns (value, nodes evaluated, nodes skipped)."""
    lo, hi = chunk["lo"], chunk["hi"]
    total = Decimal(0)
    used = 0
    skipped = 0
    with localcontext() as ctx:
        ctx.prec = digits + 40
        for k in range(lo, hi):
            t = h * k
            em = (-t).exp()
            if t > 0:
                uf = math.exp(float(t) - math.exp(-float(t)))
            else:
                uf = math.exp(float(t)) if float(t) > -700 else 0.0
            if uf > float(u_cut):
                skipped += 1
                continue
            u = (t - em).exp()
            total += g_series(u, t, digits)
            used += 1
    return total, used, skipped


# ---------------------------------------------------------------------------
# Constants: discretisation error and tails (decimal estimate, no Arb)
# ---------------------------------------------------------------------------

def pi_decimal(prec: int) -> Decimal:
    """pi by the Gauss arithmetic-geometric mean."""
    with localcontext() as ctx:
        ctx.prec = prec + 20
        a = Decimal(1)
        b = Decimal(1) / Decimal(2).sqrt()
        t = Decimal(1) / 4
        p = Decimal(1)
        for _ in range(int(math.log2(prec)) + 6):
            an = (a + b) / 2
            b = (a * b).sqrt()
            t -= p * (a - an) ** 2
            a = an
            p *= 2
        ctx.prec = prec
        return +(a + b) ** 2 / (4 * t)


def cos_decimal(x: Decimal, prec: int) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = prec + 20
        term = Decimal(1)
        total = Decimal(1)
        n = 0
        eps = Decimal(10) ** (-(prec + 10))
        while True:
            n += 2
            term = -term * x * x / (n * (n - 1))
            total += term
            if abs(term) < eps:
                break
        ctx.prec = prec
        return +total


def sin_decimal(x: Decimal, prec: int) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = prec + 20
        term = x
        total = x
        n = 1
        eps = Decimal(10) ** (-(prec + 10))
        while True:
            n += 2
            term = -term * x * x / (n * (n - 1))
            total += term
            if abs(term) < eps:
                break
        ctx.prec = prec
        return +total


def strip_constants_decimal(prec: int) -> tuple[Decimal, Decimal]:
    """Lemma 4 constant M, mirroring ``strip_constants`` in golomb_dickman.py."""
    with localcontext() as ctx:
        ctx.prec = prec
        e = Decimal(1).exp()
        d = Decimal(7) / 5
        c = cos_decimal(d, prec)
        if c <= 0:
            raise ArithmeticError("cos(7/5) is not positive")
        a = 2 * c
        C1 = (GAMMA + 2 + (e - 2) / 2).exp()
        M1 = C1 * (-a).exp() * (2 / a + 1 / (a * a))
        e2 = Decimal(2).exp()
        em2 = (-Decimal(2)).exp()
        B2 = 2 * GAMMA.exp() * (e2 * e2) * (2 * e2 + e2.exp() / 2).exp()
        M2 = 2 * B2
        th = d + em2 * sin_decimal(d, prec)
        kap = cos_decimal(th, prec)
        if kap <= 0:
            raise ArithmeticError("kappa is not positive")
        kp = kap * (-em2).exp()
        facE = ((em2 - 2).exp()).exp()
        M3 = facE * (1 + em2) * (-kp * e2).exp() / kp
        return M1 + M2 + M3, th


def plan_check(meta: dict[str, Any], cert: dict[str, Any], decimals: int) -> tuple[bool, list[str]]:
    """Recompute B_disc, tau_R, tau_L and compare with the certificate."""
    notes: list[str] = []
    ok = True
    S = meta["S"]
    with localcontext() as ctx:
        ctx.prec = decimals + 80
        h = Decimal(meta["hn"]) / Decimal(meta["hd"])
        M, _th = strip_constants_decimal(decimals + 40)
        pi = pi_decimal(decimals + 40)
        d = Decimal(7) / 5
        Bdisc = 2 * M / ((2 * pi * d / h).exp() - 1)

        t_r = h * meta["K_R"]
        u_r = (t_r - (-t_r).exp()).exp()
        tail_r = (-u_r).exp()

        t_l = h * (-meta["K_L"])
        u_l = (t_l - (-t_l).exp()).exp()
        tail_l = GAMMA.exp() * u_l * u_l / 2

        scale = Decimal(1 << S)

        def to_upper(value: Decimal) -> int:
            return int((value * scale).to_integral_value(rounding="ROUND_FLOOR")) + 1

        B = to_upper(Bdisc)
        TR = to_upper(tail_r)
        TL = to_upper(tail_l)

        log10 = lambda v: (v.ln() / LN10) if v > 0 else Decimal("-Infinity")

    print(f"  log10 M      = {log10(M):.4f}")
    print(f"  log10 B_disc = {log10(Bdisc):.4f}   (requirement 10^-{meta['D'] + meta['G'] + 6})")
    print(f"  log10 tau_R  = {log10(tail_r):.4f}")
    print(f"  log10 tau_L  = {log10(tail_l):.4f}")
    print(f"  B  = {B}")
    print(f"  TR = {TR}")
    print(f"  TL = {TL}")
    print(f"  TR + TL + B = {TR + TL + B}")

    if Bdisc * (Decimal(10) ** (meta["D"] + meta["G"] + 6)) >= 1:
        notes.append("B_disc is not below 10^-(D+G+6)")
        ok = False
    if tail_r * (Decimal(10) ** (meta["D"] + meta["G"] + 8)) >= 1:
        notes.append("tau_R is not below 10^-(D+G+8)")
        ok = False
    if tail_l * (Decimal(10) ** (meta["D"] + meta["G"] + 8)) >= 1:
        notes.append("tau_L is not below 10^-(D+G+8)")
        ok = False

    implied_B = cert["implied_B"]
    implied_C = cert["implied_C"]
    print(f"  certificate implies B = {implied_B}, TR + TL + B = {implied_C}")
    if B != implied_B:
        notes.append(f"recomputed B={B} differs from certificate B={implied_B}")
        ok = False
    if TR + TL + B != implied_C:
        notes.append(f"recomputed TR+TL+B={TR + TL + B} differs from certificate {implied_C}")
        ok = False
    return ok, notes


def compare(chunk: dict[str, Any], value: Decimal, scale_bits: int, decimals: int) -> tuple[bool, str]:
    """Compare a recomputed value with the saved interval, in integers."""
    scale = 10 ** decimals
    l_int = (chunk["even"][0] + chunk["odd"][0]) * scale >> scale_bits
    h_int = (chunk["even"][1] + chunk["odd"][1]) * scale >> scale_bits
    with localcontext() as ctx:
        ctx.prec = decimals + 60
        v_int = int((value * scale).to_integral_value(rounding="ROUND_FLOOR"))
    if l_int == h_int:
        ok = v_int == l_int
        detail = f"saved pins {decimals} digits; recomputed {'agrees' if ok else 'differs'}"
    else:
        ok = l_int - 1 <= v_int <= h_int + 1
        detail = (
            f"saved interval spans {h_int - l_int} units at 10^-{decimals}; "
            f"recomputed {'inside' if ok else 'outside'}"
        )
    return ok, detail


def select_chunks(chunks: list[dict[str, Any]], h: Decimal, u_cut: Decimal, count: int) -> list[dict[str, Any]]:
    """Chunks whose nodes can all matter at the requested precision."""
    # t - exp(-t) = ln(u_cut)   ->  largest |k| with u <= u_cut
    target = float(u_cut.ln())
    lo_t, hi_t = 0.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo_t + hi_t)
        if mid - math.exp(-mid) < target:
            lo_t = mid
        else:
            hi_t = mid
    k_max = int((0.5 * (lo_t + hi_t)) / float(h))
    useful = [c for c in chunks if abs(c["lo"]) <= k_max or abs(c["hi"]) <= k_max]
    if not useful:
        return chunks[:count]
    if count >= len(useful):
        return useful
    step = (len(useful) - 1) / (count - 1) if count > 1 else 0.0
    picked = sorted({int(round(i * step)) for i in range(count)})
    return [useful[i] for i in picked]


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--ckpt", default=vc.CKPT_NAME)
    parser.add_argument("--chunks", type=int, default=6, help="number of spread-out chunks (default 6)")
    parser.add_argument(
        "--lo",
        type=int,
        action="append",
        default=None,
        help="explicit chunk start index; may be repeated",
    )
    parser.add_argument("--decimals", type=int, default=100, help="digits to compare (default 100)")
    parser.add_argument("--u-cut", type=float, default=350.0, help="skip nodes with u above this")
    parser.add_argument("--methods", action="store_true", help="run only the formula self-test")
    parser.add_argument("--plan", action="store_true", help="recompute B_disc and the two tails only")
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args(argv)

    try:
        meta = vc.parse_meta(args.root / args.ckpt / vc.META_NAME)
        chunks, errors = vc.scan_chunks(args.root / args.ckpt, meta)
        if errors:
            raise vc.AuditError("\n".join(errors))
    except vc.AuditError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    if args.methods:
        errors = methods_selftest(max(30, args.decimals // 2))
        if errors:
            for message in errors:
                print(f"FAIL: {message}", file=sys.stderr)
            return 1
        print("OK: formula self-test passed")
        return 0

    with localcontext() as ctx:
        ctx.prec = args.decimals + 80
        h = Decimal(meta["hn"]) / Decimal(meta["hd"])
    u_cut = Decimal(str(args.u_cut))

    if args.plan:
        try:
            cert = vc.parse_certificate(args.root / vc.CERT_NAME, meta)
        except vc.AuditError as exc:
            print(f"FAIL: {exc}", file=sys.stderr)
            return 1
        Lsum = sum(c["even"][0] + c["odd"][0] for c in chunks)
        Hsum = sum(c["even"][1] + c["odd"][1] for c in chunks)
        lo_fin = (meta["hn"] * Lsum) // meta["hd"]
        hi_fin = -((-meta["hn"] * Hsum) // meta["hd"])
        cert["implied_B"] = lo_fin - cert["lam_lo"]
        cert["implied_C"] = cert["lam_hi"] - hi_fin
        ok, notes = plan_check(meta, cert, args.decimals)
        if not ok:
            for message in notes:
                print(f"FAIL: {message}", file=sys.stderr)
            return 1
        print("OK: discretisation and tail bounds reproduce the certificate corrections")
        return 0

    if args.lo:
        by_lo = {c["lo"]: c for c in chunks}
        selected = []
        for lo in args.lo:
            chunk = by_lo.get(lo)
            if chunk is None:
                print(f"FAIL: no chunk starts at {lo}", file=sys.stderr)
                return 1
            selected.append(chunk)
    else:
        selected = select_chunks(chunks, h, u_cut, args.chunks)

    print(f"  precision: {args.decimals} decimals; h = {meta['hn']}/{meta['hd']}")
    print(f"  chunks: {len(selected)} of {len(chunks)}")

    failures = 0
    for chunk in selected:
        value, used, skipped = chunk_value(chunk, h, args.decimals, u_cut)
        ok, detail = compare(chunk, value, meta["S"], args.decimals)
        status = "OK  " if ok else "FAIL"
        print(
            f"  {status} chunk [{chunk['lo']}, {chunk['hi']}): "
            f"nodes used {used}, skipped {skipped}; sum ~ {value:.15e}"
        )
        print(f"         {detail}")
        if not ok:
            failures += 1

    if args.json is not None:
        import json

        args.json.write_text(
            json.dumps(
                {
                    "chunks": [c["lo"] for c in selected],
                    "decimals": args.decimals,
                    "failures": failures,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    if failures:
        print(f"FAIL: {failures} chunk(s) disagreed", file=sys.stderr)
        return 1
    print("")
    print("OK: recomputed chunk sums agree with the saved intervals")
    print("    (numerical consistency only; use verify_checkpoints_arb.py for the rigorous check)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
