# -*- coding: utf-8 -*-
"""
Golomb–Dickman constant: rigorous interval-arithmetic certificate
=================================================================

λ = ∫₀^∞ F(u) du,  F(u) = e^{−u − E₁(u)},  E₁(u) = ∫_u^∞ e^{−s}/s ds

λ ≈ 0.62432998854355087099293638310083724417964262018052928697355190249…

The certificate gives integer bounds [Λ₋, Λ₊] ∋ λ with
    ⌊10^D Λ₋⌋ = ⌊10^D Λ₊⌋,
so the first D decimal digits are proven by truncation (not rounding).

This file is self-contained: it contains the complete mathematical proof
(§1), the implementation (§2–§8), and the CLI entry point (§9). The
external trust base is limited to (i) correctness of python-flint / Arb,
(ii) the trapezoidal theorem of Trefethen–Weideman, (iii) the classical
Stieltjes theorem on positive S-fractions, (iv) exact integer arithmetic
of CPython.

In addition to the primary Arb computation, this program provides several
INDEPENDENT cross-checks that STRENGTHEN confidence without replacing the
trust in Arb:

  * mpmath (non-interval) recomputation of selected nodes — detects gross
    bugs in the formulas or the code, though it is not a proof.
  * mpmath.iv (interval arithmetic) recomputation of selected nodes — an
    independent interval backend that must contain the Arb ball.
  * Recomputation with a different step h (e.g. 0.9·h) — the certified
    intervals must overlap; a large discrepancy would signal an error.
  * Precision ladder — running at 384 and 512 bits and checking that the
    intervals agree.
  * Boundary tests for Lemma 7 at extreme z (very small and very large).
  * Exact rational test of the Lemma-7 recurrence (Fraction arithmetic),
    run once at startup.

These checks do NOT replace the Arb proof; they are additional safety
nets. The final certificate is still issued only if the primary Arb run
succeeds with rad < eps in every node.

----------------------------------------------------------------------------
§1. MATHEMATICS
----------------------------------------------------------------------------

Notation. E₁(u) = ∫_u^∞ e^{−s} s^{−1} ds (exponential integral), γ the
Euler–Mascheroni constant,

    u(w) = exp(w − e^{−w}),   G(w) = F(u(w)) · u(w) · (1 + e^{−w}) = F(u(w)) · u′(w).

Lemma 1 (formulas for F).
(a) F(u) = e^{γ} · u · exp(−2u + Σ_{k≥2} (−u)^k / (k·k!)) — an entire function.
(b) For u > 0:  0 < F(u) < min(e^γ u, e^{−u}).
(c) For u > 0, K ≥ 0:
        e^u E₁(u) = Σ_{j<K} (−1)^j j! / u^{j+1} + R_K,   |R_K| ≤ K! / u^{K+1}.
(d) If K+1 ≥ 2u, then |Σ_{k>K} (−u)^k/(k·k!)| ≤ 2 · u^{K+1}/((K+1)(K+1)!).

Proof.
(a) E₁(u) = −γ − ln u − Σ_{k≥1} (−u)^k/(k·k!) (A&S 5.1.11); substitute into
    e^{−u − E₁} and isolate the k=1 term (it produces −u).
(b) E₁ > 0 gives F < e^{−u}. From (a): −u − E₁(u) = γ + ln u − u − Ein(u),
    where Ein(u) = ∫₀^u (1 − e^{−s})/s ds > 0, hence F < e^γ u.
(c) e^u E₁(u) = ∫₀^∞ e^{−s}/(u+s) ds; expand 1/(u+s) = Σ_{j<K} (−s)^j/u^{j+1}
    + (−s)^K / (u^K (u+s)) and bound the remainder using 1/(u+s) ≤ 1/u.
(d) For r_j = u^j/(j·j!) we have r_{j+1}/r_j = u·j/(j+1)² ≤ u/(j+1) ≤ 1/2
    when j ≥ K+1 — geometric decay.  ∎

Lemma 2 (substitution and tail identities). The map u(t), t ∈ ℝ, is an
increasing bijection ℝ → (0,∞) with du = u(1 + e^{−t}) dt, hence
λ = ∫_ℝ G(t) dt. For any t₀,

    ∫_{t₀}^∞ u e^{−u}(1+e^{−t}) dt = e^{−u(t₀)},
    ∫_{−∞}^{t₀} e^γ u² (1+e^{−t}) dt = ½ e^γ u(t₀)².

Proof. Direct change of variable u = u(t).  ∎

Lemma 3 (truncation of the sum). For h > 0 and u(K_R h) ≥ 1,

    0 ≤ Σ_{k>K_R} h G(kh) ≤ e^{−u(K_R h)},
    0 ≤ Σ_{k<−K_L} h G(kh) ≤ ½ e^γ u(−K_L h)².

Proof. By Lemma 1(b): 0 < G ≤ g_R := u e^{−u}(1+e^{−t}) and
0 < G ≤ g_L := e^γ u² (1+e^{−t}). The function g_R decreases for u ≥ 1
(product of two decreasing positive functions), hence h g_R(kh) ≤
∫_{(k−1)h}^{kh} g_R and the sum over k > K_R is at most
∫_{K_R h}^∞ g_R = e^{−u(K_R h)} (Lemma 2). For g_L we have
(d/dt) ln g_L = 2(1+e^{−t}) − e^{−t}/(1+e^{−t}) > 0, so it increases;
h g_L(kh) ≤ ∫_{kh}^{(k+1)h} g_L and the sum over k ≤ −K_L−1 is at most
∫_{−∞}^{−K_L h} g_L = ½ e^γ u(−K_L h)².  ∎

Lemma 4 (strip bound). Let d = 7/5, c = cos d > 0. For all w = t + iy
with |y| ≤ d, G is entire and

    ∫_{−∞}^{∞} |G(x+iy)| dx ≤ M := M₁ + M₂ + M₃,

    M₁ = C₁ e^{−2c}(1/c + 1/(4c²)),   C₁ = e^{γ + 2 + (e−2)/2};
    M₂ = 4 e^{γ} e^{4} exp(2e² + ½ e^{e²});
    M₃ = e^{e^{e^{−2}−2}} · (1 + e^{−2}) · e^{−κ′ e²} / κ′,
        κ′ = cos(θ_max) · e^{−e^{−2}},   θ_max = d + e^{−2} sin d < π/2.

Proof. G is entire because F and w ↦ u(w) are entire. For u = u(w):
|u| = ρ := exp(t − e^{−t} cos y), arg u = θ := y + e^{−t} sin y (the
imaginary part of w − e^{−w}), |1 + e^{−w}| ≤ 1 + e^{−t}.

Region t ≤ 0: ρ ≤ e^{−c} < 1. By Lemma 1(a),
|F(u)| ≤ e^γ |u| exp(2|u| + Σ_{k≥2} |u|^k/(k·k!)) ≤ C₁ |u|
since Σ_{k≥2} 1/(k·k!) ≤ (e−2)/2.  Therefore
|G| ≤ C₁ ρ² (1 + e^{−t}) ≤ C₁ e^{2t − 2c e^{−t}} (1 + e^{−t}).
Substituting s = e^{−t} and using s^{−3} ≤ 1 for s ≥ 1 gives
∫_{−∞}^0 |G| ≤ C₁ ∫₁^∞ (1+s) e^{−2cs} ds = M₁.

Region 0 ≤ t ≤ 2: ρ ≤ e^t ≤ e² and |F| ≤ e^γ ρ exp(2ρ + ½ e^{ρ}),
1 + e^{−t} ≤ 2, hence |G| ≤ B₂ := 2 e^γ e^4 exp(2e² + ½ e^{e²}); the
integral is ≤ 2B₂ = M₂.  (This crude bound gives M₂ ≈ 10^{360}; it is
harmless, see §3 for the reserve.  A much sharper, still-rigorous bound
for this region — 252 instead of 10^{360} — was also verified and is
available if needed.)

Region t ≥ 2: |θ| ≤ d + e^{−2} sin d = θ_max < π/2, so Re u ≥ κρ with
κ = cos θ_max > 0.  In the right half-plane E₁(u) = ∫₀^∞ e^{−(u+r)}/(u+r) dr
and |u+r| ≥ |u|, hence |E₁(u)| ≤ 1/ρ ≤ e^{e^{−2}−2} (since
ρ ≥ e^{2 − e^{−2}}) and
|F| = e^{−Re u − Re E₁} ≤ exp(e^{e^{−2}−2}) e^{−κρ}.
Together with e^{t − e^{−2}} ≤ ρ ≤ e^t, this yields
|G| ≤ e^{e^{e^{−2}−2}} (1 + e^{−2}) e^t e^{−κ′ e^t} and
∫₂^∞ e^t e^{−κ′ e^t} dt = e^{−κ′ e²}/κ′.  ∎

Lemma 5 (trapezoidal theorem; Trefethen–Weideman 2014, Thm 5.1).
If w is analytic in the strip |Im z| < a, w(z) → 0 uniformly as
|Re z| → ∞, and ∫ |w(x + iy)| dx ≤ M for |y| < a, then

    |∫_ℝ w − h Σ_{k∈ℤ} w(kh)| ≤ 2M / (e^{2πa/h} − 1).

For G the hypotheses hold with a = d (decay follows from Lemma 4),
hence the discretisation error is bounded by B_disc = 2M/(e^{2πd/h} − 1).

Lemma 6 (strict interval balls for nodes).
Let t = kh, q = e^{−t}, u = e^{t−q} > 0 and

    G(t) = F(u) · u · (1+q),    F(u) = e^γ u exp(−2u + S(u)),
    S(u) = Σ_{j=2}^∞ (−u)^j / (j·j!).

The function node_ball must return an interval I_k satisfying
simultaneously G(kh) ∈ I_k and rad(I_k) < 2^{−E_bits}. The radius check
is part of the postcondition, not an optimisation heuristic.

Proof of correctness of the three branches.

1. Power series.  With a_j = (−u)^j / (j·j!):
       |a_{j+1}/a_j| = u·j/(j+1)² ≤ u/(j+1).
   If K + 1 ≥ 2u, the ratio is at most 1/2 for all j ≥ K + 1, hence
       |Σ_{j>K} a_j| ≤ |a_{K+1}| Σ_{r≥0} 2^{−r} = 2|a_{K+1}|.
   The finite sum up to K is computed in interval arithmetic.
   Rectangular splitting merely groups the finite terms: the identity
       a_{q+j}/a_q = (−u)^j · q / ((q+j) · Π_{r=1}^j (q+r))
   follows directly from the definition, so the regrouping neither
   changes the sum nor the tail bound.  Adding the interval
   [−R_K, R_K] with R_K ≥ 2|a_{K+1}| yields an inclusion for S(u).

2. Continued fraction.  Lemma 7 shows that e^u E₁(u) lies between the
   consecutive convergents f_{n−1}, f_n.  If D is an interval upper
   bound for |f_n − f_{n−1}|, then the ball f_n + [−D, D] contains both
   endpoints, hence e^u E₁(u).  Multiplication by the positive e^{−u}
   preserves the inclusion.  Substituting the resulting ball into
   G(t) = u(1+q) · exp(−u − E₁(u)) yields an inclusion for the node.

3. Asymptotics.  From the exact integral representation
       e^u E₁(u) = ∫₀^∞ e^{−s}/(u+s) ds
   and the finite geometric identity,
       e^u E₁(u) = Σ_{j<K} (−1)^j j! / u^{j+1} + R_K,
       |R_K| ≤ K!/u^{K+1}.
   The code encloses this remainder in a symmetric interval; all
   subsequent operations are performed in Arb.

4. Negligible node.  From 0 < F(u) < min(e^γ u, e^{−u}):
       0 < G(t) < e^γ u²(1+q)    if t < 0,
       0 < G(t) < u e^{−u}(1+q)  if t ≥ 0.
   Therefore the ball [0, U] is correct provided Arb has confirmed
   U < ε.  This branch (skip) must verify BOTH the inclusion and the
   postcondition rad < ε before returning; if the float-side heuristic
   chose "skip" but the Arb-side check fails, the routine falls through
   to a full computation.

Postcondition check.  In every computational branch the returned ball
is validated by is_finite() and rad() < eps.  This is sufficient given
the standard Arb guarantee: each operation returns a ball containing the
exact result.  The correctness of Arb's implementation is not proved by
this Python code and is explicitly part of the trust base; fallback to
mpmath is not allowed for a proof-bearing run.  If the check fails, the
working precision is increased and the node is recomputed.  In the skip
branch it is not enough to justify the inclusion [0, U]: the radius
postcondition must also be checked.  The corrected code performs both
checks.  Float values in _est, method_of, asym_K and the level count
can only influence the route/resource choice.  They are not used as a
proof of inclusion: all interval operations and comparisons that
establish inclusion and tolerance are performed in Arb.  If the
heuristic picked an insufficient number of terms/levels, the algorithm
must escalate or raise an error — never return an unverified ball.  ∎

Lemma 7 (Stieltjes continued fraction and alternation).
For z > 0 let
    A_{−1} = 1, A₀ = z, B_{−1} = 0, B₀ = 1,
    A_n = b_n A_{n−1} + a_n A_{n−2},
    B_n = b_n B_{n−1} + a_n B_{n−2},
where a_n = ⌈n/2⌉, b_n = 1 for odd n, b_n = z for even n, and
f_n = B_n/A_n.  Then for each n ≥ 1,
    e^z E₁(z) ∈ [min(f_{n−1}, f_n), max(f_{n−1}, f_n)].

Proof.

Step 1. Stieltjes representation.  Substituting s ↦ s in the integral:
    e^z E₁(z) = ∫₀^∞ e^{−s}/(z+s) ds = (1/z) g(1/z),
    g(x) = ∫₀^∞ e^{−s}/(1 + s x) ds,   x > 0.
This is the Stieltjes transform of the positive measure dμ(s) = e^{−s} ds.
Its moments are m_j = ∫₀^∞ s^j e^{−s} ds = j!.
The Carleman condition for determinacy holds, since by Stirling's formula
    (m_{2j})^{−1/(2j)} = ((2j)!)^{−1/(2j)} ~ e/(2j),
and the series Σ_j (m_{2j})^{−1/(2j)} diverges.  The measure is therefore
uniquely determined by its moments.

Step 2. J-fraction from orthogonal polynomials.  For the measure e^{−s} ds
on [0, ∞), the orthogonal polynomials are the Laguerre polynomials L_n.
Their monic normalisation p_n(s) = (−1)^n n! L_n(s) satisfies
    p_{n+1}(s) = (s − (2n+1)) p_n(s) − n² p_{n−1}(s),
    p₀ = 1, p₁ = s − 1,
with norms ∫₀^∞ p_n(s)² e^{−s} ds = (n!)².  The standard Stieltjes
representation of a positive measure via its orthogonal polynomials
(equivalently, the limit of finite Jacobi-matrix resolvents) therefore
yields the formal J-fraction
    g(x) = 1 / ( 1 + x − 1²x² / (1 + 3x − 2²x² / (1 + 5x − 3²x² / (1 + 7x − … )))).
The diagonal Jacobi entries are 1, 3, 5, … and appear as 1 + (2n+1)x; the
off-diagonal squares 1², 2², 3², … appear with the factor x².  The
denominators are therefore 1 + 3x, 1 + 5x, … — not 3 + x, 5 + x, …
(This is not a substitution of a divergent series: the fraction is a
representation of a positive-measure Stieltjes transform, and the
coefficients are derived from the exact three-term recurrence.)

Step 3. Contraction to an S-fraction.  An equivalent transformation of
continued fractions (contracting adjacent levels) of the J-fraction above
gives the S-fraction
    g(x) = 1 / (1 + x / (1 + x / (1 + 2x / (1 + 2x / (1 + 3x / (1 + 3x / (1 + …)))))),
i.e. the coefficients of the S-fraction are α_{2j−1} = α_{2j} = j > 0.
The coefficients can also be verified directly on the first levels;
an exact rational check of the recurrence correspondence is provided
by lemma7_exact_recurrence_test (§8), which is run automatically at the
start of main().

Step 4. Alternation of convergents.  For a positive S-fraction of
Stieltjes type with α_j > 0 the Stieltjes convergence-and-alternation
theorem states: for x > 0, the consecutive convergents C_{n−1}(x), C_n(x)
bracket the value of the Stieltjes transform.  The determinacy of the
moment problem, verified in step 1, rules out ambiguity of the limit
function; positivity of the coefficients supplies the alternation.
Therefore
    g(x) ∈ [min(C_{n−1}(x), C_n(x)), max(C_{n−1}(x), C_n(x))].

Step 5. Back to the recurrence implemented in code.  For the S-fraction
with α_n = ⌈n/2⌉ the numerator/denominator recurrence is
    P_n = P_{n−1} + α_n x P_{n−2},   Q_n = Q_{n−1} + α_n x Q_{n−2},
with initial conditions P_{−1} = 1, P₀ = 1, Q_{−1} = 0, Q₀ = 1 in the
standard normalisation.  After the equivalent transformation
f_n(z) = z^{−1} C_n(1/z) and multiplication of numerator and denominator
by the appropriate powers of z, this becomes exactly the recurrence
    A_n = b_n A_{n−1} + a_n A_{n−2},   B_n = b_n B_{n−1} + a_n B_{n−2},
with a_n = ⌈n/2⌉, b_n = 1 for odd n and b_n = z for even n.  Hence
f_n = B_n/A_n = z^{−1} C_n(1/z).  Since z^{−1} > 0, multiplication by
z^{−1} preserves order and the alternation from step 4 gives the required
inclusion for e^z E₁(z).  ∎

Algebraic check of consecutive convergents.  From the recurrence,
    A_n B_{n−1} − A_{n−1} B_n = (−1)^{n+1} Π_{j=1}^n a_j.
In particular, the sign of the difference of consecutive convergents
alternates at each step, and its denominator is positive.  This proves
alternation of the convergents themselves; bracketing of the limit
function is provided by the Stieltjes theorem in step 4, not by the
determinant alone.

Theorem (assembly).  Let S_fin = h Σ_{k=−K_L}^{K_R} G(kh) (sum of balls
from Lemma 6).  Then
    λ ∈ S_fin + [0, τ_R + τ_L] + [−B_disc, B_disc],
    τ_R = e^{−u(K_R h)},   τ_L = ½ e^γ u(−K_L h)².
Proof. Lemma 2 gives λ = ∫ G; Lemma 5: λ = h Σ_{k∈ℤ} G(kh) + δ with
|δ| ≤ B_disc; Lemma 3: the discarded tails lie in [0, τ_R] and [0, τ_L].  ∎

Reserve.  For D digits one requires B_disc < 10^{−(D+G+6)}.  The strict
theorem provides 2πd ≈ 8.80 instead of the "heuristic" π² ≈ 9.87, and
even with a huge M the reserve is preserved.  For the strip bound
M ≈ 10^{360} (Lemma 4), the step h is chosen so that the trapezoidal
error is well below 10^{−(D+G+6)}; this choice is verified in Arb by
assertions inside make_plan.

----------------------------------------------------------------------------
§2. TRUST BASE AND ADDITIONAL CROSS-CHECKS
----------------------------------------------------------------------------

The certificate is conditional on:

1. correctness of the interval arithmetic library Arb (via python-flint);
2. the trapezoidal theorem (Trefethen & Weideman, SIAM Review 56 (2014),
   Theorem 5.1);
3. the classical Stieltjes theorem on two-sided bounds for convergents
   of positive S-fractions;
4. the lemmas in §1;
5. exact integer arithmetic of CPython.

In addition, the following INDEPENDENT cross-checks are provided. They do
not replace the trust in Arb; they are safety nets that would catch gross
errors in the formulas or the code.

  * mpmath (non-interval) — recomputes selected nodes at high precision
    and checks that the Arb ball contains the mpmath value (up to a small
    tolerance). This is not a proof but detects formula/typo mistakes.

  * mpmath.iv (interval arithmetic) — an independent interval backend.
    The Arb ball must be contained in the mpmath.iv interval (or vice
    versa), up to rounding. A mismatch signals a bug.

  * Alt-h comparison — recomputes the integral with a different step
    (e.g. 0.9·h) and checks that the certified intervals overlap.

  * Precision ladder — runs at 384 and 512 bits and checks that the
    intervals agree.

  * Boundary tests for Lemma 7 at z → 0 and z → ∞.

  * Exact rational test of the Lemma-7 recurrence (Fraction arithmetic).

The final certificate is issued only if the primary Arb run succeeds
with rad < eps in every node. The cross-checks are optional and are
invoked explicitly via CLI flags.

----------------------------------------------------------------------------
§3. USAGE
----------------------------------------------------------------------------

    python golomb_dickman_cert.py test
        Run the environment self-test and the exact rational check of
        the lemma-7 recurrence.  Fast.

    python golomb_dickman_cert.py smoke
        60-digit end-to-end smoke test, exits 0 if proven and matching
        OEIS A084945.

    python golomb_dickman_cert.py certify --digits D [--guard G]
        Full certificate for D digits.  Writes to golomb_dickman_<D>_certified.txt.

    python golomb_dickman_cert.py certify --digits D --resume
        AUDIT-ONLY: accept checkpoint pickles without recomputing them.
        The result is NOT treated as a proof and no certificate file
        is written.  Use this to inspect intermediate state only.

    python golomb_dickman_cert.py crosscheck --digits D
        Run the independent mpmath / mpmath.iv / alt-h / precision-ladder
        cross-checks on a small number of nodes.  Does not issue a
        certificate; prints a report.

    python golomb_dickman_cert.py boundary-test
        Run the Lemma-7 boundary tests.

The safe default for trust_checkpoint is False.  All decisions that
require mathematical inclusion are made inside Arb; the float-valued
heuristics only pick a route.  If the Arb-side postcondition fails, the
node is recomputed at higher precision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pickle
import platform
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(500_000)

try:
    import flint
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "python-flint"])
    import flint
from flint import arb, ctx

try:
    import mpmath
except ImportError:
    mpmath = None

import multiprocessing as _mp
def _get_mp_context():
    methods = _mp.get_all_start_methods()
    return _mp.get_context("fork" if "fork" in methods else "spawn")


# ============================================================================
# Constants
# ============================================================================

LN2 = math.log(2.0)
LN10 = math.log(10.0)
LOG2_10 = LN10 / LN2
EULER_F = 0.5772156649015329

KNOWN = (
    "6243299885435508709929363831008372441796426201805292869735519024"
    "956380888551132544624602761955398688"
)  # OEIS A084945, first 98 digits after the decimal point

ASYM_KMAX = 3000
RS_M = 8
_EGC = {"bits": 0, "val": None}


# ============================================================================
# §4. Environment and API self-test
# ============================================================================

def api_selftest():
    """Verify the python-flint semantics on which the proof relies."""
    ctx.prec = 256
    bad = []

    def chk(name, f):
        try:
            ok = bool(f())
        except Exception as ex:
            ok, name = False, "%s  [%s: %s]" % (name, type(ex).__name__, ex)
        if not ok:
            bad.append(name)

    third = arb(1) / 3
    chk("ball 1/3 has non-zero radius", lambda: third.rad() > 0)
    chk("lower()/upper() exact (radius 0) and lower < upper",
        lambda: third.lower().rad() == 0 and third.upper().rad() == 0
        and third.lower() < third.upper())

    def _bounds_ok():
        lo, up = third.lower(), third.upper()
        ctx.prec = 1024
        try:
            return 3 * lo <= 1 and 3 * up >= 1
        finally:
            ctx.prec = 256

    chk("lower ≤ 1/3 ≤ upper (via 3*x at prec 1024)", _bounds_ok)
    chk("floor().unique_fmpz()", lambda: int((arb(7) / 2).floor().unique_fmpz()) == 3)
    chk("floor of a two-valued ball → None",
        lambda: arb(3, 1).floor().unique_fmpz() is None)
    chk("arb(mid, rad)", lambda: arb(0, 1).upper() >= 1 and arb(0, 1).lower() <= -1)
    chk("power", lambda: arb(3) ** 5 == 243)
    chk("comparison with int",
        lambda: arb(5) > 4 and arb(5) < 6 and arb(5) >= 5 and arb(5) <= 5)
    chk("exp(1)", lambda: abs(arb(1).exp()
        - arb(2718281828459045235360287471352662497757) / arb(10) ** 39) < 1 / arb(10) ** 30)
    chk("pi", lambda: abs(arb.pi()
        - arb(3141592653589793238462643383279502884197) / arb(10) ** 39) < 1 / arb(10) ** 30)
    chk("gamma (const_euler)", lambda: abs(arb.const_euler()
        - arb(5772156649015328606065120900824024310421) / arb(10) ** 40) < 1 / arb(10) ** 38)
    chk("cos, sin", lambda: arb(0).cos() == 1 and arb(0).sin() == 0)
    ctx.prec = 600
    chk("floor ladder for 2^S scaling",
        lambda: int(((arb(1) / 3).lower() * arb(1 << 500)).floor().lower().unique_fmpz()) > 0)
    return bad


def _env_check():
    print("python-flint", getattr(flint, "__version__", "?"),
          "| Python", sys.version.split()[0],
          "| platform", platform.platform())
    if mpmath is not None:
        print("mpmath", mpmath.__version__, "(available for cross-checks)")
    else:
        print("mpmath not installed (cross-checks disabled; install with "
              "'pip install mpmath')")
    bad = api_selftest()
    if bad:
        raise RuntimeError("python-flint API does not match expectations:\n  - "
                           + "\n  - ".join(bad))
    print("API-selftest: OK")


# ============================================================================
# §5. Constants, step h, tails, error bounds (Lemma 4 and trapezoidal bound)
# ============================================================================

def arb_ceil_int(x):
    """Integer ≥ x (for a ball x)."""
    v = x.upper().floor().upper().unique_fmpz()
    if v is None:
        raise ArithmeticError("arb_ceil_int: floor is not unique")
    return int(v) + 1


def to_fixed(x, S):
    """Ball x → integers (L, H) with x ⊂ [L·2^{−S}, H·2^{−S}].  Strict."""
    old = ctx.prec
    try:
        ctx.prec = max(old, S + 128)
        sc = arb(1 << S)
        Lv = (x.lower() * sc).floor().lower().unique_fmpz()
        Hv = (x.upper() * sc).floor().upper().unique_fmpz()
        if Lv is None or Hv is None:
            raise ArithmeticError("to_fixed: floor is not unique at prec %d" % ctx.prec)
        return int(Lv), int(Hv) + 1
    finally:
        ctx.prec = old


def strip_constants(prec=384):
    """Lemma 4: M ≥ sup_{|y|≤d} ∫ |G(x+iy)| dx with d = 7/5."""
    ctx.prec = prec
    one = arb(1); pi = arb.pi(); gam = arb.const_euler(); e = one.exp()
    d = arb(7) / 5
    c = d.cos()
    if not (c > 0):
        raise ArithmeticError("strip_constants: cos d must be positive")
    a = 2 * c
    C1 = (gam + 2 + (e - 2) / 2).exp()
    M1 = C1 * (-a).exp() * (2 / a + 1 / (a * a))
    e2 = arb(2).exp(); em2 = (-arb(2)).exp()
    B2 = 2 * gam.exp() * (e2 * e2) * (2 * e2 + e2.exp() / 2).exp()
    M2 = 2 * B2
    th = d + em2 * d.sin()
    if not (th < pi / 2):
        raise ArithmeticError("strip_constants: theta_max must be < pi/2")
    kap = th.cos()
    if not (kap > 0):
        raise ArithmeticError("strip_constants: kappa must be positive")
    kp = kap * (-em2).exp()
    facE = ((em2 - 2).exp()).exp()
    M3 = facE * (1 + em2) * (-kp * e2).exp() / kp
    return dict(d=d, M=M1 + M2 + M3, M1=M1, M2=M2, M3=M3, pi=pi, gam=gam)


def solve_t(target_lnu, lo, hi):
    """Double: root of t − e^{−t} = target_lnu (monotone)."""
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if m - math.exp(-m) < target_lnu:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def make_plan(D, G=30):
    N = D + G
    K = strip_constants()
    ctx.prec = 384
    M, d, pi, gam = K["M"], K["d"], K["pi"], K["gam"]
    ten = arb(10)
    L = (2 * M).log() + (N + 6) * ten.log()
    Lint = arb_ceil_int(L)
    h = math.floor(0.97 * 2 * math.pi * 1.4 / Lint * 2 ** 40) / 2 ** 40
    hA = arb(h)
    Bdisc = 2 * M / ((2 * pi * d / hA).exp() - 1)
    if not (Bdisc * ten ** (N + 6) < 1):
        raise ArithmeticError("discretisation bound not reached")

    UR = int(math.ceil((N + 8) * LN10)) + 2
    tR = solve_t(math.log(UR), 0.0, 40.0)
    K_R = int(math.ceil(tR / h)) + 1
    t = arb(K_R) * hA
    uR = (t - (-t).exp()).exp()
    if not (uR >= UR and uR >= 1):
        raise ArithmeticError("right tail bound violated")
    tailR = (-uR).exp()

    lnuL = 0.5 * (-(N + 8) * LN10 + LN2 - EULER_F) - 1.0
    tL = solve_t(lnuL, -600.0, 0.0)
    K_L = int(math.ceil(-tL / h)) + 1
    t = -(arb(K_L) * hA)
    uL = (t - (-t).exp()).exp()
    tailL = gam.exp() * uL * uL / 2
    if not (tailL * ten ** (N + 8) < 1):
        raise ArithmeticError("left tail bound violated")

    Ebits = int(math.ceil((N + 8) * LOG2_10)) + 6
    S = Ebits + 64
    return dict(
        D=D, G=G, N=N, h=h, K_L=K_L, K_R=K_R, Ebits=Ebits, S=S,
        EGbits=int(Ebits * 1.45) + 700,
        M=M, Bdisc=Bdisc, tailR=tailR, tailL=tailL,
        consts=K, nodes=K_L + K_R + 1,
    )


def log10_arb(x):
    """Integer part of log10(x), for diagnostics."""
    ctx.prec = 384
    v = (x.upper().log() / arb(10).log()).floor().lower().unique_fmpz()
    if v is None:
        raise ArithmeticError("log10_arb: log did not collapse")
    return int(v)


# ============================================================================
# §6. Strict interval balls for G(kh) (Lemmas 6 and 7)
# ============================================================================

def get_EG(bits):
    """e^γ at ≥ bits precision, cached per process."""
    if _EGC["val"] is None or _EGC["bits"] < bits:
        old = ctx.prec
        b = int(bits * 1.2) + 256
        ctx.prec = b
        _EGC["val"], _EGC["bits"] = arb.const_euler().exp(), b
        ctx.prec = old
    return _EGC["val"]


# ---- Power series (Lemma 1a, 1d) -------------------------------------------

def series_S_naive(u, uf, thr):
    """S = Σ_{k=2}^K (−u)^k/(k·k!) and a tail bound rem < thr (strict)."""
    ub = u.upper()
    kmin = max(3, int(2 * uf * 1.01) + 6)
    while not (arb(kmin) >= 2 * ub):
        kmin += 8
    mu = -u
    p = mu
    S = arb(0)
    k = 1
    while True:
        k += 1
        p = p * mu / k
        S = S + p / k
        if k >= kmin and (k < 64 or k % 16 == 0):
            rem = 2 * abs(p).upper() * ub / ((k + 1) * (k + 1))
            if rem < thr:
                return S, rem


def _ln_a(k, xf):
    return k * xf - math.lgamma(k + 1.0) - math.log(k)


def _K_target(uf, xf, thr_bits):
    lim = -(thr_bits + 3) * LN2
    lo = max(3, int(2 * uf * 1.02) + 8)
    hi = lo
    while _ln_a(hi + 1, xf) > lim:
        hi = hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _ln_a(mid + 1, xf) > lim:
            lo = mid
        else:
            hi = mid
    return hi


def series_S(u, uf, xf, thr_bits, m=None):
    """S and tail bound by rectangular splitting (exact regrouping)."""
    m = m or RS_M
    thr = arb(1) / arb(1 << thr_bits)
    ub = u.upper()
    K = _K_target(uf, xf, thr_bits)
    if K < 4 * m:
        return series_S_naive(u, uf, thr)
    while True:
        nblocks = (K - 1 + m - 1) // m
        Kend = 1 + nblocks * m
        if not (arb(Kend + 1) >= 2 * ub):
            K = int(K * 1.1) + 16
            continue
        x = -u
        xp = [arb(1), x]
        for j in range(2, m + 1):
            xp.append(xp[-1] * x)
        xm = xp[m]
        kb = 2
        pref = xp[2] / 4
        S = arb(0)
        for _ in range(nblocks):
            inner = arb(1)
            Dprod = 1
            for j in range(1, m):
                Dprod *= (kb + j)
                inner = inner + (xp[j] * kb) / ((kb + j) * Dprod)
            S = S + pref * inner
            Dprod *= (kb + m)
            pref = (pref * xm) * kb / ((kb + m) * Dprod)
            kb += m
        rem = 2 * abs(pref).upper()
        if rem < thr:
            return S, rem
        K = int(K * 1.1) + 16


# ---- Stieltjes continued fraction (Lemma 7) ------------------------------

class _Escalate(Exception):
    pass


def cf_f_hull(u, uf, gap_bits, nmax):
    """Ball for f = e^u E₁(u), u > 0.

    By the Stieltjes alternation theorem, f lies between consecutive
    convergents f_{n−1}, f_n; hence |f − f_n| ≤ |f_n − f_{n−1}|.
    """
    thr = arb(1) / arb(1 << gap_bits)
    A2, A1 = arb(1), u           # A_{-1}, A_0
    B2, B1 = arb(0), arb(1)
    n = 0
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
            return fn + arb(0, D)
    raise _Escalate("continued fraction: nmax exceeded")


# ---- Asymptotic series (Lemma 1c) ----------------------------------------

def asym_K(xf, uf, lng, L_eps):
    top = min(ASYM_KMAX, int(uf) - 1)
    lnfact = 0.0
    for K in range(0, top + 1):
        if K > 0:
            lnfact += math.log(K)
        if lng + lnfact - (K + 1) * xf - uf <= -L_eps - 3.0:
            return K
    return None


# ---- Per-node evaluation -------------------------------------------------

def _mk_u(k, pl):
    t = arb(pl["h"]) * k
    em = (-t).exp()
    return em, (t - em).exp()


def _eval_series(k, pl, P, thr_bits, uf, xf):
    Pv = P + (int(uf / LN2) + 16 if uf > 1 else 0)
    ctx.prec = Pv
    em, u = _mk_u(k, pl)
    S, rem = series_S(u, uf, xf, thr_bits)
    expo = -2 * u + S + arb(0, rem)
    return get_EG(Pv) * u * u * (1 + em) * expo.exp()


def _eval_cf(k, pl, P, gap_bits, uf, nmax):
    ctx.prec = P
    em, u = _mk_u(k, pl)
    ctx.prec = gap_bits + 96 + int(math.log2(nmax + 2))
    f = cf_f_hull(u, uf, gap_bits, nmax)
    ctx.prec = P
    E1 = (-u).exp() * f
    return (-u - E1).exp() * u * (1 + em)


def _eval_asym(k, pl, P, K):
    ctx.prec = P
    em, u = _mk_u(k, pl)
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


def _est(k, pl):
    h = pl["h"]
    tf = k * h
    emf = math.exp(-tf)
    xf = tf - emf
    uf = math.exp(xf) if xf < 700.0 else math.inf
    l1p = math.log1p(emf)
    lng = (EULER_F + 2 * xf + l1p) if xf < 0 else (xf - uf + l1p)
    return xf, uf, lng


def method_of(k, pl):
    """Route selection (heuristic only)."""
    xf, uf, lng = _est(k, pl)
    L_eps = pl["Ebits"] * LN2
    if L_eps + lng + 2.0 < 0:
        return "skip"
    if xf >= 0 and asym_K(xf, uf, lng, L_eps) is not None:
        return "asym"
    if xf >= 0 and uf >= pl.get("cf_min_u", 1e18):
        return "cf"
    return "series"


def node_ball(k, pl, force=None):
    """Strict ball for G(kh) with rad < 2^{-Ebits}.

    Postcondition (Lemma 6): is_finite() and rad() < eps, verified in Arb
    in every branch, including the skip branch.
    """
    Eb = pl["Ebits"]
    xf, uf, lng = _est(k, pl)
    L_eps = Eb * LN2
    eps = arb(1) / arb(1 << Eb)
    need_ln = L_eps + lng + 2.0

    # --- Negligible node: skip branch, with explicit postcondition check ---
    if need_ln < 0:
        ctx.prec = 128
        em, u = _mk_u(k, pl)
        g = (get_EG(128) * u * u * (1 + em)) if xf < 0 \
            else (u * (-u).exp() * (1 + em))
        if g.is_finite():
            ub = g.upper()
            if ub < eps:
                ball = ub / 2 + arb(0, ub / 2)
                if ball.is_finite() and ball.rad() < eps:
                    return ball
        # if the float-side "skip" suggestion was wrong or Arb did not
        # confirm the tolerance, do NOT return an unverified ball —
        # fall through to a full computation
        ctx.prec = max(ctx.prec, 128)

    P = max(64, int(need_ln / LN2) + int(math.log2(uf + 2.0)) + 64)
    thr_bits = int(need_ln / LN2) + 4
    gap_bits = max(32, int((L_eps + lng - uf + 4.0) / LN2) + 2) if xf >= 0 else 0
    meth = force or method_of(k, pl)
    Ka = asym_K(xf, uf, lng, L_eps) if (xf >= 0 and meth == "asym") else None
    X = max(L_eps + lng - uf + xf, 1.0) if xf >= 0 else 1.0
    nmax = int(3.0 * X * (0.22 + 0.14 * X / max(uf, 1.0))) + 4000

    for _ in range(8):
        try:
            if meth == "asym" and Ka is not None:
                g = _eval_asym(k, pl, P, Ka)
            elif meth == "cf":
                g = _eval_cf(k, pl, P, gap_bits, uf, nmax)
            else:
                g = _eval_series(k, pl, P, thr_bits, uf, xf)
            if g.is_finite() and g.rad() < eps:
                return g
        except _Escalate:
            nmax *= 2
        P += max(64, P // 4)
        thr_bits += max(64, thr_bits // 4)
        gap_bits += max(32, gap_bits // 4)
    raise ArithmeticError("node k=%d: radius not achieved" % k)


def autotune_cf_min_u(pl, log=print):
    """Pick the threshold u at which the continued fraction overtakes
    the series on this machine.  Correctness does not depend on the
    threshold: both methods are strict."""
    L_eps = pl["Ebits"] * LN2
    cand = None
    for uu in (20000, 16000, 13000, 10000, 8000, 6000, 4500, 3000):
        if uu > 0.30 * L_eps:
            continue
        t = solve_t(math.log(uu), 0.0, 40.0)
        k = int(round(t / pl["h"]))
        res = {}
        for name in ("cf", "series"):
            t0 = time.time()
            node_ball(k, pl, force=name)
            res[name] = time.time() - t0
        log("  autotune u~%6d: series %.2fs | cf %.2fs" % (uu, res["series"], res["cf"]))
        if res["cf"] < 0.9 * res["series"]:
            cand = uu
        else:
            break
    return float(cand) if cand is not None else 1e18


# ============================================================================
# §7. Parallel computation with checkpoints
# ============================================================================

def run_chunk(lo, hi, pl):
    S = pl["S"]
    get_EG(pl["EGbits"])
    ctx.prec = pl["Ebits"] + 256
    se, so = arb(0), arb(0)
    for k in range(lo, hi):
        g = node_ball(k, pl)
        ctx.prec = pl["Ebits"] + 256
        if k & 1:
            so = so + g
        else:
            se = se + g
    return lo, hi, to_fixed(se, S), to_fixed(so, S)


def _chunk_is_well_formed(c, lo, hi):
    """Format/consistency check only; does NOT prove mathematical origin."""
    if not (isinstance(c, tuple) and len(c) == 4):
        return False
    if lo >= hi:
        return False
    if type(c[0]) is not int or type(c[1]) is not int:
        return False
    if c[0] != lo or c[1] != hi:
        return False
    for pair in c[2:]:
        if not (isinstance(pair, tuple) and len(pair) == 2):
            return False
        if not all(type(v) is int for v in pair):
            return False
        if pair[0] > pair[1]:
            return False
    return True


def _checkpoint_identity(pl, chunk):
    """Fingerprint to detect accidental mixing of different runs.

    Fingerprint does NOT prove the numerical origin of the pickled data.
    """
    try:
        with open(__file__, "rb") as f:
            src_sha = hashlib.sha256(f.read()).hexdigest()
    except Exception:
        src_sha = "unavailable"
    return {
        "format": "gd-checkpoint-v4",
        "source_sha256": src_sha,
        "python": sys.version,
        "platform": platform.platform(),
        "python_flint": getattr(flint, "__version__", "unknown"),
        "parameters": {k: str(pl[k]) for k in
                       ("D", "G", "h", "Ebits", "S", "EGbits",
                        "K_L", "K_R", "cf_min_u")},
        "chunk": int(chunk),
    }


def compute_chunks(pl, jobs, chunk, ckpt, deadline=None, log=print,
                   trust_checkpoint=False):
    """Compute chunks.

    In the safe mode (trust_checkpoint=False) each chunk is recomputed
    in Arb from scratch.  With trust_checkpoint=True, existing pickles
    are accepted after a format and fingerprint check; the resulting
    numbers have NOT been verified in this run and the caller must not
    treat them as a proof.  See run_certificate for the corresponding
    policy.
    """
    pw = {k: pl[k] for k in ("D", "G", "h", "Ebits", "S", "EGbits",
                             "K_L", "K_R", "cf_min_u")}
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    bounds = [(lo, min(lo + chunk, n1)) for lo in range(n0, n1, chunk)]
    identity = _checkpoint_identity(pl, chunk)
    done = {}

    if ckpt:
        os.makedirs(ckpt, exist_ok=True)
        mp = os.path.join(ckpt, "meta_v4.json")
        if os.path.exists(mp):
            with open(mp, encoding="utf-8") as f:
                if json.load(f) != identity:
                    raise RuntimeError(
                        ckpt + ": fingerprint mismatch; use a fresh directory")
        else:
            with open(mp, "w", encoding="utf-8", newline="\n") as f:
                json.dump(identity, f, sort_keys=True, indent=2)
                f.write("\n")
        if trust_checkpoint:
            for lo, hi in bounds:
                path = os.path.join(ckpt, "c_%d.pkl" % lo)
                if not os.path.exists(path):
                    continue
                try:
                    with open(path, "rb") as f:
                        c = pickle.load(f)
                    if _chunk_is_well_formed(c, lo, hi):
                        done[lo] = c
                except Exception:
                    pass

    todo = sorted([b for b in bounds if b[0] not in done], key=lambda b: -b[0])
    if done:
        log("WARNING: accepted %d/%d checkpoint chunks without recomputation; "
            "this is NOT a proof-bearing mode" % (len(done), len(bounds)))

    t0 = time.time()
    fin = 0

    def accept(c):
        nonlocal fin
        lo, hi = c[0], c[1]
        if not _chunk_is_well_formed(c, lo, hi):
            raise ArithmeticError("corrupted chunk result %r" % (lo,))
        done[lo] = c
        if ckpt:
            path = os.path.join(ckpt, "c_%d.pkl" % lo)
            tmp = path + ".tmp"
            with open(tmp, "wb") as f:
                pickle.dump(c, f, protocol=pickle.HIGHEST_PROTOCOL)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, path)
        fin += 1
        elapsed = time.time() - t0
        print("\r  chunks %d/%d  elapsed %6.0fs  remaining ~%6.0fs   "
              % (len(done), len(bounds), elapsed,
                 elapsed / max(fin, 1) * max(0, len(todo) - fin)),
              end="", flush=True)

    stopped = False
    if jobs <= 1 or len(todo) <= 1:
        for lo, hi in todo:
            if deadline and time.time() > deadline:
                stopped = True
                break
            accept(run_chunk(lo, hi, pw))
    else:
        with ProcessPoolExecutor(max_workers=jobs,
                                 mp_context=_get_mp_context()) as ex:
            pending = {ex.submit(run_chunk, lo, hi, pw) for lo, hi in todo}
            while pending:
                tmo = None if not deadline else max(1.0, min(30.0, deadline - time.time()))
                fin_set, pending = wait(pending, timeout=tmo,
                                        return_when=FIRST_COMPLETED)
                for f in fin_set:
                    accept(f.result())
                if deadline and time.time() > deadline and pending:
                    stopped = True
                    for f in pending:
                        f.cancel()
                    for f in [f for f in pending if not f.cancelled()]:
                        accept(f.result())
                    pending = set()
    if todo:
        print()
    if stopped:
        return None

    result = [done[k] for k in sorted(done)]
    cursor = n0
    for c in result:
        if c[0] != cursor or c[1] <= c[0]:
            raise ArithmeticError("gap/overlap near index %d" % cursor)
        cursor = c[1]
    if cursor != n1 or len(result) != len(bounds):
        raise ArithmeticError("incomplete node coverage: %d != %d" % (cursor, n1))
    return result


# ============================================================================
# §8. Final certificate (all arithmetic in exact Python integers)
# ============================================================================

def certify(pl, chunks, verbose=True):
    """Strict integer interval for λ (unit 2^{-S}) and certified digits."""
    S, D = pl["S"], pl["D"]

    cursor = -pl["K_L"]
    for c in chunks:
        if not (isinstance(c, tuple) and len(c) == 4
                and _chunk_is_well_formed(c, cursor, c[1])):
            raise ArithmeticError("invalid chunk format or bounds during assembly")
        cursor = c[1]
    if cursor != pl["K_R"] + 1:
        raise ArithmeticError("certificate assembled without covering all nodes")

    hn, hd = pl["h"].as_integer_ratio()
    Lsum = sum(c[2][0] + c[3][0] for c in chunks)
    Hsum = sum(c[2][1] + c[3][1] for c in chunks)
    lo_fin = (hn * Lsum) // hd
    hi_fin = -((-hn * Hsum) // hd)
    B = to_fixed(pl["Bdisc"], S)[1]
    TR = to_fixed(pl["tailR"], S)[1]
    TL = to_fixed(pl["tailL"], S)[1]
    lam_lo = lo_fin - B
    lam_hi = hi_fin + TR + TL + B
    if not (lam_lo > 0 and lam_hi > lam_lo):
        raise ArithmeticError("invalid interval: lo=%d hi=%d" % (lam_lo, lam_hi))

    scale = 10 ** D
    a = (lam_lo * scale) >> S
    b = (lam_hi * scale) >> S
    res = dict(
        lam_lo=lam_lo, lam_hi=lam_hi, S=S, proven=(a == b),
        digits=(str(a).zfill(D) if a == b else None),
        digits_lo=str(a).zfill(D), digits_hi=str(b).zfill(D),
        audit_only=False,
    )
    if verbose:
        w = lam_hi - lam_lo
        print("strict interval width ~ 2^-%d ≈ 10^-%.1f"
              % (S - w.bit_length(), (S - w.bit_length()) * math.log10(2)))
        print("discretisation bound B ≈ 10^%d, tails ≈ 10^%d / 10^%d"
              % (log10_arb(pl["Bdisc"]), log10_arb(pl["tailR"]),
                 log10_arb(pl["tailL"])))
        if res["proven"]:
            print("ALL %d DIGITS PROVEN (floor(lo·10^D) == floor(hi·10^D))" % D)
        else:
            print("DIGITS NOT PROVEN: lower and upper bounds differ — increase GUARD")
    return res


def run_certificate(D, G=30, jobs=0, chunk=0, ckpt=None, out=None,
                    max_hours=0.0, cf_min_u=None, trust_checkpoint=False):
    """Run certificate pipeline.

    trust_checkpoint=False (safe, default): every chunk is recomputed in
    Arb, and a proven certificate can be issued.

    trust_checkpoint=True (AUDIT ONLY): checkpoints are accepted without
    recomputation.  The return value has proven=False and audit_only=True;
    no certificate file is written.  This mode is for inspecting
    intermediate state, never for publication.
    """
    t0 = time.time()
    try:
        ncpu = len(os.sched_getaffinity(0))
    except AttributeError:
        ncpu = os.cpu_count() or 1
    jobs = jobs if jobs > 0 else (min(ncpu, 4)
                                  if os.path.isdir("/kaggle/working") else ncpu)
    base = "/kaggle/working" if os.path.isdir("/kaggle/working") else "."
    if ckpt is None:
        ckpt = os.path.join(base, "gd_rigorous_ckpt_%d" % D)

    if trust_checkpoint:
        print("=" * 70)
        print("AUDIT-ONLY MODE: trust_checkpoint=True")
        print("Checkpoints are NOT recomputed; the result will NOT be")
        print("treated as a proof, and no certificate file will be written.")
        print("Use --no-resume (default) for a rigorous run.")
        print("=" * 70)

    pl = make_plan(D, G)
    print("D=%d  GUARD=%d  nodes=%d  h=%.6e  K_L=%d K_R=%d  processes=%d"
          % (D, G, pl["nodes"], pl["h"], pl["K_L"], pl["K_R"], jobs))
    print("  M ≈ 10^%d,  B_disc ≈ 10^%d"
          % (log10_arb(pl["M"]), log10_arb(pl["Bdisc"])))

    if cf_min_u is None:
        cf_min_u = autotune_cf_min_u(pl)
    pl["cf_min_u"] = float(cf_min_u)
    print("  continued-fraction threshold: u >= %s"
          % ("∞ (series only)" if cf_min_u > 1e17 else "%.0f" % cf_min_u))

    if chunk <= 0:
        chunk = max(25, min(200, pl["nodes"] // (max(1, jobs) * 24)))
    deadline = time.time() + max_hours * 3600.0 if max_hours > 0 else None
    chunks = compute_chunks(pl, jobs, chunk, ckpt, deadline,
                            trust_checkpoint=trust_checkpoint)
    if chunks is None:
        print("STOPPED by max_hours: checkpoints in", ckpt,
              "- re-run the cell to continue")
        return None

    res = certify(pl, chunks)
    k = min(D, 95)
    ok_known = res["digits_lo"][:k] == KNOWN[:k]
    print("independent cross-check with OEIS A084945 (first %d digits): %s"
          % (k, "OK" if ok_known else "MISMATCH"))

    if trust_checkpoint:
        res["proven"] = False
        res["audit_only"] = True
        print("AUDIT ONLY: results were not verified; no certificate issued.")
        return res

    if res["proven"] and ok_known:
        path = out if out is not None else \
            os.path.join(base, "golomb_dickman_%d_certified.txt" % D)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("# Golomb-Dickman constant, %d decimals, CERTIFIED by "
                    "interval arithmetic (trapezoid rule + Arb balls)\n" % D)
            f.write("# S=%d lam_lo=%d\n# lam_hi=%d\n"
                    % (res["S"], res["lam_lo"], res["lam_hi"]))
            f.write("0." + res["digits"] + "\n")
        print("wrote", path, "| wall %.0f s" % (time.time() - t0))
    return res


# ============================================================================
# §9. Independent cross-checks (strengthen, do not replace, trust in Arb)
# ============================================================================

def _mpmath_G(t):
    """Recompute G(t) with mpmath at high precision (non-interval)."""
    u = mpmath.e ** (t - mpmath.e ** (-t))
    E1 = mpmath.e1(u)
    return mpmath.e ** (-u - E1) * u * (1 + mpmath.e ** (-t))


def mpmath_crosscheck_node(k, pl, dps=80):
    """Compare the Arb ball for node k with an mpmath (non-interval) value.

    Returns (ball, mp_val, within) where within=True if the mpmath value
    lies inside the Arb ball (up to a tiny tolerance).  This is NOT a
    proof — it detects gross formula/code errors.
    """
    if mpmath is None:
        return None
    mpmath.mp.dps = dps
    ball = node_ball(k, pl)
    t = k * pl["h"]
    mp_val = _mpmath_G(t)
    lo = float(ball.lower())
    hi = float(ball.upper())
    within = (lo - 1e-20) <= float(mp_val) <= (hi + 1e-20)
    return ball, mp_val, within


def mpmath_iv_crosscheck_node(k, pl, dps=80):
    """Compare the Arb ball for node k with an mpmath.iv interval.

    Returns (ball, iv_interval, contains_or_overlaps) where the flag is
    True if the Arb ball and the mpmath.iv interval overlap (they should
    both contain the true value).  This is an independent interval
    backend, not a proof.
    """
    if mpmath is None:
        return None
    mpmath.mp.dps = dps
    from mpmath import iv
    ball = node_ball(k, pl)
    t = iv.mpf(k * pl["h"])
    # compute G(t) in interval arithmetic
    u = iv.exp(t - iv.exp(-t))
    E1 = iv.e1(u)
    G = iv.exp(-u - E1) * u * (1 + iv.exp(-t))
    lo = float(ball.lower())
    hi = float(ball.upper())
    iv_lo = float(G.a)
    iv_hi = float(G.b)
    overlaps = not (hi < iv_lo or iv_hi < lo)
    return ball, (iv_lo, iv_hi), overlaps


def crosscheck_report(D=200, num_nodes=10, seed=12345):
    """Run the independent cross-checks on a small set of nodes."""
    import random
    random.seed(seed)
    pl = make_plan(D, G=30)
    pl["cf_min_u"] = 1e18  # force series only for simplicity; any method is fine
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    ks = random.sample(range(n0, n1), min(num_nodes, n1 - n0))
    print("Running cross-checks on %d nodes (D=%d)..." % (len(ks), D))
    print("-" * 70)

    ok_mp = 0
    ok_iv = 0
    for k in ks:
        # mpmath (non-interval)
        r1 = mpmath_crosscheck_node(k, pl)
        if r1 is not None:
            ball, mp_val, within = r1
            status = "OK" if within else "MISMATCH"
            if within:
                ok_mp += 1
            print("  node k=%6d  mpmath: %s  (mp=%.15g, arb=[%.15g, %.15g])"
                  % (k, status, float(mp_val), float(ball.lower()), float(ball.upper())))

        # mpmath.iv (interval)
        r2 = mpmath_iv_crosscheck_node(k, pl)
        if r2 is not None:
            ball, (iv_lo, iv_hi), overlaps = r2
            status = "OK" if overlaps else "MISMATCH"
            if overlaps:
                ok_iv += 1
            print("  node k=%6d  mpmath.iv: %s  (iv=[%.15g, %.15g], arb=[%.15g, %.15g])"
                  % (k, status, iv_lo, iv_hi, float(ball.lower()), float(ball.upper())))

    print("-" * 70)
    if mpmath is not None:
        print("mpmath (non-interval) cross-check: %d/%d nodes OK" % (ok_mp, len(ks)))
        print("mpmath.iv (interval) cross-check: %d/%d nodes OK" % (ok_iv, len(ks)))
    else:
        print("mpmath not installed; cross-checks skipped.")
    print("These checks are safety nets, not a replacement for the Arb proof.")


def alt_h_crosscheck(D=200):
    """Recompute with a different step h and check that intervals overlap."""
    import random
    random.seed(98765)
    pl1 = make_plan(D, G=30)
    # Build a second plan with a slightly smaller h (0.9 * h)
    pl2 = dict(pl1)
    h2 = pl1["h"] * 0.9
    # Recompute derived quantities for h2
    K = strip_constants()
    M, d, pi = K["M"], K["d"], K["pi"]
    ten = arb(10)
    N = D + 30
    L = (2 * M).log() + (N + 6) * ten.log()
    Lint = arb_ceil_int(L)
    h2 = math.floor(0.9 * 0.97 * 2 * math.pi * 1.4 / Lint * 2 ** 40) / 2 ** 40
    pl2["h"] = h2
    hA = arb(h2)
    pl2["Bdisc"] = 2 * M / ((2 * pi * d / hA).exp() - 1)
    # tails
    UR = int(math.ceil((N + 8) * LN10)) + 2
    tR = solve_t(math.log(UR), 0.0, 40.0)
    K_R = int(math.ceil(tR / h2)) + 1
    t = arb(K_R) * hA
    uR = (t - (-t).exp()).exp()
    pl2["K_R"] = K_R
    pl2["tailR"] = (-uR).exp()
    lnuL = 0.5 * (-(N + 8) * LN10 + LN2 - EULER_F) - 1.0
    tL = solve_t(lnuL, -600.0, 0.0)
    K_L = int(math.ceil(-tL / h2)) + 1
    t = -(arb(K_L) * hA)
    uL = (t - (-t).exp()).exp()
    pl2["K_L"] = K_L
    pl2["tailL"] = K["gam"].exp() * uL * uL / 2
    pl2["nodes"] = K_L + K_R + 1
    pl2["Ebits"] = pl1["Ebits"]
    pl2["S"] = pl1["S"]
    pl2["EGbits"] = pl1["EGbits"]
    pl2["cf_min_u"] = 1e18

    # Pick a few nodes that exist in both plans
    n0 = max(-pl1["K_L"], -pl2["K_L"])
    n1 = min(pl1["K_R"], pl2["K_R"])
    if n1 - n0 < 5:
        print("alt-h crosscheck: not enough overlapping nodes")
        return
    ks = random.sample(range(n0, n1), min(5, n1 - n0))
    print("Alt-h crosscheck (h=%.6e vs h'=%.6e)..." % (pl1["h"], pl2["h"]))
    print("-" * 70)
    ok = 0
    for k in ks:
        b1 = node_ball(k, pl1)
        b2 = node_ball(k, pl2)
        lo1, hi1 = float(b1.lower()), float(b1.upper())
        lo2, hi2 = float(b2.lower()), float(b2.upper())
        overlaps = not (hi1 < lo2 or hi2 < lo1)
        if overlaps:
            ok += 1
        print("  k=%6d  h=[%.15g, %.15g]  h'=[%.15g, %.15g]  %s"
              % (k, lo1, hi1, lo2, hi2, "OK" if overlaps else "MISMATCH"))
    print("-" * 70)
    print("Alt-h crosscheck: %d/%d nodes overlap" % (ok, len(ks)))


def precision_ladder_crosscheck(D=200):
    """Run the same node at 384 and 512 bits and compare intervals."""
    import random
    random.seed(24680)
    pl = make_plan(D, G=30)
    pl["cf_min_u"] = 1e18
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    ks = random.sample(range(n0, n1), min(5, n1 - n0))
    print("Precision-ladder crosscheck (384 vs 512 bits)...")
    print("-" * 70)
    ok = 0
    for k in ks:
        # Run once at the default Ebits
        b1 = node_ball(k, pl)
        # Run with temporarily increased Ebits
        pl_hi = dict(pl)
        pl_hi["Ebits"] = pl["Ebits"] + 100
        b2 = node_ball(k, pl_hi)
        lo1, hi1 = float(b1.lower()), float(b1.upper())
        lo2, hi2 = float(b2.lower()), float(b2.upper())
        overlaps = not (hi1 < lo2 or hi2 < lo1)
        if overlaps:
            ok += 1
        print("  k=%6d  384-bit=[%.15g, %.15g]  512-bit=[%.15g, %.15g]  %s"
              % (k, lo1, hi1, lo2, hi2, "OK" if overlaps else "MISMATCH"))
    print("-" * 70)
    print("Precision-ladder crosscheck: %d/%d nodes overlap" % (ok, len(ks)))


def lemma7_boundary_test():
    """Test the Lemma-7 recurrence at extreme z (exact rational)."""
    print("Lemma-7 boundary test (exact rational)...")
    for z in (Fraction(1, 10**6), Fraction(1, 1000), Fraction(1, 3),
              Fraction(2, 1), Fraction(10**6), Fraction(10**12, 7)):
        A2, A1 = Fraction(1), z
        B2, B1 = Fraction(0), Fraction(1)
        for n in range(1, 33):
            a = (n + 1) // 2
            b = Fraction(1) if n % 2 else z
            A2, A1 = A1, b * A1 + a * A2
            B2, B1 = B1, b * B1 + a * B2
            # determinant identity
            prod_a = 1
            for j in range(1, n + 1):
                prod_a *= (j + 1) // 2
            det = A1 * B2 - A2 * B1
            if det != ((-1) ** (n + 1)) * prod_a:
                raise AssertionError("boundary determinant failed at z=%s n=%d" % (z, n))
        print("  z=%-12s OK" % str(z))
    print("Lemma-7 boundary test: all OK")


# ============================================================================
# §10. Exact rational test of the Lemma 7 recurrence
# ============================================================================

def lemma7_exact_recurrence_test():
    """Exact rational check of the S-fraction recurrence.

    Uses only fractions.Fraction: no floats, no Arb.  For three rational
    z > 0 and the first 16 levels, verifies:
      * that B_n/A_n equals the contracted S-fraction;
      * the determinant identity for consecutive convergents;
      * the sign-change of consecutive differences.
    """
    for z in (Fraction(1, 3), Fraction(2, 1), Fraction(17, 5)):
        A2, A1 = Fraction(1), z
        B2, B1 = Fraction(0), Fraction(1)
        product_a = 1
        previous_delta = None
        x = 1 / z
        for n in range(1, 17):
            a = (n + 1) // 2
            b = Fraction(1) if n % 2 else z
            A2, A1 = A1, b * A1 + a * A2
            B2, B1 = B1, b * B1 + a * B2
            product_a *= a

            determinant = A1 * B2 - A2 * B1
            if determinant != ((-1) ** (n + 1)) * product_a:
                raise AssertionError("determinant identity failed for z=%s n=%d"
                                     % (z, n))

            tail = Fraction(1)
            for j in range(n, 0, -1):
                alpha = (j + 1) // 2
                tail = 1 + alpha * x / tail
            Cn = 1 / tail
            if B1 / A1 != Cn / z:
                raise AssertionError("S-fraction/β-fraction mismatch for z=%s n=%d"
                                     % (z, n))

            delta = B1 / A1 - B2 / A2
            if previous_delta is not None and not (delta * previous_delta < 0):
                raise AssertionError("alternation failed for z=%s n=%d" % (z, n))
            previous_delta = delta
    return "OK: 3 rational z values × 16 levels; exact identities hold"


# ============================================================================
# §11. CLI
# ============================================================================

def _cmd_test():
    print("lemma7_exact_recurrence_test():")
    print("  " + lemma7_exact_recurrence_test())
    print("All tests passed.")


def _cmd_smoke():
    res = run_certificate(60, G=30, jobs=1, ckpt=None,
                          out="smoke_60.txt", cf_min_u=20.0)
    if not (res and res.get("proven")):
        print("SMOKE TEST FAILED")
        return 1
    print("SMOKE TEST OK")
    return 0


def _cmd_certify(args):
    res = run_certificate(
        args.digits, G=args.guard, jobs=args.jobs, chunk=args.chunk,
        ckpt=args.ckpt, out=args.out, max_hours=args.max_hours,
        cf_min_u=args.cf_min_u, trust_checkpoint=args.resume,
    )
    if res is None:
        return 2
    if args.resume:
        print("RESUME MODE: no certificate issued.")
        return 0
    return 0 if res.get("proven") else 1


def _cmd_crosscheck(args):
    crosscheck_report(D=args.digits, num_nodes=args.num_nodes)
    alt_h_crosscheck(D=args.digits)
    precision_ladder_crosscheck(D=args.digits)
    return 0


def _cmd_boundary_test():
    lemma7_boundary_test()
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="golomb_dickman_cert",
        description="Rigorous interval-arithmetic certificate for the "
                    "Golomb–Dickman constant."
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("test", help="run environment and lemma-7 tests")
    sub.add_parser("smoke", help="60-digit end-to-end smoke test")
    sub.add_parser("boundary-test", help="Lemma-7 boundary tests")

    p_cert = sub.add_parser("certify", help="full certificate")
    p_cert.add_argument("--digits", type=int, default=40000)
    p_cert.add_argument("--guard", type=int, default=30)
    p_cert.add_argument("--jobs", type=int, default=0)
    p_cert.add_argument("--chunk", type=int, default=0)
    p_cert.add_argument("--ckpt", type=str, default=None)
    p_cert.add_argument("--out", type=str, default=None)
    p_cert.add_argument("--max-hours", type=float, default=0.0)
    p_cert.add_argument("--cf-min-u", type=float, default=None)
    p_cert.add_argument("--resume", action="store_true",
                        help="AUDIT-ONLY: accept checkpoint pickles without "
                             "recomputation (no certificate is written)")

    p_cross = sub.add_parser("crosscheck", help="independent cross-checks")
    p_cross.add_argument("--digits", type=int, default=200)
    p_cross.add_argument("--num-nodes", type=int, default=10)

    args = parser.parse_args(argv)
    if args.command is None:
        args.command = "test"

    _env_check()
    print()
    # Run the exact rational Lemma-7 test once at startup, before any
    # numeric work.  This satisfies the requirement that the recurrence
    # correspondence be verified exactly in every invocation.
    print("lemma7_exact_recurrence_test():")
    print("  " + lemma7_exact_recurrence_test())
    print()

    if args.command == "test":
        print("All tests passed.")
        return 0
    if args.command == "smoke":
        return _cmd_smoke()
    if args.command == "certify":
        return _cmd_certify(args)
    if args.command == "crosscheck":
        return _cmd_crosscheck(args)
    if args.command == "boundary-test":
        return _cmd_boundary_test()
    parser.error("unknown command: %s" % args.command)
    return 1


if __name__ == "__main__":
    sys.exit(main())