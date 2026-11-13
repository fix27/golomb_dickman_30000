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
(§1), the implementation (§2–§8), and the CLI entry point (§9).  The
external trust base is limited to

  (i)   correctness of python-flint / Arb,
  (ii)  the trapezoidal theorem of Trefethen–Weideman,
  (iii) the classical Stieltjes theorem on positive S-fractions
        (convergence and alternation),
  (iv)  the classical Stieltjes expansion of a determinate moment
        problem into a J-fraction (Jacobi-matrix resolvent),
  (v)   exact integer arithmetic of CPython.

The S→J contraction used in Lemma 7, step 3 is NOT part of the trust
base: it is proved from first principles in Lemma 7, step 3 below.

Independent cross-checks.  In addition to the primary Arb computation,
the program provides several INDEPENDENT checks that STRENGTHEN
confidence without replacing the trust in Arb.  Every cross-check raises
ArithmeticError on mismatch; the CLI aggregates statuses and exits
non-zero on any failure.

  * exact rational tests of Lemma 7 (Fraction arithmetic): the
    recurrence, the determinant identity, alternation of consecutive
    convergents, the S→J contraction at sampled x, and the S→J
    contraction as polynomial identities.  These are run automatically
    at startup.
  * mpmath (non-interval) recomputation of selected nodes.
  * mpmath.iv (interval arithmetic) recomputation of selected nodes.
  * direct numerical check of cf_f_hull against mpmath.e1 on a grid of u.
  * forced-branch exercise (series / cf / asym) on appropriate nodes.
  * full-quadrature alt-h comparison — computes the FULL λ interval twice
    (once with step h, once with 0.9·h) and checks overlap.
  * precision ladder — runs the same node at two precisions and compares
    Arb intervals directly (no float conversion).
  * boundary tests for Lemma 7 at extreme z.

The final certificate is issued only if the primary Arb run succeeds
with rad < eps in every node AND the certified digits agree with the
known prefix of OEIS A084945 (up to min(D, len(KNOWN)) digits).

Parity conventions (defensive).  The tail bounds of Lemma 1(c,d) have a
sign that depends on the truncation index.  The code uses a symmetric
ball [−rem, +rem] for the remainder (arb(0, rem) = {x : |x| ≤ rem}), so
the sign is irrelevant for inclusion.  In addition, and as a defensive
measure against any future change of convention, the code forces the
truncation index into the "safe" parity:

  * in the power-series branch the naive sum exits only on odd k, and
    the rectangular splitting uses an even block size (asserted);
  * in the asymptotic branch the index K is forced to be even.

Under either convention (symmetric [−rem, +rem] or one-sided [0, rem])
the returned ball contains the true value.

Version note (v18).  Relative to v17:

  1. Lemma 7, step 3 is now a complete self-contained proof of the S→J
     contraction, with explicit base cases (m = 0 and m = 1) and a
     rigorous induction step that derives
         P_{2m+1} = (1 + b_m x) P_{2m-1} − a_m x² P_{2m-3}
     from the S-recurrence.  The two failed attempts and the "which is
     wrong / hmm / doesn't match" phrases from v17 are removed.  Wall's
     theorem is now used ONLY for the convergence of the S-fraction (as
     part of Stieltjes' theorem), not for the contraction itself.  The
     contraction is therefore removed from the external trust base.

  2. Lemma 7, step 5 is rewritten: the B_n derivation is now explicit
     for both parities, matching the A_n derivation.  The garbled
     intermediate expression in v17 is replaced.

  3. A new startup test, contraction_symbolic_test, verifies the S→J
     contraction at the level of polynomial identities in x (not just
     at specific numerical x).  It checks P_{2m+1} = R_m and Q_{2m+1} = S_m
     for m = 0..6.

  4. run_certificate compares the OEIS prefix at length min(D, len(KNOWN))
     rather than the previous min(D, 95).  Since KNOWN has 100 digits,
     the audit is now stronger by 5 digits.

  All other behaviour is unchanged from v17.

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
integral is ≤ 2B₂ = M₂.

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

For G the hypotheses hold with a = d, hence the discretisation error
is bounded by B_disc = 2M/(e^{2πd/h} − 1).

Lemma 6 (strict interval balls for nodes).
Let t = kh, q = e^{−t}, u = e^{t−q} > 0 and

    G(t) = F(u) · u · (1+q),    F(u) = e^γ u exp(−2u + S(u)),
    S(u) = Σ_{j=2}^∞ (−u)^j / (j·j!).

The function node_ball must return an interval I_k satisfying
simultaneously G(kh) ∈ I_k and rad(I_k) < 2^{−E_bits}.

Proof of correctness of the three branches.

1. Power series.  With a_j = (−u)^j / (j·j!):
       |a_{j+1}/a_j| = u·j/(j+1)² ≤ u/(j+1).
   If K + 1 ≥ 2u, the ratio is at most 1/2 for all j ≥ K + 1, hence
       |Σ_{j>K} a_j| ≤ |a_{K+1}| Σ_{r≥0} 2^{−r} = 2|a_{K+1}|.
   The finite sum up to K is computed in interval arithmetic; the tail
   is enclosed in the symmetric interval [−R_K, R_K] with R_K ≥ 2|a_{K+1}|.
   The code additionally forces the exit index to be odd so that even
   under the one-sided convention [0, R_K] the interval still contains
   the true value.  Rectangular splitting merely groups the finite
   terms; the block size is asserted to be even, which keeps
   Kend = 1 + nblocks·m odd (safe parity).

2. Continued fraction.  Lemma 7 shows that e^u E₁(u) lies between the
   consecutive convergents f_{n−1}, f_n.  If D is an interval upper
   bound for |f_n − f_{n−1}|, then the ball f_n + [−D, D] contains both
   endpoints, hence e^u E₁(u).

3. Asymptotics.  From the exact integral representation
       e^u E₁(u) = Σ_{j<K} (−1)^j j! / u^{j+1} + R_K,   |R_K| ≤ K!/u^{K+1},
   where R_K = (−1)^K ∫₀^∞ e^{−s} s^K / (u^K (u+s)) ds has sign (−1)^K.
   The code encloses the remainder in the symmetric interval [−rem, rem]
   with rem = K!/u^{K+1}; in addition, K is forced to be even so that
   under the one-sided convention [0, rem] the interval also contains
   R_K (positive).

4. Negligible node.  From 0 < F(u) < min(e^γ u, e^{−u}):
       0 < G(t) < e^γ u²(1+q)    if t < 0,
       0 < G(t) < u e^{−u}(1+q)  if t ≥ 0.
   Therefore the ball [0, U] is correct provided Arb has confirmed
   U < ε.  This branch (skip) must verify BOTH the inclusion and the
   postcondition rad < ε before returning.

Postcondition check.  In every computational branch the returned ball
is validated by is_finite() and rad() < eps.  Float values in _est,
method_of, asym_K and the level count only influence the route choice;
all operations that establish inclusion and tolerance are performed in
Arb.  If the heuristic picks an insufficient number of terms/levels,
the algorithm must escalate or raise an error — never return an
unverified ball.  ∎

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
and Σ_j (m_{2j})^{−1/(2j)} diverges.

Step 2. J-fraction from the Laguerre Jacobi matrix.  For the measure
e^{−s} ds on [0, ∞), the monic Laguerre polynomials p_n(s) = (−1)^n n! L_n(s)
satisfy
    p_{n+1}(s) = (s − (2n+1)) p_n(s) − n² p_{n−1}(s),
with p₀ = 1, p₁ = s − 1 and norms ∫₀^∞ p_n(s)² e^{−s} ds = (n!)².
The Jacobi matrix has diagonal b_n = 2n+1 = 1, 3, 5, … and off-diagonal
squares a_n = n² (n ≥ 1).  The Stieltjes transform therefore has the
formal J-fraction expansion (standard theory of the Jacobi-matrix
resolvent for a determinate moment problem):
    g(x) = 1/(1+x − 1²x²/(1+3x − 2²x²/(1+5x − 3²x²/(1+7x − …)))).

Step 3. S-fraction expansion and S→J contraction.  Every positive-
definite determinate Stieltjes transform has an S-fraction expansion
    g(x) = 1/(1 + α_1 x/(1 + α_2 x/(1 + α_3 x/(…)))),
with α_n > 0 (Wall 1948, §51).  We identify the coefficients via the
following self-contained lemma.

Lemma (S→J contraction).  Let α_n > 0 for n ≥ 1, and define
    b_0 = α_1,   b_n = α_{2n} + α_{2n+1},   a_n = α_{2n-1} α_{2n}   (n ≥ 1).
Let C_n = P_n/Q_n be the n-th convergent of the S-fraction, defined by
    P_{−1} = 0, P_0 = 1;   Q_{−1} = 1, Q_0 = 1;
    P_n = P_{n−1} + α_n x P_{n−2},   Q_n = Q_{n−1} + α_n x Q_{n−2}   (n ≥ 1).
Let J^{(m)} = R_m/S_m be the m-th approximant of the J-fraction, defined by
    R_{−1} = 0, R_0 = 1;   S_{−1} = 1, S_0 = 1 + b_0 x;
    R_m = (1 + b_m x) R_{m−1} − a_m x² R_{m−2}   (m ≥ 1);
    S_m = (1 + b_m x) S_{m−1} − a_m x² S_{m−2}   (m ≥ 1).
Then for every m ≥ 0, P_{2m+1} = R_m and Q_{2m+1} = S_m as polynomials
in x.

Proof of the lemma.

Base cases.
  m = 0:  P_1 = P_0 + α_1 x P_{−1} = 1 − 0 = 1 = R_0,
          Q_1 = Q_0 + α_1 x Q_{−1} = 1 + α_1 x = 1 + b_0 x = S_0.  ✓
  m = 1:  Compute P_2 = P_1 + α_2 x P_0 = 1 + α_2 x, hence
          P_3 = P_2 + α_3 x P_1 = 1 + (α_2 + α_3) x = 1 + b_1 x,
          and R_1 = (1 + b_1 x) R_0 − a_1 x² R_{−1} = 1 + b_1 x.  So P_3 = R_1.  ✓
          Compute Q_2 = Q_1 + α_2 x Q_0 = 1 + (α_1 + α_2) x, hence
          Q_3 = Q_2 + α_3 x Q_1 = 1 + (α_1 + α_2 + α_3) x + α_1 α_3 x².
          S_1 = (1 + b_1 x) S_0 − a_1 x² S_{−1}
              = (1 + b_1 x)(1 + b_0 x) − a_1 x²
              = 1 + (b_0 + b_1) x + (b_0 b_1 − a_1) x².
          With b_0 = α_1, b_1 = α_2 + α_3, a_1 = α_1 α_2:
              b_0 + b_1 = α_1 + α_2 + α_3,
              b_0 b_1 − a_1 = α_1(α_2 + α_3) − α_1 α_2 = α_1 α_3.
          Hence S_1 = Q_3.  ✓

Induction step for P.  Assume P_{2m−1} = R_{m−1} and P_{2m−3} = R_{m−2},
with the convention P_{−1} = R_{−1} = 0 for m = 1.  From the S-recurrence:
    P_{2m+1} = P_{2m} + α_{2m+1} x P_{2m−1}            (a)
    P_{2m}   = P_{2m−1} + α_{2m} x P_{2m−2}            (b)
    P_{2m−1} = P_{2m−2} + α_{2m−1} x P_{2m−3}          (c)
From (c):  P_{2m−2} = P_{2m−1} − α_{2m−1} x P_{2m−3}.
Substituting into (b):
    P_{2m} = P_{2m−1} + α_{2m} x [P_{2m−1} − α_{2m−1} x P_{2m−3}]
          = (1 + α_{2m} x) P_{2m−1} − α_{2m−1} α_{2m} x² P_{2m−3}.
Substituting this into (a):
    P_{2m+1} = (1 + α_{2m} x) P_{2m−1} − α_{2m−1} α_{2m} x² P_{2m−3}
             + α_{2m+1} x P_{2m−1}
            = (1 + (α_{2m} + α_{2m+1}) x) P_{2m−1} − (α_{2m−1} α_{2m}) x² P_{2m−3}
            = (1 + b_m x) P_{2m−1} − a_m x² P_{2m−3}.
By induction, P_{2m−1} = R_{m−1} and P_{2m−3} = R_{m−2}, hence
    P_{2m+1} = (1 + b_m x) R_{m−1} − a_m x² R_{m−2} = R_m.  ✓
Identical algebra gives Q_{2m+1} = S_m.  ∎

Application to α_n = ⌈n/2⌉.  The contraction formulas give
    b_0 = ⌈1/2⌉ = 1,
    b_n = ⌈2n/2⌉ + ⌈(2n+1)/2⌉ = n + (n+1) = 2n+1,
    a_n = ⌈(2n−1)/2⌉ · ⌈2n/2⌉ = n · n = n².
These are precisely the coefficients of the Laguerre J-fraction from
Step 2.  Hence P_{2m+1}/Q_{2m+1} = J^{(m)} for every m.

Passage to the limit.  Since α_n > 0, the S-fraction converges
(Stieltjes; Wall 1948, Theorem 51.1), and its even and odd subsequences
converge to the same limit.  The J-fraction converges to the Stieltjes
transform g(x) (standard theory of the Jacobi-matrix resolvent).  Since
C_{2m+1} = J^{(m)}, the S-fraction converges to g(x).  Thus
    g(x) = 1/(1 + x/(1 + x/(1 + 2x/(1 + 2x/(1 + 3x/(1 + 3x/(…))))))),
with α_n = ⌈n/2⌉.

Step 4. Alternation of convergents.  For a positive S-fraction of
Stieltjes type with α_j > 0, the Stieltjes convergence-and-alternation
theorem states: for x > 0, the consecutive convergents C_{n−1}(x), C_n(x)
bracket the value of the Stieltjes transform.  Therefore
    g(x) ∈ [min(C_{n−1}(x), C_n(x)), max(C_{n−1}(x), C_n(x))].

Step 5. Correspondence between the S-fraction convergents and the
code's recurrence.  Let C_n(x) = P_n(x)/Q_n(x) be the classical
convergents of the S-fraction with α_n = ⌈n/2⌉, satisfying
    P_{−1} = 0, P_0 = 1;   Q_{−1} = 1, Q_0 = 1;
    P_n = P_{n−1} + α_n x P_{n−2},   Q_n = Q_{n−1} + α_n x Q_{n−2}.

Claim.  For every n ≥ −1:
    A_n(z) = z^{⌈(n+1)/2⌉} Q_n(1/z),      B_n(z) = z^{⌈(n−1)/2⌉} P_n(1/z),
with the convention z^{−1} · 0 = 0 for n = −1.

Proof by induction on n.  Write x = 1/z throughout.

Base cases.
  n = −1:  ⌈(n+1)/2⌉ = ⌈0/2⌉ = 0, ⌈(n−1)/2⌉ = ⌈−1⌉ = −1.  Thus
      A_{−1} = z⁰ Q_{−1}(x) = 1,
      B_{−1} = z^{−1} P_{−1}(x) = z^{−1} · 0 = 0,
  matching the code's A_{−1} = 1, B_{−1} = 0.  ✓
  n = 0:  ⌈(n+1)/2⌉ = ⌈1/2⌉ = 1, ⌈(n−1)/2⌉ = ⌈−1/2⌉ = 0.  Thus
      A_0 = z¹ Q_0(x) = z,
      B_0 = z⁰ P_0(x) = 1,
  matching the code's A_0 = z, B_0 = 1.  ✓

Induction step for A_n.  Two parities.

Case n = 2m even.  Then ⌈(n+1)/2⌉ = m+1, ⌈n/2⌉ = m, ⌈(n−1)/2⌉ = m.
The code's recurrence is A_n = z A_{n−1} + a_n A_{n−2} with a_n = ⌈n/2⌉ = m.
By induction hypothesis A_{n−1} = z^m Q_{n−1}(x), A_{n−2} = z^m Q_{n−2}(x).
Hence
    A_n = z · z^m Q_{n−1}(x) + m · z^m Q_{n−2}(x)
        = z^m [z Q_{n−1}(x) + m Q_{n−2}(x)]
        = z^m (1/x) [Q_{n−1}(x) + m x Q_{n−2}(x)]        (since z = 1/x)
        = z^m (1/x) Q_n(x)                                (by S-recurrence)
        = z^{m+1} Q_n(x)
        = z^{⌈(n+1)/2⌉} Q_n(1/z).  ✓

Case n = 2m+1 odd.  Then ⌈(n+1)/2⌉ = m+1, ⌈n/2⌉ = m+1, ⌈(n−1)/2⌉ = m.
The code's recurrence is A_n = 1 · A_{n−1} + a_n A_{n−2} with a_n = m+1.
By induction A_{n−1} = z^{m+1} Q_{n−1}(x), A_{n−2} = z^m Q_{n−2}(x).  Hence
    A_n = z^{m+1} Q_{n−1}(x) + (m+1) z^m Q_{n−2}(x)
        = z^m [z Q_{n−1}(x) + (m+1) Q_{n−2}(x)]
        = z^m (1/x) [Q_{n−1}(x) + (m+1) x Q_{n−2}(x)]
        = z^m (1/x) Q_n(x)
        = z^{m+1} Q_n(x) = z^{⌈(n+1)/2⌉} Q_n(1/z).  ✓

Induction step for B_n.  Two parities.

Case n = 2m even.  Then ⌈(n−1)/2⌉ = m.  The induction hypothesis applied
to indices n−1 and n−2 gives
    B_{n−1} = z^{⌈(n−2)/2⌉} P_{n−1}(x) = z^{m−1} P_{n−1}(x),
    B_{n−2} = z^{⌈(n−3)/2⌉} P_{n−2}(x) = z^{m−1} P_{n−2}(x).
The code's recurrence is B_n = z B_{n−1} + a_n B_{n−2} with a_n = m.  Hence
    B_n = z · z^{m−1} P_{n−1}(x) + m · z^{m−1} P_{n−2}(x)
        = z^{m−1} [z P_{n−1}(x) + m P_{n−2}(x)]
        = z^{m−1} (1/x) [P_{n−1}(x) + m x P_{n−2}(x)]
        = z^{m−1} (1/x) P_n(x)
        = z^m P_n(x)
        = z^{⌈(n−1)/2⌉} P_n(1/z).  ✓

Case n = 2m+1 odd.  Then ⌈(n−1)/2⌉ = m.  The induction hypothesis gives
    B_{n−1} = z^{⌈(n−2)/2⌉} P_{n−1}(x) = z^m P_{n−1}(x),
    B_{n−2} = z^{⌈(n−3)/2⌉} P_{n−2}(x) = z^{m−1} P_{n−2}(x).
The code's recurrence is B_n = 1 · B_{n−1} + a_n B_{n−2} with a_n = m+1.
Hence
    B_n = z^m P_{n−1}(x) + (m+1) z^{m−1} P_{n−2}(x)
        = z^{m−1} [z P_{n−1}(x) + (m+1) P_{n−2}(x)]
        = z^{m−1} (1/x) [P_{n−1}(x) + (m+1) x P_{n−2}(x)]
        = z^{m−1} (1/x) P_n(x)
        = z^m P_n(x)
        = z^{⌈(n−1)/2⌉} P_n(1/z).  ✓

Therefore A_n = z^{⌈(n+1)/2⌉} Q_n(1/z) and B_n = z^{⌈(n−1)/2⌉} P_n(1/z) for
every n ≥ −1.  ∎

Step 6. Back to f_n.  Since ⌈(n−1)/2⌉ − ⌈(n+1)/2⌉ = −1 for every integer n,
    f_n(z) = B_n(z)/A_n(z) = z^{−1} P_n(1/z) / Q_n(1/z) = z^{−1} C_n(1/z).
Since z^{−1} > 0 for z > 0, multiplication preserves order and
    e^z E₁(z) = z^{−1} g(1/z) ∈ [min(f_{n−1}, f_n), max(f_{n−1}, f_n)].
∎

Algebraic check of consecutive convergents.  From the recurrence,
    A_n B_{n−1} − A_{n−1} B_n = (−1)^{n+1} Π_{j=1}^n a_j.
In particular, the sign of the difference of consecutive convergents
alternates at each step; the bracketing of the limit function is provided
by the Stieltjes theorem in step 4, not by the determinant alone.

Theorem (assembly).  Let S_fin = h Σ_{k=−K_L}^{K_R} G(kh).  Then
    λ ∈ S_fin + [0, τ_R + τ_L] + [−B_disc, B_disc],
    τ_R = e^{−u(K_R h)},   τ_L = ½ e^γ u(−K_L h)².
Proof. Lemma 2 gives λ = ∫ G; Lemma 5: λ = h Σ_{k∈ℤ} G(kh) + δ with
|δ| ≤ B_disc; Lemma 3: discarded tails lie in [0, τ_R] and [0, τ_L].  ∎

----------------------------------------------------------------------------
§2. TRUST BASE AND ADDITIONAL CROSS-CHECKS
----------------------------------------------------------------------------

The certificate is conditional on:

1. correctness of the interval arithmetic library Arb (via python-flint);
2. the trapezoidal theorem (Trefethen & Weideman, SIAM Review 56 (2014),
   Theorem 5.1);
3. the classical Stieltjes theorem on positive S-fractions, for
   convergence and alternation of the S-fraction;
4. the classical Stieltjes expansion of a determinate moment problem
   into a J-fraction (Jacobi-matrix resolvent);
5. the lemmas in §1;
6. exact integer arithmetic of CPython.

NOTE: the S→J contraction used in Lemma 7, step 3 is NOT part of the
trust base.  It is proved from first principles (Lemma 7, step 3).  The
Wall contraction theorem is not cited as a black box.

The cross-checks (§9) are safety nets, not part of the trust base.  Each
either returns True or raises ArithmeticError on mismatch; the CLI
aggregates the statuses and exits non-zero on any failure.

----------------------------------------------------------------------------
§3. USAGE
----------------------------------------------------------------------------

    python golomb_dickman_cert.py test
        Environment self-test and exact rational Lemma-7 checks
        (recurrence, contraction at sampled x, contraction as
        polynomial identities, boundary z).  Fast.

    python golomb_dickman_cert.py smoke
        60-digit end-to-end smoke test.

    python golomb_dickman_cert.py certify --digits D [--guard G]
        Full certificate.  Writes golomb_dickman_<D>_certified.txt.
        Raises ArithmeticError on (a) known-digit mismatch, or
        (b) requested digits not proven.

    python golomb_dickman_cert.py certify --digits D --resume
        AUDIT-ONLY: accept checkpoint pickles without recomputing them.
        The result is NOT a proof; no certificate file is written.

    python golomb_dickman_cert.py certify --digits D --no-resume
        Safe mode (the default): recompute every chunk in Arb.

    python golomb_dickman_cert.py crosscheck --digits D
        Independent mpmath / mpmath.iv / cf / forced-branch /
        full-quadrature alt-h / precision-ladder cross-checks.
        Requires mpmath.  Exits non-zero if any check fails.

    python golomb_dickman_cert.py boundary-test
        Lemma-7 boundary tests (exact rational).

Requirements: python-flint and (for cross-checks) mpmath.  No pip
install is attempted by this module.  See requirements.txt.
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
from contextlib import contextmanager
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(500_000)

try:
    import flint
except ImportError:
    print("ERROR: python-flint is required.", file=sys.stderr)
    print("Install it with: pip install python-flint", file=sys.stderr)
    print("(see requirements.txt for the pinned version)", file=sys.stderr)
    raise
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
)  # OEIS A084945, first 100 digits after the decimal point

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
    chk("arb(mid, rad) is symmetric",
        lambda: arb(0, 1).upper() >= 1 and arb(0, 1).lower() <= -1)
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


def _heuristic_cf_min_u(D):
    """Pre-measured heuristic threshold for cf_min_u."""
    if D < 500:
        return 3000.0
    elif D < 5000:
        return 8000.0
    else:
        return 13000.0


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


def series_S_naive(u, uf, thr):
    """S = Σ_{k=2}^K (−u)^k/(k·k!) and a tail bound rem < thr (strict).

    Exit is only permitted on odd K, so the first discarded term
    (−u)^{K+1} is positive and the one-sided enclosure [0, rem] contains
    the true tail.
    """
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
        if k >= kmin and (k < 64 or k % 16 == 1):
            rem = 2 * abs(p).upper() * ub / ((k + 1) * (k + 1))
            if rem < thr:
                if k % 2 == 1:
                    return S, rem
                else:
                    continue


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
    if m % 2 != 0:
        raise ArithmeticError("series_S: block size m must be even, got %d" % m)
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


class _Escalate(Exception):
    pass


def cf_f_hull(u, uf, gap_bits, nmax):
    """Ball for f = e^u E₁(u), u > 0."""
    thr = arb(1) / arb(1 << gap_bits)
    A2, A1 = arb(1), u
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


def asym_K(xf, uf, lng, L_eps):
    top = min(ASYM_KMAX, int(uf) - 1)
    lnfact = 0.0
    for K in range(0, top + 1):
        if K > 0:
            lnfact += math.log(K)
        if lng + lnfact - (K + 1) * xf - uf <= -L_eps - 3.0:
            return K
    return None


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
    """Strict ball for G(kh) with rad < 2^{-Ebits}."""
    Eb = pl["Ebits"]
    xf, uf, lng = _est(k, pl)
    L_eps = Eb * LN2
    eps = arb(1) / arb(1 << Eb)
    need_ln = L_eps + lng + 2.0

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
        ctx.prec = max(ctx.prec, 128)

    P = max(64, int(need_ln / LN2) + int(math.log2(uf + 2.0)) + 64)
    thr_bits = int(need_ln / LN2) + 4
    gap_bits = max(32, int((L_eps + lng - uf + 4.0) / LN2) + 2) if xf >= 0 else 0
    meth = force or method_of(k, pl)
    Ka = asym_K(xf, uf, lng, L_eps) if (xf >= 0 and meth == "asym") else None
    if Ka is not None and (Ka % 2):
        Ka += 1
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
    """Benchmark cf_min_u on this machine.  Slow; only for tuning."""
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
    """Fingerprint to detect accidental mixing of different runs."""
    try:
        with open(__file__, "rb") as f:
            src_sha = hashlib.sha256(f.read()).hexdigest()
    except Exception:
        src_sha = "unavailable"
    return {
        "format": "gd-checkpoint-v8",
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
    """Compute chunks.  Bounded memory: at most 2·jobs in flight."""
    pw = {k: pl[k] for k in ("D", "G", "h", "Ebits", "S", "EGbits",
                             "K_L", "K_R", "cf_min_u")}
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    bounds = [(lo, min(lo + chunk, n1)) for lo in range(n0, n1, chunk)]
    identity = _checkpoint_identity(pl, chunk)
    done = {}

    if ckpt:
        os.makedirs(ckpt, exist_ok=True)
        mp = os.path.join(ckpt, "meta_v8.json")
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
        max_in_flight = max(2, jobs * 2)
        todo_iter = iter(todo)
        pending = set()

        def submit_next():
            try:
                lo, hi = next(todo_iter)
            except StopIteration:
                return False
            pending.add(ex.submit(run_chunk, lo, hi, pw))
            return True

        with ProcessPoolExecutor(max_workers=jobs,
                                 mp_context=_get_mp_context()) as ex:
            while len(pending) < max_in_flight and submit_next():
                pass
            while pending:
                tmo = None if not deadline else max(1.0, min(30.0, deadline - time.time()))
                fin_set, pending = wait(pending, timeout=tmo,
                                        return_when=FIRST_COMPLETED)
                for f in fin_set:
                    accept(f.result())
                while len(pending) < max_in_flight:
                    if not submit_next():
                        break
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
                    max_hours=0.0, cf_min_u=None, trust_checkpoint=False,
                    autotune=False):
    """Run certificate pipeline.

    On successful completion with proven=True and matching OEIS prefix,
    writes the certificate file.  Raises ArithmeticError if (a) the
    known-digit cross-check fails, or (b) the digits are not proven.
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
        print("Pass --no-resume (the default) for a rigorous run.")
        print("=" * 70)

    pl = make_plan(D, G)
    print("D=%d  GUARD=%d  nodes=%d  h=%.6e  K_L=%d K_R=%d  processes=%d"
          % (D, G, pl["nodes"], pl["h"], pl["K_L"], pl["K_R"], jobs))
    print("  M ≈ 10^%d,  B_disc ≈ 10^%d"
          % (log10_arb(pl["M"]), log10_arb(pl["Bdisc"])))

    if cf_min_u is None:
        if autotune:
            print("  autotune: benchmarking cf_min_u (this may take minutes)...")
            cf_min_u = autotune_cf_min_u(pl)
        else:
            cf_min_u = _heuristic_cf_min_u(D)
            print("  using heuristic cf_min_u = %.0f "
                  "(pass --cf-min-u for a measured value, or --autotune to benchmark)"
                  % cf_min_u)
    pl["cf_min_u"] = float(cf_min_u)
    print("  continued-fraction threshold: u >= %s"
          % ("∞ (series only)" if cf_min_u > 1e17 else "%.0f" % cf_min_u))

    if chunk <= 0:
        chunk = max(25, min(200, pl["nodes"] // (max(1, jobs) * 24)))
    deadline = time.time() + max_hours * 3600.0 if max_hours > 0 else None
    if deadline:
        print("  soft deadline: %.1f hours (running tasks may exceed it; "
              "use `timeout` for a hard stop)" % max_hours)
    chunks = compute_chunks(pl, jobs, chunk, ckpt, deadline,
                            trust_checkpoint=trust_checkpoint)
    if chunks is None:
        print("STOPPED by max_hours: checkpoints in", ckpt,
              "- re-run the cell to continue")
        return None

    res = certify(pl, chunks)
    # Fix v18: compare against the FULL known prefix (len(KNOWN)), not 95.
    k = min(D, len(KNOWN))
    ok_known = res["digits_lo"][:k] == KNOWN[:k]
    print("independent cross-check with OEIS A084945 (first %d digits): %s"
          % (k, "OK" if ok_known else "MISMATCH"))

    if trust_checkpoint:
        res["proven"] = False
        res["audit_only"] = True
        print("AUDIT ONLY: results were not verified; no certificate issued.")
        return res

    if not ok_known:
        raise ArithmeticError(
            "known-digit cross-check with OEIS A084945 failed; "
            "the computation is INVALID — do not use this result")
    if not res["proven"]:
        raise ArithmeticError(
            "requested %d digits were not proven; "
            "lower and upper bounds differ — increase GUARD" % D)

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
# §9. Independent cross-checks
# ============================================================================

def _mpmath_G(t):
    """Recompute G(t) with mpmath at high precision (non-interval)."""
    u = mpmath.e ** (t - mpmath.e ** (-t))
    E1 = mpmath.e1(u)
    return mpmath.e ** (-u - E1) * u * (1 + mpmath.e ** (-t))


def _exact_t(k, pl):
    """Return t = k·h as an exact Fraction."""
    return Fraction(*pl["h"].as_integer_ratio()) * k


def _mpf_of_fraction(frac):
    """Convert a Fraction to mpmath.mpf at current precision."""
    return mpmath.mpf(frac.numerator) / mpmath.mpf(frac.denominator)


def _relative_contains(ball, mp_val, rel_tol="1e-30"):
    """Check mp_val ∈ ball up to a relative tolerance (no float)."""
    if ball is None or mp_val is None:
        return False
    try:
        lo = mpmath.mpf(str(ball.lower()))
        hi = mpmath.mpf(str(ball.upper()))
    except Exception:
        return False
    tol = mpmath.mpf(rel_tol)
    if mp_val == 0:
        return lo <= tol and hi >= -tol
    scale = abs(mp_val)
    return lo - tol * scale <= mp_val <= hi + tol * scale


@contextmanager
def _mpmath_prec(dps, iv_dps=None):
    """Temporarily set global mpmath (and optionally mpmath.iv) precision."""
    old_dps = mpmath.mp.dps
    old_iv_dps = None
    if iv_dps is not None:
        from mpmath import iv
        old_iv_dps = iv.dps
        iv.dps = iv_dps
    mpmath.mp.dps = dps
    try:
        yield
    finally:
        mpmath.mp.dps = old_dps
        if old_iv_dps is not None:
            from mpmath import iv
            iv.dps = old_iv_dps


def mpmath_crosscheck_node(k, pl, dps=100):
    """Compare the Arb ball for node k with an mpmath value."""
    if mpmath is None:
        raise ArithmeticError("mpmath is required for mpmath_crosscheck_node")
    with _mpmath_prec(dps):
        ball = node_ball(k, pl)
        t = _mpf_of_fraction(_exact_t(k, pl))
        mp_val = _mpmath_G(t)
        within = _relative_contains(ball, mp_val, rel_tol="1e-30")
        return ball, mp_val, within


def mpmath_iv_crosscheck_node(k, pl, dps=100):
    """Compare the Arb ball for node k with an mpmath.iv interval.

    Overlap is checked at full working precision; no float conversion.
    iv.dps is set explicitly.
    """
    if mpmath is None:
        raise ArithmeticError("mpmath is required for mpmath_iv_crosscheck_node")
    with _mpmath_prec(dps, iv_dps=dps):
        from mpmath import iv
        ball = node_ball(k, pl)
        t_exact = _exact_t(k, pl)
        t_iv = iv.mpf(t_exact.numerator) / iv.mpf(t_exact.denominator)
        u = iv.exp(t_iv - iv.exp(-t_iv))
        E1 = iv.e1(u)
        G = iv.exp(-u - E1) * u * (1 + iv.exp(-t_iv))
        lo_arb = mpmath.mpf(str(ball.lower()))
        hi_arb = mpmath.mpf(str(ball.upper()))
        overlaps = not (hi_arb < G.a or G.b < lo_arb)
        return ball, (G.a, G.b), overlaps


def cf_numerical_crosscheck():
    """Verify cf_f_hull against mpmath.e1 on a grid of 5 u values."""
    if mpmath is None:
        raise ArithmeticError("mpmath is required for cf_numerical_crosscheck")
    old_dps = mpmath.mp.dps
    old_prec = ctx.prec
    try:
        mpmath.mp.dps = 50
        print("CF numerical cross-check vs mpmath.e1...")
        print("-" * 70)
        ok = 0
        total = 0
        failures = []
        for uf in (10.0, 100.0, 1000.0, 10000.0, 100000.0):
            u = arb(uf)
            ctx.prec = 500
            try:
                ball = cf_f_hull(u, uf, gap_bits=100, nmax=200000)
            except _Escalate:
                failures.append("u=%g: cf did not converge" % uf)
                continue
            total += 1
            mp_val = mpmath.e ** mpmath.mpf(uf) * mpmath.e1(mpmath.mpf(uf))
            within = _relative_contains(ball, mp_val, rel_tol="1e-20")
            if within:
                ok += 1
            else:
                failures.append("u=%g mismatch" % uf)
            print("  u=%8g  mpmath=%.15g  arb_lo=%.15g  arb_hi=%.15g  %s"
                  % (uf, float(mp_val), float(ball.lower()), float(ball.upper()),
                     "OK" if within else "MISMATCH"))
        print("-" * 70)
        print("CF numerical crosscheck: %d/%d OK" % (ok, total))
        if total != 5:
            raise ArithmeticError(
                "CF numerical crosscheck: only %d/5 u values tested" % total)
        if failures:
            raise ArithmeticError("CF numerical crosscheck failed: "
                                  + "; ".join(failures))
        return True
    finally:
        mpmath.mp.dps = old_dps
        ctx.prec = old_prec


def forced_branch_crosscheck(D=200):
    """Exercise series/cf/asym branches."""
    if mpmath is None:
        raise ArithmeticError("mpmath is required for forced_branch_crosscheck")
    old_dps = mpmath.mp.dps
    try:
        mpmath.mp.dps = 100
        pl = make_plan(D, G=30)

        print("Forced-branch cross-check (series, cf, asym)...")
        print("-" * 70)

        candidates = {"series": [], "cf": [], "asym": []}
        n1 = min(pl["K_R"] + 1, 20000)
        for k in range(0, n1):
            xf, uf, lng = _est(k, pl)
            if not math.isfinite(uf):
                continue
            if 1.0 < uf < 100.0 and len(candidates["series"]) < 3:
                candidates["series"].append(k)
            elif 1000.0 < uf < 1e15 and len(candidates["cf"]) < 3:
                candidates["cf"].append(k)
            elif uf > 1e15 and len(candidates["asym"]) < 3:
                candidates["asym"].append(k)
            if all(len(v) >= 3 for v in candidates.values()):
                break

        for method, ks in candidates.items():
            if len(ks) < 3:
                raise ArithmeticError(
                    "forced-branch: only %d nodes found for %s (need 3); "
                    "the plan D=%d may be too small"
                    % (len(ks), method, D))

        failures = []
        for method, ks in candidates.items():
            for k in ks:
                try:
                    ball = node_ball(k, pl, force=method)
                except Exception as e:
                    failures.append("%s k=%d: %s" % (method, k, e))
                    print("  k=%6d method=%-6s  ERROR: %s" % (k, method, e))
                    continue
                t = _mpf_of_fraction(_exact_t(k, pl))
                mp_val = _mpmath_G(t)
                within = _relative_contains(ball, mp_val, rel_tol="1e-25")
                if not within:
                    failures.append("%s k=%d mismatch" % (method, k))
                print("  k=%6d method=%-6s  arb_lo=%.10g  arb_hi=%.10g  %s"
                      % (k, method, float(ball.lower()), float(ball.upper()),
                         "OK" if within else "MISMATCH"))
        print("-" * 70)
        if failures:
            raise ArithmeticError("forced-branch crosscheck failed: "
                                  + "; ".join(failures))
        return True
    finally:
        mpmath.mp.dps = old_dps


def _make_alt_plan(pl, step_factor=0.9):
    """Return a copy of pl with h scaled by step_factor."""
    pl2 = dict(pl)
    K = strip_constants()
    M, d, pi, gam = K["M"], K["d"], K["pi"], K["gam"]
    ten = arb(10)
    N = pl["N"]
    L = (2 * M).log() + (N + 6) * ten.log()
    Lint = arb_ceil_int(L)
    h2 = math.floor(step_factor * 0.97 * 2 * math.pi * 1.4 / Lint * 2 ** 40) / 2 ** 40
    pl2["h"] = h2
    hA = arb(h2)
    pl2["Bdisc"] = 2 * M / ((2 * pi * d / hA).exp() - 1)
    if not (pl2["Bdisc"] * ten ** (N + 6) < 1):
        raise ArithmeticError("alt-h plan: discretisation bound not reached")
    UR = int(math.ceil((N + 8) * LN10)) + 2
    tR = solve_t(math.log(UR), 0.0, 40.0)
    K_R = int(math.ceil(tR / h2)) + 1
    t = arb(K_R) * hA
    uR = (t - (-t).exp()).exp()
    if not (uR >= UR and uR >= 1):
        raise ArithmeticError("alt-h plan: right tail bound violated")
    pl2["K_R"] = K_R
    pl2["tailR"] = (-uR).exp()
    lnuL = 0.5 * (-(N + 8) * LN10 + LN2 - EULER_F) - 1.0
    tL = solve_t(lnuL, -600.0, 0.0)
    K_L = int(math.ceil(-tL / h2)) + 1
    t = -(arb(K_L) * hA)
    uL = (t - (-t).exp()).exp()
    pl2["K_L"] = K_L
    pl2["tailL"] = gam.exp() * uL * uL / 2
    pl2["nodes"] = K_L + K_R + 1
    return pl2


def _compute_full_lambda_interval(pl):
    """Compute a rigorous λ interval [λ_lo, λ_hi] in Arb for the given plan."""
    old_prec = ctx.prec
    try:
        get_EG(pl["EGbits"])
        ctx.prec = pl["Ebits"] + 256
        total = arb(0)
        for k in range(-pl["K_L"], pl["K_R"] + 1):
            total = total + node_ball(k, pl)
            ctx.prec = pl["Ebits"] + 256
        hn, hd = pl["h"].as_integer_ratio()
        hA = arb(hn) / arb(hd)
        lam_raw = total * hA
        lam_lo = lam_raw.lower() - pl["Bdisc"].upper()
        lam_hi = (lam_raw.upper() + pl["tailR"].upper()
                  + pl["tailL"].upper() + pl["Bdisc"].upper())
        return lam_lo, lam_hi
    finally:
        ctx.prec = old_prec


def alt_h_crosscheck(D=100):
    """Recompute the FULL λ integral with h and 0.9·h; must overlap."""
    print("Alt-h full-quadrature crosscheck (D=%d)..." % D)
    print("-" * 70)
    pl1 = make_plan(D, G=30)
    pl1["cf_min_u"] = 1e18
    pl2 = _make_alt_plan(pl1, step_factor=0.9)
    pl2["cf_min_u"] = 1e18
    print("  plan A: h=%.6e  nodes=%d  K_L=%d  K_R=%d"
          % (pl1["h"], pl1["nodes"], pl1["K_L"], pl1["K_R"]))
    print("  plan B: h=%.6e  nodes=%d  K_L=%d  K_R=%d"
          % (pl2["h"], pl2["nodes"], pl2["K_L"], pl2["K_R"]))

    lam_lo_a, lam_hi_a = _compute_full_lambda_interval(pl1)
    lam_lo_b, lam_hi_b = _compute_full_lambda_interval(pl2)

    print("  plan A λ interval: [%s, %s]" % (lam_lo_a, lam_hi_a))
    print("  plan B λ interval: [%s, %s]" % (lam_lo_b, lam_hi_b))
    overlaps = not (lam_hi_a < lam_lo_b or lam_hi_b < lam_lo_a)
    print("  overlap: %s" % overlaps)
    print("-" * 70)
    if not overlaps:
        raise ArithmeticError(
            "alt-h crosscheck: full-quadrature intervals do not overlap")
    return True


def precision_ladder_crosscheck(D=200):
    """Run the same node at two precisions and compare Arb intervals."""
    import random
    random.seed(24680)
    pl = make_plan(D, G=30)
    pl["cf_min_u"] = 1e18
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    ks = random.sample(range(n0, n1), min(5, n1 - n0))
    print("Precision-ladder crosscheck (Ebits vs Ebits+100)...")
    print("-" * 70)
    failures = []
    for k in ks:
        b1 = node_ball(k, pl)
        pl_hi = dict(pl)
        pl_hi["Ebits"] = pl["Ebits"] + 100
        b2 = node_ball(k, pl_hi)
        overlaps = not (b1.upper() < b2.lower() or b2.upper() < b1.lower())
        if not overlaps:
            failures.append("k=%d" % k)
        print("  k=%6d  low-prec=[%s, %s]" % (k, b1.lower(), b1.upper()))
        print("  k=%6d  high-prec=[%s, %s]  %s"
              % (k, b2.lower(), b2.upper(),
                 "OK" if overlaps else "MISMATCH"))
    print("-" * 70)
    if failures:
        raise ArithmeticError("precision-ladder crosscheck failed at "
                              + ", ".join(failures))
    print("Precision-ladder crosscheck: %d nodes overlap" % len(ks))
    return True


def crosscheck_report(D=200, num_nodes=10, seed=12345):
    """Node-level cross-checks.  Raises on any mismatch or if mpmath is
    unavailable."""
    if mpmath is None:
        raise ArithmeticError("mpmath is required for crosscheck_report")
    import random
    random.seed(seed)
    pl = make_plan(D, G=30)
    n0, n1 = -pl["K_L"], pl["K_R"] + 1
    ks = random.sample(range(n0, n1), min(num_nodes, n1 - n0))
    print("Running node-level cross-checks on %d nodes (D=%d)..." % (len(ks), D))
    print("-" * 70)

    ok_mp = 0
    ok_iv = 0
    total_mp = 0
    total_iv = 0
    failures = []
    for k in ks:
        r1 = mpmath_crosscheck_node(k, pl)
        total_mp += 1
        ball, mp_val, within = r1
        if within:
            ok_mp += 1
        else:
            failures.append("mpmath k=%d" % k)
        print("  node k=%6d  mpmath: %s"
              % (k, "OK" if within else "MISMATCH"))

        r2 = mpmath_iv_crosscheck_node(k, pl)
        total_iv += 1
        ball, (iv_lo, iv_hi), overlaps = r2
        if overlaps:
            ok_iv += 1
        else:
            failures.append("mpmath.iv k=%d" % k)
        print("  node k=%6d  mpmath.iv: %s"
              % (k, "OK" if overlaps else "MISMATCH"))

    print("-" * 70)
    print("mpmath (non-interval):  %d/%d OK" % (ok_mp, total_mp))
    print("mpmath.iv (interval):   %d/%d OK" % (ok_iv, total_iv))
    if failures:
        raise ArithmeticError("crosscheck_report failed: "
                              + "; ".join(failures))
    return True


# ============================================================================
# §10. Exact rational tests (run at startup)
# ============================================================================

def lemma7_exact_recurrence_test():
    """Exact rational check of the S-fraction recurrence (Fraction only).

    For three rational z > 0 and the first 16 levels, verifies:
      * B_n/A_n = z^{-1} C_n(1/z) (S-fraction match);
      * the determinant identity;
      * alternation of consecutive differences.

    This tests the algebraic recurrence; the bracketing of the limit
    function is supplied by the Stieltjes theorem (Lemma 7, step 4), not
    by this test.
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


def contraction_identity_test():
    """Exact rational check of the S→J contraction at sampled x values.

    For each rational sample x > 0, verifies that the odd S-fraction
    convergents C_1, C_3, C_5, C_7 (with α_n = ⌈n/2⌉) equal the J-fraction
    approximants J_0, J_1, J_2, J_3 (with b_0 = 1, b_n = 2n+1, a_n = n²).
    """
    def alpha(n):
        return (n + 1) // 2

    def S_conv(n_levels, x):
        P_prev, P_cur = Fraction(0), Fraction(1)
        Q_prev, Q_cur = Fraction(1), Fraction(1)
        for m in range(1, n_levels + 1):
            a = alpha(m)
            P_prev, P_cur = P_cur, P_cur + a * x * P_prev
            Q_prev, Q_cur = Q_cur, Q_cur + a * x * Q_prev
        return P_cur / Q_cur

    def J_conv(n_levels, x):
        b = lambda k: 2 * k + 1
        a_ = lambda k: k * k
        val = Fraction(1) / (1 + b(n_levels) * x)
        for k in range(n_levels - 1, -1, -1):
            val = Fraction(1) / (1 + b(k) * x - a_(k + 1) * x * x * val)
        return val

    for x in (Fraction(1, 7), Fraction(1, 3), Fraction(2, 1), Fraction(17, 5)):
        for n in range(4):
            c_val = S_conv(2 * n + 1, x)
            j_val = J_conv(n, x)
            if c_val != j_val:
                raise AssertionError(
                    "contraction mismatch: x=%s, n=%d, C_%d=%s, J_%d=%s"
                    % (x, n, 2 * n + 1, c_val, n, j_val))
    return "OK: S→J contraction verified for x ∈ {1/7, 1/3, 2, 17/5}, n=0..3"


def contraction_symbolic_test():
    """Verify the S→J contraction at the level of polynomial identities in x.

    Computes P_1, P_3, …, P_{2M+1} from the S-fraction recurrence and
    R_0, R_1, …, R_M from the J-fraction recurrence, as polynomials with
    Fraction coefficients.  Verifies P_{2m+1} == R_m and Q_{2m+1} == S_m
    for m = 0..M (default M = 6).

    This is stronger than contraction_identity_test: it verifies the
    identity as polynomial identities, not just as evaluations at
    specific x.  Together with Lemma 7, step 3 (which proves the identity
    for all m by induction), this gives a strong sanity check on the
    implementation.
    """
    def alpha(n):
        return (n + 1) // 2

    def poly_trim(p):
        p = list(p)
        while len(p) > 1 and p[-1] == 0:
            p = p[:-1]
        return p

    def poly_add(a, b):
        n = max(len(a), len(b))
        r = [(a[i] if i < len(a) else Fraction(0)) +
             (b[i] if i < len(b) else Fraction(0))
             for i in range(n)]
        return poly_trim(r)

    def poly_mul_x_scalar(a, k):
        # a * k * x
        return poly_trim([Fraction(0)] + [c * k for c in a])

    def poly_mul_x2_scalar(a, k):
        # a * k * x²
        return poly_trim([Fraction(0), Fraction(0)] + [c * k for c in a])

    def poly_eq(a, b):
        return poly_trim(a) == poly_trim(b)

    M = 6
    n_max = 2 * M + 1

    # S-fraction polynomials P_n, Q_n
    P = {-1: [Fraction(0)], 0: [Fraction(1)]}
    Q = {-1: [Fraction(1)], 0: [Fraction(1)]}
    for n in range(1, n_max + 1):
        a = alpha(n)
        P[n] = poly_add(P[n-1], poly_mul_x_scalar(P[n-2], a))
        Q[n] = poly_add(Q[n-1], poly_mul_x_scalar(Q[n-2], a))

    # J-fraction polynomials R_m, S_m
    b = lambda k: alpha(2*k) + alpha(2*k+1) if k >= 1 else alpha(1)
    a_ = lambda k: alpha(2*k-1) * alpha(2*k) if k >= 1 else Fraction(1)
    R = {-1: [Fraction(0)], 0: [Fraction(1)]}
    Sj = {-1: [Fraction(1)],
          0: poly_add([Fraction(1)], poly_mul_x_scalar([Fraction(1)], b(0)))}
    for m in range(1, M + 1):
        bm, am = b(m), a_(m)
        R[m] = poly_add(poly_add(R[m-1], poly_mul_x_scalar(R[m-1], bm)),
                        poly_mul_x2_scalar(R[m-2], -am))
        Sj[m] = poly_add(poly_add(Sj[m-1], poly_mul_x_scalar(Sj[m-1], bm)),
                         poly_mul_x2_scalar(Sj[m-2], -am))

    for m in range(M + 1):
        if not poly_eq(P[2*m+1], R[m]):
            raise AssertionError("P_{2m+1} != R_m for m=%d" % m)
        if not poly_eq(Q[2*m+1], Sj[m]):
            raise AssertionError("Q_{2m+1} != S_m for m=%d" % m)

    return ("OK: P_{2m+1} = R_m and Q_{2m+1} = S_m as polynomial identities "
            "for m = 0..%d" % M)


def lemma7_boundary_test():
    """Exact rational test of Lemma-7 recurrence at extreme z."""
    print("Lemma-7 boundary test (exact rational)...")
    for z in (Fraction(1, 10**6), Fraction(1, 1000), Fraction(1, 3),
              Fraction(2, 1), Fraction(10**6), Fraction(10**12, 7)):
        A2, A1 = Fraction(1), z
        B2, B1 = Fraction(0), Fraction(1)
        product_a = 1
        x = 1 / z
        for n in range(1, 33):
            a = (n + 1) // 2
            b = Fraction(1) if n % 2 else z
            A2, A1 = A1, b * A1 + a * A2
            B2, B1 = B1, b * B1 + a * B2
            product_a *= a

            det = A1 * B2 - A2 * B1
            if det != ((-1) ** (n + 1)) * product_a:
                raise AssertionError("boundary determinant failed at z=%s n=%d"
                                     % (z, n))

            tail = Fraction(1)
            for j in range(n, 0, -1):
                alpha = (j + 1) // 2
                tail = 1 + alpha * x / tail
            Cn = 1 / tail
            if B1 / A1 != Cn / z:
                raise AssertionError("boundary S-fraction mismatch at z=%s n=%d"
                                     % (z, n))
        print("  z=%-22s OK" % str(z))
    print("Lemma-7 boundary test: all OK")


# ============================================================================
# §11. CLI
# ============================================================================

def _run_all_rational_tests():
    print("lemma7_exact_recurrence_test():")
    print("  " + lemma7_exact_recurrence_test())
    print("contraction_identity_test():")
    print("  " + contraction_identity_test())
    print("contraction_symbolic_test():")
    print("  " + contraction_symbolic_test())


def _cmd_test():
    _run_all_rational_tests()
    print("All tests passed.")
    return 0


def _cmd_smoke():
    try:
        res = run_certificate(60, G=30, jobs=1, ckpt=None,
                              out="smoke_60.txt", cf_min_u=20.0)
    except ArithmeticError as e:
        print("SMOKE TEST FAILED:", e, file=sys.stderr)
        return 1
    if not (res and res.get("proven")):
        print("SMOKE TEST FAILED (not proven)")
        return 1
    print("SMOKE TEST OK")
    return 0


def _cmd_certify(args):
    try:
        res = run_certificate(
            args.digits, G=args.guard, jobs=args.jobs, chunk=args.chunk,
            ckpt=args.ckpt, out=args.out, max_hours=args.max_hours,
            cf_min_u=args.cf_min_u, trust_checkpoint=args.resume,
            autotune=args.autotune,
        )
    except ArithmeticError as e:
        print("CERTIFICATE FAILED:", e, file=sys.stderr)
        return 1
    if res is None:
        return 2
    if args.resume:
        print("RESUME MODE: no certificate issued.")
        return 0
    return 0 if res.get("proven") else 1


def _cmd_crosscheck(args):
    """Run all cross-checks.  Requires mpmath.  Returns 0 iff all pass."""
    if mpmath is None:
        print("ERROR: crosscheck requires mpmath.", file=sys.stderr)
        print("Install with: pip install mpmath", file=sys.stderr)
        return 1

    failures = []

    def run(name, fn):
        try:
            ok = fn()
            if ok is not True:
                failures.append(name + ": returned non-True")
                print("[FAIL] %s: returned non-True" % name)
        except Exception as e:
            print("[FAIL] %s: %s" % (name, e))
            failures.append(name + ": " + str(e))

    run("crosscheck_report",
        lambda: crosscheck_report(D=args.digits, num_nodes=args.num_nodes))
    run("cf_numerical_crosscheck", cf_numerical_crosscheck)
    run("forced_branch_crosscheck",
        lambda: forced_branch_crosscheck(D=args.digits))
    run("alt_h_crosscheck", lambda: alt_h_crosscheck(D=min(args.digits, 100)))
    run("precision_ladder_crosscheck",
        lambda: precision_ladder_crosscheck(D=args.digits))

    print()
    if failures:
        print("CROSSCHECK FAILURES:")
        for f in failures:
            print("  - " + f)
        return 1
    print("ALL CROSSCHECKS PASSED")
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

    sub.add_parser("test", help="environment + rational Lemma-7 tests")
    sub.add_parser("smoke", help="60-digit end-to-end smoke test")
    sub.add_parser("boundary-test", help="Lemma-7 boundary tests")

    p_cert = sub.add_parser("certify", help="full certificate")
    p_cert.add_argument("--digits", type=int, default=40000)
    p_cert.add_argument("--guard", type=int, default=30)
    p_cert.add_argument("--jobs", type=int, default=0)
    p_cert.add_argument("--chunk", type=int, default=0)
    p_cert.add_argument("--ckpt", type=str, default=None)
    p_cert.add_argument("--out", type=str, default=None)
    p_cert.add_argument("--max-hours", type=float, default=0.0,
                        help="SOFT time limit (running tasks may exceed it)")
    p_cert.add_argument("--cf-min-u", type=float, default=None,
                        help="pre-measured cf_min_u threshold")
    p_cert.add_argument("--autotune", action="store_true",
                        help="run cf_min_u benchmark before the run (slow)")
    p_cert.add_argument("--resume", dest="resume", action="store_true",
                        default=False,
                        help="AUDIT-ONLY: accept checkpoints (no cert)")
    p_cert.add_argument("--no-resume", dest="resume", action="store_false",
                        help="safe mode (default): recompute every chunk")

    p_cross = sub.add_parser("crosscheck", help="independent cross-checks")
    p_cross.add_argument("--digits", type=int, default=200)
    p_cross.add_argument("--num-nodes", type=int, default=10)

    args = parser.parse_args(argv)
    if args.command is None:
        args.command = "test"

    _env_check()
    print()
    _run_all_rational_tests()
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