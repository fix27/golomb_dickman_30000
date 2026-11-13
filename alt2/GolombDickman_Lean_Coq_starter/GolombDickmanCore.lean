import Mathlib

/-!
# Golomb–Dickman: formally checked algebraic core

This file formalizes the local algebraic contraction step used in Lemma 7
of `golomb_dickman_24.py`.  It does NOT yet formalize the analytic input
(Stieltjes convergence/alternation, Laguerre resolvents, strip estimates,
or Arb interval arithmetic).
-/

namespace GolombDickman

/--
One contraction step for the numerator or denominator of an S-fraction.

The hypotheses are the three consecutive S-recurrence equations:
  q = t + c*x*r,
  s = q + b*x*t,
  n = s + a*x*q.

The conclusion is the corresponding second-order J-recurrence.
-/
theorem s_to_j_contraction_step
    {R : Type*} [CommRing R]
    (x a b c r t q s n : R)
    (hq : q = t + c * x * r)
    (hs : s = q + b * x * t)
    (hn : n = s + a * x * q) :
    n = (1 + (b + a) * x) * q - (b * c) * x^2 * r := by
  rw [hn, hs, hq]
  ring

/-- The coefficient identities for alpha_n = ceil(n/2), in nat arithmetic. -/
theorem contracted_diagonal (m : Nat) :
    (m + (m + 1) : Nat) = 2 * m + 1 := by omega

theorem contracted_off_diagonal (m : Nat) :
    (m * m : Nat) = m ^ 2 := by ring

/--
Local determinant identity for the continuant recurrence.
If A_n = b*A_(n-1)+a*A_(n-2) and
   B_n = b*B_(n-1)+a*B_(n-2),
then the adjacent determinant is multiplied by -a.
-/
theorem continuant_determinant_step
    {R : Type*} [CommRing R]
    (a b Aprev Aprevprev Bprev Bprevprev : R) :
    (b * Aprev + a * Aprevprev) * Bprev
      - Aprev * (b * Bprev + a * Bprevprev)
    = -a * (Aprev * Bprevprev - Aprevprev * Bprev) := by
  ring

end GolombDickman
