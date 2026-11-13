import Mathlib.Data.Real.Basic
import Mathlib.Data.Int.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.IntervalCases

namespace GolombDickman

open Polynomial

/-! ### Коэффициенты β-дроби из леммы 7

  a_n = ⌈n/2⌉,  b_n = 1 при нечётном n,  b_n = z при чётном n.
  Работаем в ℤ[z], т.е. X обозначает z. -/

/-- a_n = ⌈n/2⌉. -/
def aCoeff (n : ℕ) : ℤ := ((n + 1) / 2 : ℕ)

/-- b_n(z) = 1 (n нечётное) или z = X (n чётное). -/
noncomputable def bCoeff (n : ℕ) : ℤ[X] :=
  if n % 2 = 1 then (1 : ℤ[X]) else X

/-- A_n, индексированные так: An 0 = A_{-1} = 1, An 1 = A_0 = z,
    An (n+2) = A_{n+1}. -/
noncomputable def An : ℕ → ℤ[X]
  | 0     => 1
  | 1     => X
  | n + 2 => bCoeff (n + 1) * An (n + 1) + C (aCoeff (n + 1)) * An n

noncomputable def Bn : ℕ → ℤ[X]
  | 0     => 0
  | 1     => 1
  | n + 2 => bCoeff (n + 1) * Bn (n + 1) + C (aCoeff (n + 1)) * Bn n

/-- "Определитель" D(k) = A_{k-1} B_{k-2} - A_{k-2} B_{k-1} в исходной
    индексации, переписанный через An/Bn: D(k) = An(k) Bn(k-1) - An(k-1) Bn(k). -/
noncomputable def Ddet : ℕ → ℤ[X]
  | 0     => 0
  | k + 1 => An (k + 1) * Bn k - An k * Bn (k + 1)

/-- Произведение a_1 a_2 ⋯ a_n. -/
def aProd : ℕ → ℤ
  | 0     => 1
  | n + 1 => aProd n * aCoeff (n + 1)

/-! ### Лемма 7: детерминантное тождество

  D(n+1) = (-1)^{n+1} · a_1 ⋯ a_n,   n ≥ 0.  (индексация: D(1) = -1) -/

lemma Ddet_succ (n : ℕ) :
    Ddet (n + 2) = - C (aCoeff (n + 1)) * Ddet (n + 1) := by
  -- Раскрываем рекурсию An, Bn на шаге n+1, затем группируем.
  simp only [Ddet, An, Bn]
  -- После simp: bCoeff(n+1)*An(n+1)+C(a_{n+1})*An n  и аналогично для B.
  ring

theorem Ddet_eq (n : ℕ) :
    Ddet (n + 1) = C (((-1 : ℤ) ^ (n + 1)) * aProd n) := by
  induction n with
  | zero =>
      -- D(1) = A_0 B_{-1} - A_{-1} B_0 = z·0 - 1·1 = -1.
      simp [Ddet, An, Bn, aProd]
  | succ n ih =>
      rw [Ddet_succ, ih, aProd]
      push_cast
      -- (-a_{n+1}) · (-1)^{n+1} a_1⋯a_n = (-1)^{n+2} a_1⋯a_{n+1}.
      have hpow : (-1 : ℤ) ^ (n + 2) = -((-1 : ℤ) ^ (n + 1)) := by ring
      rw [hpow]; ring

end GolombDickman