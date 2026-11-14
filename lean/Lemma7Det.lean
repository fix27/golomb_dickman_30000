import Mathlib.Data.Int.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.IntervalCases

/-!
# Лемма 7, часть 1: детерминантное тождество для β-дроби

Коэффициенты `a_n = ⌈n/2⌉`, `b_n = 1` при нечётном `n` и `b_n = z`
при чётном.  Определитель `Ddet` подчиняется
`D(n+1) = (-1)^{n+1}·a_1⋯a_n`.

Источник: `archive/Формальное доказательство.lean`, часть 2; совпадает с
`archive/Лемма 7 тождество детерминанта (полное доказательство).lean`.
Текст доказательств перенесён без изменений.
-/

namespace GolombDickman

open Polynomial


def aCoeff (n : ℕ) : ℤ := ((n + 1) / 2 : ℕ)

noncomputable def bCoeff (n : ℕ) : ℤ[X] :=
  if n % 2 = 1 then (1 : ℤ[X]) else X

noncomputable def An : ℕ → ℤ[X]
  | 0     => 1
  | 1     => X
  | n + 2 => bCoeff (n + 1) * An (n + 1) + C (aCoeff (n + 1)) * An n

noncomputable def Bn : ℕ → ℤ[X]
  | 0     => 0
  | 1     => 1
  | n + 2 => bCoeff (n + 1) * Bn (n + 1) + C (aCoeff (n + 1)) * Bn n

noncomputable def Ddet : ℕ → ℤ[X]
  | 0     => 0
  | k + 1 => An (k + 1) * Bn k - An k * Bn (k + 1)

def aProd : ℕ → ℤ
  | 0     => 1
  | n + 1 => aProd n * aCoeff (n + 1)

lemma Ddet_succ (n : ℕ) :
    Ddet (n + 2) = - C (aCoeff (n + 1)) * Ddet (n + 1) := by
  simp only [Ddet, An, Bn]
  ring

theorem Ddet_eq (n : ℕ) :
    Ddet (n + 1) = C (((-1 : ℤ) ^ (n + 1)) * aProd n) := by
  induction n with
  | zero =>
      simp [Ddet, An, Bn, aProd]
  | succ n ih =>
      rw [Ddet_succ, ih, aProd]
      push_cast
      have hpow : (-1 : ℤ) ^ (n + 2) = -((-1 : ℤ) ^ (n + 1)) := by ring
      rw [hpow]; ring

end GolombDickman
