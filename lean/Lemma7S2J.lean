import Mathlib.Data.Rat.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-!
# Лемма 7, часть 2: стягивание S-дроби в J-дробь

Рекурсии числителей и знаменателей S-дроби (`Ps`, `Qs`) и J-дроби
(`Rs`, `Ss`), трёхшаговые тождества `Ps_step3`, `Qs_step3` и главная
теорема `S_to_J_contraction`: `Ps (2m) = Rs m`, `Qs (2m) = Ss m`.

Источник: `archive/Формальное доказательство.lean`, часть 3.
Текст доказательств перенесён без изменений.
-/

namespace GolombDickman

open Polynomial


section S2J

variable (α : ℕ → ℚ)

/-- Числители S-дроби (сдвиг индекса: Ps 0 = P_{-1} = 0, Ps 1 = P_0 = 1). -/
noncomputable def Ps : ℕ → ℚ[X]
  | 0     => 0
  | 1     => 1
  | n + 2 => Ps (n + 1) + C (α (n + 1)) * X * Ps n

/-- Знаменатели S-дроби (Qs 0 = Q_{-1} = 1, Qs 1 = Q_0 = 1). -/
noncomputable def Qs : ℕ → ℚ[X]
  | 0     => 1
  | 1     => 1
  | n + 2 => Qs (n + 1) + C (α (n + 1)) * X * Qs n

/-- Коэффициенты J-дроби: b_0 = α_1, b_m = α_{2m} + α_{2m+1}, a_m = α_{2m-1} α_{2m}. -/
def jb : ℕ → ℚ
  | 0     => α 1
  | m + 1 => α (2 * (m + 1)) + α (2 * (m + 1) + 1)

def ja : ℕ → ℚ
  | 0     => 1
  | m + 1 => α (2 * m + 1) * α (2 * (m + 1))

/-- Числители J-дроби (Rs 0 = R_{-1} = 0, Rs 1 = R_0 = 1). -/
noncomputable def Rs : ℕ → ℚ[X]
  | 0     => 0
  | 1     => 1
  | m + 2 => (1 + C (jb α (m + 1)) * X) * Rs (m + 1) - C (ja α (m + 1)) * X ^ 2 * Rs m

/-- Знаменатели J-дроби (Ss 0 = S_{-1} = 1, Ss 1 = S_0 = 1 + b_0 X). -/
noncomputable def Ss : ℕ → ℚ[X]
  | 0     => 1
  | 1     => 1 + C (jb α 0) * X
  | m + 2 => (1 + C (jb α (m + 1)) * X) * Ss (m + 1) - C (ja α (m + 1)) * X ^ 2 * Ss m

/-- Трехшаговое рекуррентное соотношение для P_{2m+3} (шаги (a), (b), (c)). -/
lemma Ps_step3 (m : ℕ) :
    Ps α (2 * (m + 2)) =
      (1 + C (jb α (m + 1)) * X) * Ps α (2 * (m + 1)) - C (ja α (m + 1)) * X ^ 2 * Ps α (2 * m) := by
  have h4 : 2 * (m + 2) = (2 * m + 2) + 2 := by ring
  have h3 : 2 * m + 3 = (2 * m + 1) + 2 := by ring
  have h2 : 2 * (m + 1) = (2 * m) + 2 := by ring
  rw [h4, Ps]
  have h3' : 2 * m + 2 + 1 = (2 * m + 1) + 2 := by ring
  rw [h3', Ps, h2, Ps]
  simp only [jb, ja, map_add, map_mul]
  have idx1 : 2 * m + 2 + 1 = 2 * (m + 1) + 1 := by ring
  have idx2 : 2 * m + 1 + 1 = 2 * (m + 1) := by ring
  rw [idx1, idx2]
  ring

/-- Трехшаговое рекуррентное соотношение для Q_{2m+3}. -/
lemma Qs_step3 (m : ℕ) :
    Qs α (2 * (m + 2)) =
      (1 + C (jb α (m + 1)) * X) * Qs α (2 * (m + 1)) - C (ja α (m + 1)) * X ^ 2 * Qs α (2 * m) := by
  have h4 : 2 * (m + 2) = (2 * m + 2) + 2 := by ring
  have h2 : 2 * (m + 1) = (2 * m) + 2 := by ring
  rw [h4, Qs]
  have h3' : 2 * m + 2 + 1 = (2 * m + 1) + 2 := by ring
  rw [h3', Qs, h2, Qs]
  simp only [jb, ja, map_add, map_mul]
  have idx1 : 2 * m + 2 + 1 = 2 * (m + 1) + 1 := by ring
  have idx2 : 2 * m + 1 + 1 = 2 * (m + 1) := by ring
  rw [idx1, idx2]
  ring

/-- Теорема о стягивании S → J: P_{2m+1} = R_m и Q_{2m+1} = S_m как многочлены от X. -/
theorem S_to_J_contraction (m : ℕ) :
    Ps α (2 * m) = Rs α m ∧ Qs α (2 * m) = Ss α m := by
  induction m using Nat.strong_induction_on with
  | h m ih =>
    match m with
    | 0 =>
        simp [Ps, Qs, Rs, Ss]
    | 1 =>
        simp [Ps, Qs, Rs, Ss, jb]
    | m + 2 =>
        have ⟨ihP1, ihQ1⟩ := ih (m + 1) (Nat.lt_succ_self _)
        have ⟨ihP0, ihQ0⟩ := ih m (Nat.lt_trans (Nat.lt_succ_self m) (Nat.lt_succ_self (m + 1)))
        constructor
        · rw [Ps_step3, ihP1, ihP0, Rs]
        · rw [Qs_step3, ihQ1, ihQ0, Ss]

end S2J

end GolombDickman
