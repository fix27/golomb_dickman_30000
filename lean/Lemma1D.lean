import Mathlib.Data.Real.Basic
import Mathlib.Data.Int.Basic
import Mathlib.Data.Rat.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Topology.Algebra.InfiniteSum.Real
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# Лемма 1(d): геометрическая оценка хвоста ряда

Оценка остатка экспоненциального ряда, нужная для Леммы 1.

* `tail_ratio`      — отношение соседних членов хвоста не превосходит `1/2`;
* `tail_term_bound` — член хвоста мажорируется геометрической прогрессией;
* `lemma_1_d`       — `|Σ_{k>K} (-u)^k/(k·k!)| ≤ 2·u^{K+1}/((K+1)·(K+1)!)`.

Источник: `archive/Формальное доказательство.lean`, часть 1.
Текст доказательств перенесён без изменений.
-/

namespace GolombDickman

open Real Polynomial


/-- Отношение соседних членов ряда: |a_{k+1}/a_k| ≤ 1/2 при k ≥ K + 1 ≥ 2u. -/
lemma tail_ratio {u : ℝ} (hu : 0 < u) {K k : ℕ}
    (hK : 2 * u ≤ K + 1) (hk : K + 1 ≤ k) :
    u * k / ((k + 1 : ℝ) ^ 2) ≤ 1 / 2 := by
  have hkj : (K + 1 : ℝ) ≤ k := by exact_mod_cast hk
  have hu' : u ≤ (K + 1 : ℝ) / 2 := by linarith
  have h1 : u * k ≤ ((K + 1 : ℝ) / 2) * k :=
    mul_le_mul_of_nonneg_right hu' (Nat.cast_nonneg k)
  have h2 : ((K + 1 : ℝ) / 2) * k ≤ ((k : ℝ) / 2) * (k + 1) := by
    nlinarith
  calc u * k / ((k + 1 : ℝ) ^ 2)
      ≤ (((k : ℝ) / 2) * (k + 1)) / ((k + 1 : ℝ) ^ 2) := by
        apply div_le_div_of_nonneg_right _ (by positivity)
        linarith
    _ = (k : ℝ) / (2 * (k + 1)) := by field_simp; ring
    _ ≤ 1 / 2 := by
        rw [div_le_div_iff₀ (by positivity) (by norm_num)]
        nlinarith

/-- Модуль r-го члена хвоста ряда мажорируется геометрической прогрессией (1/2)^r. -/
lemma tail_term_bound {u : ℝ} (hu : 0 < u) {K : ℕ} (hK : 2 * u ≤ K + 1) (r : ℕ) :

    |(-u) ^ (K + 1 + r) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|
    ≤ (u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial)) * (1 / 2 : ℝ) ^ r := by
  induction r with
  | zero =>
      simp only [add_zero, pow_zero, mul_one, abs_div, abs_pow, abs_neg, abs_of_pos hu]
      have hpos : 0 < ((K + 1 : ℕ) : ℝ) * (K + 1).factorial := by positivity
      rw [abs_of_pos hpos]
  | succ r ih =>
      have hk : K + 1 ≤ K + 1 + r := Nat.le_add_right _ _
      have hratio := tail_ratio hu hK hk
      have hpos_r : 0 < ((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial := by positivity
      have hpos_r1 : 0 < ((K + 1 + (r + 1) : ℕ) : ℝ) * (K + 1 + (r + 1)).factorial := by positivity
      have hstep :

          |(-u) ^ (K + 1 + (r + 1)) / (((K + 1 + (r + 1) : ℕ) : ℝ) * (K + 1 + (r + 1)).factorial)|
          = |(-u) ^ (K + 1 + r) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|
            * (u * (K + 1 + r : ℕ) / (((K + 1 + r : ℕ) : ℝ) + 1) ^ 2) := by
        rw [abs_div, abs_div, abs_of_pos hpos_r1, abs_of_pos hpos_r]
        rw [abs_pow, abs_pow, abs_neg, abs_of_pos hu]
        have hidx : K + 1 + (r + 1) = (K + 1 + r) + 1 := by ring
        rw [hidx, pow_succ, Nat.factorial_succ]
        push_cast
        field_simp
        ring
      rw [hstep, pow_succ, ← mul_assoc]
      calc |(-u) ^ (K + 1 + r) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|
              * (u * (K + 1 + r : ℕ) / (((K + 1 + r : ℕ) : ℝ) + 1) ^ 2)
          ≤ ((u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial)) * (1 / 2 : ℝ) ^ r)
              * (1 / 2) := by
            apply mul_le_mul ih hratio (by positivity) (by positivity)
        _ = (u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial)) * ((1 / 2 : ℝ) ^ r * (1 / 2)) := by
            ring

/-- Лемма 1(d): |Σ_{k>K} (-u)^k/(k·k!)| ≤ 2·u^{K+1}/((K+1)·(K+1)!). -/
theorem lemma_1_d {u : ℝ} (hu : 0 < u) {K : ℕ} (hK : 2 * u ≤ K + 1) :
    |∑' (r : ℕ), ((-u) ^ (K + 1 + r)) /
        (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|
    ≤ 2 * u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial) := by
  set C := u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial)
  have hC : 0 ≤ C := by positivity
  have hsum_geom : Summable (fun r : ℕ => C * (1 / 2 : ℝ) ^ r) :=
    Summable.mul_left C (summable_geometric_of_lt_one (by norm_num) (by norm_num))
  have hsum_abs : Summable (fun r : ℕ =>

      |(-u) ^ (K + 1 + r) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|) :=
    Summable.of_nonneg_of_le (fun _ => abs_nonneg _) (tail_term_bound hu hK) hsum_geom
  calc |∑' (r : ℕ), ((-u) ^ (K + 1 + r)) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)|
      ≤ ∑' (r : ℕ), |((-u) ^ (K + 1 + r)) / (((K + 1 + r : ℕ) : ℝ) * (K + 1 + r).factorial)| :=
        abs_tsum_le_tsum_abs hsum_abs
    _ ≤ ∑' (r : ℕ), C * (1 / 2 : ℝ) ^ r :=
        tsum_le_tsum (tail_term_bound hu hK) hsum_abs hsum_geom
    _ = C * ∑' (r : ℕ), (1 / 2 : ℝ) ^ r := tsum_mul_left
    _ = C * (1 - 1 / 2 : ℝ)⁻¹ := by rw [tsum_geometric_of_lt_one (by norm_num) (by norm_num)]
    _ = 2 * u ^ (K + 1) / (((K + 1 : ℕ) : ℝ) * (K + 1).factorial) := by
        dsimp [C]; norm_num; ring

end GolombDickman
