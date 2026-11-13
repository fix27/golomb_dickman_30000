import Mathlib.Data.Real.Basic
import Mathlib.Data.Int.Basic
import Mathlib.Data.Rat.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Topology.Algebra.InfiniteSum.Real
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

namespace GolombDickman

open Real Polynomial

/-! ### Часть 1. Лемма 1(d): оценка остатка степенного ряда для S(u) -/

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


/-! ### Часть 2. Лемма 7: детерминантное тождество -/

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


/-! ### Часть 3. Лемма 7 (шаг 3): стягивание S-дроби в J-дробь -/

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
