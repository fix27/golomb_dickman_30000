namespace GolombDickman
open Real

/-- |a_{k+1}/a_k| ≤ 1/2 при k ≥ K+1 ≥ 2u. -/
lemma tail_ratio {u : ℝ} (hu : 0 < u) {K k : ℕ}
    (hK : 2 * u ≤ K + 1) (hk : K + 1 ≤ k) :
    u * k / ((k + 1 : ℝ) ^ 2) ≤ 1 / 2 := by
  have hk1 : (0 : ℝ) < k + 1 := by positivity
  have hkj : (k : ℝ) ≥ K + 1 := by exact_mod_cast hk
  have hu' : u ≤ (K + 1 : ℝ) / 2 := by linarith
  have h1 : u * k ≤ ((K + 1 : ℝ) / 2) * k := by
    exact mul_le_mul_of_nonneg_right hu' (Nat.cast_nonneg k)
  have h2 : ((K + 1 : ℝ) / 2) * k ≤ ((k : ℝ) / 2) * (k + 1) := by
    have : (K + 1 : ℝ) ≤ k := hkj
    nlinarith
  calc u * k / ((k + 1 : ℝ) ^ 2)
      ≤ (((k : ℝ) / 2) * (k + 1)) / ((k + 1 : ℝ) ^ 2) := by
        apply div_le_div_of_nonneg_right _ (by positivity)
        linarith
    _ = (k : ℝ) / (2 * (k + 1)) := by field_simp; ring
    _ ≤ 1 / 2 := by
        rw [div_le_div_iff (by positivity) (by norm_num)]
        nlinarith

/-- Лемма 1(d): |Σ_{k>K} (-u)^k/(k·k!)| ≤ 2·u^{K+1}/((K+1)·(K+1)!). -/
theorem lemma_1_d {u : ℝ} (hu : 0 < u) {K : ℕ}
    (hK : 2 * u ≤ K + 1) :
    |∑' (r : ℕ), ((-u) ^ (K + 1 + r)) /
        ((K + 1 + r : ℝ) * (K + 1 + r).factorial)|
    ≤ 2 * u ^ (K + 1) / ((K + 1 : ℝ) * (K + 1).factorial) := by
  -- |a_{K+1+r}| ≤ |a_{K+1}| · (1/2)^r по tail_ratio → геометрическая сумма.
  sorry

end GolombDickman