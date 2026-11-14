namespace GolombDickman
open Real MeasureTheory

/-- Лемма 3, правая часть:
    для h > 0 и u(K_R h) ≥ 1,  Σ_{k>K_R} h·G(kh) ≤ e^{−u(K_R h)}. -/
theorem lemma_3_right {h : ℝ} (hh : 0 < h)
    {KR : ℤ} (hu : 1 ≤ uR (KR * h)) :
    ∑' (n : ℕ), h * G ((KR + 1 + n) * h) ≤ exp (-(uR (KR * h))) := by
  -- g_R монотонно убывает при u ≥ 1, поэтому h·g_R(kh) ≤ ∫_{(k-1)h}^{kh} g_R.
  -- Сумма мажорируется ∫_{KR h}^∞ g_R = e^{−u(K_R h)} (лемма 2).
  sorry

end GolombDickman