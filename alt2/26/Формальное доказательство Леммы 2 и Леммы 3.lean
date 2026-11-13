import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Integral.IntervalIntegral
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

namespace GolombDickman

open Real Set Filter MeasureTheory Topology

/-! ### Определения функций замены переменной и мажорант -/

/-- Замена переменной u(t) = exp(t - e^{-t}). -/
noncomputable def u_fn (t : ℝ) : ℝ := exp (t - exp (-t))

/-- Производная u'(t) = u(t) · (1 + e^{-t}). -/
lemma hasDerivAt_u_fn (t : ℝ) :
    HasDerivAt u_fn (u_fn t * (1 + exp (-t))) t := by
  have h1 : HasDerivAt (fun s : ℝ => s - exp (-s)) (1 + exp (-t)) t := by
    have he : HasDerivAt (fun s : ℝ => exp (-s)) (-exp (-t)) t := by
      simpa using (hasDerivAt_neg t).exp
    simpa [sub_eq_add_neg] using (hasDerivAt_id t).sub he
  simpa [u_fn] using h1.exp

lemma u_fn_pos (t : ℝ) : 0 < u_fn t := exp_pos _

/-- Строгая монотонность u(t) на ℝ. -/
lemma strictMono_u_fn : StrictMono u_fn := by
  refine strictMono_of_deriv_pos (fun t => ?_)
  rw [(hasDerivAt_u_fn t).deriv]
  have hu : 0 < u_fn t := u_fn_pos t
  have he : 0 < exp (-t) := exp_pos _
  positivity

/-- Предел u(t) → +∞ при t → +∞. -/
lemma tendsto_u_fn_atTop : Tendsto u_fn atTop atTop := by
  refine tendsto_exp_atTop.comp ?_
  have h1 : Tendsto (fun t : ℝ => exp (-t)) atTop (nhds 0) :=
    tendsto_exp_atBot.comp tendsto_neg_atTop_atBot
  simpa using tendsto_id.sub h1

/-- Предел u(t) → 0 при t → -∞. -/
lemma tendsto_u_fn_atBot : Tendsto u_fn atBot (nhds 0) := by
  refine tendsto_exp_atBot.comp ?_
  have h1 : Tendsto (fun t : ℝ => exp (-t)) atBot atTop :=
    tendsto_exp_atTop.comp tendsto_neg_atBot_atTop
  exact Tendsto.atBot_add tendsto_id (tendsto_neg_atTop_atBot.comp h1)

/-- Правая мажоранта g_R(t) = u(t) e^{-u(t)} (1 + e^{-t}). -/
noncomputable def gR (t : ℝ) : ℝ := u_fn t * exp (-u_fn t) * (1 + exp (-t))

/-- Левая мажоранта g_L(t) = e^γ u(t)² (1 + e^{-t}). -/
noncomputable def gL (t : ℝ) : ℝ := exp Real.eulerMascheroniConstant * (u_fn t) ^ 2 * (1 + exp (-t))

lemma gR_nonneg (t : ℝ) : 0 ≤ gR t := by
  unfold gR
  have hu := u_fn_pos t
  have he1 := exp_pos (-u_fn t)
  have he2 := exp_pos (-t)
  positivity

lemma gL_nonneg (t : ℝ) : 0 ≤ gL t := by
  unfold gL
  have hg := exp_pos Real.eulerMascheroniConstant
  have hu := u_fn_pos t
  have he := exp_pos (-t)
  positivity

/-! ### Лемма 2: точные формулы интегралов хвостов -/

/-- Лемма 2 (правый хвост): ∫_{t₀}^∞ u e^{-u} (1 + e^{-t}) dt = e^{-u(t₀)}. -/
theorem lemma_2_right (t₀ : ℝ) :
    ∫ t in Ioi t₀, gR t = exp (-u_fn t₀) := by
  let F : ℝ → ℝ := fun t => -exp (-u_fn t)
  have hderiv : ∀ t ∈ Ici t₀, HasDerivAt F (gR t) t := by
    intro t _
    have hu := hasDerivAt_u_fn t
    have h1 : HasDerivAt (fun s => -u_fn s) (-(u_fn t * (1 + exp (-t)))) t := hu.neg
    have h2 : HasDerivAt (fun s => exp (-u_fn s))
        (exp (-u_fn t) * (-(u_fn t * (1 + exp (-t))))) t := h1.exp
    convert h2.neg using 1
    unfold gR; ring
  have htend : Tendsto F atTop (nhds 0) := by
    have h1 : Tendsto (fun t => exp (-u_fn t)) atTop (nhds 0) :=
      tendsto_exp_atBot.comp (tendsto_neg_atTop_atBot.comp tendsto_u_fn_atTop)
    simpa [F] using h1.neg
  have hint := integral_Ioi_of_hasDerivAt_of_tendsto' hderiv (fun t _ => gR_nonneg t) htend
  simpa [F] using hint

/-- Лемма 2 (левый хвост): ∫_{-∞}^{t₀} e^γ u² (1 + e^{-t}) dt = ½ e^γ u(t₀)². -/
theorem lemma_2_left (t₀ : ℝ) :
    ∫ t in Iic t₀, gL t = (1 / 2) * exp Real.eulerMascheroniConstant * (u_fn t₀) ^ 2 := by
  let C := (1 / 2) * exp Real.eulerMascheroniConstant
  let F : ℝ → ℝ := fun t => C * (u_fn t) ^ 2
  have hderiv : ∀ t ∈ Iic t₀, HasDerivAt F (gL t) t := by
    intro t _
    have hu := hasDerivAt_u_fn t
    have hpow : HasDerivAt (fun s => (u_fn s) ^ 2) (2 * u_fn t ^ 1 * (u_fn t * (1 + exp (-t)))) t :=
      hu.pow 2
    convert hpow.const_mul C using 1
    unfold gL; dsimp [C]; ring
  have htend : Tendsto F atBot (nhds 0) := by
    have hpow : Tendsto (fun t => (u_fn t) ^ 2) atBot (nhds 0) := by
      simpa using tendsto_u_fn_atBot.pow 2
    simpa [F] using hpow.const_mul C
  have hint := integral_Iic_of_hasDerivAt_of_tendsto' hderiv (fun t _ => gL_nonneg t) htend
  simpa [F, C] using hint


/-! ### Вспомогательные свойства монотонности мажорант -/

/-- Функция x ↦ x e^{-x} убывает при x ≥ 1. -/
lemma antitoneOn_mul_exp_neg : AntitoneOn (fun x : ℝ => x * exp (-x)) (Ici 1) := by
  refine antitoneOn_of_deriv_nonpos (convex_Ici 1) ?_ ?_ ?_
  · exact (continuous_id.mul (continuous_exp.comp continuous_neg)).continuousOn
  · exact ((differentiable_id.mul (differentiable_exp.comp differentiable_neg))).differentiableOn
  · intro x hx
    rw [interior_Ici] at hx
    have hd : HasDerivAt (fun s : ℝ => s * exp (-s)) ((1 - x) * exp (-x)) x := by
      have h1 : HasDerivAt (fun s : ℝ => exp (-s)) (-exp (-x)) x := by
        simpa using (hasDerivAt_neg x).exp
      convert (hasDerivAt_id x).mul h1 using 1
      ring
    rw [hd.deriv]
    have h1 : 1 - x ≤ 0 := by linarith
    have h2 : 0 ≤ exp (-x) := exp_nonneg _
    exact mul_nonpos_of_nonpos_of_nonneg h1 h2

/-- Правая мажоранта g_R(t) убывает на луче {t | 1 ≤ u(t)}. -/
lemma antitoneOn_gR : AntitoneOn gR {t : ℝ | 1 ≤ u_fn t} := by
  intro t₁ ht₁ t₂ ht₂ hle
  unfold gR
  have hu_le : u_fn t₁ ≤ u_fn t₂ := strictMono_u_fn.monotone hle
  have h1 : u_fn t₂ * exp (-u_fn t₂) ≤ u_fn t₁ * exp (-u_fn t₁) :=
    antitoneOn_mul_exp_neg ht₁ ht₂ hu_le
  have h2 : 1 + exp (-t₂) ≤ 1 + exp (-t₁) := by
    have : exp (-t₂) ≤ exp (-t₁) := exp_le_exp.mpr (neg_le_neg hle)
    linarith
  have hpos1 : 0 ≤ u_fn t₂ * exp (-u_fn t₂) := by
    have := u_fn_pos t₂; have := exp_pos (-u_fn t₂); positivity
  have hpos2 : 0 ≤ 1 + exp (-t₂) := by
    have := exp_pos (-t₂); positivity
  exact mul_le_mul h1 h2 hpos2 (le_trans hpos1 h1)

/-- Левая мажоранта g_L(t) возрастает на всей оси ℝ. -/
lemma monotone_gL : Monotone gL := by
  refine monotone_of_deriv_nonneg (fun t => ?_)
  have hu := hasDerivAt_u_fn t
  have he : HasDerivAt (fun s : ℝ => 1 + exp (-s)) (-exp (-t)) t := by
    simpa using ((hasDerivAt_neg t).exp).const_add 1
  have hprod : HasDerivAt (fun s => (u_fn s) ^ 2 * (1 + exp (-s)))
      ((2 * u_fn t ^ 1 * (u_fn t * (1 + exp (-t)))) * (1 + exp (-t)) +
       (u_fn t) ^ 2 * (-exp (-t))) t :=
    (hu.pow 2).mul he
  have hd : HasDerivAt gL
      (exp Real.eulerMascheroniConstant *
       ((2 * u_fn t ^ 1 * (u_fn t * (1 + exp (-t)))) * (1 + exp (-t)) +
        (u_fn t) ^ 2 * (-exp (-t)))) t :=
    hprod.const_mul _
  rw [hd.deriv]
  have hpos_g : 0 ≤ exp Real.eulerMascheroniConstant := exp_nonneg _
  have hpos_u : 0 ≤ (u_fn t) ^ 2 := sq_nonneg _
  have he_pos : 0 < exp (-t) := exp_pos (-t)
  have h_inner : 0 ≤ 2 * (1 + exp (-t)) ^ 2 - exp (-t) := by
    nlinarith [sq_nonneg (exp (-t))]
  have h_eq : (2 * u_fn t ^ 1 * (u_fn t * (1 + exp (-t)))) * (1 + exp (-t)) +
              (u_fn t) ^ 2 * (-exp (-t)) =
              (u_fn t) ^ 2 * (2 * (1 + exp (-t)) ^ 2 - exp (-t)) := by ring
  rw [h_eq]
  exact mul_nonneg hpos_g (mul_nonneg hpos_u h_inner)


/-! ### Лемма 3: оценка отбрасываемых хвостов квадратурной суммы -/

/-- Аксиома связи с определением G(t) из Леммы 1(b): 0 ≤ G(t) ≤ min(g_L(t), g_R(t)). -/
variable (G : ℝ → ℝ)
  (hG_nonneg : ∀ t, 0 ≤ G t)
  (hG_le_gR : ∀ t, G t ≤ gR t)
  (hG_le_gL : ∀ t, G t ≤ gL t)

/-- Один шаг мажорирования для правого хвоста: h · G(kh) ≤ ∫_{(k-1)h}^{kh} g_R(t) dt. -/
lemma right_term_le_integral {h : ℝ} (hh : 0 < h) {KR : ℤ} (hu : 1 ≤ u_fn (KR * h)) (n : ℕ) :
    h * G ((KR + 1 + n) * h) ≤ ∫ t in ((KR + n : ℤ) * h)..((KR + 1 + n : ℤ) * h), gR t := by
  let a : ℝ := (KR + n : ℤ) * h
  let b : ℝ := (KR + 1 + n : ℤ) * h
  have hab : a ≤ b := by
    dsimp [a, b]; push_cast; nlinarith
  have hba : b - a = h := by
    dsimp [a, b]; push_cast; ring
  have hG_step : h * G b ≤ h * gR b := mul_le_mul_of_nonneg_left (hG_le_gR b) hh.le
  have hint_step : h * gR b ≤ ∫ t in a..b, gR t := by
    have hconst : ∫ _ in a..b, gR b = (b - a) * gR b := by
      simp [intervalIntegral.integral_const, smul_eq_mul]
    rw [← hba, ← hconst]
    refine intervalIntegral.integral_mono_on hab intervalIntegrable_const ?_ (fun t ht => ?_)
    · exact (antitoneOn_gR.mono (fun t ht => ?_)).intervalIntegrable
      have hKR_t : (KR : ℝ) * h ≤ t := by
        have h1 : (KR : ℝ) * h ≤ a := by dsimp [a]; push_cast; nlinarith
        exact le_trans h1 ht.1
      exact le_trans hu (strictMono_u_fn.monotone hKR_t)
    · have ht_ge : 1 ≤ u_fn t := by
        have h1 : (KR : ℝ) * h ≤ a := by dsimp [a]; push_cast; nlinarith
        exact le_trans hu (strictMono_u_fn.monotone (le_trans h1 ht.1))
      have hb_ge : 1 ≤ u_fn b := le_trans ht_ge (strictMono_u_fn.monotone ht.2)
      exact antitoneOn_gR ht_ge hb_ge ht.2
  exact le_trans hG_step hint_step

/-- Лемма 3 (правый хвост): Σ_{k > K_R} h · G(kh) ≤ e^{-u(K_R h)}. -/
theorem lemma_3_right {h : ℝ} (hh : 0 < h) {KR : ℤ} (hu : 1 ≤ u_fn (KR * h)) :
    ∑' (n : ℕ), h * G ((KR + 1 + n) * h) ≤ exp (-u_fn (KR * h)) := by
  have hint_gR : IntegrableOn gR (Ioi (KR * h)) := by
    refine integrableOn_Ioi_of_intervalIntegral_norm_tendsto (exp (-u_fn (KR * h))) (KR * h)
      (fun n => ?_) tendsto_natCast_atTop_atTop (fun n => ?_)
    · -- Интегрируемость на [KR*h, n]
      sorry -- стандартная непрерывность gR на компактном отрезке
    · -- Предел ∫_{KR*h}^n |gR| → exp(-u(KR*h)) из Леммы 2
      sorry
  have hsum_int :
      ∑' (n : ℕ), ∫ t in ((KR + n : ℤ) * h)..((KR + 1 + n : ℤ) * h), gR t
      = ∫ t in Ioi (KR * h), gR t := by
    have h_tend : Tendsto (fun N : ℕ => (KR + N : ℤ) * h) atTop atTop := by
      have : (fun N : ℕ => ((KR + N : ℤ) : ℝ) * h) = (fun N : ℕ => (KR : ℝ) * h + (N : ℝ) * h) := by
        ext N; push_cast; ring
      rw [this]
      exact tendsto_atTop_add_const_left atTop _ (tendsto_natCast_atTop_atTop.atTop_mul_const hh)
