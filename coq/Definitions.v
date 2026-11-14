From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.

Module GolombDickman.

Local Open Scope ring_scope.

(* ========================================================================= *)
(* Общие определения: коэффициенты, числители/знаменатели и определитель.    *)
(* ========================================================================= *)

(* a_n = ceil(n/2) *)
Definition aCoeff (n : nat) : int := (n.+1 %/ 2)%:Z.

(* b_n = 1 при нечётном n и z = 'X при чётном n *)
Definition bCoeff (n : nat) : {poly int} :=
  if (n %% 2 == 1)%N then 1 else 'X.

(* A_n: An 0 = A_{-1} = 1, An 1 = A_0 = z *)
Fixpoint An (n : nat) : {poly int} :=
  match n with
  | 0 => 1
  | 1 => 'X
  | n'.+2 => bCoeff n'.+1 * An n'.+1 + (aCoeff n'.+1)%:P * An n'
  end.

(* B_n: Bn 0 = B_{-1} = 0, Bn 1 = B_0 = 1 *)
Fixpoint Bn (n : nat) : {poly int} :=
  match n with
  | 0 => 0
  | 1 => 1
  | n'.+2 => bCoeff n'.+1 * Bn n'.+1 + (aCoeff n'.+1)%:P * Bn n'
  end.

(* D(k) = A_{k-1} B_{k-2} - A_{k-2} B_{k-1} *)
Definition Ddet (k : nat) : {poly int} :=
  match k with
  | 0 => 0
  | k'.+1 => An k'.+1 * Bn k' - An k' * Bn k'.+1
  end.

(* a_1 * ... * a_n *)
Fixpoint aProd (n : nat) : int :=
  match n with
  | 0 => 1
  | n'.+1 => aProd n' * aCoeff n'.+1
  end.

(* ========================================================================= *)
(* S- и J-дроби (Лемма 7, часть 2).                                          *)
(* ========================================================================= *)

Section S2JDefs.

Variable alpha : nat -> rat.

(* Ps 0 = P_{-1} = 0, Ps 1 = P_0 = 1 *)
Fixpoint Ps (n : nat) : {poly rat} :=
  match n with
  | 0 => 0
  | 1 => 1
  | n'.+2 => Ps n'.+1 + (alpha n'.+1)%:P * 'X * Ps n'
  end.

(* Qs 0 = Q_{-1} = 1, Qs 1 = Q_0 = 1 *)
Fixpoint Qs (n : nat) : {poly rat} :=
  match n with
  | 0 => 1
  | 1 => 1
  | n'.+2 => Qs n'.+1 + (alpha n'.+1)%:P * 'X * Qs n'
  end.

(* b_0 = alpha_1; b_m = alpha_{2m} + alpha_{2m+1} для m >= 1.
   Замечание: в исходном "Формальное доказательство.coq" здесь была
   ошибка сдвига (alpha_{2m+1} + alpha_{2m+2}); исправлено по
   lean/Lemma7S2J.lean. *)
Definition jb (m : nat) : rat :=
  if m is m'.+1
  then alpha (2 * (m'.+1))%N + alpha (2 * (m'.+1) + 1)%N
  else alpha 1%N.

(* a_0 = 1; a_m = alpha_{2m-1} * alpha_{2m} для m >= 1 *)
Definition ja (m : nat) : rat :=
  if m is m'.+1
  then alpha (2 * m' + 1)%N * alpha (2 * (m'.+1))%N
  else 1.

(* Rs 0 = R_{-1} = 0, Rs 1 = R_0 = 1 *)
Fixpoint Rs (m : nat) : {poly rat} :=
  match m with
  | 0 => 0
  | 1 => 1
  | m'.+2 =>
      (1 + (jb m'.+1)%:P * 'X) * Rs m'.+1 -
      (ja m'.+1)%:P * 'X^2 * Rs m'
  end.

(* Ss 0 = S_{-1} = 1, Ss 1 = S_0 = 1 + b_0 X *)
Fixpoint Ss (m : nat) : {poly rat} :=
  match m with
  | 0 => 1
  | 1 => 1 + (jb 0)%:P * 'X
  | m'.+2 =>
      (1 + (jb m'.+1)%:P * 'X) * Ss m'.+1 -
      (ja m'.+1)%:P * 'X^2 * Ss m'
  end.

End S2JDefs.

End GolombDickman.
