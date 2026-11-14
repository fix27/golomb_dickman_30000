From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.

Module GolombDickman.

Local Open Scope ring_scope.

(* ========================================================================= *)
(* 1. Детерминантное тождество для β-дроби (Лемма 7, алгебраическая проверка) *)
(* ========================================================================= *)

(* a_n = ⌈n/2⌉ *)
Definition aCoeff (n : nat) : int := (n.+1 %/ 2)%:Z.

(* b_n(z) = 1 (n нечётное) или z = 'X (n чётное) *)
Definition bCoeff (n : nat) : {poly int} :=
  if (n %% 2 == 1)%N then 1 else 'X.

(* Индексация со сдвигом на 1: An 0 = A_{-1} = 1, An 1 = A_0 = z *)
Fixpoint An (n : nat) : {poly int} :=
  match n with
  | 0 => 1
  | 1 => 'X
  | n'.+2 => bCoeff n'.+1 * An n'.+1 + (aCoeff n'.+1)%:P * An n'
  end.

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

Fixpoint aProd (n : nat) : int :=
  match n with
  | 0 => 1
  | n'.+1 => aProd n' * aCoeff n'.+1
  end.

Lemma Ddet_succ n :
  Ddet n.+2 = - (aCoeff n.+1)%:P * Ddet n.+1.
Proof.
  rewrite /Ddet /An /Bn /=.
  by ring.
Qed.

Theorem Ddet_eq n :
  Ddet n.+1 = ((-1) ^+ n.+1 * aProd n)%:P.
Proof.
  elim: n => [| n ih].
  - by rewrite /Ddet /An /Bn /=; ring.
  - rewrite Ddet_succ ih /aProd /=.
    by rewrite exprS; ring.
Qed.

(* ========================================================================= *)
(* 2. Стягивание S-дроби в J-дробь (Лемма 7, шаг 3)                          *)
(* ========================================================================= *)

Section S2J.

Variable alpha : nat -> rat.

(* Числители и знаменатели S-дроби: Ps 0 = P_{-1} = 0, Ps 1 = P_0 = 1 *)
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

(* Коэффициенты J-дроби: b_0 = alpha_1, b_m = alpha_{2m} + alpha_{2m+1}, a_m = alpha_{2m-1} * alpha_{2m} *)
Definition jb (m : nat) : rat :=
  if m is m'.+1 then alpha (2 * m'.+1)%N + alpha (2 * m'.+1).+1 else alpha 1%N.

Definition ja (m : nat) : rat :=
  if m is m'.+1 then alpha (2 * m' + 1)%N * alpha (2 * m'.+1)%N else 1.

(* Числители и знаменатели J-дроби: Rs 0 = R_{-1} = 0, Rs 1 = R_0 = 1 *)
Fixpoint Rs (m : nat) : {poly rat} :=
  match m with
  | 0 => 0
  | 1 => 1
  | m'.+2 => (1 + (jb m'.+1)%:P * 'X) * Rs m'.+1 - (ja m'.+1)%:P * 'X^2 * Rs m'
  end.

(* Ss 0 = S_{-1} = 1, Ss 1 = S_0 = 1 + b_0 X *)
Fixpoint Ss (m : nat) : {poly rat} :=
  match m with
  | 0 => 1
  | 1 => 1 + (jb 0)%:P * 'X
  | m'.+2 => (1 + (jb m'.+1)%:P * 'X) * Ss m'.+1 - (ja m'.+1)%:P * 'X^2 * Ss m'
  end.

(* Трехшаговое тождество для P_n (аналог шагов (a), (b), (c) из доказательства) *)
Lemma Ps_step3 (m : nat) :
  Ps (2 * m + 4) =
    (1 + (jb m.+1)%:P * 'X) * Ps (2 * m + 2) - (ja m.+1)%:P * 'X^2 * Ps (2 * m).
Proof.
  have h1 : (2 * m + 4 = (2 * m + 2).+2)%N by lia.
  have h2 : (2 * m + 3 = (2 * m + 1).+2)%N by lia.
  have h3 : (2 * m + 2 = (2 * m).+2)%N by lia.
  rewrite h1 /= h2 /= h3 /=.
  rewrite /jb /ja !mulnS !addnS !addn0 !rmorphD !rmorphM /=.
  by ring.
Qed.

Lemma Qs_step3 (m : nat) :
  Qs (2 * m + 4) =
    (1 + (jb m.+1)%:P * 'X) * Qs (2 * m + 2) - (ja m.+1)%:P * 'X^2 * Qs (2 * m).
Proof.
  have h1 : (2 * m + 4 = (2 * m + 2).+2)%N by lia.
  have h2 : (2 * m + 3 = (2 * m + 1).+2)%N by lia.
  have h3 : (2 * m + 2 = (2 * m).+2)%N by lia.
  rewrite h1 /= h2 /= h3 /=.
  rewrite /jb /ja !mulnS !addnS !addn0 !rmorphD !rmorphM /=.
  by ring.
Qed.

(* Основная теорема: P_{2m+1} = R_m и Q_{2m+1} = S_m для всех m >= -1
   (в нашей сдвинутой на 1 индексации: Ps (2*m) = Rs m и Qs (2*m) = Ss m) *)
Theorem S_to_J_contraction (m : nat) :
  Ps (2 * m) = Rs m /\ Qs (2 * m) = Ss m.
Proof.
  elim/ltn_ind: m => [[| [| m]] ih].
  - (* m = 0 (индекс -1) *)
    by split.
  - (* m = 1 (индекс 0) *)
    split; rewrite /Ps /Qs /Rs /Ss /jb /=; ring.
  - (* m + 2 (индекс m + 1 >= 1) *)
    have [ihP1 ihQ1] := ih m.+1 (ltnSn _).
    have [ihP0 ihQ0] := ih m (ltnW (ltnSn _)).
    have hidx : (2 * m.+2 = 2 * m + 4)%N by lia.
    have hidx1 : (2 * m.+1 = 2 * m + 2)%N by lia.
    split.
    + rewrite hidx Ps_step3 -hidx1 ihP1 ihP0 /=.
      by ring.
    + rewrite hidx Qs_step3 -hidx1 ihQ1 ihQ0 /=.
      by ring.
Qed.

End S2J.

End GolombDickman.
