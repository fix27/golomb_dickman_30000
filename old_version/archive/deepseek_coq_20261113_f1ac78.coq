From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.

Module GolombDickman.

Local Open Scope ring_scope.

(* a_n = ⌈n/2⌉ *)
Definition aCoeff (n : nat) : int := (n.+1 %/ 2)%:Z.

(* b_n = 1 (n odd) или z (n even); работаем в {poly int}, X обозначает z *)
Definition bCoeff (n : nat) : {poly int} :=
  if (n %% 2 == 1)%N then 1 else 'X.

(* A_n(z), индексация: An 0 = A_{-1} = 1, An 1 = A_0 = z. *)
Fixpoint An (n : nat) : {poly int} :=
  match n with
  | 0 => 1
  | 1 => 'X
  | n'.+2 =>
      bCoeff n'.+1 * An n'.+1 + (aCoeff n'.+1)%:P * An n'
  end.

Fixpoint Bn (n : nat) : {poly int} :=
  match n with
  | 0 => 0
  | 1 => 1
  | n'.+2 =>
      bCoeff n'.+1 * Bn n'.+1 + (aCoeff n'.+1)%:P * Bn n'
  end.

(* D(k) = An(k)·Bn(k-1) − An(k-1)·Bn(k) *)
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

(* Шаг рекурсии для D: D(k+2) = -a_{k+1} · D(k+1). *)
Lemma Ddet_succ n :
  Ddet n.+2 = - (aCoeff n.+1)%:P * Ddet n.+1.
Proof.
  (* Раскрываем: bCoeff(n+1)*An(n+1) + a*An(n), аналогично Bn. *)
  rewrite /Ddet /An /Bn /=.
  (* Разбор чётности bCoeff(n+1); в каждом случае — полиномиальное тождество. *)
  case: (bCoeff n.+1 =P 1) => _; ring.
  case: (bCoeff n.+1 =P 'X) => _; ring.
Qed.

(* Детерминантное тождество: D(n+1) = (-1)^{n+1} · a_1 ⋯ a_n. *)
Theorem Ddet_eq n :
  Ddet n.+1 = ((-1) ^+ n.+1 * aProd n)%:P.
Proof.
  elim: n => [| n ih].
  - (* D(1) = X·0 − 1·1 = −1. *) by rewrite /Ddet /An /Bn /=; ring.
  - rewrite Ddet_succ ih /aProd /=.
    (* (−a)·((−1)^{n+1} aProd n) = (−1)^{n+2} (aProd n · a). *)
    by rewrite exprS; ring.
Qed.

End GolombDickman.