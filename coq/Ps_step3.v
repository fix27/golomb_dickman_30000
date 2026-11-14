From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.
From Coq Require Import Ring.
From Coq Require Import Lia.
From GolombDickmanCoq Require Import Definitions.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.
Import GolombDickman.
Local Open Scope ring_scope.

(* ========================================================================= *)
(* Лемма 7, часть 2: трёхшаговое тождество для числителей S-дроби.          *)
(* P_{2m+4} = (1 + b_{m+1} X) P_{2m+2} - a_{m+1} X^2 P_{2m}.               *)
(* ========================================================================= *)

Lemma Ps_step3 (alpha : nat -> rat) (m : nat) :
  Ps alpha (2 * m + 4) =
    (1 + (jb alpha m.+1)%:P * 'X) * Ps alpha (2 * m + 2) -
    (ja alpha m.+1)%:P * 'X^2 * Ps alpha (2 * m).
Proof.
  have h1 : (2 * m + 4 = (2 * m + 2).+2)%N by lia.
  have h2 : (2 * m + 3 = (2 * m + 1).+2)%N by lia.
  have h3 : (2 * m + 2 = (2 * m).+2)%N by lia.
  rewrite h1 /= h2 /= h3 /=.
  rewrite /jb /ja !mulnS !addnS !addn0 !rmorphD !rmorphM /=.
  by ring.
Qed.
