From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.
From Coq Require Import Ring.
From GolombDickmanCoq Require Import Definitions.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.
Import GolombDickman.
Local Open Scope ring_scope.

(* ========================================================================= *)
(* Лемма 7, часть 1: шаг определителя.                                      *)
(* D(k+2) = -a_{k+1} * D(k+1).                                              *)
(* ========================================================================= *)

Lemma Ddet_succ n :
  Ddet n.+2 = - (aCoeff n.+1)%:P * Ddet n.+1.
Proof.
  rewrite /Ddet /An /Bn /=.
  by ring.
Qed.
