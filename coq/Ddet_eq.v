From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.
From GolombDickmanCoq Require Import Definitions Ddet_succ.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.
Import GolombDickman.
Local Open Scope ring_scope.

(* ========================================================================= *)
(* Лемма 7, часть 1: детерминантное тождество.                              *)
(* D(n+1) = (-1)^{n+1} * a_1 * ... * a_n.                                   *)
(* ========================================================================= *)

Theorem Ddet_eq n :
  Ddet n.+1 = ((-1) ^+ n.+1 * aProd n)%:P.
Proof.
  elim: n => [| n ih].
  - (* D(1) = X*0 - 1*1 = -1 *)
    by rewrite /Ddet /An /Bn /=; ring.
  - rewrite Ddet_succ ih /aProd /=.
    (* (-a) * ((-1)^{n+1} aProd n) = (-1)^{n+2} (aProd n * a) *)
    by rewrite exprS; ring.
Qed.
