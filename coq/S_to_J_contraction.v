From mathcomp Require Import all_ssreflect all_algebra.
From mathcomp Require Import poly.
From Coq Require Import Ring.
From Coq Require Import Lia.
From GolombDickmanCoq Require Import Definitions Ps_step3 Qs_step3.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory Num.Theory.
Import GolombDickman.
Local Open Scope ring_scope.

(* ========================================================================= *)
(* Лемма 7, часть 2: стягивание S-дроби в J-дробь.                          *)
(* P_{2m} = R_m и Q_{2m} = S_m.                                             *)
(* ========================================================================= *)

Theorem S_to_J_contraction (alpha : nat -> rat) (m : nat) :
  Ps alpha (2 * m) = Rs alpha m /\ Qs alpha (2 * m) = Ss alpha m.
Proof.
  elim/ltn_ind: m => [[| [| m]] ih].
  - (* m = 0: Ps 0 = Rs 0 = 0, Qs 0 = Ss 0 = 1 *)
    by split.
  - (* m = 1: Ps 2 = Rs 1 = 1, Qs 2 = Ss 1 = 1 + alpha_1 X *)
    split; rewrite /Ps /Qs /Rs /Ss /jb /=; ring.
  - (* m + 2: индукционный шаг *)
    have [ihP1 ihQ1] := ih m.+1 (ltnSn _).
    have [ihP0 ihQ0] := ih m (ltnW (ltnSn _)).
    have hidx  : (2 * m.+2 = 2 * m + 4)%N by lia.
    have hidx1 : (2 * m.+1 = 2 * m + 2)%N by lia.
    split.
    + rewrite hidx Ps_step3 -hidx1 ihP1 ihP0 /=.
      by ring.
    + rewrite hidx Qs_step3 -hidx1 ihQ1 ihQ0 /=.
      by ring.
Qed.
