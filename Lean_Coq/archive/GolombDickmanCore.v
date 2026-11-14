From Coq Require Import Reals Ring Arith Lia.
Open Scope R_scope.

(** Golomb--Dickman: formally checked algebraic core.

    This file formalizes the local algebraic contraction step used in
    Lemma 7 of golomb_dickman_24.py.  It does not yet formalize the
    analytic input (Stieltjes convergence/alternation, Laguerre
    resolvents, strip estimates, or Arb interval arithmetic).
*)

Module GolombDickman.

(** One contraction step for the numerator or denominator of an S-fraction.

    The assumptions are three consecutive S-recurrence equations:
      q = t + c*x*r,
      s = q + b*x*t,
      n = s + a*x*q.

    The conclusion is the corresponding second-order J-recurrence.
*)
Theorem s_to_j_contraction_step :
  forall (x a b c r t q s n : R),
    q = t + c*x*r ->
    s = q + b*x*t ->
    n = s + a*x*q ->
    n = (1 + (b+a)*x)*q - (b*c)*x^2*r.
Proof.
  intros x a b c r t q s n Hq Hs Hn.
  rewrite Hn, Hs, Hq.
  ring.
Qed.

(** Coefficient identities for alpha_n = ceil(n/2), expressed
    in the even/odd indexing used in Lemma 7. *)
Theorem contracted_diagonal :
  forall m : nat, (m + (m+1))%nat = (2*m+1)%nat.
Proof. intros; lia. Qed.

Theorem contracted_off_diagonal :
  forall m : nat, (m*m)%nat = (m^2)%nat.
Proof. intros; nia. Qed.

(** Local determinant identity for the continuant recurrence.
    If A_n = b*A_(n-1)+a*A_(n-2) and likewise for B, then
    the adjacent determinant is multiplied by -a. *)
Theorem continuant_determinant_step :
  forall (a b Aprev Aprevprev Bprev Bprevprev : R),
    (b*Aprev + a*Aprevprev)*Bprev
      - Aprev*(b*Bprev + a*Bprevprev)
    = -a*(Aprev*Bprevprev - Aprevprev*Bprev).
Proof. intros; ring. Qed.

End GolombDickman.
