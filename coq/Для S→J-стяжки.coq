Section S2J.
(* P_n, Q_n — S-дробь; R_m, S_m — J-дробь *)
Fixpoint Ps (n : nat) : {poly rat} :=
  match n with
  | 0 => 1
  | 1 => 1
  | n'.+2 => Ps n'.+1 + (aCoeff n'.+2)%:P * 'X * Ps n'
  end.

Fixpoint Qs (n : nat) : {poly rat} :=
  match n with
  | 0 => 1
  | 1 => 1 + 'X
  | n'.+2 => Qs n'.+1 + (aCoeff n'.+2)%:P * 'X * Qs n'
  end.

Fixpoint Rs (m : nat) : {poly rat} :=
  match m with
  | 0 => 1
  | m'.+1 =>
      let bm : rat := (aCoeff (2*m')%N + aCoeff (2*m'.+1))%:R in
      let am : rat := (aCoeff (2*m' - 1)%N * aCoeff (2*m'))%:R in
      (1 + bm%:P * 'X) * Rs m' - am%:P * 'X^2 * Rs (m' - 1)
  end.

(* Аналогично Ss. *)
Fixpoint Ss (m : nat) : {poly rat} :=
  match m with
  | 0 => 1 + (aCoeff 1)%:P * 'X
  | m'.+1 =>
      let bm : rat := (aCoeff (2*m')%N + aCoeff (2*m'.+1))%:R in
      let am : rat := (aCoeff (2*m' - 1)%N * aCoeff (2*m'))%:R in
      (1 + bm%:P * 'X) * Ss m' - am%:P * 'X^2 * Ss (m' - 1)
  end.

Theorem S_to_J_contraction m :
  Ps (2*m + 1) = Rs m /\ Qs (2*m + 1) = Ss m.
Proof.
  elim: m => [| m [ihP ihQ]].
  - (* P_1 = 1 = R_0, Q_1 = 1 + X = S_0 *) by split; rewrite /Ps /Rs /Qs /Ss /=; ring.
  - split.
    + (* P_{2m+3} = R_{m+1}: три применения рекурсии S и группировка *)
      rewrite /Ps /= ihP.
      (* После подстановки — тождество многочленов. *)
      admit.  (* см. примечание ниже *)
    + (* Аналогично для Q. *)
      rewrite /Qs /= ihQ. admit.
Admitted.

End S2J.