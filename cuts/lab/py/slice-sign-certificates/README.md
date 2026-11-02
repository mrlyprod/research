# slice-sign-certificates

- Certifies every numerical step of the even-half theorems at bases 3 and 5 of the paper The Slice Sign Law, Theorems 6.1 and 7.1, and regenerates the exact certificates, determinants and anchors the topic's claim rows print.
- The carry machine is built from the digit polynomial `P_D = A_q^(D-1)(A_q + D t^m)` alone; the column test `q beta_(K+1) < fill beta_K` and the row `(M^T)^t v` run on the reflection-folded core in exact integers.
- Every real constant is an `mpmath.iv` enclosure at 80 digits with outward rounding, asserted strictly below the decimal the paper prints; every check raises on failure.

## VERBS

- `base3`: `R`, the seven constants of Lemmas 7.5 and 7.6, the tails of Lemma 7.4 at `K = 1..7`; the 182 inequalities `(*)` at even `D = 38..400` with `K_1` built from enclosures of `log R`, `s_D` and `A` and asserted `>= K* + 1`; the Section 8 bounds `2E < 8D/(D+2)^2`, `T(D) < 1e-70`, `T < E/1000`, `wA < 1.1E` and the window margin; the majorant `F(400)` that closes every `D >= 400`; exact certificates at even `D = 2..36`.
- `base5`: Lemma 6.4 (`C_2`, `C_3`, the three branches, and the exact maxima of `Q_K(n)/Q_K(1)` over all residues at `K = 2, 3, 4` with each branch's maximum); Proposition 6.5 at `D = 18..32`, failing at `16`, at every even `D = 34..600` and on every depth range `K = 2..9` to `D = 10^6`; the constants of the proof of Theorem 6.1; the exact staircase to `D = 64`, depth 2 at every even `D = 16..320`, depth 3 at every even `D = 66..316`, the row-certificate depth `t(D)` at every even `D <= 400`, and the Perron masses at `D = 8`.
- `bases`: exact certificates at bases 7, 9, 11 (Fact 8.5), `V(L) > 0` towers off base 3, the sign law by exact Sturm counts in 82 cases, `min g_7`.
- `odd`: `sgn Delta_D = (-1)^D` at base 3, `D = 2..40`; `Delta_D != 0` in the odd classes `D = 1 mod q` to `D = 80` at bases 3 and 5 (Bareiss with a mod-prime cross-check); `v_2(det M_even)` by 2-adic elimination at base 5, even `D <= 156`, and base 3, `D = 4 mod 6`, `D <= 118`, against Bareiss at `D = 2..30` (base 3) and `D = 120, 130` (base 5).
- `anchors`: the `D = 3` rung, the column-sum orientation at `D = 3..50`, the threshold `3^(D-1)`, `V` at `D = 4, 6`, the first breaks of the `W_k` sign alternation, and the product formula of Theorem 5.1 at 60 digits.

## RESULTS

- `(*)` holds at all 182 even `D = 38..400`: slack `1.11528` at `D = 38`, `966.7` at `D = 60`, `4.357e42` at `D = 400`; `K_1 = 82, 90, 203, 9075` at `D = 38, 40, 60, 400`; the widest enclosure of `K*` is `1.8e-66` relative; `E(38) <= 0.0896632`.
- Past `D = 400`, with `K_1 = floor(K*) + 2`: `K_1 + 1 <= D^2/16` and `K_1 >= D` give `(D+2)(E + delta) <= F(D) + G(D)`, both decreasing for `D >= 18`, and `F(400) + G(400) < 1.02e-42`.
- Every printed constant of Lemmas 6.4, 7.4, 7.5, 7.6 is an upper bound: true values `0.76795411`, `0.32647214`, `0.26514803`, `0.7528156995`, `0.7052517912`, `0.6137009811`, `0.2266815969`, `0.8900158975`, `0.6696596527`, `0.6418791704`; `9^K LT(K) <= 0.20774`, `9^K (D-2) ST(D, K) <= 0.68539`.
- Base 5 `v_2(det M_even)` is `27` at `D = 120, 132` and `28` at `D = 130, 144`, so the bound `<= 26` to `D = 156` once printed is false; every value stays below `2(D-1) + v_2(D+4)`, so `Delta_D != 0` stands.
- Base 3, `D = 4 mod 6`, `D <= 118`: `v_2(det M_even) <= 17`, reading `1, 2, 3, 1, 4` from `D = 4`; base 5 reads `2, 1, 2, 3, 2` at `D = 4..12`.

## RUN

```
uv run python -I research/cuts/lab/py/slice-sign-certificates/certify.py
uv run python -I research/cuts/lab/py/slice-sign-certificates/certify.py base3
```

- Verbs run in the order given, all five without arguments; about 25 s on one core, `base5` and `odd` about 10 s each.

## WITNESSES

- The paper: the proofs of Theorems 6.1 and 7.1, Lemmas 6.4 and 7.4 to 7.6, Proposition 6.5, the Section 8 bounds, Facts 5.5, 6.7, 7.8 and 8.5.
- `claims/conjecture-s-even-half.md`: both theorem-entire rows, the certificate rows at bases 3, 5, 7, 9, 11, the 2-adic rows, and the refuted depth-2 and dip rows.
- `claims/conjecture-s-odd-half.md`, `claims/slice-sign-law-in-every-dimension.md`, `claims/diagonal-slice-ladder.md`, `claims/conjecture-s-the-even-half-transient.md`, `claims/slice-ladder-rate-and-dead-routes.md`: the rows that name this study.
