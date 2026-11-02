# gaussian-franel

- Franel's identity one field up: an exact formula for the Fourier `L^2` discrepancy of the Gaussian Farey set on the torus `C/Z[i]`, written as a gcd-weighted double sum over Gaussian Mertens sums.
- The node set is the one [spun-stack](../spun-stack/) lights: `w = u/d` in lowest terms in `Z[i]`, denominators one per associate class with `N(d) <= N`, numerators the reduced residues mod `d`, read modulo `Z[i]` in the unit square.
- Three results: the Gaussian Ramanujan sum in Mobius form (Theorem 1), the identity itself (Theorem 2), and the two-sided equivalence with the Riemann hypothesis for `zeta_{Q(i)}(s) = zeta(s) L(s, chi_-4)` (Theorem 3).
- Everything exact where the claim is exact: Gaussian Euclidean gcd, complete residue systems on the gcd box, sums of roots of unity reduced modulo the cyclotomic polynomial, and the identity itself an exact rational.

## THE OBJECT

- Pairing: `<lambda, w> = Re(lambda w)`, which is `p x - q y` for `lambda = p + qi` and `w = x + iy`. It is `Z[i]`-periodic in `w`, so `e_lambda(w) = e(<lambda, w>)` is a character of `C/Z[i]`, and `lambda -> e_lambda` is an isomorphism of `Z[i]` onto the dual group.
- The choice matters once. Under this pairing the character `e_lambda` is trivial on `(1/d) Z[i] / Z[i]` exactly when `d | lambda`; under the Euclidean pairing `Re(conj(lambda) w)` the same condition reads `conj(d) | lambda`. The two differ by relabelling the dual by conjugation, which fixes `N(lambda)`, so no printed number below depends on the choice.
- The node set: `G_N = {u/d mod Z[i] : [d] an associate class with N(d) <= N, u mod d, gcd(u, d) = 1}`, of size `m = sum_{[d], N(d) <= N} Phi(d)`, the Gaussian totient sum. This is the lit set of the exact spun stack.
- Against the complex Farey set of [arXiv:2407.04380](https://arxiv.org/abs/2407.04380), `G_T = {pr(p/q) : p, q in Z[i], 0 < |q| <= T}`: the two sets are equal, with `N = T^2`, the factor being 1. Letting `q` run over all four associates and `p` over all of `Z[i]` produces each torus point exactly once anyway, since a reduced denominator class `[d]` contributes the `Phi(d)` points `u/d` however its associates are listed. **Verified** by building `G_T` literally, `q` over every element of norm at most `T^2` and `p` over the box `[0, N(q))^2`, which covers a complete residue system because `N(q)` and `i N(q)` both lie in the ideal: the two sets agree at `T = 2, 3, 4, 5, 6` with `4, 24, 64, 176, 320` points (`check_sayous`).

## THE FUNCTIONAL

- There is no rank. Franel's `sum_j delta_j^2` needs the linear order of `F_Q` on `[0,1]` and `C/Z[i]` has none, so the twin is a choice and is named as one.
- The choice: `D_2(N)^2 = (1/m^2) sum_{lambda != 0} |S_N(lambda)|^2 / N(lambda)^2` with `S_N(lambda) = sum_{w in G_N} e(<lambda, w>)`, where `N(lambda) = |lambda|^2`.
- Why this one. Writing `nu_N` for the deviation of the empirical measure from Haar, the sum is `|| nu_N * K ||_2^2` for the kernel with `hat K(lambda) = N(lambda)^-1`, so `D_2(N) = 4 pi^2 || U_N ||_2` with `U_N` the periodic Newtonian potential of `nu_N`, `-Laplacian U_N = nu_N`. It is translation-invariant, invariant under the unit group of `Z[i]`, needs no fundamental domain and no anchored box, and it metrises weak-* convergence to Haar, so it is a Weyl criterion in `L^2`. The anchored `L^2` star-discrepancy fails all four.
- The classical twin of the same functional is Edwards' `I`. With `S(k) = sum_{rho in F_Q} e(k rho)` and `B1bar` the sawtooth, `int_0^1 (sum_rho B1bar(u + rho))^2 du = (1/(4 pi^2)) sum_{k != 0} |S(k)|^2 / k^2`, so the classical `m^2 D_2^2` is `4 pi^2 I`, and Edwards section 12.2 evaluates `I` as `sum_{a,b} M(Q/a) M(Q/b) gcd(a,b)^2 / (12 a b)`. Theorem 2 is that evaluation one field up.

## THEOREM 1, THE EXPONENTIAL SUM

- Partition `G_N` by reduced denominator class: `S_N(lambda) = sum_{[d], N(d) <= N} c_d(lambda)` with the Gaussian Ramanujan sum `c_d(lambda) = sum_{u mod d, gcd(u,d) = 1} e(<lambda, u/d>)`. Well defined: `u -> u + dt` moves `u/d` by `t in Z[i]` and `<lambda, t>` is an integer; replacing `d` by an associate permutes the reduced residues, so `c_d` depends only on the ideal.
- **Kluyver in `Z[i]`.** `c_d(lambda) = sum_{e | gcd(d, lambda)} mu_G(d/e) N(e)`, the sum over ideal divisors, `mu_G` the Mobius function on ideals of `Z[i]`, with the convention `gcd(d, 0) = d` so that `c_d(0) = Phi(d)`. **Proved.**
- Proof. Every residue `u mod d` has `gcd(u, d) = d/e` for a unique ideal `e | d`, and then `u/d = u'/e` with `u'` reduced mod `e`, so `sum_{e | d} c_e(lambda) = sum_{u mod d} e(<lambda, u/d>)`. The right side is the sum of a character over the finite abelian group `Z[i]/(d)`; that character is trivial iff `lambda/d in Z[i]`, tested on `u = 1` and `u = i`, so the sum is `N(d)` when `d | lambda` and `0` otherwise. Mobius inversion over the divisor lattice of the ideal `(d)`, which is a lattice of ideals because `Z[i]` is a principal ideal domain, gives the stated formula.
- **The Mertens form.** Writing `d = e f` and `M_G(x) = sum_{[f], N(f) <= x} mu_G(f)` for the Gaussian Mertens function over associate classes, `S_N(lambda) = sum_{[e] | lambda, N(e) <= N} N(e) M_G(N / N(e))`. **Proved**, by exchanging the two sums. At `lambda = 0` this is the zero mode `m = sum_{[e], N(e) <= N} N(e) M_G(N/N(e))`, and at any unit `lambda` it is `S_N(lambda) = M_G(N)`.
- **Verified** exactly at `N = 50`. All 2720 Gaussian Ramanujan sums, one for each of the 40 associate classes and each of the 68 nonzero `lambda` with `N(lambda) <= 20`, are computed by literal summation over the reduced residues as integer vectors of roots of unity reduced modulo `Phi_{N(d)}`, and every one is a rational integer equal to the Mobius formula: 0 mismatches. The 68 class sums equal the Mertens form: 0 mismatches. The 68 literal sums over all 672 nodes agree with both in floating point, worst error `1.281e-13` (`check_theorem_1`).
- **Verified**, the zero mode: `120`, `672`, `10608` at norm bounds `20`, `50`, `200`, read both as `sum Phi(d)` and as `sum N(e) M_G(N/N(e))` (`main`). The middle two are the node counts [spun-stack](../spun-stack/) prints.

## THEOREM 2, FRANEL ONE FIELD UP

- **The identity.** With `F(N) = sum_{[a],[b] : N(a), N(b) <= N} (N(gcd(a,b))^2 / (N(a) N(b))) M_G(N/N(a)) M_G(N/N(b))`,

```
m^2 D_2(N)^2 = 4 zeta_K(2) F(N),    zeta_K(2) = zeta(2) L(2, chi_-4) = zeta(2) * Catalan = 1.506703
```

- **Proved.** Substitute the Mertens form of Theorem 1, expand the square, and exchange: for fixed ideals `a, b` the inner sum is `sum_{lambda != 0, lcm(a,b) | lambda} N(lambda)^-2 = 4 zeta_K(2) / N(lcm(a,b))^2`, the 4 because every nonzero ideal has four generators. Then `N(a) N(b) / N(lcm(a,b))^2 = N(gcd(a,b))^2 / (N(a) N(b))`. The double sum is finite because `M_G(N/N(a))` vanishes for `N(a) > N`, so the right side is exact, and `F(N)` is an exact rational.
- The kernel is the gcd matrix `N(gcd(a,b))^2 / (N(a) N(b))`, the Gaussian twin of the `gcd(m,n)^2/(mn)` that is the layer Gram matrix of the parity carpet stack on [the stack page](../../../notes/levels-dilations-nodes.md). The same Smith kernel carries the moire correlation law and the Franel identity.
- **Verified** at `N = 20` and `N = 50` against the Fourier side truncated at `N(lambda) <= 200000`, the tail bounded by `m^2 (4 zeta_K(2) - sum_{0 < N(lambda) <= 200000} N(lambda)^-2)` since `|S_N(lambda)| <= m`. At `N = 20`, `m = 120`, `F(N) = 114917096/3663075 = 31.371756`, identity `189.071678` against truncated Fourier side `189.069795`, gap `0.001883` inside the tail bound `0.226192`. At `N = 50`, `m = 672`, `F(N) = 17870882021826065419/177236423132266875 = 100.830753`, identity `607.687997` against `607.677451`, gap `0.010546` inside the tail bound `7.093392` (`check_theorem_2`).
- The classical control is the same code path on the rational data. At `Q = 40` the direct Farey enumeration gives `m = 490` and `sum delta_v^2 = 0.0104270117`, and `C(Q) - 1 = 12 m sum delta_v^2` holds as an identity of exact rationals, which is Edwards section 12.2 regenerated here; the sieve route reads `C(40) = 62.310828579` against the exact `62.310828579` (`classical_exact`, `farey_delta_square`, `franel_form`). **Verified.**

## THEOREM 3, THE EQUIVALENCE

- **Statement.** `F(N) = O(N^{1+eps})` for every `eps > 0`, equivalently `m^2 D_2(N)^2 = O(N^{1+eps})`, equivalently `D_2(N) = O(N^{-3/2+eps})`, is equivalent to the Riemann hypothesis for `zeta_{Q(i)}(s) = zeta(s) L(s, chi_-4)`. That hypothesis is RH and GRH for `chi_-4` together and is strictly stronger than RH; the product is named every time it appears and is never shortened.
- The load-bearing input is `M_G(x) = O(x^{1/2+eps})` for every `eps > 0` if and only if `zeta_K` has no zero with `Re s > 1/2`. It is not called standard here; both directions are stated with their hypotheses.
- **The input, backward.** The bound gives the hypothesis. **Proved**, elementarily. Write `m(n) = sum_{N(a) = n} mu_G(a)`, so `M_G(x) = sum_{n <= x} m(n)` and `sum_n m(n) n^{-s} = 1/zeta_K(s)` on `Re s > 1`. Partial summation gives `sum_{n <= x} m(n) n^{-s} = M_G(x) x^{-s} + s int_1^x M_G(t) t^{-s-1} dt`, and under the bound both pieces converge as `x -> infinity` whenever `Re s > 1/2`. A Dirichlet series is analytic in its half plane of convergence, so `1/zeta_K` continues analytically to `Re s > 1/2`, and a zero of `zeta_K` there would be a pole of it. Nothing about `zeta_K` beyond its Euler product is used.
- **The input, forward.** The hypothesis gives the bound. **Proved**, by Littlewood's argument transcribed to `zeta_K`, with the two analytic inputs cited rather than assumed. Perron at half-integer `x` writes `M_G(x) = (1/(2 pi i)) int_{Re s = 1 + 1/log x} x^s / (s zeta_K(s)) ds` up to the standard truncation error, and the contour is pushed left to `Re s = 1/2 + eps` and closed on horizontal segments at heights where `zeta_K` is not small. Both inputs are proved for a Dedekind zeta in [arXiv:2109.06665](https://arxiv.org/abs/2109.06665), *On a Mertens-type conjecture for number fields*, Math. Proc. Cambridge Philos. Soc., read at source. Its Lemma 5.4 gives, under the hypothesis and for `|t| >= 1`, `|log zeta_K(s)| <= n_K (log(1/(1 - sigma)) + O((log tau)^{2 - 2 sigma} / ((1 - sigma) log log tau)))` on `1/2 + 1/log log tau <= sigma <= 1 - 1/log log tau` with `tau = |t| + 4`, which exponentiates to `1/zeta_K(sigma + it) << tau^delta` for every `delta > 0`, uniformly on `sigma >= 1/2 + eps`. Its Lemma 2.4 gives, under the same hypothesis, heights `T_n in [n, n+1)` with `|zeta_K(sigma + i T_n)| >= exp(-C log n / log log n)` for `-1 <= sigma <= 2`, which are the horizontal segments the contour closes on. Taking `T = x^2` makes each of the three pieces `O(x^{1/2 + eps + delta})`, so the bound holds at half-integers and hence everywhere, `M_G` moving between consecutive half-integers by the number of ideals of one norm, which is `O(x^eps)`.
- The biconditional is assembled here, not quoted. That paper states neither half as an equivalence; what it supplies is its `(1.1)`, `1/zeta_K(s) = s int_1^infinity M_K(x) x^{-s-1} dx`, which is the backward direction's Mellin identity one field up, its Lemma 5.4 and its Lemma 2.4, and the two paragraphs above are the assembly. The argument's shape is Edwards section 12.1, read at source, which supplies the completion Littlewood omitted; that section cites Titchmarsh Theorem 14.25(c) territory for the rational biconditional and is `K = Q` only.
- **Theorem 3, forward.** The hypothesis gives the decay. **Proved**, through the forward input above. Under `|M_G(x)| <= C_eps x^{1/2+eps}`, each term of Theorem 1 is at most `C_eps N^{1/2+eps} N(e)^{1/2}`, and only `e` with `N(e) <= min(N, N(lambda))` contribute, so `|S_N(lambda)| <= C_eps N^{1/2+eps} d(lambda) min(N, N(lambda))^{1/2}` with `d` the number of ideal divisors. Splitting the `lambda` sum at `N(lambda) = N` gives `sum_{N(lambda) <= N} d(lambda)^2 / N(lambda) = O(log^4 N)` and `N sum_{N(lambda) > N} d(lambda)^2 / N(lambda)^2 = O(log^3 N)`, so `m^2 D_2(N)^2 = O(N^{1+3eps})`.
- **Theorem 3, backward.** The decay gives the hypothesis. **Proved** outright, the backward input being elementary. In one line: the four units `lambda` have `N(lambda) = 1` and `S_N(lambda) = M_G(N)`, so dropping every other term of the Fourier sum gives `M_G(N)^2 <= zeta_K(2) F(N)`. Hence `F(N) = O(N^{1+eps})` forces `M_G(N) = O(N^{1/2+eps})`, and the backward input turns that into the absence of zeros right of the critical line.
- The two halves are not equally cheap. Backward needs nothing but the Euler product and partial summation; forward needs Perron and two conditional estimates for `zeta_K`, both read at source in the paper cited above. The backward direction is easier than Franel's because the functional is a Fourier sum by construction while `sum_j delta_j^2` is a rank statistic. What is not free is the identity: the content is Theorem 2, which turns a two-dimensional discrepancy into a finite gcd-weighted quadratic form in Gaussian Mertens sums, and the threshold `N^{1+eps}`, which is derived rather than inherited.
- The inequality of the backward direction is **Verified** at every rung of the meter, `M_G(N)^2 / (zeta_K(2) F(N))` never rising above `0.017` (`main`).
- What the equivalence is not. The Gaussian Farey set on `C/Z[i]` is not the rational Farey sequence carrying a character weight. Huxley's theorem for `lambda(q) = chi(q)` reaches the zeros of one Dirichlet `L`-function through weighted rational Farey points; this object is unweighted, two-dimensional, and invariant under the unit group, and it reaches `zeta` and `L(s, chi_-4)` together because `zeta_K` factors. Whether a number-field Farey sequence in the announced second part of that paper already carries a discrepancy statement is not settled here, that part being unread.

## WHAT IT PRINTS

- `check_sayous` at `T = 2, 3, 4, 5, 6`: the literal complex Farey set `G_T` and the node set at `N = T^2` are equal, with `4, 24, 64, 176, 320` points.
- `check_residues` to norm bound 200: 158 associate classes, `0` residue systems of the wrong size, `0` collisions mod `d`, `0` disagreements with the Gaussian totient. The residues of `d` are the box `[0, g) x [0, N(d)/g)` with `g = gcd(Re d, Im d)`, a fundamental domain of the ideal because its area is the index and its two sides are the Hermite normal form of the lattice.
- `check_theorem_1` at `N = 50` and `N(lambda) <= 20`: 672 nodes, 68 characters, 2720 exact Gaussian Ramanujan sums, `0` Mobius-formula mismatches, `0` Mertens-form mismatches, `0` literal node-sum mismatches at `1e-9`, worst literal error `1.281e-13`.
- `main`, the zero mode: `sum Phi(d)` and `sum N(e) M_G(N/N(e))` both read `120`, `672`, `10608` at norm bounds `20`, `50`, `200`.
- `readout`: `sum_{N(a) <= x} M_G(x/N(a)) = 1` at every `x` from 1 to 2000, the Gaussian twin of the collapsing global readout of the Mertens meter.
- `check_theorem_2` at `N = 20, 50`: the exact rational `F(N)`, the identity, the Fourier side truncated at `N(lambda) <= 200000`, the gap and the printed tail bound, with the gap inside the bound at both.
- `farey_delta_square` and `classical_exact` at `Q = 40`: `m = 490`, `sum delta_v^2 = 0.0104270117`, and `C(Q) - 1 = 12 m sum delta_v^2` as exact rationals.
- `franel_form` on the rational data at `Q = 125` to `8000`: `S2 * Q = (C(Q) - 1) Q / (12 Phi(Q))` reads `0.5395, 0.5848, 0.6241, 0.6387, 0.6560, 0.6538, 0.6564`, regenerating the Farey discrepancy table of [the Farey page](../../../../farey/notes/farey.md) with no Farey enumeration anywhere.
- `franel_form` on the Gaussian data at `N = 100` to `64000`: the meter table, `m`, `F(N)`, `F(N)/N`, `D_2(N)^2 N^3`, the local log-log slope of `F`, `M_G(N)`, the ratio `M_G(N)^2 / (zeta_K(2) F(N))`, and the classical `C(N)/N` beside it.

| `N` | `m` | `F(N)` | `F(N)/N` | `D_2(N)^2 N^3` | slope | `M_G(N)` | classical `C(N)/N` |
|---|---|---|---|---|---|---|---|
| 100 | 2600 | 212.1213 | 2.121213 | 189.1147 | - | -2 | 1.877964 |
| 250 | 16424 | 619.1589 | 2.476636 | 216.1484 | 1.1691 | 0 | 2.140099 |
| 500 | 65784 | 1278.4558 | 2.556912 | 222.5578 | 1.0460 | -3 | 2.282045 |
| 1000 | 260944 | 2725.4800 | 2.725480 | 241.2326 | 1.0921 | -1 | 2.332467 |
| 2000 | 1045088 | 5556.4599 | 2.778230 | 245.2845 | 1.0277 | -8 | 2.394800 |
| 4000 | 4176032 | 11394.1497 | 2.848537 | 252.0124 | 1.0361 | -1 | 2.385130 |
| 8000 | 16680488 | 23773.2448 | 2.971656 | 263.6505 | 1.0610 | -21 | 2.394592 |
| 16000 | 66694240 | 47603.5196 | 2.975220 | 264.1861 | 1.0017 | 18 | 2.442443 |
| 32000 | 266670328 | 96361.4734 | 3.011296 | 267.6034 | 1.0174 | 38 | 2.499486 |
| 64000 | 1067245288 | 193284.7209 | 3.020074 | 268.0999 | 1.0042 | -70 | 2.478573 |

- Reading the table: the local slope of `F` walks to `1.0042` and `F(N)/N` climbs slowly from `2.12` to `3.02`, which is the shape a bounded power of a logarithm has and is what the equivalence predicts; the classical column climbs the same way. No exponent is claimed beyond this window, which is 10 nested points and cannot separate `N^{1+eps}` from `N^{1.02}`. **Verified**, as a window.

## RUN

- `uv run python research/stack/lab/py/gaussian-franel/gaussian_franel.py`
- From the repository root. One core, under one second.
- Domain: the identification with the complex Farey set at `T = 2` to `6`, residue systems to norm bound 200, exact checks at norm bound 50 with `N(lambda) <= 20`, the identity at norm bounds 20 and 50 against the Fourier side to `N(lambda) <= 200000`, the classical control at `Q = 40` and the classical meter to `Q = 8000`, the Gaussian meter to `N = 64000`.
- Nothing is written to disk.

## WITNESSES

- Theorem 1: `ramanujan_exact` against `ramanujan_formula` and `sum_formula`, with `literal_sum` over the node set as the third route.
- Theorem 2: `franel_exact` against `lambda_side`, with the printed tail bound; `farey_delta_square` against `classical_exact` as the published control.
- Theorem 3: the meter columns of `main`, and `M_G(N)^2 / (zeta_K(2) F(N)) <= 1` at every rung.
- The node set and its count: `node_set` and `totient_class`, agreeing with the Gaussian totient sums of [spun-stack](../spun-stack/); `check_residues` for the residue systems and `check_sayous` for the identification with the complex Farey set of the literature.
