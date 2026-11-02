---
title: A Power Saving for Mobius Sums on Missing-Digit Integers, under GRH
lead: Under the generalized Riemann hypothesis the Mobius function summed over the integers that avoid one digit in base `q` is at most `A_F(x)^(1 - delta')` for every fixed `delta' < delta_q`, at every base `q >= 3690`, with `delta_q` explicit and tending to `1/4`; the proof is five elementary steps around one quoted exponential-sum bound, every bound proved on the page and the wall read off two closed forms, the same steps run under any zero-free half plane for Dirichlet `L`-functions at a base computed for each, and the `l^1` norm they pay is proved never to fall below one, so the method yields nothing unconditional and nothing below `q^(3/4)` digits. Sharper constants and certified exponents elsewhere on this tree have since moved the wall to `1499` and to `34`, and the unconditional dissection removes the hypothesis from `584`.
date: 2026-10-17
figure: paper-sparse-mertens-under-grh
---

Take the whole numbers whose base-`q` expansion never uses one chosen digit, and add up the Mobius function over the ones below a bound `x`. The trivial estimate for that sum is the number of terms, `A_F(x)`, of size `x^(alpha_q)` with `alpha_q = log(q-1)/log q`, and among the sources read here nothing better is in print. Assuming the generalized Riemann hypothesis for Dirichlet `L`-functions, this paper beats it by a fixed power at every base `q >= 3690` with one digit removed: `|M_F(x)| <<_(q,delta') A_F(x)^(1 - delta')` for every fixed `delta' < delta_q`, where `delta_q` is explicit and tends to `1/4` as the base grows. The proof expands the digit indicator into additive frequencies, pays one `l^1` norm for the expansion, and quotes one published exponential-sum bound at each frequency; the content is that past a computable base the norm costs less than the set's own mass. The same five steps run on any zero-free half plane for Dirichlet `L`-functions, at a base computed for each, and the floor of the method is proved rather than guessed: the `l^1` norm can never be smaller than one, so the argument yields nothing unconditional and nothing at all about digit sets of fewer than `q^(3/4)` digits. The wall `3690` is the wall of the digit-blind kernel bound proved here, the sign of one closed-form difference at that base; the Mobius page of this tree and two sibling papers have since moved it, and Section 8 says where.

## Introduction

![The 343 frequencies a/343 of base 7 with the digit 3 removed, at level 3, one vertical bar each across the square, read left to right with the frequency 0 at the centre; each bar's height is the cube root of |hat F_3(a/343)|/216, so the centre bar, where the transform is 216 = 6^3, fills the height, and the bars repeat one comb at the scales 1/7 and 1/49; the ink runs from blue at the short bars to yellow at the tall. The figure asserts that the 343 absolute values add to a total between the floor 343 of Theorem 1.4 and the proved bound (7 P_7(1))^3 = 18964.4 of Lemma 4.2, and that their squares add to 343 times 216.](paper-sparse-mertens-under-grh)

Fix a base `q` and throw away one digit. In base ten with the digit `7` removed, what is left is `1, 2, 3, 4, 5, 6, 8, 9, 10, ..., 16, 18, ..., 69, 80, 81, ...`, and the pattern continues forever, thinning as the numbers get longer: below `10^n` there are about `9^n` of them instead of `10^n`. Sets like this are called missing-digit or ellipsephic sets; on this tree a digit set at a base is a design, and the set of integers it spells is `S_F`. They are sparse, they are easy to describe, and they are hard to do arithmetic on, because knowing the digits of `n` tells you almost nothing about the primes dividing `n`.

Now put the Mobius function on them. `mu(n)` is `0` when `n` has a repeated prime factor, and otherwise `+1` or `-1` according to whether `n` has an even or an odd number of prime factors. Summing `mu` over all integers up to `x` gives Mertens' function, and how much it cancels is the Riemann hypothesis. Summing it over the missing-digit integers up to `x` gives `M_F(x) = sum_(n in S_F, n <= x) mu(n)`, and the question is the same one: does it cancel? The number of terms is `A_F(x)`, of size `x^(alpha_q)`, and `|M_F(x)| <= A_F(x)` is free. The figure is the object the proof pays for: the digit transform of base `7` missing `3` at level `3`, all `343` frequencies drawn once, a self-similar comb whose `l^1` mass is the whole cost of the method.

Write `F` for the surviving digits, `S_F` for the integers all of whose base-`q` digits lie in `F`, `A_F(x) = #{n in S_F : n <= x}` and `M_F(x) = sum_(n in S_F, n <= x) mu(n)`; `alpha_q = log_q(q-1)` is the mass exponent, so that `A_F(x)` has order `x^(alpha_q)`. The one-step constant `B_q(F)`, the proved constants `Phi_q` and `P_q(m)`, the cost exponent `c_q`, the saving `delta_q`, the certificate `gap_q`, the hypothesis `Z(a)` and its exponent `b(a)` are all defined in Section 2; `Z(1/2)` is the generalized Riemann hypothesis in its Dirichlet form.

**Theorem 1.1 (the power saving under the generalized Riemann hypothesis).** Assume `Z(1/2)`: no Dirichlet `L`-function has a zero in `Re s > 1/2`. Let `q >= 3690`, let `F` omit exactly one digit, and let `eps > 0`. Then for all `x >= 2`,

```
|M_F(x)| <<_(q,eps) x^(3/4 + c_q + eps) ,
```

with `c_q = log_q P_q(1)` the cost exponent of Definition 2.3. Moreover `3/4 + c_q < alpha_q` at every such `q`, so with `delta_q = (alpha_q - 3/4 - c_q)/alpha_q > 0` the bound reads: for every fixed `delta' < delta_q`,

```
|M_F(x)| <<_(q,delta') A_F(x)^(1 - delta') .
```

The implied constants depend on `q` and on `eps`, respectively on `delta'`, and on nothing else; no uniformity in `q` is claimed. The endpoint `A_F(x)^(1 - delta_q)` is not proved: `delta_q` is fixed before `eps` is chosen.

**Corollary 1.2 (several avoided digits).** Let `E` be the set of avoided digits, `|E| = m >= 1`, and suppose `gap_q(m) = (q - m) q^(-3/4) - P_q(m) > 0`. Then under `Z(1/2)` the same five steps give, for every `eps > 0` and all `x >= 2`, `|M_F(x)| <<_(q,eps) x^(3/4 + c_q(m) + eps)`, and hence `|M_F(x)| <<_(q,delta') A_F(x)^(1 - delta')` for every fixed `delta' < delta_q`, with `c_q(m) = log_q P_q(m)` and `delta_q` formed from it and from `alpha_q = log_q(q - m)` as in Definition 2.4.

**Theorem 1.3 (the same saving on any zero-free half plane).** Fix `a` with `1/2 <= a < 1` and assume `Z(a)`: no Dirichlet `L`-function has a zero in `Re s > a`. Let `F` omit exactly one digit and let `q >= q_0(a)`, the least base with `gap_q(b(a), 1) > 0`, where `b(a)` is the exponent of Definition 2.5. Then for every `eps > 0` and all `x >= 2`, `|M_F(x)| <<_(q,eps) x^(b(a) + c_q + eps)`, and `b(a) + c_q < alpha_q`, so `|M_F(x)| <<_(q,delta') A_F(x)^(1 - delta')` for every fixed `delta' < delta_q(a) = (alpha_q - b(a) - c_q)/alpha_q > 0`. The least base `q_0(a)` exists at every `a`, and `gap_q(b(a), 1) > 0` at every `q >= q_0(a)`, so the hypothesis is a half line and not a window (Proposition 5.2). With `m` avoided digits the same chain runs whenever `gap_q(b(a), m) > 0`. All the caveats of Theorem 1.1 are unchanged.

**Theorem 1.4 (the `l^1` floor).** For every base `q >= 2`, every digit set `F` with `k = |F| >= 1` and every real `t`,

```
sum_(r mod q) |hat F((t + r)/q)|^2 = qk     and hence     sum_(r mod q) |hat F((t + r)/q)| >= q .
```

Consequently `B_q(F) >= q`, the recursion of Lemma 4.2 never contracts, and the normalised `l^1` mass `q^(-L) Sigma_L` is at least `1` at every level `L`.

The proof of Theorem 1.1 in one paragraph. The indicator of a length-`L` digit string over `F` expands into additive characters modulo `q^L`, and because digits are independent the Fourier coefficient factors over digit positions: `hat F_L(t) = prod_(j < L) hat F(q^j t)` with `hat F(t) = sum_(d in F) e(dt)`. Expanding `M_F(x)` this way turns it into `q^(-L) sum_(a mod q^L) hat F_L(a/q^L) sum_(n <= x) mu(n) e(-na/q^L)`, so the whole sum is bounded by the `l^1` norm of the transform times the largest Mobius exponential sum. Under the generalized Riemann hypothesis that largest sum is `x^(3/4 + eps)`, uniformly in the frequency, by a theorem of Baker and Harman (1991). The `l^1` norm obeys a one-step recursion: peeling off the last digit costs a single constant `B_q(F) = sup_t sum_(r mod q) |hat F((t + r)/q)|` per digit, so the normalised norm is at most `(B_q(F)/q)^L <= P_q(1)^L = q^(L c_q)`. What remains is to bound `B_q(F)` explicitly, and that is elementary trigonometry on a grid of `q` equally spaced points: the Dirichlet kernel contributes at most `Phi_q`, about `(2/pi) q log q`, and the one missing digit contributes at most `q` by an exact Parseval identity. Then `c_q = log_q(1 + Phi_q/q)` is of size `log log q/log q`, and the arithmetic of Section 6 does the rest: the cost `3/4 + c_q` falls below the mass `alpha_q` once, at `q = 3690`, and never comes back.

Nothing in that chain except the quoted Mobius bound depends on the exponent `3/4`. So the same five steps run on any hypothesis of the form "no Dirichlet `L`-function vanishes in `Re s > a`", with the exponent `b(a)` that Baker and Harman (1991) and Zhang (2024) attach to it, and the price of a weaker hypothesis is paid entirely in the base: Theorem 1.3 and the rungs of Fact 6.4. A zero-free half plane at `Re s > 19/20` still buys the saving, at a base with `84` decimal digits.

Where the wall now stands. The base `3690` is the least base at which the digit-blind kernel bound `B_q(F) <= q P_q(1)` of Lemma 4.4 closes the gap, and nothing in this paper is quoted below it. [The Mobius page](../notes/grh-power-saving.md) of this tree proves a sharper one-step bound at one avoided digit, exact on the shifted grid and read through a chord of the cosecant, which closes the same chain from `1499` at every avoided digit and from `1032` at the end digits `0` and `q-1`, each wall an up-set of an exhaustive scan and a half line by the floor proved here; those constants are proved there and not here, and none of them touches the `m >= 2` corollary, the ladder or the floor. [The First Base Below a Quarter](first-base-below-a-quarter.md) replaces the one-step constant by the certified `l^1` exponent of the transform itself and closes the chain, under the same hypothesis, at base `21` missing `0` and at every base from `34` at every digit, from `125` by a proof. [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md) removes the hypothesis altogether from base `584`, and from `115` by certificate, at the price the floor of Section 7 predicts: the saving falls from a power of `A_F(x)` to `exp(-c sqrt(log x))`. [An Unconditional Mertens Bound on the Half Interval](half-interval-mobius.md) runs the certified chain on a different design, half the digits kept as an interval, closing under the hypothesis at every odd base from `3789`. What stands in this paper is the chain itself with every bound proved on the page and its wall read off two closed forms at one base, at every number `m` of avoided digits, the ladder over every zero-free half plane, which none of the three carries, the proved monotone floor that makes every wall a half line without a scan, and the `l^1` floor with its three consequences.

What this does not do, and cannot. The method has a floor, and Theorem 1.4 proves it: for every digit set, the shifted-grid energy is exactly `qk`, so `B_q(F) >= q`, while `P_q(m) >= 1` gives `c_q(m) >= 0`. The `l^1` recursion therefore never contracts, and three consequences follow at once (Corollary 7.1). The decomposition needs `alpha_q > 3/4`, that is more than `q^(3/4)` surviving digits, so it says nothing whatever about a set like `{0, 1}` in base `3`. Put Davenport's unconditional bound into the quoted step instead of the conditional one and the answer exceeds `A_F(x)`: no unconditional theorem hides inside this argument. And rearranging the same two factors by Cauchy-Schwarz is worse than trivial. That is the honest picture: a dense-column theorem, conditional, with a proved wall on both sides of it; Section 9 lists what it leaves undone.

## Definitions

Throughout, `e(t) = exp(2 pi i t)`, `||t||` is the distance from `t` to the nearest integer, `gamma = 0.5772156649...` is Euler's constant, `mu` is the Mobius function with `mu(0) = 0`, and `log` is the natural logarithm, `log_q` the logarithm to base `q`. Implied constants are made explicit where they matter and are otherwise absolute unless a subscript says otherwise.

**Definition 2.1 (the set and its mass).** Fix an integer base `q >= 3` and a set `E` of avoided digits, a subset of `{0, 1, ..., q-1}` with `m = |E| >= 1`. Put `F = {0, 1, ..., q-1}` less `E` and `k = |F| = q - m`; `q` is the base and `k` the fill, and the two letters keep these meanings throughout. `S_F` is the set of integers `n >= 1` whose every base-`q` digit lies in `F`, `A_F(x) = #{n in S_F : n <= x}` its counting function and `M_F(x) = sum_(n in S_F, n <= x) mu(n)` the Mobius sum over it. The mass exponent is `alpha_q = log_q k`, the dimension of the design. Counting digit strings gives `A_F(q^L) >= k^L - 1` at every `L >= 1`, and hence `A_F(x) >>_q x^(alpha_q)` for all `x >= q`: with `L = floor(log_q x) >= 1`, `A_F(x) >= k^L - 1 >= k^L/2 = q^(L alpha_q)/2 >= x^(alpha_q)/(2q)`. So `alpha_q` is the exponent of the set.

**Definition 2.2 (the digit transform and its level form).** The digit transform of `F` is `hat F(t) = sum_(d in F) e(dt)`, that of the avoided digits is `hat E(t) = sum_(d in E) e(dt)`, and `D_q(t) = sum_(d < q) e(dt)` is the full Dirichlet kernel, so that `hat F = D_q - hat E`, `|D_q(t)| = |sin(pi q t)/sin(pi t)|` for `t` not an integer and `D_q(t) = q` at the integers. Let `D_L` be the set of integers `sum_(j < L) f_j q^j` with every `f_j in F`, the values of the length-`L` digit strings over `F`, leading zeros allowed; it has `k^L` members, all in `[0, q^L)`. The level-`L` transform is

```
hat F_L(t) = sum_(n in D_L) e(nt) = prod_(j < L) hat F(q^j t) ,
```

the product formula holding because the digits of a string are chosen independently (Lemma 4.1). Its `l^1` mass on the grid of level `L` is `Sigma_L = sum_(0 <= a < q^L) |hat F_L(a/q^L)|`, with `Sigma_0 = 1`.

**Definition 2.3 (the one-step constant and the proved bound).** The one-step constant of the shifted-grid recursion is

```
B_q(F) = sup_(t in R) sum_(r mod q) |hat F((t + r)/q)| ,
```

the largest `l^1` mass of one digit over a shifted grid of `q` points. Write `H(n) = log n + gamma + 1/(2n)`, an upper bound for the `n`-th harmonic number `H_n`, since `H_n - log n - 1/(2n)` is nondecreasing by the trapezoid rule on the convex `1/z` and tends to `gamma`. Set

```
Phi_q = (4/pi) q + (2q/pi) H(ceil((q-2)/2)) + (1 - 2/pi)(q - 2) + 0.727 ,     P_q(m) = sqrt(m) + Phi_q/q .
```

Lemma 4.4 proves `B_q(F) <= q P_q(m)` for every `F` with `|E| = m`. The cost exponent is `c_q(m) = log_q P_q(m)`, written `c_q` when `m = 1`, defined from the proved bound `P_q(m)` and never from the exact supremum `B_q(F)`. Note `Phi_q = (2/pi) q log q (1 + o(1))`, so `c_q = (log log q + log(2/pi) + o(1))/log q`.

**Definition 2.4 (the saving and the certificate).** For `m` avoided digits put

```
delta_q = (alpha_q - 3/4 - c_q(m))/alpha_q ,     gap_q(m) = k q^(-3/4) - P_q(m) .
```

Since `alpha_q - 3/4 - c_q(m) = log_q(k q^(-3/4)/P_q(m))`, the two tests agree: `3/4 + c_q(m) < alpha_q` if and only if `gap_q(m) > 0`, and

```
delta_q = log(1 + gap_q(m)/P_q(m)) / (alpha_q log q) .     (2.1)
```

Form (2.1) never differences two numbers of size `1` to produce one of size `10^(-6)`, and it is the form every printed `delta_q` margin of Section 6 is computed from; the `gap_q` margins are the difference itself, which is what they are.

**Definition 2.5 (the zero-free hypothesis and the exponent it buys).** For `1/2 <= a < 1`, let `Z(a)` be the hypothesis that `L(s, chi)` has no zero in the half plane `Re s > a`, for every Dirichlet character `chi` of every modulus. Thus `Z(1/2)` is the generalized Riemann hypothesis in its Dirichlet form. Let

```
b_BH(a) = a + 1/4          on 1/2 <= a < 11/20 ,
b_BH(a) = 4/5              on 11/20 <= a < 3/5 ,
b_BH(a) = (a + 1)/2        on 3/5 <= a < 1 ,
b_Z(a) = (8a - 7a^2)/(4 - 2a)   on 1/2 <= a <= 4/7 ,   b_Z(a) = +infinity outside that range ,
b(a) = min(b_BH(a), b_Z(a)) .
```

These are the exponents of Lemmas 3.1 and 3.2; `b(1/2) = 3/4`, and `b(a) >= 3/4` throughout (Lemma 5.1). For `b in [3/4, 1)` and `m` avoided digits set `gap_q(b, m) = k q^(-b) - P_q(m)`, so that `gap_q(m) = gap_q(3/4, m)`; let `q_0(a)` be the least `q >= 3` with `gap_q(b(a), 1) > 0`, and let `Q(b)` be the least `q >= 40` with `(1 - b)(q - 2)(q + 1)^(-b) >= 1.291`. Proposition 5.2 shows that `q_0(a)` exists and exceeds `Q(b(a))` at every `a`.

## The quoted inputs

Exactly one analytic input enters the proof of Theorem 1.1, and it enters once. It is stated here, with its own hypothesis, in the form in which it is used and no stronger; nothing below re-proves it or enlarges it. Write

```
S(x, theta) = sum_(n <= x) mu(n) e(n theta) ,     T(x) = max_theta |S(x, theta)| ,
```

the maximum over all real `theta`.

**Lemma 3.1 (Baker and Harman 1991, quoted).** Let `1/2 <= a < 1` and assume `Z(a)`. Then for every `eps > 0`, `T(x) <<_(a,eps) x^(b_BH(a) + eps)`, with `b_BH` as in Definition 2.5 and the implied constant depending at most on `a` and `eps`. In particular, under `Z(1/2)` the exponent is `3/4 + eps`. This is the Theorem of pp. 193-194 of that paper, read at source; it is restated at source in Porritt (2018) and Zhang (2024), and the frequencies this paper feeds it are the `a/q^L`, well inside its uniformity.

**Lemma 3.2 (Zhang 2024, Theorem 1.1, quoted).** Let `1/2 <= a <= 4/7` and assume `Z(a)`. Then for every `eps > 0`, `T(x) <<_(a,eps) x^(b_Z(a) + eps)` with `b_Z(a) = (8a - 7a^2)/(4 - 2a)`. At `a = 1/2` this is again `3/4`, so the endpoint of the generalized Riemann hypothesis is not moved. Zhang states the hypothesis for the characters of one fixed modulus `l`; `Z(a)`, the hypothesis carried through this paper, asks the half plane of every character of every modulus and so implies it at every `l`.

**Lemma 3.3 (Davenport 1937, quoted as recorded in Porritt 2018, pp. 298-299).** Unconditionally, for every `A > 0`, `T(x) <<_A x (log x)^(-A)`.

Lemma 3.3 is used only in Section 7, to prove that this decomposition yields nothing unconditional. The conjectured truth in place of `3/4` is `1/2`, also recorded in Porritt (2018); Section 7 shows that even that would not reach the sparse range.

## The five steps

The proof of Theorem 1.1 is five steps. Four of them are proved here from the definitions; the fifth is Lemma 3.1, used once and exactly as stated. None of the four names an exponent, which is why Theorem 1.3 costs nothing extra.

**Lemma 4.1 (orthogonality and the product formula).** For every `L >= 1` and every integer `n` with `0 <= n < q^L`,

```
1_(D_L)(n) = q^(-L) sum_(0 <= a < q^L) hat F_L(a/q^L) e(-na/q^L) ,
```

and `hat F_L(t) = prod_(j < L) hat F(q^j t)` for every real `t`.

**Proof.** The additive characters `n -> e(na/q^L)`, `a` running modulo `q^L`, satisfy `sum_(a mod q^L) e(na/q^L) = q^L 1_(q^L | n)`. Hence for `0 <= n < q^L`,

```
q^(-L) sum_(a mod q^L) hat F_L(a/q^L) e(-na/q^L) = q^(-L) sum_(v in D_L) sum_(a mod q^L) e((v - n) a/q^L) = #{v in D_L : v = n mod q^L} ,
```

and `D_L` lies in `[0, q^L)`, so the count is `1_(D_L)(n)`. For the product formula, every `v in D_L` is uniquely `sum_(j < L) f_j q^j` with `f_j in F`, the choices being independent, so

```
hat F_L(t) = sum_(f_0 in F) ... sum_(f_(L-1) in F) prod_(j < L) e(f_j q^j t) = prod_(j < L) sum_(f in F) e(f q^j t) = prod_(j < L) hat F(q^j t) .  □
```

**Lemma 4.2 (one constant per digit).** `Sigma_L <= B_q(F) Sigma_(L-1)` for every `L >= 1`, hence `Sigma_L <= B_q(F)^L` and

```
q^(-L) Sigma_L <= (B_q(F)/q)^L .
```

**Proof.** Peel the digit at position `0`: the product formula gives `hat F_L(t) = hat F(t) hat F_(L-1)(qt)`. Write `a = a' + s q^(L-1)` with `0 <= a' < q^(L-1)` and `0 <= s < q`, and put `t' = a'/q^(L-1)`. Then `q (a/q^L) = t' + s`, and `hat F_(L-1)` has period `1`, so `hat F_(L-1)(q a/q^L) = hat F_(L-1)(t')`; and `a/q^L = (t' + s)/q`. Therefore

```
Sigma_L = sum_(0 <= a' < q^(L-1)) |hat F_(L-1)(t')| sum_(s mod q) |hat F((t' + s)/q)| <= B_q(F) sum_(0 <= a' < q^(L-1)) |hat F_(L-1)(t')| = B_q(F) Sigma_(L-1) ,
```

by Definition 2.3. With `Sigma_0 = 1` this gives `Sigma_L <= B_q(F)^L`. □

**Lemma 4.3 (two inequalities).** For `0 <= v <= 1`, `sin(pi v) <= 4v(1 - v)`. For `0 < z <= pi/2`, `1/sin z <= 1/z + 1 - 2/pi`.

**Proof.** Let `f(v) = 4v(1 - v) - sin(pi v)`, so `f(0) = f(1/2) = f(1) = 0`, `f(v) = f(1 - v)` and `f'(1/2) = 4 - 4 - pi cos(pi/2) = 0`. Since `f''(v) = pi^2 sin(pi v) - 8`, `f` is concave where `sin(pi v) < 8/pi^2` and convex where `sin(pi v) > 8/pi^2`; on `[0, 1/2]` that is one concave piece `[0, v_0]` followed by one convex piece `[v_0, 1/2]`, with `sin(pi v_0) = 8/pi^2`. On the convex piece `f'` is nondecreasing and `f'(1/2) = 0`, so `f' <= 0` there, `f` decreases to `f(1/2) = 0`, and `f >= 0` on `[v_0, 1/2]`, in particular `f(v_0) >= 0`. On the concave piece `f` lies above the chord from `(0, 0)` to `(v_0, f(v_0))`, which is nonnegative, so `f >= 0` on `[0, v_0]`. Symmetry gives `[1/2, 1]`.

For the second, let `h(z) = 1/sin z - 1/z` on `(0, pi/2]`. Then `h'(z) = 1/z^2 - cos z/sin^2 z`, so `h' >= 0` exactly when `sin^2 z >= z^2 cos z`. Put `s = sin(z/2)`, so `cos z = 1 - 2s^2 >= 0` and `sin^2 z = 4s^2(1 - s^2)` on this range. From `z/2 <= tan(z/2) = s/sqrt(1 - s^2)` we get `z^2 <= 4s^2/(1 - s^2)`, whence

```
sin^2 z - z^2 cos z >= 4s^2(1 - s^2) - 4s^2(1 - 2s^2)/(1 - s^2) = 4s^2 [(1 - s^2)^2 - (1 - 2s^2)]/(1 - s^2) = 4s^6/(1 - s^2) >= 0 .
```

So `h` increases on `(0, pi/2]` and `h(z) <= h(pi/2) = 1 - 2/pi`. □

**Lemma 4.4 (the kernel bound).** For every `q >= 3` and every `E` with `|E| = m`,

```
sup_t sum_(r mod q) |D_q((t + r)/q)| <= Phi_q ,     sum_(r mod q) |hat E((t + r)/q)| <= q sqrt(m)   for every t ,
```

and hence `B_q(F) <= Phi_q + q sqrt(m) = q P_q(m)`.

**Proof.** Fix `t` and put `v = {t}/q in [0, 1/q)`, so that the `q` points `u_r = (t + r)/q` are, modulo `1`, exactly `v + r/q` for `r = 0, ..., q-1`: a grid of spacing `1/q`. Because `q u_r = t + r`, `|sin(pi q u_r)| = |sin(pi t)|` for every `r`; write `S = |sin(pi t)| = sin(pi w)` with `w = qv in [0, 1)`. Writing `d_r = ||u_r||`, `|D_q(u_r)| = S/sin(pi d_r)` when `d_r > 0`, and `|D_q(u_r)| = q` when `d_r = 0`, which happens only if `w = 0`; then `S = 0`, every other term vanishes, and the whole sum is `q <= Phi_q`. Assume `w in (0, 1)` from here.

The points `v + r/q <= 1/2` sit at distance `(w + r)/q` from `0` and the rest at distance `(q - r - w)/q` from `1`, so the multiset of distances `d_r` is `{(w + j)/q}` and `{(j + 1 - w)/q}` for `j = 0, 1, 2, ...`, truncated to `q` values in all. The two values at `j = 0` are `w/q` and `(1 - w)/q`, both below `1/q <= 1/3`. By Lemma 4.3 with `z = pi w/q` and `z = pi(1 - w)/q`, those two contribute

```
S (1/sin(pi w/q) + 1/sin(pi(1 - w)/q)) <= sin(pi w) (q/(pi w) + q/(pi(1 - w)) + 2 - 4/pi) <= (4/pi) q + 2 - 4/pi ,
```

using `sin(pi w) <= 4w(1 - w)` on the first two terms and `sin(pi w) <= 1` on the third; and `2 - 4/pi = 0.7267... < 0.727`.

Each of the remaining `q - 2` points has `d_r >= j/q` for its own index `j >= 1`, and `d_r <= 1/2` always, so `pi d_r <= pi/2` and Lemma 4.3 gives `S/sin(pi d_r) <= 1/(pi d_r) + 1 - 2/pi <= q/(pi j) + 1 - 2/pi`. On each side the indices run `1, 2, ...` without a gap, and at most `n = ceil((q-2)/2) = ceil(q/2) - 1` of them occur on either side: the `r` with `v + r/q <= 1/2` other than `r = 0` satisfy `1 <= r < q/2`, and the `r` with `v + r/q > 1/2` other than `r = q - 1` satisfy `q/2 - 1 < r <= q - 2`, each a range of at most `ceil(q/2) - 1` integers. Summing and using `sum_(j <= n) 1/j <= H(n)`,

```
sum_(remaining) |D_q(u_r)| <= (2q/pi) H(n) + (q - 2)(1 - 2/pi) .
```

Adding the two bounds gives `Phi_q`.

For the avoided part, expand the square and sum over the grid: since the elements of `E` are distinct residues modulo `q`,

```
sum_(r mod q) |hat E((t + r)/q)|^2 = sum_(d, d' in E) e((d - d') t/q) sum_(r mod q) e((d - d') r/q) = q #{(d, d') in E x E : d = d'} = qm ,
```

an exact identity, not an estimate. Cauchy-Schwarz over the `q` points turns it into `sum_r |hat E((t + r)/q)| <= sqrt(q) sqrt(qm) = q sqrt(m)`. Finally `|hat F| = |D_q - hat E| <= |D_q| + |hat E|` pointwise, so `B_q(F) <= Phi_q + q sqrt(m)`. □

**Proof of Theorem 1.1 and Corollary 1.2.** Let `beta` denote the exponent supplied by the quoted input, so `beta = 3/4` here and `beta = b(a)` in Theorem 1.3; nothing in this proof uses the value of `beta` beyond `beta < 1`. Put `delta = (alpha_q - beta - c_q(m))/alpha_q`, which is the `delta_q` of Definition 2.4 at `beta = 3/4`. Fix `eps > 0`, let `x >= 2`, and set `L = ceil(log_q(x + 1))`, so that `x < q^L <= q(x + 1)`. By Lemma 3.1, `T(y) <<_eps y^(beta + eps)`.

Suppose first `0 in F`. Then an integer `n` with `1 <= n < q^L` lies in `S_F` if and only if `n in D_L`, since padding with leading zeros is allowed. Hence by Lemma 4.1,

```
M_F(x) = sum_(1 <= n <= x) mu(n) 1_(D_L)(n) = q^(-L) sum_(0 <= a < q^L) hat F_L(a/q^L) S(x, -a/q^L) ,
```

the interchange being a finite sum. Taking absolute values and applying Lemma 4.2 and then Lemma 4.4,

```
|M_F(x)| <= q^(-L) Sigma_L T(x) <= (B_q(F)/q)^L T(x) <= P_q(m)^L T(x) = q^(L c_q(m)) T(x) .
```

Since `q^L <= q(x + 1) <= 3qx` and `c_q(m) >= 0`, because `P_q(m) >= sqrt(m) >= 1`, we have `q^(L c_q(m)) <= (3qx)^(c_q(m)) <<_q x^(c_q(m))`, so `|M_F(x)| <<_(q,eps) x^(beta + c_q(m) + eps)`.

Suppose instead `0` is not in `F`. Then every digit string over `F` has a nonzero leading digit, so `D_l` is exactly the set of members of `S_F` with precisely `l` digits, and these blocks are disjoint. For `l < L` the block `D_l` lies in `[1, q^l)`, and the display above run at level `l` with `q^l - 1` in place of `x` gives `|sum_(n in D_l) mu(n)| <= q^(l c_q(m)) T(q^l - 1) <<_eps q^(l(beta + c_q(m) + eps))`; the top block, cut at `x`, gives `q^(L c_q(m)) T(x)`. Summing over `l <= L` is a geometric series of ratio `q^(beta + c_q(m) + eps) > 1`, so the total is `<<_(q,eps) q^(L(beta + c_q(m) + eps)) <<_(q,eps) x^(beta + c_q(m) + eps)` again. In both cases the finitely many `x` in any bounded range are absorbed by the implied constant.

It remains to convert to the `A_F` yardstick, and that is where the hypothesis on `q` enters. By Definition 2.1, `A_F(x) >>_q x^(alpha_q)` for `x >= q`, so

```
|M_F(x)| <<_(q,eps) x^(beta + c_q(m) + eps) <<_q A_F(x)^((beta + c_q(m) + eps)/alpha_q) = A_F(x)^(1 - delta + eps/alpha_q) ,
```

valid for `x >= q`, which is a saving exactly when `beta + c_q(m) < alpha_q`, that is when `gap_q(beta, m) > 0` by Definition 2.4. For `2 <= x < q` the sum `M_F(x)` runs over fewer than `q` terms and `A_F(x)` is either `0`, when `M_F(x) = 0` too, or at least `1`, so those `x` are absorbed by the implied constant as well. Given a fixed `delta' < delta`, choosing `eps = alpha_q(delta - delta')` yields `|M_F(x)| <<_(q,delta') A_F(x)^(1 - delta')`. Because `eps` is chosen after `delta'`, and `eps > 0` is forced, the endpoint `delta' = delta` is not attained by this argument.

For Theorem 1.1, `beta = 3/4` and `m = 1`, and `gap_q(1) > 0` for every `q >= 3690`: at `q = 3690` by Fact 6.2, one evaluation of two closed forms, and from there on by Proposition 5.2(i) at `b = 3/4`, whose floor is `Q(3/4) = 723 < 3690`. For Corollary 1.2 the hypothesis `gap_q(m) > 0` is assumed outright. □

**Proposition 4.5 (the shape at large base).** As `q -> infinity` with `m = 1` fixed,

```
c_q = (log log q + log(2/pi) + o(1))/log q -> 0 ,     alpha_q -> 1 ,     delta_q -> 1/4 .
```

So the exponent `1 - delta'` that Theorem 1.1 delivers can be pushed as close to `3/4` as one likes by taking `q` large: the exponent that the generalized Riemann hypothesis gives for `T(x)` on the whole line, transplanted to the missing-digit set and measured against that set's own mass. The convergence is logarithmic and nothing better.

**Proof.** `H(n) = log n + O(1)` and `n = ceil((q-2)/2)`, so `Phi_q = (2q/pi) log q + O(q)` and `P_q(1) = 1 + Phi_q/q = (2/pi) log q (1 + o(1))`. Hence `c_q = log_q P_q(1) = (log log q + log(2/pi) + o(1))/log q -> 0`, while `alpha_q = log_q(q - 1) = 1 + log_q(1 - 1/q) -> 1`. Therefore `delta_q = (alpha_q - 3/4 - c_q)/alpha_q -> 1/4`, and `1 - delta_q + eps = 3/4 + o(1)` as `q -> infinity` and `eps -> 0`. □

**Remark 4.6 (no hypothesis on the digits).** No hypothesis on `gcd(F)` is imposed, and none is needed: the proof never uses multiplicativity or any structure of `S_F` beyond the digit factorisation. In the range of Corollary 1.2 one holds anyway: `gap_q(m) > 0` forces `sqrt(m) <= P_q(m) < (q - m) q^(-3/4) < q^(1/4)`, so `m < sqrt(q)`, and `gap_q(m) <= gap_q(1)` forces `q >= 3690` by Fact 6.2; fewer than `floor(q/2)` avoided digits leave two consecutive digits in `F`, hence `gcd(F) = 1`. At small `k` a coprimality hypothesis does become necessary for other statements about `S_F`, and it must then be stated; here it would be decoration.

## The ladder

**Lemma 5.1 (the two exponent tables).** On `(1/2, 4/7)` one has `b_Z(a) < b_BH(a)`, with equality at `a = 1/2` (both `3/4`) and at `a = 4/7` (both `4/5`); and `b(a) >= 3/4` on all of `[1/2, 1)`.

**Proof.** Three factorisations, each an identity of rational functions:

```
b_Z(a) - (a + 1/4) = -5 (a - 1/2)(a - 2/5) / (4 - 2a) ,
b_Z(a) - 4/5       = -7 (a - 4/7)(a - 4/5) / (4 - 2a) ,
b_Z(a) - 3/4       = -7 (a - 1/2)(a - 6/7) / (4 - 2a) .
```

On `4 - 2a > 0` the signs read off: the first is negative on `(1/2, 11/20)`, where `b_BH(a) = a + 1/4`, the second on `[11/20, 4/7)`, where `b_BH(a) = 4/5`, and the third is positive on `(1/2, 4/7]`, with the first and second vanishing at `a = 1/2` and `a = 4/7` respectively. For `b_BH`, `a + 1/4 >= 3/4` and `(a + 1)/2 >= 3/4` when `a >= 1/2`, and `4/5 >= 3/4`. □

**Proposition 5.2 (the floor, and that every rung has a wall).** Let `b in [3/4, 1)` and `Q = Q(b)` as in Definition 2.5. Then (i) `gap_(q+1)(b, 1) > gap_q(b, 1)` for every `q >= Q`; (ii) `gap_q(b, 1) < 0` for every `3 <= q <= Q`; (iii) consequently `q_0(a)` exists and satisfies `q_0(a) > Q(b(a))`, with `gap_q(b(a), 1) > 0` for every `q >= q_0(a)`, so the hypothesis `q >= q_0(a)` of Theorem 1.3 describes a half line, not a window.

**Proof.** Write `u = 1 - b in (0, 1/4]` and

```
P_q^-(1) = 1 + 4/pi + (2/pi)(log((q-2)/2) + gamma) + (1 - 2/pi)(q-2)/q ,
```

which is at most `P_q(1)`, because `n = ceil((q-2)/2) >= (q-2)/2` and the terms `1/(2n)` and `0.727/q` are dropped.

(i) The step. First the cost side: `P_(q+1)(1) - P_q(1) < 1.291/(q-2)` for `q >= 40`. Passing from `q` to `q + 1`, the index `n` rises by `0` or `1`. If it does not rise the whole change is negative: the term `(1 - 2/pi)(q-2)/q` rises by `2(1 - 2/pi)/(q(q+1))`, the term `0.727/q` falls by `0.727/(q(q+1))`, and `2(1 - 2/pi) = 0.72676... < 0.727`. If `n` rises by one, then `H(n+1) - H(n) = log((n+1)/n) + 1/(2(n+1)) - 1/(2n) < 1/n <= 2/(q-2)`, contributing under `(4/pi)/(q-2)` through its factor `2/pi`; the term `(1 - 2/pi)(1 - 2/q)` rises by `2(1 - 2/pi)/(q(q+1))`, which is at most `0.017/(q-2)` once `q^2 - 41.75 q + 85.5 >= 0`, that is for `q >= 40`; and `0.727/q` falls. Adding, `P_(q+1)(1) - P_q(1) < (4/pi + 0.017)/(q-2) < 1.291/(q-2)`.

Now the mass side. With `f(q) = (q-1) q^(-b)` one has `f'(q) = q^(-b-1)((1-b) q + b) > (1-b) q^(-b)`, so

```
f(q+1) - f(q) = int_q^(q+1) f'(t) dt > (1-b)(q+1)^(-b) .
```

Hence `gap_(q+1)(b, 1) - gap_q(b, 1) > (1-b)(q+1)^(-b) - 1.291/(q-2)`, which is positive if `(1-b)(q-2)(q+1)^(-b) >= 1.291`, and that quantity is strictly increasing in `q` because `d/dq [(q-2)(q+1)^(-b)] = (q+1)^(-b-1)(q + 1 - b(q-2)) > 0`. So the step is positive at every `q >= Q(b)`.

(ii) Nothing closes at or below the floor. Since `q^(-b) <= q^(-3/4)` for `q >= 1` and `b >= 3/4`, `gap_q(b, 1) <= gap_q(3/4, 1)`, and the latter is negative on `3 <= q <= 3689` by Fact 6.2. So assume `Q = Q(b) >= 3690` and consider `3690 <= q <= Q`.

Both factors of `(1-b)(q-2)(q+1)^(-b)` fall as `b` rises at fixed `q`, so `Q(b)` rises with `b` and `Q >= Q(3/4) = 723`. Minimality of `Q` gives `(1-b)(Q-3) Q^(-b) < 1.291`, hence

```
(Q-1) Q^(-b) < 1.291/u + 2 Q^(-b) < 1.291/u + 0.015 ,
```

the slack because `2 * 723^(-3/4) < 0.015`; and `(1-b)(Q-2)(Q+1)^(-b) >= 1.291` with `Q - 2 < Q + 1` gives `(Q+1)^u > 1.291/u`, that is `log(Q+1) > log(1.291/u)/u`. Feeding both into `P_Q^-(1)` through `log((Q-2)/2) >= log(Q+1) - log(1448/721)` and `(Q-2)/Q >= 721/723`,

```
P_Q^-(1) > (2/pi) log(1.291/u)/u + kappa ,     kappa = 1 + 4/pi + (2/pi)(gamma - log(1448/721)) + (1 - 2/pi)(721/723) ,
```

and `kappa = 2.55916... > 0.015 + 2.544`, so

```
gap_Q(b, 1) < 1.291/u + 0.015 - (2/pi) log(1.291/u)/u - kappa < Lambda(u)/u ,     Lambda(u) = 1.291 - (2/pi) log(1.291/u) - 2.544 u ,
```

the coefficient `2.544` assembling the three ingredients `0.015`, `log(1448/721)` and `721/723`. Now `Lambda'(u) = 2/(pi u) - 2.544 > 0` for `u <= 1/4 < 2/(2.544 pi)`, so `Lambda(u) <= Lambda(1/4) = -0.39014... < 0` on `(0, 1/4]`, and dividing a negative number by `u <= 1/4` gives `gap_Q(b, 1) < 4 Lambda(1/4) < -1.56`.

Next, `U(q) = q^(1-b) - P_q^-(1)` dominates the gap, since `(q-1) q^(-b) < q^(1-b)` and `P_q^-(1) <= P_q(1)`; and it exceeds it by

```
U(q) - gap_q(b, 1) = q^(-b) + (2/pi)(log(2n/(q-2)) + 1/(2n)) + 0.727/q < 0.004     (q >= 3690) ,
```

because `q^(-b) <= 3690^(-3/4) < 0.0022`, `log(2n/(q-2)) <= log((q-1)/(q-2)) < 0.0003`, `1/(2n) < 0.0003` and `0.727/q < 0.0002`. So `U(Q) < -1.556`.

Finally `U` has exactly one interior minimum on `q > 2`. Indeed

```
U'(q) = (1-b) q^(-b) - (2/pi)/(q-2) - 2(1 - 2/pi)/q^2 ,
```

and `U'(q) > 0` if and only if `[(2/pi) q^b/(q-2) + 2(1 - 2/pi) q^(b-2)]/(1-b) < 1`. The bracket is strictly decreasing, since `d/dq [q^b/(q-2)] = q^(b-1)((b-1) q - 2b)/(q-2)^2 < 0` and `b - 2 < 0`, and it falls from `+infinity` to `0`; so `U'` is negative and then positive, and `U` attains its maximum on any interval at an endpoint. Hence on `[3690, Q]`,

```
gap_q(b, 1) < U(q) <= max(U(3690), U(Q)) .
```

We have `U(Q) < -1.556`. And `U(3690)` falls as `b` rises, with `U(3690) < -0.95` at `b = 1417/1850`; any `b < 1417/1850` has `Q(b) <= Q(1417/1850) = 1486 < 3690` and an empty interval, so in the case at hand `b >= 1417/1850` and `U(3690) < -0.95`. Both endpoints are therefore negative and `gap_q(b, 1) < 0` on `[3690, Q]`.

(iii) By (ii) the gap is negative up to `Q(b)`, by (i) it increases strictly from `Q(b)` on, and `gap_q(b, 1) -> infinity` as `q -> infinity` since `q^(1-b)` outgrows `P_q(1) = O(log q)`. So a least `q` with a positive gap exists, exceeds `Q(b)`, and the gap stays positive above it. The constants of this proof, `Q(b) >= 723`, `gap_(Q(b))(b, 1) < -1.56`, the slack under `0.004`, `U(3690) < -0.95` wherever `Q(b) >= 3690`, `Q(1417/1850) = 1486`, `Lambda(1/4) = -0.39014...` and `2 * 723^(-3/4) < 0.015`, are recomputed by the generator of Section 6 on the grid `b = 0.750, 0.755, ..., 0.975`. □

**Proof of Theorem 1.3.** Identical to the proof of Theorem 1.1 with `beta = b(a)` in place of `3/4`: Lemmas 4.1, 4.2 and 4.4 and the assembly never name the exponent, and the geometric series still has ratio above `1` because `b(a) < 1`. The quoted input is Lemma 3.1 at `a` when `b(a) = b_BH(a)` and Lemma 3.2 at `a` when `b(a) = b_Z(a)`; both have hypothesis exactly `Z(a)`, and `b(a)` is by definition the smaller. The condition `gap_q(b(a), 1) > 0` holds for every `q >= q_0(a)` by Proposition 5.2, and with `m` avoided digits the condition `gap_q(b(a), m) > 0` is assumed outright. □

## The numbers

Every number in this section is printed by one generator, [`mobius/lab/rs/mertens-numerology`](../lab/rs/mertens-numerology/), from the formulas of Definitions 2.3 to 2.5 alone, in the stated finite domain and to the stated direction of rounding; the run that printed them is named in Section 10. Throughout, `alpha_q` is truncated down at six decimals, `c_q` rounded up at five and `delta_q` rounded down at five, each from the unrounded value with a directional guard of `10^(-12)`, so every printed digit is a true bound in its own direction and `alpha_q` never prints as `1.000000`. Scientific rows carry a relative guard of `10^(-10)`.

**Fact 6.1 (the table, one avoided digit).** For `m = 1` and the ten bases of Table 1, the values are correct in the directions just stated, and the last column records the sign of `gap_q(1)`. Generator: `mobius/lab/rs/mertens-numerology`, the table block.

| `q` | `alpha_q` | `c_q` (up) | `delta_q` (down) | closes |
|---|---|---|---|---|
| 1000 | 0.999855 | 0.28087 | -0.03102 | no |
| 2000 | 0.999934 | 0.26335 | -0.01342 | no |
| 3000 | 0.999958 | 0.25430 | -0.00434 | no |
| 3689 | 0.999966 | 0.24997 | -0.00001 | no |
| 3690 | 0.999967 | 0.24997 | 0.00000 | yes |
| 5000 | 0.999976 | 0.24393 | 0.00605 | yes |
| 10^4 | 0.999989 | 0.23141 | 0.01858 | yes |
| 10^5 | 0.999999 | 0.19906 | 0.05094 | yes |
| 10^6 | 0.999999 | 0.17589 | 0.07411 | yes |
| 10^9 | 0.999999 | 0.13305 | 0.11695 | yes |

**Table 1.** The mass exponent, the cost exponent and the saving at ten bases. The rounded columns cannot display the sign at the wall; Fact 6.2 certifies it.

**Fact 6.2 (the wall at `3690`, and the margin).** With `m = 1`: `gap_q(1) < 0` at every `3 <= q <= 3689`, exhaustively, and `gap_3690(1) > 0`. At the two bases either side of the wall, to ten significant digits, the two `delta_q` margins through the cancellation-free form (2.1) and the two `gap_q` margins as the direct difference they are,

```
delta_3689 <= -2.395807653 * 10^(-6) ,     gap_3689(1) <= -1.533059397 * 10^(-4) ,
delta_3690 >=  5.863425182 * 10^(-6) ,     gap_3690(1) >=  3.752213034 * 10^(-4) .
```

Both margins sit below the fifth printed digit of Table 1, so that table neither displays nor certifies the sign there. Beyond the wall `gap_q(1)` rises at every one of the `96310` steps of `3690 <= q <= 10^5`, the smallest step being at least `3.172 * 10^(-5)` and occurring at `q = 99998`, at the top of the range rather than at the wall; Proposition 5.2(i) proves the rise from `723` on, so the scan is a check on the proof and not its support. The two tests `3/4 + c_q < alpha_q` and `gap_q(1) > 0` are asserted to agree at every integer `3 <= q < 20000`, and the two forms of `delta_q` to agree to `10^(-12)` at the wall. Generator: `mobius/lab/rs/mertens-numerology`, the wall and margin blocks.

**Fact 6.3 (how many digits may be removed).** The largest `m` with `gap_q(m) > 0` at `q = 10^4`, `10^5` and `10^6` is `6`, `78` and `451`, each asserted maximal by testing `m` and `m + 1`, with the exponents of Table 2. Generator: `mobius/lab/rs/mertens-numerology`, the `m`-corollary block.

| `q` | max `m` | `alpha_q` | `c_q(m)` (up) | `delta_q` (down) |
|---|---|---|---|---|
| 10^4 | 6 | 0.999934 | 0.24865 | 0.00129 |
| 10^5 | 78 | 0.999932 | 0.24972 | 0.00022 |
| 10^6 | 451 | 0.999967 | 0.24994 | 0.00003 |

**Table 2.** Corollary 1.2 at three bases, against `sqrt(q) = 100, 316, 1000`. Since `P_q(m) = sqrt(m) + Phi_q/q` and `Phi_q/q` is of size `(2/pi) log q`, the condition holds asymptotically for all `m <= q^(1/2)(1 - o(1))`: at `m = q^(1/2)(1 - eta)` the gap is `q^(1/4)(1 - sqrt(1 - eta)) - m q^(-3/4) - Phi_q/q`, positive once `eta q^(1/4)` outgrows `(2/pi) log q`.

**Fact 6.4 (the rungs of the ladder).** The values of `q_0(a)` and `Q(b)` in Table 3 are correct. A rung is meaningless quoted without its `a` and the table its exponent comes from, and two different rungs share `b = 4/5`; "both" means the two tables agree there. A row printed with `<=` gives an upper bound on the least `q` and not the least `q` itself: a wall prints as an exact integer only when it sits below `2^53` and both neighbouring gaps exceed `1024` ulps of the terms differenced. The exponents `b(a)` are carried as exact rationals and compared by cross multiplication, every wall below `4 * 10^6` is reproduced by an exhaustive scan from `q = 3` against the bisection, `Q(b) < q_0(a)` is asserted at every rung, and no `q < 3690` closes at any rung. Generator: `mobius/lab/rs/mertens-numerology`, the ladder block.

| `a` | `b(a)` | source | `q_0(a)` | `Q(b)` |
|---|---|---|---|---|
| 1/2 | 3/4 | both | 3690 | 723 |
| 13/25 | 1417/1850 | Zhang | 8578 | 1486 |
| 11/20 | 913/1160 | Zhang | 33547 | 4754 |
| 4/7 | 4/5 | both | 92317 | 11221 |
| 3/5 | 4/5 | Baker-Harman | 92317 | 11221 |
| 2/3 | 5/6 | Baker-Harman | 3107080 | 216023 |
| 3/4 | 7/8 | Baker-Harman | 6939524168 | 129458304 |
| 4/5 | 9/10 | Baker-Harman | <= 3.09358 * 10^13 | 128606353005 |
| 9/10 | 19/20 | Baker-Harman | <= 3.23663 * 10^34 | <= 1.73431 * 10^28 |
| 19/20 | 39/40 | Baker-Harman | <= 9.24614 * 10^83 | <= 3.30712 * 10^68 |

**Table 3.** Every common zero-free half plane buys the saving; the price is paid in the base. Reading down, `log_10 q_0(a)` rounded up at two decimals is `3.57, 3.94, 4.53, 4.97, 4.97, 6.50, 9.85, 13.50, 34.52, 83.97`.

**Fact 6.5 (what a weaker half plane spends first).** At `q = 10^7` the largest `m` with `gap_q(b(a), m) > 0`, for the rungs whose wall lies below `10^7`, is `1971` at `a = 1/2`, then `1002` at `a = 13/25`, `365` at `a = 11/20`, `176` at `a = 4/7` and at `a = 3/5`, and `8` at `a = 2/3`. Each is asserted maximal. Generator: `mobius/lab/rs/mertens-numerology`, the ladder `m`-corollary block.

## The wall of the method

Everything above is bought with one trade: a supremum over frequencies against an `l^1` mass. That trade has a floor, and the floor is an identity rather than an estimate. It is worth proving because it says exactly which problems this argument can never touch, and it says so unconditionally.

**Proof of Theorem 1.4.** The identity is the computation of Lemma 4.4 with `F` in place of `E`: the digits of `F` are distinct residues modulo `q`, so the off-diagonal terms are killed by `sum_(r mod q) e((d - d') r/q) = 0`, leaving `qk`. Since `|hat F| <= k` pointwise, `sum_r |hat F((t + r)/q)| >= (sum_r |hat F((t + r)/q)|^2)/max_r |hat F((t + r)/q)| >= qk/k = q`. So `B_q(F) >= q`. In the peel of Lemma 4.2 the inner sum over `s` is at least `q` at every `t'`, so the recursion runs with `>=` in place of `<=`: `Sigma_L >= q Sigma_(L-1)` and so `Sigma_L >= q^L`. □

**Corollary 7.1 (three things this decomposition cannot do).** Let `beta` be the exponent supplied to the quoted step and let `F` be any digit set.

(i) Sparse columns are out of reach. The chain closes only if `alpha_q > beta`. With `beta = 3/4` that is `k > q^(3/4)`, which `F = {0, 1}` in base `3` already fails, its mass exponent being `log_3 2 = 0.6309...`; and with the conjectured `beta = 1/2` in place of `3/4` it is still `k > q^(1/2)`, which every family of fixed digit count fails as soon as `q > k^2`. No hypothesis about `L`-functions moves either threshold.

(ii) Nothing unconditional is hiding here. Feeding Lemma 3.3 into the quoted step produces, as the bound on `|M_F(x)|`, the product of `q^(-L) Sigma_L` with `x (log x)^(-A)`, which is at least `x (log x)^(-A)` and exceeds `A_F(x)` by the factor `x^(1 - alpha_q)(log x)^(-A) -> infinity`. Moreover no unconditional uniform power saving can be substituted, because such an input is itself a zero-free half plane: if `T(x) << x^(1 - eta)` then in particular `M(x) = sum_(n <= x) mu(n) << x^(1 - eta)`, and `1/zeta(s) = s int_1^infinity M(u) u^(-s-1) du` then continues `1/zeta` analytically to `Re s > 1 - eta`, so `zeta` has no zero there.

(iii) Rearranging the same two factors is worse than trivial. Cauchy-Schwarz with Parseval on both sides, `sum_(a mod q^L) |hat F_L(a/q^L)|^2 = q^L k^L` and `sum_(a mod q^L) |S(q^L, a/q^L)|^2 = q^L sum_(n <= q^L) mu(n)^2 ~ (6/pi^2) q^(2L)`, gives at `x = q^L` the exponent `(1 + alpha_q)/2`, which exceeds `alpha_q` at every `alpha_q < 1`.

**Proof.** (i) is `alpha_q > beta + c_q(m)` together with `c_q(m) >= 0`, which holds because `P_q(m) >= 1`; Theorem 1.4 says the same for the exact constant, `B_q(F)/q >= 1`, so no sharpening of Lemma 4.4 can rescue the sparse range. (ii) and (iii) are the displayed computations, using Theorem 1.4 for the lower bound `q^(-L) Sigma_L >= 1` in (ii) and the two Parseval identities in (iii): the first is Lemma 4.1's orthogonality squared, the second is the same orthogonality on `n, n'` in `[1, q^L]`, a complete residue system, together with `sum_(n <= y) mu(n)^2 = (6/pi^2) y + O(sqrt y)`. □

**Remark 7.2 (a limit of the method, not of the problem).** The same phenomenon is named by the papers this one leans on. Maynard (2022) calls the exponent `4/5` of its minor-arc input a limit of its basic method, and its Theorem 1.3 takes fewer than `q^(1/5 - eps)` avoided digits; in Maynard (2019) the corresponding thresholds are `s <= q^(23/80)` avoided digits, and `s <= q - q^(57/80)` when they are consecutive, both measured by the `l^1` exponent of the digit transform. So the shape of the obstruction is known and belongs to the field. What Theorem 1.4 adds is the exact constant behind it in this setting: the shifted-grid energy is `qk` on the nose at every digit set and every `t`, so `B_q(F) >= q` with no room at all. [The Mobius page](../notes/grh-power-saving.md) sharpens the floor at one avoided digit to `B_q(F) >= 2(q - 1)`, by letting `t -> 0` on the grid, so that there the method needs `k > (2 - 2/q) q^(3/4)` and no exact constant can bring its base below `17`; that sharper floor is proved there and not here.

## What is in print

Every source below is listed in the references. Maynard (2019), Maynard (2022), Nath (2024), Leng and Sawhney (2025), Kim (2024), Baker and Harman (1991), Zhang (2024) and Porritt (2018) are read at source; Erdos, Mauduit and Sarkozy (1998) and Konyagin (2001) are not read at source and are carried for the statements in their titles, the distribution of missing-digit sets in residue classes; Dartyge and Mauduit (2000) is carried from its restatement in Maynard (2019), and Davenport (1937) from Porritt (2018).

Missing-digit sets carry a serious literature, and the arithmetic functions it has reached do not include `mu`. Erdos, Mauduit and Sarkozy (1998) and Konyagin (2001) study their distribution in residue classes; Dartyge and Mauduit (2000) get almost primes; Maynard (2019) proves there are infinitely many primes missing a digit in base ten, and Maynard (2022) extends that to large bases with many digits missing; Nath (2024) proves Bombieri-Vinogradov theorems for `Lambda(n) 1_(S_F)(n)`; Leng and Sawhney (2025) settle ternary Goldbach on such a set; Kim (2024) computes the divisor function over one, and frames the set's lack of multiplicative structure as what blocks the standard approaches. In every one of these the set enters through the Fourier norms of its digit transform, exactly as it does here.

What could not be found anywhere is the Mobius sum itself. In Maynard (2019) a whole-text search of the 70-page source finds "Mobius" once, in "By Mobius inversion" opening the proof of its Proposition 7.1, and "Mertens" and "Liouville" nowhere; in Nath (2024) the Mobius function is defined on p. 5 and used as `mu_(<= U)` in a Vaughan-type identity for `Lambda`, and no Mobius sum over the set is stated; in Maynard (2022) no Mobius or Mertens sum appears. So the positioning is a search and not a non-existence result: those three texts were read at their sources, the rest of the bibliography was read for its statements or carried as Section 8 says, and in none of them is there a bound on `sum_(n <= x, n in S_F) mu(n)`. The card this paper claims is precise, and each coordinate on it matters: the sequence is `mu` restricted to `S_F`; the norm is a pointwise bound in `x`, with no averaging and no modulus; the range is `q >= 3690` with one digit removed, or `m` digits under `gap_q(m) > 0`; and the saving is a fixed power of the set's own counting function, conditional on the generalized Riemann hypothesis. Change any coordinate and the claim is not this paper's.

One thing the reader should not conclude is that the shape of such a result is new. Maynard (2022) records that its asymptotic for `sum_(n < q^j) Lambda(n) 1_(S_F)(n)` has error terms that could be replaced by effective ones of size `O((q-1)^j exp(-c j^(1/2)))`, because the moduli occurring are highly composite and Siegel zeros do not intervene. A saving of that shape over a missing-digit set is therefore already in print, for `Lambda`. What is different here is the function, and what is worse here is the hypothesis.

On this tree the chain has since been run three more ways, each named in the introduction. [The Mobius page](../notes/grh-power-saving.md) proves a sharper one-step constant at one avoided digit and moves the wall of Theorem 1.1 to `1499`, and to `1032` at the end digits, the walls up-sets of exhaustive scans and half lines by the comparison `P'_q(e_0) < P_q(1)` with Proposition 5.2 here. [The First Base Below a Quarter](first-base-below-a-quarter.md), Theorem 5.1 and Corollary 5.2 there, replaces the one-step constant by the certified `l^1` exponent of the transform and moves it to `34` at every digit and to `21` at the digit `0`. [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md) removes the hypothesis from `584`, splitting the frequencies by Dirichlet approximation, which Corollary 7.1(ii) here says it must, and settling for `exp(-c sqrt(log x))`, which that corollary says it must too. On the half interval, half the digits kept as an interval, [An Unconditional Mertens Bound on the Half Interval](half-interval-mobius.md), Theorem 1.3 there, runs the certified chain under the hypothesis from every odd base `3789`. None of the three carries the ladder of Theorem 1.3, the corollary at `m >= 2` avoided digits with proved constants, or the floor of Proposition 5.2.

## Open problems

Four things are left undone, and none of them is close.

First, the unconditional problem in this frame. By Corollary 7.1(ii) no unconditional bound on `M_F(x)` follows from this decomposition. The dissection of [the large-base paper](../../coprime/papers/unconditional-mertens-at-large-base.md) gets one from base `584` by splitting the frequencies by Dirichlet approximation and using `mu` in progressions to the moduli dividing a power of `q`, with the saving `exp(-c sqrt(log x))` and a tiny `c`; whether any unconditional route returns a power of `A_F(x)` on a dense column is open, and the floor of Theorem 1.4 says the route cannot be this one.

Second, the sparse regime. Corollary 7.1(i) closes the range `k <= q^(3/4)` to this method permanently, and the natural target there is the one the census of [the Mobius page](../notes/meter-transfer-census.md) is built around.

**Conjecture 9.1 (the square-root conjecture).** For every digit set with `2 <= k <= q - 1` and squarefree digit gcd, `theta(F) = 1/2`: the running maximum `max_(y <= x) |M_F(y)|` is `A_F(x)^(1/2 + o(1))`. The 47 running-maximum exponents of the census across the bases `3`, `4`, `5` and `10` read `0.4465` to `0.5358` at finite depth, the full-set controls `0.4413` to `0.4517` at the same depths, and the finite tables decide nothing. Generator: `mobius/lab/rs/mobius-designs`.

Even the conjectural `beta = 1/2` in the quoted step leaves the floor at `k > q^(1/2)`, so the conjecture is beyond every version of this argument and not merely beyond its conditional form.

Third, the constant `0.727` and everything downstream of it. Lemma 4.4 splits `|hat F| <= |D_q| + |hat E|`, and that triangle inequality is where sharpness dies at moderate `q`. The Mobius page takes the loss back at one avoided digit, exactly on the grid and through a chord of the cosecant, which lowers the wall to `1499`; it also reads the exact supremum `B_q(F)` on a grid of `t`, and that reading would close the chain at `927` at every avoided digit and at `304` at the end digits, but it is a reading with no upper certificate on the supremum and no monotonicity in `q`, so it bounds nothing and enters no statement (Verified there; generator `mobius/lab/py/mrly-pairing`, verb `onestep`). The sharper single-arc proposition that Baker and Harman prove on the same pages as Lemma 3.1, which beats the uniform maximum at low denominators, buys no exponent here: the `l^1` mass of the transform sits at the top denominator, and the Mobius page refutes the per-denominator weighting as a lever on `c_q` or on `q_0`.

Fourth, whether `3690` is sharp for this method. It is the least base at which the proved constant `P_q(1)` closes the gap, and not the least base at which the true constant `B_q(F)/q` does. A proof that `q -> B_q(F)/q` crosses the mass exponent exactly once, which the proved bound has (Proposition 5.2) and the exact constant is not known to have, is the first thing anyone should want; the seat of the supremum is not `t = 1/2` at every base, so the proof cannot start by reading one point of the grid.

## Reproducibility

One generator prints every number of Section 6 and asserts each against its printed value, so a mismatch aborts the run: `cargo run --release --manifest-path research/mobius/lab/rs/mertens-numerology/Cargo.toml`, from the repository root, in about a second, prints Table 1 with its rounding directions and guard, the exhaustive scan of Fact 6.2 below `3690` with its four ten-digit margins and the `96310` steps to `10^5`, Tables 2 and 3 with the `log_10` trend, Fact 6.5, the sharpened walls of the Mobius page beside the step-3 wall `3690`, and a cost-out of the saving against a hypothetical defect that this paper does not use. Its tests, `cargo test` on the same manifest, pin every printed row as a string, the agreement of the two tests `3/4 + c_q < alpha_q` and `gap_q(1) > 0` at every `3 <= q < 20000`, the ten-digit margins in both directions and the agreement of the two forms of `delta_q`, the constants of Proposition 5.2 on the grid `b = 0.750, 0.755, ..., 0.975`, the walls of Table 3 below `4 * 10^6` against an exhaustive scan from `q = 3`, the `l^1` floor of Theorem 1.4, `Phi_q` against the exact shifted-grid kernel sum on a `4001`-point grid of `t` at `q = 50`, `101` and `200`, where the bound is loose by under `20` percent, and the asymptotic shape of Proposition 4.5, `c_q` within `0.01` of `(log log q + log(2/pi))/log q` at `q = 10^12`. An independent script of the standard library alone, linked from [the Mobius page](../notes/grh-power-saving.md), recomputes every number of Section 6 from the formulas of Definitions 2.3 to 2.5 and agrees to every printed digit. `uv run python research/mobius/lab/py/mrly-pairing/pairing.py onestep`, in about half a minute, prints the sharpened constants of the Mobius page, their walls `1499` and `1032`, and the readings `927` and `304` of Section 9. The figure is `figures/paper-sparse-mertens-under-grh.ts`, drawn in the browser from the level-`3` transform of base `7` with the digit `3` removed; its census `figures/paper-sparse-mertens-under-grh.json` holds the proved bound `(7 P_7(1))^3 = 18964.43...`, written by `research/shared/census/src/bin/paper-sparse-mertens-under-grh.rs` from Definition 2.3, and the figure asserts at draw time that the `343` absolute values it draws peak at `216` at the frequency `0`, have energy `343 * 216` (Theorem 1.4), and add to a total at least `343` and below the bound (Lemma 4.2). The theorems and their proofs print no number beyond these.

## References

- Baker and Harman 1991, Exponential sums formed with the Mobius function, J. London Math. Soc. 43, 193-198. [doi.org/10.1112/jlms/s2-43.2.193](https://doi.org/10.1112/jlms/s2-43.2.193)
- Zhang 2024, On an exponential sum related to the Mobius function, Proc. Amer. Math. Soc. 152, 1373-1376. [arxiv.org/abs/2204.04613](https://arxiv.org/abs/2204.04613)
- Davenport 1937, On some infinite series involving arithmetical functions (II), Quart. J. Math. os-8, 313-320, carried from Porritt 2018. [doi.org/10.1093/qmath/os-8.1.313](https://doi.org/10.1093/qmath/os-8.1.313)
- Porritt 2018, A note on exponential-Mobius sums over F_q[t], Finite Fields Appl. 51, 298-305. [doi.org/10.1016/j.ffa.2018.02.005](https://doi.org/10.1016/j.ffa.2018.02.005)
- Maynard 2019, Primes with restricted digits, Invent. Math. 217, 127-218. [link.springer.com](https://link.springer.com/article/10.1007/s00222-019-00865-6)
- Maynard 2022, Primes and polynomials with restricted digits, Int. Math. Res. Not. 2022, 10626-10648. [doi.org/10.1093/imrn/rnab002](https://doi.org/10.1093/imrn/rnab002)
- Nath 2024, Primes with a missing digit: distribution in arithmetic progressions and an application in sieve theory, J. London Math. Soc. 109, e12837. [arxiv.org/abs/2108.09212](https://arxiv.org/abs/2108.09212)
- Leng and Sawhney 2025, Vinogradov's theorem for primes with restricted digits. [arxiv.org/abs/2409.06894](https://arxiv.org/abs/2409.06894)
- Kim 2024, The divisor function over integers with a missing digit. [arxiv.org/abs/2411.09076](https://arxiv.org/abs/2411.09076)
- Erdos, Mauduit and Sarkozy 1998, On arithmetic properties of integers with missing digits I: distribution in residue classes, J. Number Theory 70, 99-120. [doi.org/10.1006/jnth.1998.2229](https://doi.org/10.1006/jnth.1998.2229)
- Konyagin 2001, Arithmetic properties of integers with missing digits: distribution in residue classes, Period. Math. Hungar. 42, 145-162. [link.springer.com](https://link.springer.com/article/10.1023/A:1015256809636)
- Dartyge and Mauduit 2000, Nombres presque premiers dont l'ecriture en base r ne comporte pas certains chiffres, J. Number Theory 81, 270-291, carried from Maynard 2019. [doi.org/10.1006/jnth.1999.2458](https://doi.org/10.1006/jnth.1999.2458)
