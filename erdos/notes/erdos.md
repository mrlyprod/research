---
title: Erdos problems
lead: The open Erdos problems that are questions about digit designs, ranked by the odds of a checked result and by who reads them, and problem 124 attacked: sums of powers of several bases, where a bitset certificate makes each complete level a finite computation, every base set up to `10` with `sum 1/(d-1) > 1` and `gcd = 1` is certified at the first four levels, and the two-part route of the literature is proved unable to reach any base set below `sum 1/(d-1) = 2`, while two window criteria fall to explicit incomplete multisets; at `sum 1/(d - 1) = 1`, Ellison's bound on `|2^x - 3^y|` and Laurent's on linear forms in two logarithms make each level of every base set with no two bases powers of one integer one finite run, and all `62` such sets of at most six bases up to `47` are complete at the first four levels, all but `{3, 4, 7}` at the fifth.
figure: research-erdos
slug: erdos
---

A digit design is a set of integers whose base-`d` digits lie in a set `F`, and a surprising share of the open problems at [erdosproblems.com](https://www.erdosproblems.com) are statements about one, about a sumset of two in different bases, or about the carries that Kummer's theorem turns into divisibility. This page ranks those problems, then attacks the top one, problem 124, with the machinery of [Two bases](../../cobham/notes/cobham.md): the sumset of designs across bases held as one bit array.

`lab/rs/mixed-powers` is the generator behind every number below.

## The objects

- `P(d, k)` is the set of sums of distinct powers `d^j` with `j >= k`, which is `d^k` times the base-`d` design `{0, 1}`.
- A base set is a finite `D = {d_1 < ... < d_r}` of integers `>= 3`; `M_k(D)` is the multiset of the powers `d_i^j`, `j >= k`, one copy per base, so `sum_i P(d_i, k)` is the set `FS(M_k(D))` of sums of distinct terms of `M_k(D)`.
- `sigma(D) = sum_i 1/(d_i - 1)` and `C_k(D) = sum_i d_i^k/(d_i - 1)`.
- `M_k(D)` is complete when `FS(M_k(D))` contains every large integer; `F(D, k)` is then the largest integer outside it.
- `a_1 <= a_2 <= ...` is `M_k(D)` sorted, `S_n = a_1 + ... + a_n`, and `P_n` the set of sums of distinct terms among the first `n`, a subset of `[0, S_n]` symmetric under `x -> S_n - x`.

## The ranked scouting table

Rank is the odds of a checked result in one sitting times the readership; status is the site's own; a problem this tree already carries is marked and not reopened.

| rank | problem | statement | site | the angle |
|---|---|---|---|---|
| 1 | [124](https://www.erdosproblems.com/124) | for `D` with `sigma(D) >= 1` and `gcd(D) = 1`, is `sum_i P(d_i, k)` cofinite for every `k >= 1`? | open | a sumset of scaled `{0, 1}` designs across bases; the surplus `sigma - 1` and a bitset certificate per level; attacked below |
| 2 | [817](https://www.erdosproblems.com/817) | least `N` with an `n`-set in `[1, N]` whose subset sums hold no 3-term progression; is it `>> 3^n`? | open | the powers of `3` give the base-3 design `{0, 1}`; admissible means the `{0, 1, 2}` digit sums are distinct, a base-3 numeration; the discussion thread reports `o(3^n)` and exact values to `n = 6`; `g_3(7) = 474` below, and at `k = 4` `g_4(n) = 1, 3, 5, 14, 40, 79` to `n = 6` |
| 3 | [406](https://www.erdosproblems.com/406) | only finitely many powers of `2` lie in the base-3 design `{0, 1}`? | open | the orbit of `x -> 2x` against a base-3 design, a 3-adic carry automaton; the site records every digit present for `16 <= n <= 5.9 10^21`, past any census here; the 3-adic relaxation at four exponents is attacked below |
| 4 | [376](https://www.erdosproblems.com/376) | infinitely many `n` with `binomial(2n, n)` prime to `105`? | open | carried: three designs at bases `3, 5, 7`, [Two bases](../../cobham/notes/cobham.md) Object T, not reopened |
| 5 | [271](https://www.erdosproblems.com/271) | growth of the greedy sequences without 3-term progressions | open | `A(1)` is the base-3 design `{0, 1}`; which starts give a design, read as digit rules |
| 6 | [168](https://www.erdosproblems.com/168) | density of the largest subset of `[1, N]` with no `{n, 2n, 3n}` | open | a problem on the lattice of `2^a 3^b`, the two-base object of [Two bases](../../cobham/notes/cobham.md) |
| 7 | [1146](https://www.erdosproblems.com/1146) | is `{2^m 3^n}` an essential component? | open | the same two-base lattice, against Schnirelmann density |
| 8 | [354](https://www.erdosproblems.com/354) | is `{floor(2^s alpha)} U {floor(2^t beta)}` complete for `alpha/beta` irrational? | open | two binary digit streams of `alpha` and `beta` summed with carries |
| 9 | [396](https://www.erdosproblems.com/396) | for every `k` an `n` with `(n - k) ... n` dividing `binomial(2n, n)`? | open | Kummer: each divisibility is a carry count of `n + n` in each base, a finite carry automaton per prime |
| 10 | [377](https://www.erdosproblems.com/377) | is the sum of `1/p` over `p <= n` not dividing `binomial(2n, n)` bounded? | open | `p` does not divide exactly when `n` lies in the base-`p` design of digits below `p/2`; a bound is a transversality of many designs |
| 11 | [1145](https://www.erdosproblems.com/1145) | `A + B` cofinite and `a_n/b_n -> 1` force an unbounded representation count? | open | the page's own example is two binary designs, even and odd digit places, the first family to test |
| 12 | [349](https://www.erdosproblems.com/349) | for which `t, alpha` is `floor(t alpha^n)` complete? | open | expansions in the non-integer base `alpha` |

Problem 125 is Object S of [Two bases](../../cobham/notes/cobham.md) and stays there.

## Chapters

- [Problem 124](problem-124.md): Problem 124 attacked with a bitset certificate: what is in print, the surplus certificate, the two-part route proved dead below two, the census of certified base sets, the missing lemma, the two routes, the windowed lemma and the deepest hole.
- [Sigma one](sigma-one.md): Problem 124 where `sum 1/(d - 1) = 1`: an exact deficit writes a failing step as close powers of two bases, so a lower bound on `|e^a - f^b|`, Ellison's for `{2, 3}`-smooth pairs and Laurent's on linear forms in two logarithms with an exact enumeration for the rest, makes each level of every such base set with no two bases powers of one integer one finite run that halts exactly when the level is complete; all `62` such sets of at most six bases up to `47` are complete at `k <= 4`, and all but `{3, 4, 7}` at `k = 5`.
- [Problems 817 and 406](problems-817-406.md): Problem 817 at n = 7, its values at k = 4 to n = 6, and problem 406 in the 3-adic integers, each settled by a finite computation.

## What stays open

- **Conjecture.** `M_k({3, 4, 5})` is complete at every `k`, the least base set outside every theorem whose proof is read here; checked at `k <= 5` with no derivation. A proof has to hold uniformly in `k`, while the only uniform tool read here, the two-part criterion, is dead below `sigma = 2`.
- **Proved (the target as one inequality).** If `1 < sigma(D) < 2`, then `M_k(D)` is complete exactly when some `n` has `2 T_n - 1 <= (sigma(D) - 1) a_(n+1) - C_k(D)`, with `T_n = 1 + max{x <= S_n/2 : x not in P_n}`. Proof: at such an `n`, `T_n < a_(n+1)` since `sigma - 1 < 1`, and the surplus clause of the certificate is empty because `a_(m+1) >= a_(n+1)` for `m >= n`; conversely, along the steps of the gap proof `T_n = F(D, k) + 1` while `a_(n+1)` grows without bound.
- So a proof for `{3, 4, 5}`, where `sigma - 1 = 1/12`, is a bound uniform in `k` that brings the deepest hole of `P_n` below about `a_(n+1)/24`. A criterion whose window clause a windowless multiset satisfies cannot supply it: the chain-excluded lemma is refuted above by a multiset with no window, and the pair route makes every `F` of the census from `k = 4` on; the windowed lemma is refuted by `V`, whose windows open at every fifth step. Such a bound cannot hold while the births run: the step after a rise has `S_(n+1) - a_(n+2) <= 2 T_(n+1) - 2`, so the one inequality fails there, every rise lands past `S_n/2`, and every certified cell with `k >= 2` has `T_n > 1.17 a_(n+1)` after its last base enters. At `sigma = 1` the surplus is `-C`, and an exact deficit and Ellison's bound make each level a finite run, complete exactly when its run closes, for every base set holding two `{2, 3}`-smooth bases and no two bases powers of one integer ([Sigma one](sigma-one.md)): `{3, 4, 7}` closes through `k = 4`, the other `51` such sets of at most six bases up to `47` through `k = 5`. The `10` sets of that range with one smooth base, from `{3, 5, 7, 13}`, wait on a lower bound for a linear form in two logarithms.
