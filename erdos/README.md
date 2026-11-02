---
title: Erdos problems
lead: The open Erdos problems that are questions about digit designs, ranked by the odds of a checked result and by who reads them, and problem 124 attacked: sums of powers of several bases, where a bitset certificate makes each complete level a finite computation, every base set up to `10` with `sum 1/(d-1) > 1` and `gcd = 1` is certified at the first four levels, and the two-part route of the literature is proved unable to reach any base set below `sum 1/(d-1) = 2`, while two window criteria fall to explicit incomplete multisets; at `sum 1/(d - 1) = 1`, Ellison's bound on `|2^x - 3^y|` and Laurent's on linear forms in two logarithms make each level of every base set with no two bases powers of one integer one finite run, and all `62` such sets of at most six bases up to `47` are complete at the first four levels, all but `{3, 4, 7}` at the fifth.
figure: research-erdos
slug: erdos
date: 2026-09-23
---

# Erdos problems

The open Erdos problems that are questions about digit designs, ranked by the odds of a checked result and by who reads them, and problem 124 attacked: sums of powers of several bases, where a bitset certificate makes each complete level a finite computation, every base set up to `10` with `sum 1/(d-1) > 1` and `gcd = 1` is certified at the first four levels, and the two-part route of the literature is proved unable to reach any base set below `sum 1/(d-1) = 2`, while two window criteria fall to explicit incomplete multisets; at `sum 1/(d - 1) = 1`, Ellison's bound on `|2^x - 3^y|` and Laurent's on linear forms in two logarithms make each level of every base set with no two bases powers of one integer one finite run, and all `62` such sets of at most six bases up to `47` are complete at the first four levels, all but `{3, 4, 7}` at the fifth.

- [notes/erdos.md](notes/erdos.md): Erdos problems
- [notes/problem-124.md](notes/problem-124.md): Problem 124
- [notes/problems-817-406.md](notes/problems-817-406.md): Problems 817 and 406
- [notes/sigma-one.md](notes/sigma-one.md): Sigma one
- [claims/erdos.md](claims/erdos.md): Erdos
- [refs.md](refs.md): the references this topic cites
- [papers/finite-certificates.md](papers/finite-certificates.md): Finite Certificates for Erdos Problems 817 and 124
- [lab/rs/mixed-powers/](lab/rs/mixed-powers/): Computes, for a finite base set `D` of integers `>= 3` and a level `k >= 1`, whether every large integer is a sum of distinct terms of the multiset `M_k(D)` of powers `d^j`, `d` in `D`, `j >= k`, one copy per base
- [figures/research-erdos.ts](figures/research-erdos.ts): research-erdos
- [figures/research-erdos.json](figures/research-erdos.json): research-erdos census
- [figures/paper-finite-certificates.ts](figures/paper-finite-certificates.ts): paper-finite-certificates

## Lost

- Problem 124 at every level by a finite certificate alone: the certificate passes completeness only downward in `k`, the deepest hole grows with `k`, and at `{3, 4, 5}` level 6 leaves the box; the levels closed are in [claims/erdos.md](claims/erdos.md).
- Problem 125's upper density by finite runs: the rows need `liminf h_k(tau_k) < infinity`, which no finite run reaches; the lower density was settled elsewhere.
- The two-part route of Fan and of Bergelson and Simmons below `sum 1/(d-1) = 2`: no two disjoint infinite parts of the powers are both syndetic there (The two-part route is dead below two).
- Linear surplus with `H_1 = {0}` as a completeness criterion: `2 F_m - 1` is incomplete; the Refuted row in [claims/erdos.md](claims/erdos.md).
- The chain-excluded lemma and its windowed repair: explicit multisets with surplus, `H_1 = {0}` and no chain are incomplete, so no window criterion decides problem 124 (The windowed lemma).
