---
title: Which levels tile
lead: Which levels of a dimension-one design tile the integers: the index lemma reads every level's cyclotomic divisors off the zero set of the digit transform, the tiling levels form an initial segment, a closed condition decides T1 at every level, the guess that the digit set must tile its residues fails at base 12, a horizon decides T2, and the classification decides every digit set whose size has at most two prime factors.
figure: research-tilings
slug: tilings
---

## The object

- `Z(F) = {m >= 2 : Phi_m | F(x)}`, the orders of the roots of unity where `hat F` vanishes; `hat F(r/m) = 0` with `gcd(r, m) = 1` exactly when `m` is in `Z(F)`.
- `Z_p(F) = {c >= 1 : p^c in Z(F)}` for a prime `p`, and `e_p = v_p(b)`, the exponent of `p` in the base.
- `S_A` is the set of prime powers `s` with `Phi_s | A(x)`. Coven and Meyerowitz read two conditions on it: **T1** `A(1) = prod_{s in S_A} Phi_s(1)`, and **T2** `Phi_{s_1 ... s_m} | A(x)` whenever `s_1, ..., s_m` in `S_A` are powers of distinct primes.
- Their theorems, read at source: T1 and T2 imply that `A` tiles (Theorem A); a tile satisfies T1 (Theorem B1); a tile whose size has at most two prime factors satisfies T2 (Theorem B2). Whether every tile satisfies T2 is open.
- The tiling depth of `F` at base `b` is the largest `n` such that `A_1, ..., A_n` all tile, `inf` when every level tiles.

## The index lemma

**Proved.** `Phi_m | A_n(x)` if and only if `m/gcd(m, b^i)` is in `Z(F)` for some `i < n`, and the multiplicity of `Phi_m` in `A_n` is `sum_{i<n}` of the multiplicity of `Phi_{m/gcd(m, b^i)}` in `F`. A primitive `m`-th root `zeta` makes `zeta^(b^i)` a primitive root of order `m/gcd(m, b^i)`, so `F(zeta^(b^i)) = 0` exactly when that order is in `Z(F)`; and `Phi_u(x^j)` is a product of distinct cyclotomic polynomials, so multiplicities add over the factors `F(x^(b^i))`. In Fourier terms the zero set of `hat F_n(t) = prod_{i<n} hat F(b^i t)` on `Q/Z` is the union over `i < n` of the preimages of the zero set of `hat F` under `t -> b^i t`; one step of it is Coven-Meyerowitz Lemma 1.1(7).

For a prime power the lemma reads off the level exponents: `p^c` is in `S_{A_n}` exactly when `c` is in `Z_p(F)` for `p` not dividing `b`, and when `c - e_p i` is in `Z_p(F)` for some `i < n` for `p | b`. Every census verdict below is computed from `Z(F)` through this lemma and never from a factorisation of the level. (**Verified**: PARI factors 752 masks, every digit set with `0` at bases 2 to 8, levels 1 to 3 and level 4 at bases up to 4, and the cyclotomic part equals the lemma with multiplicity on all 752; `index`.)

## The hand check

The Cantor level `{0, 2, 6, 8}`, level 2 of `{0, 2}` at base 3, has mask `(1 + x^2)(1 + x^6) = Phi_4^2 Phi_12`, so `S = {4}` and T1 reads `|A_2| = 4` against `Phi_4(1) = 2`: it does not tile. (**Verified**: PARI and the lemma agree on `Phi_4^2 Phi_12`, and no complement exists in `Z/N` for `4 | N <= 128`; `hand`.)

## The tiling levels form an initial segment

**Proved.** If `A_n` tiles then `A_m` tiles for every `m <= n`. Digits below the base give unique expansions, so `A_n = A_m + b^m A_(n-m)` with every element written once, and `A_n + C = Z` gives `A_m + (b^m A_(n-m) + C) = Z`. So the tiling levels are `1..N` for the depth `N`, and a digit set that does not tile `Z` has no tiling level at all: no non-tiling `F` has a tiling level.

## A digit set that tiles its residues tiles every level

**Proved.** If `F + E = Z/b` with every residue written once, then `A_n + E_n = Z/b^n` with `E_n = sum_{i<n} b^i E`, so every level tiles the integers with period `b^n`. The set `D = F + E` is a complete residue system mod `b`, and `D_n = sum_{i<n} b^i D` is one mod `b^n`: a residue `r` fixes `d_0` in `D` by `r mod b`, then `(r - d_0)/b` fixes `d_1`, and so on, the carry automaton of [cuts](../../cuts/notes/cuts.md) run on the digit polynomial `D(x) = F(x) E(x)`, which is `1 + x + ... + x^(b-1)` mod `x^b - 1`. Each `d_i` splits uniquely as `f_i + e_i`, so the sum `A_n + E_n` is direct.

The converse is the guess, and it is false (below). What survives is a condition on `Z(F)` alone.

## T1 at every level

**Proved.** T1 holds at every level `A_n` if and only if
- T1 holds for `F`, that is `k = prod_p p^|Z_p(F)|`;
- `Z_p(F)` is empty for every prime `p` not dividing `b`;
- for every `p | b` the elements of `Z_p(F)` are pairwise incongruent mod `e_p`.

By the index lemma the exponents of `p` in `S_{A_n}` are `Z_p(F)` when `p` does not divide `b` and the union of the translates `Z_p(F) + e_p i`, `i < n`, when it does. T1 at level `n` compares `v_p(k^n) = n v_p(k)` with the size of that set. Off the base the size is constant, so `v_p(k) = |Z_p(F)| = 0`. On the base the union has at most `n |Z_p(F)|` elements and level 1 gives `v_p(k) = |Z_p(F)|`, so T1 holds for every `n` exactly when the translates are pairwise disjoint, that is when no two elements of `Z_p(F)` differ by a positive multiple of `e_p`.

**Corollaries, Proved.**
- If every level tiles then `|F|` divides `b`: T1 is necessary (Theorem B1), and the condition forces `|Z_p(F)| <= e_p`.
- The T1 depth, the largest `n` with T1 at levels `1..n`, is `0` when `F` fails T1, `1` when some prime off the base has `Z_p(F)` nonempty, and otherwise the least `(c' - c)/e_p` over pairs `c < c'` in some `Z_p(F)` with `c = c' mod e_p`, `inf` when there is none.
- When `|F|` is a prime power, T2 is vacuous at every level, so the tiling depth equals the T1 depth and every level tiles exactly when the three conditions hold.
- When the base is a prime power `p^e`, every level tiles exactly when `F` tiles `Z/b`. A factor `Phi_{p^c}` with `c > e` has degree `p^(c-1)(p-1) >= p^e > deg F`, so `Z_p(F)` lies in `1..e`, the size of `F` is a power of `p`, and the Coven-Meyerowitz complement of `F` has period `lcm(S_F)`, a divisor of `b`.

The depth can stop anywhere. `{0, 1, 4, 5}` at base 6 has mask `Phi_2 Phi_8`, `Z_2 = {1, 3}` and `e_2 = 1`, so levels 1 and 2 tile and level 3 does not; `{0, 1, 8, 9}` at base 10 has `Z_2 = {1, 4}` and stops at level 3. (**Verified**: on all 65519 digit sets with `0` at bases 2 to 16 the census asserts the T1 depth formula against the T1 depth read level by level to level 16, the tiling depth against the least of the T1 depth, the T2 depth and 8, which equals the T1 depth at every prime-power size while T2 cuts below it on 771 sets of the other sizes, and the prime-power-base law, and prints depth `2` at bases 6, 10 and 14 and depth `3` at bases 10 and 14; `census`, `witness`.)

## The guess is false

The guess: every level tiles exactly when `F` tiles `Z/b`. One direction is the residue section above; the other fails.

**Refuted**, twice.
- `F = {0, 2}` at base 6. `{0, 2}` does not tile `Z/6`, since translation by `2` has odd cycles on `Z/6`, yet `Z(F) = {4}`, `Z_2 = {2}`, `e_2 = 1`, and every level tiles by the T1 condition. The cure is to divide by the gcd: tiling is invariant under scaling (Coven-Meyerowitz Lemma 1.4), and `{0, 1}` tiles `Z/6`.
- `F = {0, 1, 8, 9}` at base 12, with gcd `1`. Its mask is `(1 + x)(1 + x^8) = Phi_2 Phi_16`, so `Z_2 = {1, 4}`, `e_2 = 2`, and `1`, `4` are incongruent mod `2`: T1 holds at every level, and `|A_n| = 4^n` is a prime power, so every level tiles. `F` itself tiles `Z` with period `16`, `{0, 1, 8, 9} + {0, 2, 4, 6} = {0..15}`. But it does not tile `Z/12`: a complement `E` of size `3` would carry `Phi_3`, `Phi_4`, `Phi_6` and `Phi_12`, since `F` vanishes at none of those orders, and T1 for `E` would read `3` against at least `Phi_3(1) Phi_4(1) = 6`. No integer multiple `vF` tiles `Z/12` either. Scaling raises both exponents, `Z_2(vF) = {1 + v_2(v), 4 + v_2(v)}` by Coven-Meyerowitz Lemma 1.4(2), and their Lemma 2.1 puts every element of `S_vF` among the prime powers dividing `12` when `vF` tiles `Z/12`; whatever happens to `1 + v_2(v)`, the power `2^(4 + v_2(v))` never divides `12`. (**Verified**: an exact-cover search over `v = 1..12`, every class mod 12, finds no complement for `vF` or for `v{0, 3, 8, 11}`; `witness`.)

**Verified.** No base below 12 has a counterexample after the gcd: over the 65519 digit sets at bases 2 to 16, 609 tile at every level and 588 tile `Z/b`; of the 21 that tile at every level without tiling `Z/b`, 19 are cured by the gcd, and the remaining two are `{0, 1, 8, 9}` and `{0, 3, 8, 11}` at base 12, both mirror-symmetric, zero sets `{2, 16}` and `{2, 6, 16}` (`census`). The levels of `{0, 1, 8, 9}` at base 12 tile with complements `{0, 2, 4, 6}` mod `16`, `{0, 2, 16, 18}` mod `64` and `{0, 2, 64, 66}` mod `256` (`witness`). The reason is visible in the T1 condition: a counterexample after the gcd needs an exponent `c` in some `Z_p(F)` above `e_p`, which the degree bound forbids at a prime-power base; a prime size `p` is cured by the gcd, since `Phi_{p^c} | F` with `|F| = p` puts `F` in one class mod `p^(c-1)`; size `4` needs `4 | b`, and 12 is the first base divisible by `4` that is not a prime power.

## T2 at every level

**Proved.** If `F` satisfies T1 and T2 and `lcm(S_F)` divides `b`, then every level satisfies T1 and T2, so every level tiles; this is the case where `F` tiles `Z/b` through its Coven-Meyerowitz complement. Each exponent `c` in `Z_p(F)` is at most `e_p`, so they are incongruent mod `e_p`. A T2 instance at level `n` takes, for distinct primes `p_j`, exponents `c_j + e_j i_j` with `c_j` in `Z_{p_j}(F)` and `i_j < n`; at `l = max i_j` the index lemma reduces it to `prod p_j^((c_j - e_j(l - i_j))_+)`, which keeps only the primes with `i_j = l` and their exponents `c_j`, a product of elements of `S_F` at distinct primes, and T2 for `F` puts it in `Z(F)`.

**Proved (the horizon).** Assume the T1 condition. Let `C_p` be the largest exponent of `p` in an element of `Z(F)`, `G` the largest `ceil(C_p/e_p)` over the `omega` primes with `Z_p(F)` nonempty. Then T2 holds at every level if and only if it holds at levels `1..1 + G(2 omega - 1)`. A T2 instance is hardest at the level `1 + max i_j`, and its witnesses are the `l <= max i_j` whose reduction lies in `Z(F)`. Every witness has `l >= i_j - G` for every `j`, since otherwise the exponent of `p_j` exceeds `C_{p_j}`. So an instance with `min i_j > G` can be lowered to `min i_j = G` with its witnesses; and if the sorted `i_j` have a gap above `2G`, every witness lies more than `G` above each lower `i_j`, the lower primes vanish at every witness, and the instance stands or falls with its upper part. The reduced instances have `max i_j <= G + 2G(omega - 1)`. (**Verified**: on the census the horizon agrees with T2 checked six levels past it on every digit set that passes the T1 condition; `census`.)

**Proved.** `F = {0, 1, 2, 6, 7, 8}` at base 18 tiles the integers and satisfies T1 at every level, yet level 2 does not tile. The mask of `F = {0, 1, 2} + {0, 6}` is `Phi_3 Phi_4 Phi_12`, so `S_F = {3, 4}`, T1 and T2 hold, and `F` tiles with period 12; `Z_2 = {2}`, `Z_3 = {1}` satisfy the T1 condition with `e_2 = 1`, `e_3 = 2`. At level 2 the exponents are `{2, 3}` for `2` and `{1, 3}` for `3`, so T2 asks for `Phi_108 = Phi_{4 * 27}`; the index lemma gives `108/gcd(108, 1) = 108` and `108/gcd(108, 18) = 6`, neither in `Z(F) = {3, 4, 12}`. `|A_2| = 36` has two prime factors, so Theorem B2 says `A_2` does not tile. So `F` tiles `Z`, T1 holds at every level, and level 2 still fails, through T2.

**Verified.** No base below 18 has such a digit set. Over the digit sets of a size with two prime factors below the base, at bases up to 30 with at most 200000 sets a cell, those that tile, pass the T1 condition, and fail T2 at a level up to 6 number `0` of size 6 at base 12, `16` of size 6 (10 up to mirror) at base 18, `0` of size 10 at base 20, `4` of size 6 at base 24 and `258` of size 6 (141 up to mirror) at base 30, every one failing at level 2; size 12 at base 24 (1352078 sets), size 14 at base 28 (20058300) and sizes 10 and 15 at base 30 exceed the cap and are skipped; at bases below 18 the only such size is `6` at base 12 (`witness`). A level-2 failure of this kind needs an exponent above `e_p` together with T2 for `F`, and the degree of `Phi_{p^c} Phi_{q^d} Phi_{p^c q^d}` must fit below `b`.

## The classification

**Proved.** For a base `b` and a digit set `F` in `{0..b-1}`:
- every level tiles only if `F` passes the T1 condition, and then `|F|` divides `b`;
- if `|F|` is a prime power, every level tiles if and only if `F` passes the T1 condition;
- if `|F|` has two prime factors, every level tiles if and only if `F` passes the T1 condition and T2 holds at levels `1..1 + G(2 omega - 1)`;
- if `|F|` has three or more prime factors, the same condition is sufficient, and its necessity is the open necessity of T2; such a size divides `b`, is below it, and so needs `b >= 60`;
- if `b` is a prime power, every level tiles if and only if `F` tiles `Z/b`.

The first two lines and the prime-power base are T1 alone; the third is Theorems A, B1 and B2 applied level by level through the index lemma and the horizon.
