---
title: Windows of a design
lead: The infinite design of a plane code as a subshift: its `k x k` windows counted exactly through the substitution, the carpet's count, conjecturally piecewise quadratic with `p(k)/k^2` of lower limit `8` and upper limit `9.6`, every code at bases 2 and 3 of growth `1`, `k` or `k^2`, of finite type only when empty, full or on the boundary, and sofic always.
figure: research-windows
slug: windows
---

Fix a base `b` and a plane code at that base, a bitmask over the `b x b` cells of level 1 read in row-major order; let `F` be its filled cells. The level-`level` render `R_level` is the `b^level x b^level` picture with a 1 at `(i, j)` exactly when every base-`b` digit pair `(i_k, j_k)`, `k < level`, lies in `F`. It is the substitution that replaces a 1 by the level-1 tile and a 0 by the empty tile, run `level` times from one 1. This page reads the infinite design as a set of plane pictures closed under the shift: it counts the `k x k` windows that occur, asks whether finitely many forbidden windows cut the set out, as for a [subshift of finite type](https://en.wikipedia.org/wiki/Subshift_of_finite_type), and asks how few [Wang tiles](https://en.wikipedia.org/wiki/Wang_tile) can lay it.

## The objects

A **render window** is a `k x k` block that occurs in some render. For a nonempty code `R_level` sits inside `R_(level+1)` at every filled cell, so a window that occurs at one level occurs at every later one. The design's subshift `X_F` is the set of pictures `x` on all of `Z^2` whose every finite block is a render window; it is closed and invariant under the plane shift, and for every code but the full one it is the subshift of the substitution in the sense of [Mozes 1989](https://doi.org/10.1007/BF02793412); for the full code that subshift is `{0, 1}` and `X_F` is `{1}`. The **window count** `p(k)` is the number of distinct `k x k` blocks of the pictures of `X_F`. A **boundary code** is a code whose filled cells all lie in one boundary line of the `b x b` box: the top row, the bottom row, the left column or the right column.

**Proved.** If `F` is neither empty nor a boundary code, every render window is a block of `X_F`, so `p(k)` is the number of render windows; when `(0, 0)` is filled it is also the number of `k x k` windows of the quarter-plane picture `a(i, j)`, whose top-left `b^level` block is `R_level`. Since `F` is in no boundary line, some filled cell has a row digit other than `0`, some other than `b - 1`, and likewise for the column digit; a level-`n` cell whose digit word runs through all of `F` is then a 1 strictly inside `R_n`, at `(P, Q)` with `0 < P, Q < b^n - 1`. A window of `R_m` occurs in every `R_(m')`, `m' >= m`, and the copy of `R_(m')` at block `(P, Q)` of `R_(m' + n)` has all eight neighbouring blocks inside that render, so the window occurs with a margin of `b^(m')` on every side; letting `m'` grow and taking a limit of these patches gives a picture of `X_F` that shows it.

**Proved.** A boundary code has `X_F = {0}`, so `p(k) = 1`. If the filled cells lie in the top row, every 1 of every render sits in its top row; a picture of `X_F` with a 1 at `(i, j)` would show the `2 x 2` block at `(i - 1, j)` with a 1 in its lower row, which no render has. The other three lines are the same argument turned. The empty code also gives `{0}` and the full code `{1}`, so `p(k) = 1` for all three.

**Proved.** `p` is unchanged by the eight symmetries of the square acting on the code: transposing the box swaps the two digit words, and reflecting it sends every digit `d` to `b - 1 - d`, which reflects every render. The `count` verb of `lab/rs/design-windows` asserts the equality on every code at both bases for `k <= 12` (**Verified**), so the census runs one code per orbit: `6` orbits of the `16` codes at base 2 and `102` of the `512` at base 3.

The carpet is code `7` at base 2 drawn at side 3, which fills a cell of the `3 x 3` box unless both coordinates are odd; at base 3 that is every cell but the centre, code `511 - 16 = 495`. (**Verified**: the `scan` verb asserts the crate's two renders agree cell for cell at level 6.) The Sierpinski gasket is code `7` at base 2 at its own side, `a(i, j) = 1` exactly when the binary digits of `i` and `j` never share a 1, which is Pascal's triangle mod 2 in the square form `binom(i + j, i) mod 2`.

## Counting through the substitution

Write `sigma` for one round of the substitution and `L_2` for the set of `2 x 2` render windows.

**Proved.** `L_2` is reached by a closure that stops. A `2 x 2` window of `R_(level+1)` lies inside the image `sigma(B)` of some `2 x 2` window `B` of `R_level`, and every `2 x 2` window of such an image is a window of `R_(level+1)`. So the window sets of successive renders are `S, Phi(S), Phi(Phi(S)), ...` for one map `Phi` on subsets of the `16` blocks, they increase because each render contains the last, and the first repeat is `L_2`. The repeat comes by level 3 for every code at bases 2 and 3: at level 1 for `6`, level 2 for `75` and level 3 for `21` of the `102` orbits at base 3 (**Verified**, `count`).

**Proved.** For `k <= b^n + 1` the `k x k` render windows are exactly the `k x k` blocks of the pictures `sigma^n(B)`, `B` in `L_2`, each of side `2 b^n`. A window of a render of level above `n` meets at most two consecutive blocks of the level-`n` grid in each direction, since `k <= b^n + 1`, so it lies inside `sigma^n` of the `2 x 2` window of block letters beneath it; conversely `sigma^n(B)` is a block of `R_(m+n)` when `B` is a block of `R_m`. So off the boundary codes `p(k)` is a census of at most `16` pictures of side `2 b^n`, never of a render of side `b^level` with `level` large.

**Proved.** `p(k) <= 64 b^2 k^2` for every code and every `k >= 3`. Take `n` least with `b^n >= k - 1`; then `b^n < b k`, each picture has fewer than `(2 b^n)^2` positions, and there are at most `16` pictures. The same two-by-two supertile bound for automatic double sequences is Proposition 2.1 of [Berthe 2000](http://www.numdam.org/item/ITA_2000__34_5_403_0/).

The census counts distinct blocks without comparing them cell by cell. Give every `1 x 1` block its value as its name; a `(k+1) x (k+1)` block at `(x, y)` is fixed by the `k x k` block at `(x, y)`, the `k x k` block at `(x+1, y+1)` and its two remaining corners `(x+k, y)` and `(x, y+k)`, so its name is the index of that quadruple of integers in a table built afresh for each `k`. Equal names mean equal blocks at every size, by induction on `k`, and `p(k)` is the size of the table. One pass from `k = 1` to `k = top` costs `16 (2 b^n)^2` table lookups per size.

**Verified**, the brute-force control. The `scan` verb renders each code with the crate's own renderer, asserts it equals `sigma^level(1)` cell for cell, and counts the distinct `k x k` blocks of the render by storing every block's rows: every orbit at base 2 at level 7 for `k <= 17`, every orbit at base 3 at level 6 for `k <= 12`, and the carpet at level 7, side `2187`, for `k <= 16`. Every scan count equals the census count, `0` disagreements in `108` orbits. A render of level `n + s` holds every `k x k` window once `b^n + 1 >= k`, `s` the level at which `L_2` settles, by the census lemma; for the carpet, `s = 2`, so level 4 serves `k <= 10` and undercounts at `k = 11` and `12`, while level 5 serves `k <= 28`.

## The carpet

The carpet's window count opens `2, 10, 40, 74, 152, 242, 344, 442, 544, 650, 872, 1106` and is computed exactly to `k = 244` (**Verified**, `count 3 244 495`). Its `L_2` has `10` of the `16` blocks: every block but the two diagonals `10/01` and `01/10` and the four `L` shapes with three 0s, because the 0s of a carpet render form square holes ringed by 1s.

**Conjecture.** With `N = 3^n`, `n >= 0`,

```
p(k) = 6 k^2 + 12 (N - 1) k - 10 N^2 - 12 N + 8     for N + 1 <= k <= 2N + 1,
p(k) = 2 k^2 + (24 N - 4) k - 18 N^2 - 24 N + 4      for 2N + 1 <= k <= 3N + 1.
```

On each of the eight pieces with `n <= 3` the quadratic is the exact one through every computed value of the piece, `k` from `2` to `82`; the coefficients as polynomials in `N` are read off `n = 1, 2, 3` and then hold on both pieces of `n = 4`, every `k` from `82` to `244`, and at `n = 0` (**Verified**, `count 3 244 495`). The pieces agree where they meet: `p(N + 1) = 8 N^2 + 2`, `p(2N + 1) = 38 N^2 + 2`, and the second piece at `3N + 1` gives `72 N^2 + 2`, the value of the first piece one scale up. If the formula holds, `p(k)/k^2` does not converge: it is `8 + O(1/N)` at `k = N + 1` and reaches `9.6 + O(1/N)` near `k = 5N/3`, the top of `6 + 12/t - 10/t^2` at `t = k/N = 5/3`.

**Proved.** `p(k) = Theta(k^2)` for the carpet: the upper bound is the general one above, and the lower bound is the recognition lemma of the next section with `m = 1`, `p(k) >= k^2/36` for `k >= 6`.

## Growth

The upper bound `p(k) <= 64 b^2 k^2` holds for every code. The lower bound needs the level-`m` grid to be visible in a window.

**Proved, the recognition lemma.** Suppose `F` is neither empty nor a boundary code, and for some `m >= 1` every occurrence of the render `R_m` in a picture `sigma^m(B)`, `B` in `L_2`, is at one of the four offsets in `{0, b^m}^2`. Then `R_(jm)` occurs in the pictures `sigma^(jm)(B)` only at offsets in `{0, b^(jm)}^2` for every `j >= 1`, and `p(k) >= k^2 / (4 b^(2m))` for every `k >= 2 b^m`. By induction on `j`: `R_(jm)` is made of copies of `R_m` on the level-`m` grid, each copy lies in `sigma^m` of a `2 x 2` window of the block letters beneath it, which is in `L_2`, so the hypothesis aligns the whole occurrence to the level-`m` grid; peeling `m` levels off both sides, which `sigma^m` allows because it sends distinct letters to distinct blocks, leaves `R_((j-1)m)` inside `sigma^((j-1)m)(B)`. Now fix `n = jm` and look at the `2 b^n x 2 b^n` windows that hold a copy of `R_n` at offset `(u, v)`, `0 <= u, v < b^n`; they exist because `F` has a 1 inside some render with all eight neighbours present. Every such window lies in `sigma^n` of a `3 x 3`, hence in `sigma^(n+1)` of a `2 x 2` legal block, so by the alignment just proved two copies of `R_n` inside one window differ by multiples of `b^n` in both coordinates. Two windows with different offsets therefore differ, `p(2 b^n) >= b^(2n)`, and `p` never decreases in `k`, since every `(k+1)`-block restricts to a `k`-block and every `k`-block of a picture extends.

**Verified.** The hypothesis holds at `m = 1` for the gasket and for `91` of the `102` orbits at base 3, `476` codes; it fails at `m = 1, 2, 3` for the line codes below (`kind`, the words `recognised` and `unrecognised`). So `p(k) = Theta(k^2)` for those `480` codes (**Proved**, the lemma with that certificate).

The eleven remaining orbits at base 3 are the empty and the full code, the five boundary orbits and four line orbits; the five remaining at base 2 are the empty and the full code, the two boundary orbits and the diagonal pair. For each, `p(k) <= k^2/2` already at `k = 4` (`count`, the word `periodic`), so every picture of their `X_F` has a nonzero period by Theorem 1.2 of [Cyr and Kra 2015](https://doi.org/10.1090/S0002-9947-2015-06391-0), read at source: a picture with `P(n, k) <= nk/2` for one rectangle is periodic. Periodicity alone does not fix the order; the constant codes give `1`, and the line codes give the linear counts below.

- **Proved.** The middle row, base-3 code `56`: `X_F` is `{0}` and the pictures with one horizontal line of 1s, so `p(k) = k + 1`.
- **Proved.** The diagonal and the antidiagonal, base-2 codes `6` and `9`, base-3 codes `84` and `273`: one line of slope `-1` or `1`, so `p(k) = 2k`, the all-0 block and `2k - 1` offsets of the line.
- **Verified** to `k = 82`: the stripes of the base-3 digit set `{0, 1}`, code `63` (rows `0` and `1` filled), `p(k) = 2k`; the Cantor stripes, code `365` (columns `0` and `2` filled), `p(k) = 2k - 1` for `k >= 2`. **Proved:** both are `Theta(k)`. Every picture is constant along one axis, so `p(k)` is the window count of the one-dimensional design across it, at most `8 b k` by the supertile bound on the line, at most four two-letter words and fewer than `2 b^n` offsets each. It is at least `k + 1`: the gaps between the digit words grow without bound, so the line holds a word that ends in a 1 followed by 0s for ever, and its `k` windows that hold that last 1 at each of the `k` places differ from each other and from the all-0 window.

So every code at bases 2 and 3 has growth order `1`, `k` or `k^2`, and the order is `k^2` exactly for the `480` codes with the certificate (**Proved**, the lemma, the four line families and the certificate list).

The quadratic codes are not alike. Ten orbits at base 3 and the gasket at base 2 have a count that is one quadratic from `k = 1` or `2` to the end of the census (**Verified**, `count`, the word `quadratic`, to `k = 82` and `129`): the centre alone, code `16`, `k^2 + 1`, which is **Proved** since its `X_F` is `{0}` and the pictures with a single 1; the pairs `{(0,1), (1,0)}` and `{(0,0), (1,1)}`, codes `10` and `17`, `2k^2 - 2k + 2`; `{(0,2), (1,0)}`, code `12`, `2k^2 - 3k + 3`; `{(0,1), (1,1)}`, code `18`, `2k^2 - k + 1`; `{(1,0), (1,2)}`, code `40`, `2k^2 - 2k + 1` from `k = 2`; `{(0,2), (2,0)}`, code `68`, `2k^2 - 4k + 5` from `k = 2`; the products `27`, `45` and `325` of the next section, the centre and codes `18` and `40` being products too; and the gasket, `4k^2 - 6k + 4`. The other `81` orbits, the carpet among them, are not one quadratic from `k <= 3` on.

The gasket's count `4k^2 - 6k + 4` is Theorem 5.16 of [Allouche and Berthe 1997](https://doi.org/10.36045/bbms/1105730620), read at source, which counts the `u x v` blocks of `binom(m, n) mod 2` on `N^2` as `u^2 + v^2 + 2uv - 3u - 3v + 4`. That array is the triangle form, the square form sheared by `(i, j) -> (i + j, i)`, so a square window of one is a parallelogram of the other; the two square counts agree on every computed `k` (**Verified**, `count 2 129 7` prints the check), and this page proves no bijection between them.

## Products and lines

A code is a **product** when `F = A x C` for a set `A` of row digits and a set `C` of column digits; then every render is the outer product of the two line renders of `A` and `C` at the same level.

**Proved.** For a product code that is not a boundary code `p_F(k) = n_A(k) n_C(k) + z(k)`, where `n_A(k)` and `n_C(k)` count the `k`-windows of the two line designs that are not all 0 and `z(k)` is `1` when either line has an all-0 `k`-window and `0` otherwise. A `k x k` window of a render at `(x, y)` is the outer product of the line windows at `x` and at `y`, chosen independently, and off the boundary codes the render windows are the blocks of `X_F`; a 0-1 matrix of rank one that is not 0 determines both factors, its row support and its column support, and every product with a 0 factor is the all-0 window.

**Verified** against the census, `count`: the square of the base-3 digit set `{0, 1}`, code `27`, has `n = 2k - 1` on both axes and `p = (2k - 1)^2 + 1 = 4k^2 - 4k + 2`; the Cantor dust `{0, 2}^2`, code `325`, has `n = 2k - 2` and `p = 4k^2 - 8k + 5` for `k >= 2`; the mixed product `{0, 1} x {0, 2}`, code `45`, gives `(2k - 1)(2k - 2) + 1 = 4k^2 - 6k + 3` for `k >= 2`. The one-axis products, `A` or `C` the whole digit set, are the stripes of the last section, whose windows are the line windows.

**Proved.** If `F = A x C` and neither `A` nor `C` is empty, a single digit `0` or a single digit `b - 1`, the pictures of `X_F` are exactly the outer products `u v` of a picture `u` of the row line and a picture `v` of the column line. Every window of a picture is a product, so two 1s at `(i, j)` and `(i', j')` force a 1 at `(i, j')` through any window holding both; the support is a product `R x C'`, and the indicators of `R` and `C'` are line pictures because their windows are the factors. At base 2 every product is empty, full or a boundary code, so the reduction to the line is a base-3 matter, where the proper line sets are `{1}`, `{0, 1}`, `{1, 2}` and `{0, 2}`.

## Not of finite type

Write `X_r` for the pictures whose `r x r` blocks are all blocks of `X_F`. `X_F` is a [subshift of finite type](https://en.wikipedia.org/wiki/Subshift_of_finite_type) exactly when `X_F = X_r` for some `r`: a finite forbidden list fits in some `r x r` square. Off the boundary codes `L_2` is the set of `2 x 2` blocks of `X_F`, and `X_2` is the cut it makes alone.

**Proved, the scaling lemma.** Let `F` be neither empty nor a boundary code. If `X_F = X_r` and `b^n + 1 >= r`, then `sigma^n(Q)` lies in `X_F` for every picture `Q` of `X_2`, and `p(j b^n) >= |L_j(X_2)|` for every `j`, where `L_j(X_2)` is the set of `j x j` blocks of `X_2`. Every `r x r` block of `sigma^n(Q)` lies inside `sigma^n` of a `2 x 2` block of `Q`, which is in `L_2`, and `sigma^n` of a render window is a render window, a block of `X_F` by the first section; so `sigma^n(Q)` is in `X_r = X_F`. `sigma^n` sends distinct `j x j` blocks to distinct `j b^n x j b^n` blocks.

**Proved.** Hence if `|L_j(X_2)|` is not `O(j^2)`, `X_F` is not of finite type, since `p(k) <= 64 b^2 k^2`. The `kind` verb certifies the growth of `L_j(X_2)` from `L_2` by one of three finite witnesses.

- `blocks`: two `s x s` blocks, `s <= 2`, such that every `2 x 2` window of each of the `16` arrangements of them in a `2 x 2` square is in `L_2`. Every arrangement of the two on the whole `s`-grid is then in `X_2`, so `|L_(sj)(X_2)| >= 2^(j^2)`.
- `lines`: the pictures constant along rows, along columns, along diagonals or along antidiagonals form a window rule on the line, and its graph has a strongly connected part with more edges than vertices. That part carries two distinct cycles through one vertex, so free concatenations of the two give at least `2^(floor(j/c))` words of length `j`, `c` the longer cycle, and distinct words give distinct `j x j` blocks of `X_2`.
- `xor`: `L_2` is exactly the eight blocks in which one corner is the sum mod 2 of its two neighbours. Then the two sides of a block that meet at the opposite corner are free and fix the rest: for the gasket's rule, bottom-right equals top-right plus bottom-left, any row extends downward with one free bit per row and upward uniquely, so `|L_j(X_2)| = 2^(2j - 1)`.

**Proved, with the certificate.** At bases 2 and 3, `X_F` is of finite type exactly when the code is empty, full or a boundary code, and `X_F` is then a single fixed picture, `{0}` or `{1}`. Every other code carries a witness (**Verified**, `kind 2` and `kind 3`): at base 2, `10` codes are of finite type, the gasket's `4` codes carry `xor` and the diagonal pair `lines`; at base 3, `26` are of finite type, `468` carry `blocks` and `18` carry `lines`, the line codes of the last sections and code `87`, the top row with the centre and the bottom-left corner.

The carpet's witness is the pair `11/00` and `11/10`: lay these two `2 x 2` blocks on the even grid in any pattern and every `2 x 2` window is among the carpet's ten. So `X_2` holds at least `2^(j^2)` blocks of side `2j`, and the carpet's subshift is not of finite type: no `X_r` equals it. The same fact without the lemma: in a carpet render the 0s form square holes ringed by 1s, so no picture of `X_F` has an infinite horizontal strip of 0s of finite height between two rows that hold 1s; yet the limit `x` at the top edge of ever larger holes, carpet above and 0s below, and the limit `y` at their bottom edge, 0s above and carpet below, shifted to share a band of `r` rows of 0s, glue into a picture all of whose `r x r` blocks are blocks of `x` or of `y`, and it has exactly such a strip.

`X_F` holds the all-0 picture for every code but the full one, so it is not minimal once it holds a 1, and criteria built on minimal subshifts do not apply; the scaling lemma uses growth instead.

**Verified**, the forbidden-list control. A forbidden list of blocks up to `10 x 10` does not cut out the carpet's subshift, read off the census alone: `X_F = X_2` would force `p(6) >= 2^9`, `X_F = X_r` with `r <= 4` would force `p(24) >= 2^16`, and `r <= 10` would force `p(72) >= 2^16`, against `p(6) = 242`, `p(24) = 4570` and `p(72) = 42970` (`count 3 244 495`).

The gasket's `X_2` is a familiar shift. Its rule `a(i+1, j+1) = a(i, j+1) + a(i+1, j) mod 2` reads, along the antidiagonals `t = i + j` with place `i - j`, as one step of rule 90 of [the automata](../../automata/notes/automata.md), each cell the sum of its two neighbours. So `X_2` is the set of every two-sided history of rule 90 on one parity class of the plane, a set of finite type with `2^(2j - 1)` blocks of side `j`, while the gasket's `X_F` holds only the pictures whose blocks all occur in the forward history of a single seed.

## Sofic

A subshift is **sofic** when it is the image of a subshift of finite type under a letter-to-letter map.

**Verified**, read at source. Theorem 4.5 of [Mozes 1989](https://doi.org/10.1007/BF02793412) gives, for a two-dimensional substitution system with (1) every rule of height and width at least 2 and (2) property A, a tiling system of finite type and a factor map from it onto the system's subshift, made letter to letter by the modification printed just before it, the pictures whose every block is a sub-block of a block legally derived from one letter. The statements of Theorems 4.1 to 4.5 and the definitions of section 1 are read in an open scan; the proofs are not.

**Proved.** Every plane code at every base `b >= 2` meets both hypotheses, so its `X_F` is sofic. Both rules, `1` to the level-1 tile and `0` to the empty tile, are `b x b`. Property A asks that every derivation sequence of a `2 x 2` sub-block `L` of a derived block `U` be visible inside some derivation sequence of `U`; here each letter has one rule and all rules share one size, so every block has exactly one derivation, `sigma`, and the derivation of `U` shows `sigma^i(L)` at `b^i` times the place of `L`, at every occurrence, the all-0 quadruples included. Mozes' subshift also admits blocks derived from `0`, which are all-0 blocks and occur in a render of every code but the full one; for the full code `X_F = {1}` is of finite type anyway. [Aubrun and Sablik 2014](https://arxiv.org/abs/1103.0895), Remark 2, read at source, states the same: a single deterministic substitution has property A.

So `X_F` is strictly sofic, sofic and not of finite type, for the `6` codes at base 2 and the `486` at base 3 that carry a witness, and no code at either base is neither.

## The Wang face

A set of [Wang tiles](https://en.wikipedia.org/wiki/Wang_tile) is a subshift of finite type with the tiles as letters, and any two-dimensional subshift of finite type is a Wang set after recoding by blocks. The Wang face of a design is the least Wang set whose tilings, after a map from tiles to `{0, 1}`, are exactly the pictures of `X_F`.

**Proved.** For the `492` codes that carry a witness the map cannot be a recoding: being of finite type is kept by every conjugacy, so no Wang set is conjugate to `X_F`, and the face needs tiles that carry hidden information beyond the picture. For the `36` codes of finite type the face is one tile with four edges of one colour, mapped to `0` or to `1`.

**Proved.** For every code a finite face exists: `X_F` is sofic by the last section, the finite-type cover of Theorem 4.5 is a Wang set after recoding, and its letter map composes with the recoding.

The least face is open at every code with a witness, the carpet and the gasket included: Mozes' cover is built for every substitution at once and is far from least, and no lower bound beyond the trivial one is proved here. For the gasket a natural start is the eight tiles that lay the histories of rule 90, whose tilings read on one parity class and turned by half a right angle are its `X_2`; what must be added to forbid every history outside the closure of one seed is the open part.

## The census

At base 2 the `16` codes fall into `6` orbits (**Verified**, `count 2 129` and `kind 2`).

| code | orbit | cells | `X_F` | `p(k)` | type |
|---|---|---|---|---|---|
| `0` | 1 | none | `{0}` | `1` | finite |
| `1` | 4 | one cell | `{0}`, boundary | `1` | finite |
| `3` | 4 | one side | `{0}`, boundary | `1` | finite |
| `6` | 2 | a diagonal | `{0}` and one line | `2k` | sofic, `lines` |
| `7` | 4 | three cells, the gasket | the gasket's closure | `4k^2 - 6k + 4` | sofic, `xor` |
| `15` | 1 | all | `{1}` | `1` | finite |

At base 3 the `512` codes fall into `102` orbits and give `83` distinct count sequences to `k = 82` (**Verified**, `count 3 82`). By codes: `26` give the constant count `1` and are of finite type, `10` in `4` orbits give a linear count and are line codes, and `476` in `91` orbits give a count above `k^2/2` at every computed `k`, carry the recognition certificate and so grow as `k^2`; all `486` non-constant codes are strictly sofic. The gasket's count is twice [A084849](https://oeis.org/A084849) shifted by one, `2 (2(k-1)^2 + (k-1) + 1)`.

## Where the numbers live

Every number on this page is printed by `lab/rs/design-windows`, a study in the root workspace.

- `count <base> <top> [code]`: the census of the second section, one line per orbit with its window count to `top`, the exact quadratic when there is one, the carpet's and the gasket's formula check, and the class tallies; `count 3 244 495` is the carpet's run.
- `scan <base> <level> <top> [code]`: the brute-force control on a crate render.
- `kind <base>`: the finite-type verdict per orbit with its witness, and the recognition certificate.
- `pairs <base> <code>`: the set `L_2`.
