---
title: Structure against noise
lead: Self-similar designs raced against matched random cell sets on components and boundary.
figure: research-connectivity
slug: connectivity
---

A design's fractal occupies some number of cells of a grid. A random set can occupy exactly the same number of cells of the same grid. Race the two on measurable geometry - how many connected pieces, how much boundary per cell - and the question of what the [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) rule buys gets a number instead of an adjective.

The race is measurement at finite sizes; the diagonal's boundary and the boundary law of every design are Proved and flagged where they occur. The [graphs demo](/graphs/) draws the cell network of any design live, flat, in the cube and on the diagonal slice, with its tips, junctions, pieces and length beside it.

## The race

Five designs, each built by substituting its level-1 tile into itself and cross-checked against the digit rule that defines it - codes and corner order as in [the core](../../core/notes/core.md):

| design | base | `dim` | rule |
|---|---|---|---|
| gasket | 2 | 2 | code 7, corners `(0,0)`, `(0,1)`, `(1,0)` |
| diagonal | 2 | 2 | code 9, corners `(0,0)`, `(1,1)` |
| seven-of-eight | 2 | 3 | every corner but one |
| carpet | 3 | 2 | 8 of 9 cells - `bang dim 2, code 7` read at base 3 |
| sponge | 3 | 3 | the Menger sponge, 20 of 27 cells |

The control is matched exactly, not approximately: a uniform sample of exactly `N` distinct cells, the occupied count asserted equal to the design's on every draw - not a Bernoulli field with matching probability. Components are face-connected, 4-neighbour in 2D and 6-neighbour in 3D. Boundary per cell is `(2*dim*N - 2*E)/N` with `E` the number of face-adjacent occupied pairs and the grid exterior counting as void, so an isolated cell contributes the full `2*dim`. Every random figure below is a mean and a sample standard deviation: 400 seeds up to 729 cells, 200 at the 2187-cell `128 x 128` gasket, the 2401-cell `16^3` seven-of-eight and the 128-cell `128 x 128` diagonal, 100 at 6561 cells and at the `27^3` sponge, 20 at the `81^3` sponge (`lab/py/percolation-race`).

Two independently written passes race the designs inside one program, `lab/py/percolation-race` - [Kronecker](https://en.wikipedia.org/wiki/Kronecker_product) powers with union-find and a PCG64 control, and substitution over coordinate sets with breadth-first search and a Mersenne Twister control - and agree on every comparison to within a standard deviation; the widest gap printed is 0.2061 standard deviations. The numbers printed here are the second pass.

## Connectivity

Every self-similar design above except the diagonal - built to scatter, and raced on boundary instead - is a single component at every size tested. The matched random set is not close, and the gap widens with size in every family. (**Verified**, `lab/py/percolation-race`.)

| design | grid | cells | density | design comps | random comps |
|---|---|---|---|---|---|
| gasket | `32 x 32` | 243 | 0.2373 | 1 | `134.39 +/- 7.34` |
| gasket | `64 x 64` | 729 | 0.1780 | 1 | `478.09 +/- 12.62` |
| gasket | `128 x 128` | 2187 | 0.1335 | 1 | `1613.61 +/- 21.96` |
| gasket | `256 x 256` | 6561 | 0.1001 | 1 | `5257.23 +/- 32.94` |
| seven-of-eight | `8^3` | 343 | 0.6699 | 1 | `2.80 +/- 1.43` |
| seven-of-eight | `16^3` | 2401 | 0.5862 | 1 | `23.77 +/- 4.85` |
| carpet | `27 x 27` | 512 | 0.7023 | 1 | `9.43 +/- 2.89` |
| carpet | `81 x 81` | 4096 | 0.6243 | 1 | `146.58 +/- 12.52` |
| sponge | `27^3` | 8000 | 0.4064 | 1 | `560.80 +/- 22.75` |
| sponge | `81^3` | 160000 | 0.3011 | 1 | `31576.40 +/- 152.71` |

The largest-component fraction says the same thing from the other side: the design holds 1.0000 of its cells in one piece in every row, while random holds `0.0411 +/- 0.0098` at the `32 x 32` gasket, 0.0011 at `256 x 256`, and 0.0169 at the `81^3` sponge. (**Verified**, `lab/py/percolation-race`.) The result holds at base 2 and base 3, in two dimensions and three, and no gap narrows as the grid grows, which rules out a small-size artifact at the sizes reached.

## Boundary

At the same matched cell counts, the connected designs also expose less surface than random, and one design - built to scatter - exposes the most surface possible. (**Verified**, `lab/py/percolation-race`, except the Proved lines.)

| design | grid | design boundary | random boundary | maximum |
|---|---|---|---|---|
| gasket | `32 x 32` | 2.0082 | `3.0814 +/- 0.0639` | 4 |
| gasket | `256 x 256` | 2.0003 | `3.6006 +/- 0.0103` | 4 |
| diagonal | `32 x 32` | 4.0000 | `3.8808 +/- 0.0897` | 4 |
| diagonal | `128 x 128` | 4.0000 | `3.9677 +/- 0.0232` | 4 |
| seven-of-eight | `8^3` | 1.8542 | `2.4941 +/- 0.0589` | 6 |
| carpet | `81 x 81` | 0.8633 | `1.5330 +/- 0.0124` | 4 |
| sponge | `27^3` | 2.2560 | `3.6535 +/- 0.0141` | 6 |

The diagonal at level `level` is `2^level` isolated cells on a `2^level x 2^level` grid - no two cells of the pattern are ever face-adjacent - so its boundary per cell is exactly `2*dim = 4`, the theoretical maximum, at every level. (**Proved**, by construction; the exact 4.0000 is the measured confirmation.) The random set at the same density falls short of the maximum at every size tested.

One number in this table deserves its own flag: the gasket's 2.0082 is a level-5 figure, not a constant. It drifts 2.0082, 2.0027, 2.0009, 2.0003 at sides 32, 64, 128, 256, approaching 2, and the drift is exact: the next section proves it.

## The boundary law

Boundary per cell is `(2*dim*N - 2*E)/N`, so the drift is the count `E` of face-adjacent pairs at level `k`, and it has a closed form for every design of the race. (**Proved**; every number below **Verified** by `lab/py/boundary-drift`.)

**Contact lemma.** Take a tile of base `b` with `f` cells, `A` face-adjacent pairs of tile cells, and let `c` be the number of positions `u` on a face such that `(b-1, u)` and `(0, u)` are both tile cells, the same `c` on every axis. Then level `k` has `E_k = A (f^k - c^k)/(f - c)` and boundary per cell `2*dim - 2A/(f - c) + (2A/(f - c)) (c/f)^k`.

- Level `k+1` is `f` copies of level `k` at the tile cells, so `E_(k+1) = f E_k` plus the pairs across copies.
- Two cells of different copies are face-adjacent only when the copies are tile-adjacent along an axis, one on the high face and the other on the low face, at the same transverse position `u`; each of the `A` tile pairs contributes `c_k` such positions.
- A position is on the high face of level `k` iff every base-`b` digit pair `(b-1, u_j)` is a tile cell, and on the low face iff every `(0, u_j)` is; both iff every digit of `u` is one of the `c` shared positions, so `c_k = c^k`.
- Hence `E_(k+1) = f E_k + A c^k` with `E_0 = 0`, whose solution is `A (f^k - c^k)/(f - c)`, and `N = f^k`.

**Gasket.** `f = 3`, `A = 2`, `c = 1` (the corner copy touches each of the other two in one pair, because the high and low faces meet in one cell, and those two meet only at a point): `E_k = 3^k - 1 = N - 1`, and the set is connected by induction, three connected copies joined by two pairs, so every level is a tree and the boundary per cell is exactly `2 + 2/3^k`. This is the most a face-connected set can expose: face-connected means `E >= N - 1` (a spanning tree has `N - 1` pairs), so the boundary of face-connected sets in the plane of `N` cells is at most `4N - 2(N - 1) = 2N + 2`, with equality exactly for a tree. The gasket's `2.0082, 2.0027, 2.0009, 2.0003` are `k = 5..8`, and `3.0814` and `3.6006` for random sets of the same sizes sit above a ceiling that every face-connected set of those sizes obeys, so on the gasket the boundary win over random is connectivity again, not a second effect.

**Carpet.** `f = 8`, `A = 8` (the eight copies form a ring), `c = 3` (each side of the carpet is a full row, so neighbours touch along `3^k` pairs): `E_k = (8/5)(8^k - 3^k)`, perimeter `(4/5) 8^k + (16/5) 3^k`, and boundary per cell `4/5 + (16/5)(3/8)^k`. The limit `4/5` sits below a tree's `2` because `E/N` tends to `8/5`, not to 1: the carpet is thick and full of cycles, and every pair beyond a tree's `N - 1` removes 2 from the perimeter. The printed `0.8633` is `k = 4`.

| design | `f` | `A` | `c` | boundary per cell | limit | table value |
|---|---|---|---|---|---|---|
| gasket | 3 | 2 | 1 | `2 + 2/3^k` | 2 | 2.0082 at `k = 5` |
| carpet | 8 | 8 | 3 | `4/5 + (16/5)(3/8)^k` | `4/5` | 0.8633 at `k = 4` |
| seven-of-eight | 7 | 9 | 3 | `3/2 + (9/2)(3/7)^k` | `3/2` | 1.8542 at `k = 3` |
| sponge | 20 | 24 | 8 | `2 + 4 (2/5)^k` | 2 | 2.2560 at `k = 3` |
| diagonal | 2 | 0 | 0 | 4 | 4 | 4.0000 |

The three-dimensional rows use the same lemma: a face of the sponge tile is the ring of 8 cells (`c = 8`, 24 adjacent pairs among its 20 cells), and seven-of-eight has 9 adjacent pairs and a shared face of 3 cells, the high face being an `L` and the low face the full square. Every row is rebuilt at levels 1 to 8 (levels 1 to 5 for the sponge) from the Kronecker power and, at the levels listed in the lab, from the digit rule, with the cross-copy pair count `A c^(k-1)` checked directly and the exact rational boundary equal to the closed form at every level.

**Other trees.** Every design whose levels are all trees has `E_k = f^k - 1` (the lemma covers those with equal `c`; the straight lines have `c = (0,1)` and are checked on the render), so the boundary per cell in two dimensions is `2 + 2/f^k`. That is code 79 (`f = 5`, `2 + 2/5^k`), the `L` of base 4 (`f = 7`), the `L` of base `b` (`f = 2b - 1`) and every tree code of [the walk census](../../walks/notes/trees.md); the generator finds the same 8 tree codes at base 2 and 15 at base 3 by rendering and checks the formula on each.

**Words.** The recursion also closes the mixed product of [the next chapter](order-in-the-mixed-product.md). For a word `w = (a_1, ..., a_n)`, first letter outermost, the faces of the inner word are digit-wise products of the letters' faces, so `E(w) = f_1 E(a_2 ... a_n) + A_1 c_2 ... c_n` on each axis, and `E(w) = sum_i f_1 ... f_(i-1) A_i c_(i+1) ... c_n`, summed over the two axes, with the outer letters multiplying the copies and the inner letters multiplying the contact. (**Proved**; equal to the drawn count on all 54240 words of length at most 4 over the 15 codes of base 2, **Verified**, `lab/py/boundary-drift`.) The cell count `f_1 ... f_n` is order-blind and `E` is not: orderings of one multiset differ in `E`, and so in perimeter `4N - 2E`, for 78 of the 105 multisets of length 2 and 188 of the 210 of length 3 on the ten-code library, the counts of the Euler column there. The boundary row of that table counts border cells, a different observable.

## How decisive, honestly

Means hide ties, so the per-draw counts are part of the result (`lab/py/percolation-race`).

- Gasket, `32 x 32`: 0 of 400 random draws reach one component, and 0 of 100 at `256 x 256`. Decisive at every size.
- Seven-of-eight, `8^3`: 73 of 400 random draws are *also* a single component, and random's largest fraction is already 0.9937 - at density 0.67 a random set percolates too, so the 3D win at that size is a win in the mean, not per draw. One level up, at `16^3`, it is 0 of 200 and decisive.
- Diagonal, `32 x 32`: 61 of 400 random draws also hit boundary 4.0000, because at density 0.031 a random set is often already an independent set. The design wins the mean by 0.119, about 1.3 standard deviations; at `128 x 128` the ties are 27 of 200. A real edge, and a small one.

And one scope line the summary sentence invites getting wrong: at the dispersing extreme, random has *fewer* components than the design - 30.10 against 32 at `32 x 32`, 125.93 against 128 at `128 x 128` - because the diagonal is `N` isolated cells by construction. "Random wins at neither extreme" is true only of the metric named at each extreme: connectivity at the percolating one, boundary at the dispersing one. (**Verified**, `lab/py/percolation-race`.)

## Where the honest line falls

What is reached: grids to `256 x 256` in 2D and `81^3` - 160000 occupied cells - in 3D, with the seed count thinning from 400 to 20 at the largest size; one control ensemble, uniform over cell sets of the exact matched size. The race is no limit theorem; only the boundary law, which is exact at every level, reaches one. The claim the race supports is finite and plain: at matched grid and matched cell count, over the five designs and the sizes tested, self-similar structure holds together where noise shatters, and the one design built to scatter, scatters perfectly. One program, `lab/py/percolation-race`, runs both passes, prints every figure above in about nine seconds and exits nonzero on any failed check; no log is kept. `lab/py/boundary-drift` rebuilds the boundary law in about four seconds.

## Chapters

- [Order in the mixed product](order-in-the-mixed-product.md): Replacing the repeated tile by a word of different tiles, which observables depend on the order of the word, and the component exponent along a word.
- [Word order, closed and bounded](word-order-closure-and-spectral-radius.md): The remaining word pairs closed by a suffix recursion, what the invariant cone buys, and the joint spectral radius of the cocycle.
