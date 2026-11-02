---
title: The rounds
lead: Every design round on the torus is a bijection, its Feistel step an involution, and an even permutation once it writes two cells or more; a census of every split of the 3 x 3 and 4 x 4 torus, 7176 cells over three neighbourhoods, four Sources and six spreads, reads which rounds generate the alternating group with every stamp, 147 of 3240 nonlinear cells at 3 x 3 and 198 of 2142 at 4 x 4, five proved design laws account for every other nonlinear cell, and a single round is far from a random permutation, 136 cycles on average against the 6.8 of a uniform permutation of 512 states.
figure: research-random-rounds
---

A generator on a grid of `N` bits permutes the `2^N` grids, so two questions precede any statistical test: which permutations its rounds can make, and whether one round looks random. This chapter asks both of the rounds of the [random spine](random.md), for every design split of the `3 x 3` and `4 x 4` torus over three neighbourhoods, four Sources and six spreads, 4320 cells at `3 x 3` and 2856 at `4 x 4`. Every round is a bijection, and an even permutation once it writes two cells or more, so rounds and stamps generate a subgroup of `A_M` unless one round writes a single cell with a mess of odd weight, which happens in 9 cells of the catalog, all at `3 x 3`. The group of a round and every stamp is the alternating group in 147 of the 3240 nonlinear `3 x 3` cells and in 198 of the 2142 nonlinear `4 x 4` cells, and five proved laws, each a block system or a product action read off the design, explain every other nonlinear cell, with no cell undecided. A single round is nowhere near a uniform permutation: 136.43 cycles and 28.43 fixed points on average at `3 x 3`, against `H_512 = 6.8165` and 1.

One generator, `lab/rs/rounds`, prints every number below. The figure draws the whole cycle structure of one alternating `3 x 3` cell, the split `0x005f` with the net neighbourhood, `primes` and `slide`: its 512 grids are laid along their cycles, longest first, each cycle on a ground of its own ink, blue and orange in turn, with the live cells in the foreground ink. The round has 62 cycles, a longest cycle of 36 states and 4 fixed points, as `blocks` prints and the figure recounts, where a uniform permutation of 512 points has 6.8 cycles and 1 fixed point on average.

## The machine

The grid is the torus `Z_n x Z_n`, `n >= 3`: `N = n^2` cells, a state `x` in `V = F_2^N`, `M = 2^N` states. A **split** is a cell set `A`, the fill or written side, with complement `B`, the void or read side.

The **mess** writes `A` from `B`. A neighbourhood is a base-3 level-1 design with its centre popped, in the code of [core](../../core/notes/core.md): `carpet 495` (the eight Moore offsets), `net 186` (the plus) or `void 341` (the diagonals). A cell of `A` counts the live cells of `B` around it and writes `[count in S]`, `S` a Source of [automata](../../automata/notes/automata.md): `primes {2,3,5,7}`, `binary {1,2,4,8}`, `fibonacci {0,1,2,3,5,8}` or `odds {1,3,5,7}`. `odds` makes the mess linear and is the control.

A **round** is `R = t_c L F`:

- the Feistel step `F(x) = x + mess(x)`, which changes only `A` and reads only `B`;
- a spread `L`, a torus map `sigma` moving the bit at `p` to `sigma(p)`: `id`, `shift (i, j+1)`, `drop (i+1, j)`, `slide (i+1, j+1)`, `turn (j, n-1-i)`, `flip (j, i)`;
- a stamp `t_c(x) = x + c` by a design; the census stamps the split's own design, `c = A`.

The undo is the same steps reversed, `R^-1 = F L^-1 t_c`.

The catalog is enumerated, never picked: at `3 x 3` every split `1..510`, the nonempty proper base-3 level-1 codes; at `4 x 4` every two-letter base-2 word `k1 (x) k2`, `k1, k2` in `1..15`, the full word dropped, 224 masks. A split is kept once per orbit of the torus symmetries that commute with the spread: over three neighbourhoods, four Sources and six spreads, 4320 cells at `3 x 3` and 2856 at `4 x 4`.

## The sign law

**Every round is a bijection, its Feistel step an involution, and the undo its inverse. Proved.** The mess depends only on `x_B`, which `F` keeps, so `F(F(x)) = x`; `L` and `t_c` are bijections, `t_c` its own inverse. The census checks all three on every state: 0 faults in 7176 cells.

**Every round that writes two cells or more is an even permutation once `N >= 3`. Proved.**

- A stamp `t_c`, `c != 0`, is `2^(N-1)` transpositions: even once `N >= 2`.
- A spread is a product of coordinate swaps, each `2^(N-2)` transpositions: even once `N >= 3`.
- The Feistel step acts on each fibre `{x : x_B = y}` as the translation by `mess(y)`, so its sign is `(-1)^(w 2^(|A|-1))`, `w` the number of `y` with a nonzero mess: even once `|A| >= 2`, and `(-1)^w` when `|A| = 1`.

So rounds and stamps generate a subgroup of `A_M` unless a round writes one cell with a mess of odd weight: the parity fact of Coppersmith and Grossman 1975 and Even and Goldreich 1983 ([refs](../refs.md)) in design form, proved from scratch. One written cell reading `k` cells of `B` has weight `w = 2^(|B| - k) sum_{s in S} C(k, s)`, odd only when it reads all of `B`.

**In the catalog that happens only at `3 x 3`, with the carpet neighbourhood and `binary`: 9 odd cells, none at `4 x 4`. Proved; Verified.** At `4 x 4` a cell reads at most 8 of the 15 others. The `3 x 3` carpet reads all 8, and `sum_{s in S} C(8, s)` is odd only for `binary`, `107`, against `148`, `150`, `128` for `primes`, `fibonacci`, `odds`.

**The even null matches the uniform null in mean cycles and fixed points. Proved; Verified.** The fair null is a uniform even permutation: `E[cycles | even] = H_M + (-1)^M / (M (M - 1))`, off by `3.8e-6` at `M = 512`; for `M >= 3` the even permutations fixing a point are the `(M - 1)!/2` even permutations of the other `M - 1` points, so the mean fixed-point count is `M ((M - 1)!/2) / (M!/2) = 1`; `null` checks both over every even permutation of `M = 2..8` points.

## The cycle census

A uniform permutation of `M` points has `H_M` cycles on average, `6.8165` at `M = 512` and `11.6676` at `M = 65536`, one fixed point, and a longest cycle covering `0.62494` of the states at 512 (exact), tending to the Golomb-Dickman constant `0.62433` (A084945).

**A single design round is far from a random permutation. Verified.**

| cells | count | cycles, mean | cycles, median | longest/M | fixed points |
|---|---|---|---|---|---|
| `3 x 3` nonlinear | 3240 | 136.43 | 116 | 0.05378 | 28.43 |
| `3 x 3` nonlinear, a symmetry fixes the split | 1476 | 184.80 | 160 | 0.02587 | 53.38 |
| `3 x 3` odds | 1080 | 96.77 | 72 | 0.06240 | 12.77 |
| `4 x 4` nonlinear | 2142 | 16759.65 | 14236 | 0.01298 | 3081.50 |
| `4 x 4` nonlinear, a symmetry fixes the split | 1773 | 18556.72 | 16336 | 0.00816 | 3706.66 |
| `4 x 4` odds | 714 | 14295.61 | 10304 | 0.00685 | 1092.20 |

A torus symmetry fixing the split commutes with the round and roughly doubles the cycle count.

## The group

With every stamp allowed, rounds generate `G = <T, h>`, `T` the `M` translations and `h = L F` the round without its stamp. `T` is transitive, so every `G` is.

**A group containing every translation is imprimitive exactly when some proper nonzero subspace `U` has `h(x + u) + h(x)` in `U` for all `u` in `U`, all `x` and every generator `h`. Read at source.** This is the subspace fact of Caranti, Dalla Volta and Sala 2009, Section 2.1: a block through `0` is a subgroup and the blocks are its cosets. The smallest block through `0` and `v` is the closure of `{v}` under derivatives.

The census computes that closure for every `v != 0`, so a primitive verdict is a proof. A primitive cell is then searched for an element whose only cycle of length divisible by a prime `p <= M - 3` has length `p`; a power of it is a `p`-cycle, and by Jordan's theorem the group is `A_M` or `S_M`.

## The law in the design

Five proved laws force `G != A_M`; for nonlinear cells in both catalogs the census shows they are the only way to miss.

- **I, a read-closed set. Proved.** A nonempty proper union `O` of `sigma`-orbits holding every cell of `A` whose mess reads a cell of `O` makes `F_2^O` a block; the identity spread always has `O = A`.
- **II, the mess module. Proved.** Cells of `A` reading the same live cells of `B` write the same bit, so every derivative of the mess lies in the span `E` of those classes; the smallest `sigma`-invariant subspace containing `E` is a block if proper, and a constant mess makes the round affine.
- **III, the complement. Proved.** If each cell's rule is self-complementary, `[k in S] = [m - k in S]` for a cell reading `m` cells, the all-ones line is a block.
- **IV, the parity of a design partition. Proved.** For a partition of the cells, the grids even on every part form a block exactly when the part parities after a round depend only on those before it. The census tries the whole torus, the cosets of every torus subgroup and, at `4 x 4`, the `2 x 2` tiles.
- **V, a product. Proved.** If the cells split into `k` cosets of a torus subgroup that `sigma` permutes and each cell of `A` reads only inside its coset, `G` lies in `Sym(2^(N/k)) wr S_k`, far below `M!/2`: primitive perhaps, not `A_M`.

**At `3 x 3` the group of every nonlinear cell is read off the laws. Verified.** Of 3240 nonlinear cells, 147 generate `A_512`, 6 are primitive products and 3087 are imprimitive, the first law holding being I for 2901, II for 183 and III for 3. None is imprimitive without a law or primitive and undecided, and none generates `S_512`: the 9 odd cells are imprimitive. The 1080 `odds` cells are affine, so never `A_512`: 1068 are imprimitive and 12 primitive, those 12 meeting none of the five laws, each run on every cell. Only `shift`, `drop` and `slide` give `A_512`, 49 cells each.

**At `4 x 4` the same holds once the parity law is added. Verified.** Of 2142 nonlinear cells, 198 generate `A_65536`, 3 are primitive products and 1941 are imprimitive, the first law holding being I for 1802, II for 113, III for 8 and IV for 18; all 714 affine `odds` cells are imprimitive. `turn`, which fixes no cell, gives `A_65536` in 75 cells.

**A split and its complement written in turn reach the alternating group in a few cells. Verified.** With the identity spread, the 12 splits of `3 x 3` up to symmetry and complement give `A_512` in 12 of 108 nonlinear cells and the 23 words of `4 x 4` give `A_65536` in 10 of 207; every other cell is imprimitive.
