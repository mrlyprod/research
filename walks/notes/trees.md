---
title: Trees
lead: Two designs of the walk census are trees, so their walk dimensions are exact: log 15 / log 3 for code 79 and log 6 / log 2 for the base-2 design.
figure: research-walks
slug: walks
---

## Trees: the walk dimension exactly

Two rows of [the census](walks.md) are not measurements waiting for a limit. `bang dim 2, base 3, code 79` and `bang dim 2, code 7` (base 2, 3 of 4 cells) draw graphs with no cycle, at every level. On a tree the electrical resistance between two cells is their graph distance, so the resistance scale per level is read off the picture, and a theorem of Barlow, Coulhon and Kumagai turns volume growth plus that resistance into the walk dimension: `d_w = log 15 / log 3 = 2.464974` and `d_w = log 6 / log 2 = 2.584963`. Both sit on Barlow's ceiling `d_w <= 1 + d_f`. This chapter proves that, recomputes every number it uses, and lists the sibling designs the same argument closes. The cross (code 186) is Barlow's Vicsek tree `Gamma(2, 1)` (Barlow 2004, Lemma 4.2), so its walk dimension is classical; the L designs and the other row-plus-column trees are new instances of the same count.

Notation: `b` is the base, `f = 2b - 1` the fill, `P_n` the level-`n` graph (one node per filled cell of the `b^n` by `b^n` grid, one edge per face-adjacent pair, side `s = b^n`). Cell `(a, c)` of a tile is bit `b*a + c` of the code. Code 79 (`b = 3`, `f = 5`) fills the cells `(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)` and code 7 of base 2 (`b = 2`, `f = 3`) the cells `(0, 0), (0, 1), (1, 0)`: the tile is a full row and a full column meeting at `(0, 0)`, an L with arms of length `b`. Call it `L_b`. The corners of `P_n` that matter are `O = (0, 0)`, `A = (0, s - 1)` and `B = (s - 1, 0)`; the corner `(s - 1, s - 1)` is empty. (Definitions.)

## Two designs are trees

**Lemma 1.** `P_n` is a tree with `f^n` nodes. `P_n` contains the whole row `{0} x [0, s)` and the whole column `[0, s) x {0}`, and its last row `{s - 1} x [0, s)` and last column `[0, s) x {s - 1}` each hold one cell, `B` and `A`. (**Proved**, by induction.)

*Proof.* `P_{n+1}` is `f` copies of `P_n` at the tile cells. Two cells of different copies are face-adjacent only if the copies sit at face-adjacent tile cells and both cells lie on the shared face. Take copies at `(a, c)` and `(a + 1, c)`: the first shows its last row, the single cell `B`; the second shows its first row, which is full; they meet in exactly one pair, `B` of the first with `O` of the second. Copies at `(a, c)` and `(a, c + 1)` meet in exactly one pair, `A` of the first with `O` of the second. The tile cells themselves form a tree with `2b - 2` edges, a row path and a column path sharing `(0, 0)`. So `P_{n+1}` is `f` trees joined along a tree by one edge per tile edge, which is a tree. The rows and columns of `P_{n+1}` are the rows and columns of its copies at the tile cells `(0, c)` and `(a, 0)`.

Pinned to the render: `lab/py/tree-walks` rebuilds both designs cell by cell from the digit rule, and its digests of the cell list and the edge list equal those printed by `lab/rs/tree-walks` from the crate's own render and graph at all 15 pairs (code 79 at levels 1 to 6, code 7 at levels 1 to 9); every level has `edges = nodes - 1` and is connected. (**Verified**, `lab/py/tree-walks` verb `render`, `lab/rs/tree-walks`.)

**Lemma 2.** `d(O, A) = d(O, B) = s - 1`, `d(A, B) = 2 (s - 1)`, every cell is within `s - 1` of `O`, and the diameter of `P_n` is `2 (s - 1)`. (**Proved**.)

*Proof.* The row and the column are paths, so the first two hold. A tree has one simple path between two cells, and `A, O, B` is simple, so `d(A, B) = 2 (s - 1)`. For the radius, let `e_n` be the largest distance from `O` in `P_n`. In `P_{n+1}` the copy at `(0, c)` is reached by crossing `c` copies from their `O` to their `A`, `c (s - 1)`, plus `c` joining edges, then at most `e_n` inside; the copies at `(a, 0)` are the same. The largest is `c = b - 1`: `e_{n+1} = (b - 1) s + e_n`, `e_0 = 0`, so `e_n = b^n - 1`. The diameter is at most `2 e_n` and `A, B` attains it.

## Resistance and distance

In a tree the resistance between two cells equals their graph distance: the unique path is the only route, and the potential that is linear along it and constant on the branches has energy `1/d`, which gives `R >= d`, while the series bound gives `R <= d` (Barlow, Coulhon and Kumagai, Lemma 5.1, read at source; [refs](../refs.md)). Hence the corner-to-corner resistance:

`R(A, B) = 2 b^n - 2`, so `2 * 3^n - 2` for code 79 and `2^(n+1) - 2` for code 7, and `R(O, A) = R(O, B) = b^n - 1`. (**Proved**, Lemma 2 and Lemma 5.1.)

The generator computes `R(A, B)` without assuming a tree: exact star-mesh elimination over rationals on the whole graph, then the graph distance by breadth-first search, then a float sparse solve of the Laplacian system. Code 79, levels 1 to 6: `4, 16, 52, 160, 484, 1456`. Code 7, levels 1 to 9: `2, 6, 14, 30, 62, 126, 254, 510, 1022`. `R(O, A)` and `R(O, B)` are half of each. All three routes agree at every level, the float solve to `7.3e-8` at the 19683-node level and better below. (**Verified**, `lab/py/tree-walks` verb `resist`.)

A finite-level identity needs no limit and no outside theorem. For any finite weighted graph the commute time between two nodes is `E_x T_y + E_y T_x = mu(G) R(x, y)` with `mu(G)` the sum of the degrees (Barlow 2004, equation (1.4), read at source; the identity is credited there to Tetali, unread). On `P_n` that is `2 (f^n - 1) * 2 (b^n - 1)`: a walker crosses from corner to corner and back in `4 (f^n - 1)(b^n - 1)` steps, which is `(b f)^n` up to the constant 4, that is `side^(log(b f)/log b)`. Solved exactly from the hitting-time equations by rational Gaussian elimination, the commute times are `32, 768, 12896` for code 79 at levels 1 to 3 and `8, 96, 728, 4800` for code 7 at levels 1 to 4, each equal to `2 (nodes - 1) R(A, B)`. (**Proved**, Barlow 2004 (1.4) and Lemma 2; the hitting times **Verified**, `lab/py/tree-walks` verb `resist`.) So the crossing time of the side has exponent exactly `log 15 / log 3` and `log 6 / log 2`. What the next two sections add is the exponent of the mean square displacement and of the heat kernel, on the infinite tree.

## Volume growth in the graph metric

The theorem wants volume in the metric the walker feels, the graph distance, on an infinite graph. Two things need care: the infinite graph, and the metric.

**The infinite graph.** `P_n` sits inside `P_{n+1}` as the copy at the tile cell `(0, 0)`, so `G = P_0 U P_1 U ...` is an infinite tree in the quadrant. Cut at any scale `k`, the cells of `G` fall into blocks of side `b^k`, each a copy of `P_k` and an induced subtree, indexed by the cells of `G` itself: the blocks form the same tree one scale up, two blocks being joined by one edge exactly when their indices are adjacent in `G`. (Lemma 1 at every scale.)

**The chemical-distance exponent is 1 for corners and wrong for neighbours.** The graph distance between the corners `A` and `B` of the level-`n` block is `2 (b^n - 1)`: comparable to the side, exponent 1. But the graph metric is not comparable to the Euclidean one. Two cells at Euclidean distance `sqrt 2`, the facing cells of the copies at `(0, 1)` and `(1, 0)` (the `B` of one and the `A` of the other), sit at graph distance `4 b^n - 2`: the walker must go round through `O`. So the distortion between the two metrics is unbounded, and the volume exponent must be read in the graph metric. From the corner `O` the two metrics do agree: for a cell `y` of `G`, `max(y_0, y_1) <= d(O, y) <= b max(y_0, y_1)`, because `y` outside `P_{n-1}` lies in `P_n`, within `b^n - 1` of `O` by Lemma 2, and has a coordinate at least `b^(n-1)`. The generator's largest ratio `d(O, y) / max(y_0, y_1)` is `1.9959` for code 79 and `1.9961` for code 7, `2 - 1/b^L` at the cell `(b^L - 1, b^L)`. (**Proved**; the distance `4 b^n - 2` at `n = 1` to `5` and `1` to `8`, and the ratio, **Verified**, verbs `resist` and `volume`.)

**Proposition 3.** With `alpha = log f / log b`, every cell `x` of `G` and every `r >= 1` satisfy `(r / 2b)^alpha <= |B(x, r)| <= 7 (2 b r)^alpha`, where `B(x, r)` is the set of cells within graph distance `r`. (**Proved**.)

*Proof.* Lower: let `k` be the largest integer with `2 (b^k - 1) <= r`. The block of `x` at scale `k` has diameter `2 (b^k - 1)`, so it lies in `B(x, r)`, and `|B(x, r)| >= f^k`; and `2 (b^(k+1) - 1) > r` gives `f^k >= (r / 2b)^alpha`. Upper: an edge leaving a block of side `s = b^k` leaves from a corner cell of it. `O` has up to two outward neighbours (the `B` of the block before it along the first coordinate, the `A` of the block before it along the second), `A` and `B` one each; the three corners of a block are at mutual distance at least `s - 1`. A path of length `r < s - 1` therefore cannot cross a block from one corner to another. It can cross only through `O`, from one of its two outward neighbours to the other, and the block it then enters, at its `A` or its `B`, has no further outward neighbour on that corner, so the path stops there. The blocks met by the ball are the block of `x`, its at most four neighbours, and at most one more behind each of the at most two neighbours entered at their `O`: at most 7. Hence `|B(x, r)| <= 7 f^k` for `r <= b^k - 2`. Take `k` smallest with `b^k >= r + 2`; then `b^k <= b (r + 1) <= 2 b r`, so `f^k <= (2 b r)^alpha`.

The generator checks every inequality of the proof on the graph of level `L + 1` with the centres in its corner copy of level `L` (`L = 5` for code 79, `L = 8` for code 7), where balls of radius up to `b^L - 2` are exact in `G`: `V(x, 2 (b^k - 1)) >= f^k` at all 3125 and 6561 centres, `V(x, b^k - 2) <= 7 f^k`, and a ball of radius `b^k - 2` meets at most 3 blocks of side `b^k` (the proof allows 7). The ratio `V(x, r) / r^alpha` over all centres lies in `[0.46, 4.00]` for code 79 and `[0.41, 4.00]` for code 7 at `r >= 1`, and in `[0.46, 2.10]` and `[0.41, 2.33]` from `r = 3` and `r = 2`. (**Verified**, `lab/py/tree-walks` verb `volume`.)

## The theorem

**Theorem.** For `b >= 2` let `G` be the infinite tree of `L_b` and `alpha = log(2b - 1)/log b`. The simple random walk on `G` has a heat kernel with respect to the degree measure satisfying, for all cells `x, y` and times `n`,

`h_n(x, y) <= C n^(-alpha/(alpha+1)) exp(-(d(x, y)^(alpha+1) / (C n))^(1/alpha))`

and, for `n >= d(x, y)`, `h_n(x, y) + h_(n+1)(x, y) >= c n^(-alpha/(alpha+1)) exp(-(d(x, y)^(alpha+1) / (c n))^(1/alpha))`. Its walk dimension, the `beta` of these bounds, is `beta = alpha + 1 = log(b (2b - 1)) / log b`: `log 15 / log 3 = 2.464974` for code 79 and `log 6 / log 2 = 2.584963` for the base-2 design. (**Proved**, Barlow, Coulhon and Kumagai 2005, Proposition 5.2, hypotheses below.)

The source is Barlow, Coulhon and Kumagai, Theorem 1.3 and Proposition 5.2, read in the authors' preprint ([refs](../refs.md)). Theorem 1.3: a weighted graph with `p(x, y) >= p0 > 0` on edges and volume growth `VG(beta-)` for some `beta >= 2` satisfies the two-sided sub-Gaussian bound with exponent `beta` iff `R(x, y)` is comparable to `d(x, y)^beta / V(x, d(x, y))`. Proposition 5.2 is its corollary for trees. Each hypothesis, on `G`:

| hypothesis | on `G` |
|---|---|
| infinite, locally finite, connected | a nested union of connected trees, degree at most 4 |
| a tree, conductances bounded above and below | Lemma 1; all conductances 1 |
| `p(x, y) >= p0` | `p = 1/deg >= 1/4` |
| `V(x, r)` comparable to `r^alpha`, `alpha >= 1`, uniformly in the centre `x` and in `r >= 1` | Proposition 3, with `V(x, r)` the sum of degrees in the ball, within a factor 4 of the count; `alpha >= 1` because `G` is infinite and connected |
| `VG(beta-)` with `beta = alpha + 1`: `V(x, r) <= C (r/s)^gamma V(x, s)` for `r >= s > 0`, some `gamma < beta` | `gamma = alpha`, from the two-sided bound when `s >= 1` and from `V(x, s) = deg(x) >= 1` when `s < 1`, a ball of radius below 1 being a point |
| `R(x, y)` comparable to `d^beta / V(x, d)` | `R = d` (Lemma 5.1) and `V(x, d)` comparable to `d^alpha`, so `d^beta / V` is comparable to `d` |

The printed range in Proposition 5.2 is `r >= 0`; it is read as `r >= 1`, and Theorem 1.3 needs only the `VG` row, checked for all `r >= s > 0`. The uniformity in the centre is the hypothesis a corner measurement cannot supply, and Proposition 3 supplies it at every cell of `G`, not only at `O`. No hypothesis fails.

**Corollary (displacement).** Started at any cell `x` of `G`, `c n^(2/beta) <= E_x d(x, X_n)^2 <= C n^(2/beta)`. Started at `O`, the same holds for the Euclidean displacement, since `max(y_0, y_1) <= d(O, y) <= b max(y_0, y_1)`. (**Proved** from the Theorem.) The lower bound: by the upper heat kernel bound and Proposition 3, `P_x(d(x, X_n) <= eps n^(1/beta)) <= C eps^alpha`, at most 1/2 for small `eps`, so `E_x d^2 >= eps^2 n^(2/beta) / 2`. The upper bound integrates the same bound against the volume. This is the mean square displacement `MSD(t) ~ t^(2/d_w)` that [the census](walks.md) fits, with `d_w = beta`. Barlow 2004 gives the matching ceiling: a graph with volume exponent `alpha` and escape time exponent `beta` has `2 <= beta <= 1 + alpha` (Theorem 1), and `G` has `beta = 1 + alpha`.

The mean exit time is the other face of `beta`. From the corner, `E_O tau_(O, b^k)`, the mean time to leave the ball of radius `b^k`, which is exactly the block `P_k`, is an integer at every `k` computed, and equal to `b^k (f^k + 1) / 2`: `9, 117, 1701, 25353, 379809` for code 79 at `k = 1..5` and `4, 20, 112, 656, 3904, 23360, 140032, 839936` for code 7 at `k = 1..8`, by exact rational elimination on the tree. The ratio of successive terms is `13.00, 14.54, 14.90, 14.98` and `5.00, 5.60, 5.86, 5.95, 5.98, 5.99, 6.00`, tending to `b f = 15` and `6`. (**Verified**, `lab/py/tree-walks` verbs `exact` and `scale`. The closed form for every `k` is **Conjecture**: 13 exact values fit two parameters, no proof is given here.)

## What the readings say

The census table reads these two rows as finite-level estimates: eigenvalue ratios of the Laplacian on `P_n` and fits of walker displacement. The theorem is about the infinite tree, so the readings are re-run against its exact value. The generator rebuilds the spectral reading with its own eigensolver (dense, combinatorial Laplacian, `d_w = log(lambda_k(n) / lambda_k(n+1)) / log b`) and matches the lab's printed modes digit for digit. The lab's walker census was rerun as well.

| design | exact `d_w` | spectral, modes 1 to 4, levels 3 to 4 (5 to 6) | spectral, levels 4 to 5 (6 to 7) | walkers (lab rerun) |
|---|---|---|---|---|
| code 79 | 2.464974 | 2.4712 2.4707 2.4794 2.4694 | 2.4662 2.4663 2.4680 2.4662 | 2.4909, level 5 |
| code 7, base 2 | 2.584963 | 2.5894 2.5863 2.5979 2.5979 | 2.5865 2.5856 2.5894 2.5894 | 2.7462, level 8 |

(**Verified**, `lab/py/tree-walks` verb `scale` and `lab/rs/walk-dimension`.) Every spectral value sits above the exact one and falls toward it from the first pair to the second: code 79 by `0.0062` then `0.0012` on mode 1, code 7 by `0.0044` then `0.0015`. The walker value is `1.1%` above for code 79 and `6.2%` above for code 7, whose drift is the `0.29` the lab prints; the census note already tells the reader to quote the spectral value there. The exact value reads the walker figure as a finite-window overestimate on these two rows; the proof does not use any of these readings.

## The siblings

Nothing in the argument is special to `L_b` except two facts: the tile is a tree, and the faces of the tile meet in single cells. Both are properties of the tile.

**Proposition 4 (port trees, dimension 2).** Let a tile `T` of side `b` have `|T| >= 2` and be connected, and write `R_0 = {c : (b - 1, c) in T}`, `L_0 = {c : (0, c) in T}`, and `R_1, L_1` for the second coordinate. The level-`n` graph is a tree for every `n` iff the face-adjacency graph of `T` is a tree, `|R_0 meet L_0| = 1` when `T` has an adjacent pair along the first coordinate, and `|R_1 meet L_1| = 1` when it has one along the second. (**Proved**, below; **Verified** by the census.)

*Proof.* Two copies at adjacent tile cells are joined by `|R_0 meet L_0|^(n-1)` edges at level `n` (the face of `P_(n-1)` at the last first coordinate holds the cells whose second-coordinate digits all lie in `R_0`, the face at the first holds those with digits in `L_0`). One edge per tile edge on a tree gives a tree by induction. An overlap of 2 or more makes a cycle at level 2; an overlap of 0 disconnects level 2; a tile graph with a cycle or a gap fails at level 1 or 2.

The census agrees with the proposition on every code, by building the render and testing it for cycles and connectivity: of the 15 codes of base 2, 8 are trees with at least two cells at levels 1 to 4; of the 511 of base 3, 15 at levels 1 to 3; of the 65535 of base 4, 508 at levels 1 and 2. (**Verified**, verb `siblings`.) The base-2 trees are the four dominoes and the four L's. The 15 base-3 trees are the six straight lines (three rows, three columns), which draw paths, and nine designs of five cells: the four L's, the four T's, and the cross. Those nine are one class of the 72-element group of [the census](walks.md), the class of code 79: they are the nine unions of one full row and one full column of the 3 by 3 residue grid, codes 79, 121, 151, 186, 295, 316, 457, 466, 484, and the group (every permutation of the residues on each axis, and the swap of the axes) sends one such union to another. The six lines are the class of code 7. (**Proved**, by listing; the count **Verified**.) The class of code 79 therefore holds three drawn shapes up to the symmetries of the square, 79, 121 (the T) and 186 (the cross, the tile of the Vicsek graph that Barlow, Coulhon and Kumagai name as a tree of polynomial growth), and all nine are trees with `f = 5` and the same exact `d_w = log 15 / log 3`. The wrap caveat of Finding two, which splits 255 from 495, costs nothing in this class: every drawn representative has the exponent that the rule does.

Volume growth in a tree design needs one more fact, uniformity of the corner distances. A port is a cell where a copy meets its neighbour, the cell with all digits equal to the unique overlap row; a block has up to four, east, west, north and south. Distances between ports obey `d(n+1) = M d(n) + c` for an integer matrix `M` and constant vector `c`, read off the path in the tile tree. Call the design uniform when every nonzero port distance grows at one rate `lambda > 1`: all pairs reach the top Perron root of `M`, and no chain of two blocks of `M` carries that root, which would put a factor `n` in front. In a uniform design each nonzero port distance lies between `c lambda^n` and `C lambda^n` (Perron-Frobenius on the blocks of `M`), and the distance from any cell to a port of its block obeys `e_(n+1) <= e_n + C lambda^n`, so the diameter is also comparable to `lambda^n`. The block-counting step of Proposition 3 survives: two ports of a block share a cell only for one pair at most (`P-` and `Q-` iff both overlap rows are first, `P-` and `Q+` iff the first is last and the second first, and so on, four cases that exclude each other), opposite ports never share, so a path of length below the least port distance that crosses a block through a shared cell enters the next block at the port opposite the one it left by, and that port shares its cell with no other. Blocks met by a ball: the block of the centre, at most 4 neighbours, at most one more behind each: 9 in place of 7. Then `V(x, r)` is comparable to `r^alpha` with `alpha = log(fill) / log(lambda)`, the Theorem applies in the graph metric, and `d_w = log(fill * lambda) / log(lambda)`. (The block count and the growth rate are **Proved** for each uniform design as for `L_b`; the matrix, its blocks and its roots are computed by the generator, whose recursion equals breadth-first port distances on the render at levels 1 to 3 for every class: **Verified**.)

| base | tree codes | classes under the square's 8 symmetries | uniform | not uniform |
|---|---|---|---|---|
| 2 | 8 | 2 | 2 | 0 |
| 3 | 15 | 5 | 5 | 0 |
| 4 | 508 | 69 | 33 | 36 |

Base 2: the domino, `d_w = 2`, and `L_2`, `d_w = log 6 / log 2`. Base 3: the two kinds of straight line, `d_w = 2`, and the class of code 79 (three shapes), `d_w = log 15 / log 3`. Base 4, uniform, by `(fill, lambda)` and the value `log(fill * lambda)/log(lambda)` in the graph metric: the lines (4, 4, 2.0000) with 2 classes; (7, 4, 2.4037) with 3; (8, 4, 2.5000) with 6; (9, 4, 2.5850) with 12; (10, 4, 2.6610) with 3; (11, 4, 2.7297) with 1, the class of code 24414; and (9, 5, 2.3652) with 6, whose corner distances grow like `5^n` against the side `4^n`: a tree whose chemical distance exponent is `log 5 / log 4 = 1.161`, not 1, and whose displacement exponent in the graph metric, 2.3652, is a different number from the one `log(fill * lambda) / log(base) = 2.7459` that a resistance-scale reading by the side would print. `L_b` for every `b >= 2` is in the list and has `lambda = b`. (**Verified**, verb `siblings`, which prints every class.) Thirty-six base-4 classes are trees that are not uniform: some port distances grow slower than others or carry the extra factor, the volume is not a pure power, the Theorem does not apply, and the row stays open.

## What is not closed

- The general statement `d_w >= 2` for the census designs. For the trees it is `d_w = alpha + 1 >= 2`. Barlow 2004 Theorem 1 gives `beta >= 2` for graphs with both a volume exponent and an escape time exponent, and Nash-Williams type bounds give less; the other seven rows are not trees (each tile holds a 2 by 2 block or another cycle), have neither exponent proved, and keep the measurement-only status of "Music never beats mass". (**Conjecture**.)
- The finite-level readings of the other rows, and the hypothesis that any other census row has uniform volume growth in its graph metric, which would need the quantities this chapter proves for trees.
- Dimension 3. The sponge row is not a tree. The tripod, a corner cell with three arms, has the same tile structure, but the block-counting step is proved here only for two coordinates, because in three the ports can coincide in more than one pair.
- The 36 base-4 trees that are not uniform, and base 5 and above, where the tile test is cheap and the codes are many (`2^25`).
- The closed form `E_O tau_(O, b^k) = b^k (f^k + 1) / 2` for every `k`, and a rate for the decay of the spectral gaps in the table above.

> A tree has no hole to hide in: the walker pays for every step in distance and in volume, and the exponent is the sum.
