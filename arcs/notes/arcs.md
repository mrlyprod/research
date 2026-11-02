---
title: Loops in arcs
lead: A design drawn in quarter-circle arcs, its deleted cells turned once, and the closed loops the arcs make level by level.
figure: research-arcs
slug: arcs
---

Fix a dimension-two design: a mask of filled cells on a `base x base` grid, substituted into itself. Level `level` is the full `side x side` grid, `side = base^level`, each cell filled or deleted. Every cell carries two quarter-circle arcs, each joining the midpoints of two adjacent edges, the two orientations of a [Truchet tile](https://en.wikipedia.org/wiki/Truchet_tiles). A filled cell takes the arcs around its lower-left and upper-right corners; a deleted cell stays and takes the other pair, around its lower-right and upper-left corners. Arcs meeting at a shared edge midpoint join into curves: `L(level)` closed loops and some open strands ending on the boundary. The [arcs demo](/arcs/) draws any design at bases 2 to 5 this way, its loops in one ink and its strands in another, counted and set against the laws below.

## The object

- Coordinates: `x` is the column and `y` the row, row 0 at the bottom and `y` growing upward, and cell `(x, y)` is the unit square with lower-left corner `(x, y)`. Every cell on this page, of a mask or of a level, is named this way.
- The tree names a design `bang dim 2, base b, code c`: cell `(x, y)` of the `b x b` mask is filled when bit `b y + x` of `c` is set, and level `n` is the `n`-fold Kronecker power of the mask, built by `mrly_math::bang::factory::create`.
- A parity design `bang dim 2, code c` drawn at side number `m` fills cell `(x, y)` when bit `2 (y mod 2) + (x mod 2)` of `c` is set, so it is the base-`m` design whose code collects those cells. The carpet is `bang dim 2, code 7` at side number 3, the same picture as `bang dim 2, base 3, code 495`; at base 2 the two names agree.
- Level 0 is one filled cell. Write `N = side`, and `k` for the number of filled cells of the mask.
- Every count on this page is printed by `lab/rs/arc-loops`, which counts loops three ways: union-find over the edge midpoints, the mirror graph below, and the block recursion below. The first two share the cells built by `mrly-math` and one union-find structure; the block recursion shares neither.

## Strands

**Every curve is a path or a cycle, and there are exactly `2 side` strands. Proved.** An interior edge midpoint lies on two cells, and each cell has exactly one arc ending at each of its four edge midpoints, so the midpoint has degree 2; a boundary midpoint has degree 1. A graph of maximum degree 2 is a disjoint union of paths and cycles, and its path ends are its `4 side` vertices of degree 1. (**Verified**, `arc-loops check` on every code at base 2 to level 7 and at base 3 to level 4, and `arc-loops carpet` to level 7.)

## Loops are cycles of the mirror graph

The mirror graph `G` of a level has the `(side + 1)^2` lattice points as vertices and one edge per cell, the diagonal its two arcs do not cross: a filled cell joins `(x + 1, y)` to `(x, y + 1)`, a deleted cell joins `(x, y)` to `(x + 1, y + 1)`. Read the edge as a two-sided mirror and the curves are the light paths between mirrors.

**`L = c(G) - 2 side - 1`, where `c(G)` is the number of connected components of `G`, isolated lattice points included; so `L` is the cycle rank `E - V + c(G)` of `G`, and `L` is the number of components of `G` that hold no boundary lattice point. Proved.** The diagonal of a cell splits it into two right triangles and each triangle holds exactly one arc, the one around its right-angle corner. Two triangles that share a cell edge meet at its midpoint, where their arcs join, and no two triangles meet across a diagonal. So the triangles glued across cell edges are the faces of the plane graph `G` inside the square, and each face carries exactly one curve. A face that holds a boundary cell edge opens onto the unbounded face of `G` and its curve is a strand; every other face is a bounded face of `G` and its curve is a loop. Hence `L` is the number of bounded faces, which Euler's formula makes `E - V + c(G)` with `E = side^2` and `V = (side + 1)^2`. For the last clause, cut the square along all curves: `S` non-crossing chords leave `S + 1` regions on the boundary and each loop adds one region off it, while each region is the union of the corner pieces around the lattice points of one component of `G` and touches the boundary exactly when that component holds a boundary point. (**Verified**, `arc-loops check`: union-find over midpoints and the cycle rank agree on all 2688 levels checked.)

**The half turn, the transpose and the anti-transpose of the mask fix `L` at every level; the quarter turn and the two axis reflections need not. Proved.** The first three map the corner pair lower-left and upper-right to itself and commute with the Kronecker power, so they carry the picture of level `n` onto the picture of the image design at level `n`. The other four swap the two corner pairs, so they carry the picture to the image design drawn with the filled and the deleted arcs exchanged. Base 2 code 9 has `2^n - 1` loops and its quarter turn, code 6, has none; both counts are proved in [Base 2, complete](#base-2-complete). (**Verified**, `arc-loops check`: all three images of every code at base 2 to level 7 and at base 3 to level 4, 1584 images, zero mismatches.) So the census reduces each base to classes of this group of order 4.

## The block recursion

- Level `n + 1` is a `base x base` array of blocks of side `N = base^n`: the block under a filled mask cell is level `n`, the block under a deleted cell is the void block of `N x N` deleted cells.
- Each block side carries `N` edge midpoints, its ports, indexed from the lower-left corner along the side: `B_t`, `R_t`, `T_t`, `L_t` on the bottom, right, top and left, `0 <= t < N`. The strands of a block are a non-crossing perfect matching of its `4N` ports, written `B_t - L_t` and so on.

**The void block has strands `B_t - R_(N-1-t)` and `L_t - T_(N-1-t)` and no loop. Proved.** In a deleted cell the arc from the bottom edge runs to the right edge and the arc from the left edge to the top edge, so a strand entering at `B_t` climbs a staircase one cell right and one cell up at a time and leaves the right side at height `N - 1 - t`; the mirror graph is all main diagonals, each with both ends on the boundary. (**Verified**, `arc-loops check`: the glued void equals the formula on 27 sides at bases 2 to 5.)

**`L(n + 1) = k L(n) + J(n)`, where `J(n)` is the number of cycles in the glued matchings. Proved.** A loop of level `n + 1` either stays inside one block, where it is a loop of that block, or crosses a block side, where it alternates block strands and shared ports and so is a cycle of the graph that glues the block matchings along their shared sides; a void block holds no loop. The glued paths between outer ports are the matching of level `n + 1`, so the recursion carries the matching along and never needs the cells.

- `arc-loops` runs the recursion with the void block by formula and the matching of level `n` as one array of `4N` ports, so memory is linear in `side` and the census reaches `side = 3^17`. (**Verified**, `arc-loops check`: the recursion agrees with union-find over midpoints and with the cycle rank on every code at base 2 to level 7 and at base 3 to level 4, 2688 levels, zero mismatches.)
- Loops are easier than components here because the gluing sees only the matching of the block boundary, never the inside. The matching is a `4N`-point object, so finite state is not automatic; what the proofs below find is that for the carpet and at base 2 the matching is a fixed finite list of families, each family indexed by a set whose size is C-finite, and every family glues to families with the same index. Then `J(n)` is a fixed linear combination of the family sizes and `L` is C-finite.

## The carpet law

**The carpet has `L(n) = (8^n - 1)/7 - 3^n + n + 1` loops at level `n`, and its gluing adds `J(n) = 5 * 3^n - 7n - 5`. Proved.** The values run 0, 0, 3, 50, 509, 4444, 36727. (**Verified**, `arc-loops carpet`: union-find over midpoints and the cycle rank to level 7, the block recursion to level 15, zero mismatches against the law.)

The proof reads the strand matching of every level off one lemma. Write `a_j = (3^j - 1)/2`, so the base-3 digits of `a_j` are `j` ones; let `A_n = {a_0, ..., a_n}` and `A'_n = {a_0, ..., a_(n-1)}`, and note that `a_n = (N - 1)/2` is the middle port. Let `P_0` be empty and let `P_(n+1)` hold `P_n`, `N + P_n` and `2N + P_n` together with the pairs `(jN - 1 - a, jN + a)` for `j = 1, 2` and `a` in `A'_n`.

**Lemma. Proved.** The strands of carpet level `n` are exactly

- the lower-left family `B_a - L_a` and the upper-right family `R_(N-1-a) - T_(N-1-a)`, for `a` in `A_n`;
- the lower-right family `B_(N-1-a) - R_a` and the upper-left family `L_(N-1-a) - T_a`, for `a` in `A'_n`;
- the turns `X_t - X_u` on each of the four sides `X`, for every pair `(t, u)` in `P_n`.

So each side carries `2n + 1` corner ports and `|P_n| = (3^n - 2n - 1)/2` turns, from `|P_(n+1)| = 3 |P_n| + 2n`. (**Verified**, `arc-loops carpet`: the stated matching equals the glued matching at every level 0 to 10.)

Proof, by induction on `n`. Level 0 is one filled cell with strands `B_0 - L_0` and `R_0 - T_0`, which is the lemma with `A_0 = {0}` and `A'_0`, `P_0` empty. Assume the lemma at level `n` and glue the eight filled blocks around the void; name a block by its column and row, `00` the lower-left and `11` the void. A loop's length counts the block strands it uses, strands of the void block included.

- The turn pairs are symmetric under `t -> N - 1 - t`: the map permutes the three copies in `P_(n+1)` and exchanges its two new families. So on every side of a filled block the corner ports form the same set `A'_n`, `{a_n}`, `N - 1 - A'_n`, and the turns sit on the same set `P_n`.
- Turns meet turns. At each of the 8 sides shared by two filled blocks a turn `(t, u)` meets the turn `(t, u)` and closes a loop of 2 strands. The void joins the top of `10` to the left of `21` by `x -> N - 1 - x`, and the right of `01` to the bottom of `12` the same way, so by the symmetry each of those two corners closes `|P_n|` loops of 4 strands. That is `10 |P_n|` loops.
- Corner strands meet corner strands and keep their index `a`. Across a vertical seam, left block before right: lower-right `a` meets lower-left `a`, upper-right `a` meets upper-left `a`, and upper-right `a_n` meets lower-left `a_n`. Across a horizontal seam, lower block before upper: upper-left `a` meets lower-left `a`, upper-right `a` meets lower-right `a`, and upper-right `a_n` meets lower-left `a_n`. Through the void, `10` upper-left `a` meets `21` upper-left `a`, `10` upper-right meets `21` lower-left, `01` lower-right meets `12` lower-right, and `01` upper-right meets `12` lower-left, all at equal `a`.
- So the corner strands split into one network for each `a` in `A'_n`, all `n` of them the same, and one network for `a_n`.

Each network for `a` in `A'_n` uses 32 strands of the filled blocks and closes 3 loops, of 8, 4 and 4 strands: `00` upper-right, `01` lower-right, `12` lower-right, `22` lower-left, `21` upper-left, `10` upper-left; then `10` upper-right, `21` lower-left, `20` upper-left; then `01` upper-right, `02` lower-right, `12` lower-left. Its other strands make 12 paths between outer ports: `00` lower-left, `20` lower-right, `22` upper-right and `02` upper-left are the four corner strands of index `a` at level `n + 1`, and the eight pairs `00` lower-right with `10` lower-left, `10` lower-right with `20` lower-left, `00` upper-left with `01` lower-left, `01` upper-left with `02` lower-left, `20` upper-right with `21` lower-right, `21` upper-right with `22` lower-right, `02` upper-right with `12` upper-left and `12` upper-right with `22` upper-left are the new turns `(jN - 1 - a, jN + a)`, two on each side.

The network for `a_n` uses 16 strands of the filled blocks and closes no loop. With `m = a_n`, so that `N + m = a_(n+1)` and `2N + m = 3N - 1 - m`, its six outer paths are `B_m - L_m` from `00`; `B_(N+m) - L_(N+m)` through `10`, `00`, `01`; `B_(2N+m) - R_m` through `20`, `10`, the void, `21`, `20`; `L_(2N+m) - T_m` through `02`, `01`, the void, `12`, `02`; `R_(N+m) - T_(N+m)` through `21`, `22`, `12`; and `R_(2N+m) - T_(2N+m)` from `22`. These are the corner strands of index `a_n` in all four families and of index `a_(n+1)` in the lower-left and upper-right families, which completes `A_(n+1)` and `A'_(n+1)` and the lemma at level `n + 1`.

Counting the loops closed at this gluing, `J(n) = 10 |P_n| + 3n = 5 * 3^n - 7n - 5`, and `L(n + 1) = 8 L(n) + J(n)` with `L(0) = 0` solves to the law. (**Verified**, `arc-loops carpet`: `J(n) = 10 |P_n| + 3n` at all 15 gluings to level 15.)

- Read on the picture: every loop of the carpet closes at some gluing, and there it is a 2-strand loop across a seam, a 4-strand loop around one of two corners of the hole, two of its strands in the void, or one of the three loops per index `a` around the hole. No loop ever uses more than 8 block strands. (**Verified**, `arc-loops lengths`: 8 is the longest at gluings 0 to 11.)
- The four roots `8, 3, 1, 1` are the number of filled cells, the base, and the double root from the `n` corner indices.

## Base 2, complete

**At base 2 the loop count is one of four laws. Proved.**

| codes | `L(n)`, `n >= 1` | `J(n)`, `n >= 1` | OEIS |
|---|---|---|---|
| 7, 14 | `3^(n-1) - 2^n + 1` | `2^n - 2` | [A028243](https://oeis.org/A028243) |
| 11, 13 | `3^(n-1) - 2^(n-1)` | `2^(n-1)` | [A001047](https://oeis.org/A001047) |
| 9 | `2^n - 1` | 1 | [A000225](https://oeis.org/A000225) |
| the other 11 | 0 | 0 | |

(**Verified**, `arc-loops census` to level 24 and `arc-loops two`: each lemma below equals the glued matching at every level 1 to 14 and each `J` holds at all 13 gluings.) Codes 14 and 13 are the half turn and the transpose of 7 and 11, so the symmetry above carries their laws. Each law follows from `L(n + 1) = k L(n) + J(n)` and `L(1) = 0`, or `L(0) = 0` for code 9, once its lemma is proved by induction from level 1, read off the `2 x 2` picture, or from level 0 for code 9.

- **Code 7**, the void in the upper-right block. Lemma: `B_0 - L_0` and `B_1 - L_1`; the turns `(2i, 2i + 1)` for `i >= 1` on the bottom and on the left; the turns `(2i, 2i + 1)` for all `i` on the right and on the top. Gluing: across the seam of `00` and `10`, the right turns of `00` meet the left turns of `10` for `i >= 1`, `N/2 - 1` loops, and its turn `(0, 1)` meets `10` at `L_0`, `L_1`, which run to the outer `B_N`, `B_(N+1)` and make the new turn there; the seam of `00` and `01` is the transpose; the void routes the top of `10` and the right of `01` to the outer right and top sides, where its reversal keeps each turn a turn. So `J(n) = N - 2`.
- **Code 11**, the void in the upper-left block. Lemma: `B_0 - L_0`, `B_(N-1) - R_0`, `R_(N-1) - T_(N-1)`, `L_y - T_(N-1-y)` for `1 <= y < N`, and the turns `(2i - 1, 2i)` for `1 <= i < N/2` on the bottom and on the right. Gluing: a right turn `(2i - 1, 2i)` of `00` crosses `10` from left to top as `T_(N-2i)`, `T_(N-1-2i)` and meets the bottom turn of `11` at that pair, `N/2 - 1` loops; the strand `R_(N-1) - T_(N-1)` of `00` runs through `10`, `11` and the void back to itself, one loop. So `J(n) = N/2`.
- **Code 9**, the diagonal. Lemma, from level 0: the void matching with its two outer corner strands exchanged, `B_0 - L_0` and `R_(N-1) - T_(N-1)`, then `B_x - R_(N-1-x)` and `L_x - T_(N-1-x)` for `1 <= x < N`. Gluing: the strand `R_(N-1) - T_(N-1)` of `00` runs through the void `10`, `11` and the void `01` back to itself, one loop, and every other strand runs straight through. So `J(n) = 1`.
- **The other 11 codes never loop**, by the last clause of the cycle-rank theorem: every component of the mirror graph touches the boundary. Codes 0 and 15 are uniform, so every component is a full diagonal. Codes 1, 2, 4, 8 fill one corner cell of the square at every level, whose diagonal touches the boundary, and every main diagonal it cuts keeps one end on the boundary. Codes 3, 5, 10, 12 fill one full outer row or column: there each diagonal touches the boundary, and every main diagonal of the void reaches the boundary or ends on that row or column. Code 6 fills the anti-diagonal of the square, whose diagonals form one path from `(N, 0)` to `(0, N)`, and each main diagonal crosses it at most once, so each piece keeps an end on the boundary.
- Code 6 and code 9 are quarter turns of each other: the pair that shows `L` is not invariant under the full symmetry group of the square.

## The census at base 3

**Base 3 has 512 codes in 168 classes under the symmetry; 48 classes, 149 codes, never loop to level 14, and the other 120 classes give 74 distinct nonzero sequences. Verified**, `arc-loops census`, every class representative to level 14 by the block recursion.

**When `k > base`, `L(n) / k^n` converges to `C = sum_m J(m) / k^(m+1)`, positive unless `L` vanishes. Proved.** A loop closed at a gluing crosses a seam and comes back, so it uses at least two of the `2 base (base - 1) N` seam ports, and `J(n) <= base (base - 1) base^n`. Unrolling `L(n) = sum_(m < n) k^(n-1-m) J(m)` then gives a series in `(base/k)^m` that converges. The carpet has `C = 1/7`.

The fits below come from Berlekamp-Massey over a prime, lifted to the integers and checked exactly on every term, allowing a transient of up to 7 levels, and kept only with at least two terms past those that determine it. A fit is a statement about the computed range only. (**Verified**, `arc-loops census`, for every fit on its range; **Conjecture** for every `n`.)

- 61 of the 74 sequences fit a recurrence with integer roots only, of order at most 6. On all 68 sequences that fit by level 14 the largest root is `k`, and 55 of them carry the root 3. Six need a transient of at least two levels: the classes of code 30 for five levels, 173 and 189 for four, 71 for three, 177 and 379 for two. (**Conjecture**.)
- 8 sequences, 66 codes, fit with the factor `x^2 - 3x + 1`, whose roots are `phi^2` and `phi^-2` with `phi` the golden ratio. The class of code 13 is the cleanest: its gluings add the odd-indexed Fibonacci numbers, `J(n) = F_(2n-3)` for `n >= 1`, and `L` is [A104487](https://oeis.org/A104487) shifted by two. (**Verified**, `arc-loops gains`, at gluings 1 to 13; **Conjecture** at every gluing.) Every code of the 8 classes holds the lower-right corner cell and the left middle cell, up to the symmetry.
- Code 287 fits with the factor `x^2 + x + 1`, a part of period 3, with only two terms of margin at level 17. (**Conjecture**.)
- Codes 43, 171, 175 and 181 fit no recurrence of order at most 8 with two terms of margin through level 17. (**Verified**, `arc-loops census`, no such fit; whether they are C-finite is open.) Their gluings close loops of ever more block strands: at the gluing from level 11 to 12 code 43 closes a loop of 6032 block strands, void strands counted, against at most 8 for the carpet on the same count. (**Verified**, `arc-loops lengths`.)

The parity designs at side number 3 add clean gains, each verified at gluings 0 to 13 and open past them.

| parity code | base-3 code | `J(n)` | `L(n)` | OEIS |
|---|---|---|---|---|
| 14 | 186 | `2 * 3^(n-1)`, `J(0) = 1` | `2 * 5^(n-1) - 3^(n-1)`, `n >= 1` | [A081625](https://oeis.org/A081625) |
| 9 | 341 | `2 * 3^n` | `5^n - 3^n` | [A005058](https://oeis.org/A005058) |
| 6 | 170 | `2^(n+1)` | `4^n - 2^n` | [A020522](https://oeis.org/A020522) |
| 11, 13 | 471, 381 | `2^(n+1) - 2` | `(7^n - 6 * 2^n + 5)/15` | none |

(**Verified**, `arc-loops gains`, which tests each `J` at gluings 0 to 13 and each `L` at levels 0 to 14; **Conjecture** at every level. The closed forms follow from the gains by `L(n + 1) = k L(n) + J(n)`.)

## Bases 4 and 5

The first codes at bases 4 and 5 are the fifteen parity designs at side number 4 and 5 and every design that deletes one cell, 14 classes at base 4 to level 12 and 20 at base 5 to level 10. (**Verified**, `arc-loops bases`, every count; **Conjecture**, every fit, with the same rules as at base 3.)

- No loop to level 12 at side number 4 for parity codes 3 and 15, and none to level 10 at side number 5 for parity codes 5, 12 and 15; for code 15, the full square, none at any level.
- Parity codes 1 and 2 at side number 4 and code 8 at side number 5 keep 4 cells and fit `4^n - 3^n`, [A005061](https://oeis.org/A005061). Parity code 6 fits roots `8, 4` at side number 4 and `12, 4` at side number 5; parity code 9 fits `8, 3` and `13, 5`; parity code 1 at side number 5 fits `9, 5`, code 4 fits `6, 4`, code 13 fits `19, 3, 1` and code 14 fits `16, 5, 3, 1`.
- The carpet's own parity code 7 fits roots `12, 4, 1` at side number 4 and `21, 5, 1, 1` at side number 5, the second of the carpet's shape `k, base, 1, 1`.
- Deleting one cell: at base 4 every class fits roots `15, 4, 2, 1`, except the cell `(1, 3)`, `bang dim 2, base 4, code 57343`, which fits `15, 4, 3, 2, 1`. At base 5 the corner, the centre and the inner ring all fit `24, 5, 3, 1`, so the centre deletion at base 5 does not repeat the carpet's double root 1. The three cells `(x, 4)` with `0 < x < 4`, codes 25165823, 29360127 and 31457279, fit nothing to level 10.
- No fit at bases 4 and 5 has a root outside the integers.

## The OEIS

Every nonzero sequence is looked up in a local copy of [the OEIS stripped dump](https://oeis.org/stripped.gz) by its first seven terms from the first nonzero one. A hit means the window agrees; it is not a proof that the sequences agree. (**Verified**, `arc-loops census`, `arc-loops bases` and `arc-loops gains`, each row printing its hit or its absence.)

- Base 2: all 3 hit, [A028243](https://oeis.org/A028243), [A001047](https://oeis.org/A001047) and [A000225](https://oeis.org/A000225), and the laws above prove the agreement.
- Base 3: 20 of the 74 distinct nonzero sequences hit and 54 are absent. Among the absent are the carpet, 0, 3, 50, 509, 4444, whose law is proved above, the parity designs 11 and 13 at side number 3, 2, 20, 154, 1108, code 287 and the four codes with no fit.
- Bases 4 and 5: the 14 and 20 classes give 10 and 15 distinct nonzero sequences, of which 3 and 3 hit.

## Generators

- `lab/rs/arc-loops` prints every count on this page; its README lists each verb, its runtime and the lines it witnesses.
- `mrly_math::bang::factory::create` builds every level the union-find and cycle-rank counts read; the block recursion reads only the mask and computes the void block by its formula.
