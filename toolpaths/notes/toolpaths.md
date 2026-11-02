---
title: The toolpath
lead: Self-similar curves on the triangular lattice as one-stroke print paths, and the sharp turn a printer slows for: the edge-covering kind never turns 60 degrees, every point-avoiding replacement curve of order 3 to 13 (order 12 without mirrors) turns 120 degrees at some level, the junction law that makes that a finite check, the exact sharp-turn count of the Gosper curve, the bead axis, which on the Gosper curve fades as `7^(-k/2)`, slower than Hilbert's, while Peano's holds at `1/2`, and the radix census, where unit steps at every level are one identity at level 1 and a short periodic check and no plane-filling replacement curve with every flag `F` to a chord among the seven radix bases avoids points.
figure: research-toolpaths
slug: toolpaths
---

Write `w = e^(i pi/3)`. A replacement curve on the triangular lattice [`Z[w]`](https://en.wikipedia.org/wiki/Eisenstein_integer) is a generator, `n` unit steps in directions `g_0, ..., g_(n-1)` (multiples of 60 degrees) from `0` to a chord `c` with `N(c) = n`, and a flag per step saying which copy of the generator replaces that step at the next level: `F` as it is, `R` run backwards, and, when `c/conj(c)` is a unit, `M` mirrored across the chord and `MR` both; these four are every symmetry of a segment, and the mirror stays on the lattice only at the orders `n = r^2` and `3 r^2`, so `3, 4, 9, 12, 16` among those below. Level `k` is `n^k` unit steps from `0` to `c^k`, a curve of [similarity dimension](https://en.wikipedia.org/wiki/Fractal_dimension) 2. A turn is the change of heading between two steps: 0, 60 or 120 degrees, since 180 retraces a step. A curve is gentle when no level turns 120 degrees; a printer head slows at a sharp turn, so the goal is a gentle self-similar infill.

The generator is the lab study `lab/py/toolpath-turns` (`uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py <verb>`), verbs `census`, `plain`, `gosper`, `bead` and `rate`. [Arndt](https://arxiv.org/abs/1607.02433) searches only the edge-covering kind, as L-systems of turns, and gives the Gosper curve, which is not of that kind, as an L-system in his section 1.3; the flowsnake read as a radix tile is in the [radix dial](../../beneath/notes/beneath.md).

## Two kinds of curve

Two families share the name plane-filling. An edge-covering curve runs every edge of a region once and never crosses itself; it passes an interior point three times, once per pair of its six edges. A point-avoiding curve visits every point at most once; the Gosper curve is one.

**An edge-covering curve on the triangular lattice never turns 60 degrees at an interior point. Proved.** Around a point its six edges are paired by three passes, and a non-crossing pairing of six points on a circle never pairs two edges with exactly one edge between them, since that edge is then fenced in with no partner; a 60-degree turn is exactly such a pair, a straight pass pairs opposite edges and a 120-degree turn adjacent ones. So every turn of the kind is 0 or 120 degrees, which is the turn alphabet Arndt's triangular-grid search uses (his section 2.3), and a curve of that kind that turns at all turns sharply. The question with content is the point-avoiding kind.

## The junction law

Write a turn type as a triple `(t, f, f')`: the turn `t` between two steps and their flags. At the next level a step with flag `f` becomes a copy whose first and last steps sit at fixed offsets from its heading, so a turn of type `(t, f, f')` becomes one junction turn `t + a(f') - b(f)` with a new flag pair, and each copy adds the generator's own turns with flags carried along.

**The turn types of every level lie in the closure of the level-1 types under that map and the copy rule, a set of at most 96 triples; a curve is gentle at every level if and only if no triple in the closure turns 120 degrees. Proved.** Level `k+1` holds copies of level `k`, so a sharp turn or a revisited point at one level persists at every later level. Gentleness at every level is therefore a finite check (`closure` in the study), and self-avoidance is the only part that needs levels. Where a forward copy meets a reversed one the junction turn equals the parent turn; two forward copies shift it by `g_0 - g_(n-1)`.

## The census

**At orders 3, 4, 7, 9, 12 and 13, the norms from 3 to 13, every point-avoiding replacement curve turns 120 degrees at some level; mirrors are included at 3, 4 and 9, are off the lattice at 7 and 13, and are left out at 12. Verified** (`census`, `plain 12`). Each of these norms has one Eisenstein integer up to units and conjugation, so one chord per order covers every curve up to congruence. The search runs every gentle generator walk to the chord, every flag assignment that is gentle by the closure and point-avoiding at level 2, and expands the rest. A pair is a generator with its flags; a class is a pair with its reversal and, where the chord allows it, its mirror image.

| order | chord | flags | gentle walks | pairs gentle at every level, avoiding at level 2 | classes | first crossing |
|---|---|---|---|---|---|---|
| 3 | `1+w` | `F R M MR` | 0 | 0 | 0 | none |
| 4 | `2` | `F R M MR` | 0 | 0 | 0 | none |
| 7 | `2+w` | `F R` | 20 | 10 | 5 | level 3, all |
| 9 | `3` | `F R M MR` | 94 | 11696 | 2972 | level 3, all |
| 12 | `2+2w` | `F R` | 2646 | 3554 | 889 | level 3, all |
| 13 | `3+w` | `F R` | 6794 | 9680 | 4840 | level 3, all |

At order 12 with mirrors the census stops after walk 379 of 2646, with 477775 pairs, all crossing at level 3; at order 16 with `F, R` it stops after walk 97426 of 122964, with 93057 pairs, all crossing at level 3; order 16 with mirrors is not run. Every failure is a crossing: no pair at these orders is gentle by the closure and still avoiding at level 3. **The same holds at orders 12 with mirrors and 16 with or without them. Conjecture.**

## The Gosper count

The Gosper curve is the generator `0, 5, 3, 4, 0, 0, 1` to the chord `3-w` with flags `F R R F F F R`. **At level `k` it has `7^(k-1)` straight joints, `4*7^(k-1) - 1` turns of 60 degrees and `2*7^(k-1)` of 120 degrees. Verified** at levels 1 to 6 by `gosper`, with the curve point-avoiding through level 5. Its mirror image `0, 1, 3, 2, 0, 0, 5` to the chord `2+w`, with the same flags, is the flowsnake direction sequence [A229214](https://oeis.org/A229214); `gosper` checks its first 49 terms.

## The bead axis

A printed bead has an axis but no arrow, so its orientation lives in the second harmonic `h2 = sum e^(2 i theta)` over the steps; the fourth harmonic is `1` on every square-lattice step and the sixth on every triangular one, so neither compares anything. With `F, R` flags the level-`k` harmonic is `G^k` for the generator's own sum `G`, since reversal keeps directions.

**On the Gosper curve `abs(h2)^2 = 7^k` exactly, since `G = 3 + w^2` has norm 7, so its axis anisotropy `abs(h2)/7^k` is `7^(-k/2)`. Proved**, and Verified at levels 1 to 6 by `gosper`. **Hilbert's curve has `h2 = -1` at levels 1 to 6 and Peano's has `h2 = (9^k - 1)/2` at levels 1 to 4. Verified** (`bead`). Per step that is Hilbert near `1/4^k`, Gosper `7^(-k/2)` and Peano a constant `1/2`: Hilbert balances its axis fastest.

## The printer's line

**If a replacement curve of order `n` first turns 120 degrees at level `L`, then level `k >= L` has at least `n^(k-L)` such turns. Proved**: level `k` holds `n^(k-L)` disjoint copies of level `L`, each moved by one of the four segment symmetries, which keep every turn size. So sharp turns never thin out; the only lever is their density.

**Up to traversal direction, two order-7 replacement curves to the chord `2+w` avoid points through level 4: the Gosper curve's mirror image, with `2*7^(k-1)` turns of 120 degrees per `7^k` unit steps, `2/7` per step, and `0, 2, 2, 0, 0, 0, 4` with flags `R R F F R F F`, with `3*7^(k-1)`. Verified** (`rate 7`, levels 1 to 5). Every other chord of norm 7 is a rotation of `2+w` or of its mirror image `3-w`, so the Gosper curve has the fewest sharp turns at order 7. So at the orders the census exhausts, with `M, MR` searched only at 3, 4 and 9, no one-generator replacement curve is a point-avoiding path with no sharp turn; curves built from several motifs lie outside the census.

## The radix census

A radix design of [the radix dial](../../beneath/notes/beneath.md) is a base `b` in `Z[i]` or `Z[w]`, digits `d_0, ..., d_(k-1)` pairwise incongruent modulo `b`, a unit twist `u_a` per digit and place maps `phi_a(x) = (u_a x + d_a)/b`. Its level-`L` points are the scaled words `P(a_1 ... a_L) = sum_j (prod_(i<j) u_(a_i)) d_(a_j) b^(L-j)`, listed with the first digit slowest, and that list is a path through the level. The census runs at the seven bases of the dial: `2`, `1+i`, `2+i` in `Z[i]` and `2`, `1+w`, `3`, `2+w` in `Z[w]`, which are the dial's `2`, `2+omega`, `3` and `3+omega` since `omega = w^2 = w - 1`. A pair is a code over the canonical residues, read in canonical order, with a twist vector. A design has unit steps at a level when consecutive points differ by a unit; it is a curve when it has unit steps and no two words land on one point; it is one piece when the level's points are connected under unit adjacency; it is plane-filling when `k = N(b)` and no two words of one length share a place map. Each is asked at every level. The generator is the lab study `lab/rs/radix-census` (`cargo run --release --manifest-path research/toolpaths/lab/rs/radix-census/Cargo.toml`), verbs `census`, `curves`, `junction` and `named`.

### The step law

Write `z_0 = d_0/(b - u_0)` and `z_e = d_(k-1)/(b - u_(k-1))`, the fixed points of the first and last place maps. Two consecutive words at level `L` share a prefix `p` of length `L - l`, then read `a` against `a+1`, then `k-1` repeated against `0` repeated, and their step is `U(p) J_a(l)`, with `U(p)` the product of the prefix's twists and

```
J_a(l) = A_a b^(l-1) - u_(a+1) z_0 u_0^(l-1) + u_a z_e u_(k-1)^(l-1),    A_a = b (phi_(a+1)(z_0) - phi_a(z_e)).
```

**A design has unit steps at every level if and only if `phi_a(z_e) = phi_(a+1)(z_0)` for every `a` from `0` to `k-2` and `J_a(l)` is a unit for `l` from 1 to the order of the unit group, 4 on `Z[i]` and 6 on `Z[w]`. Proved.** A run of `m` first digits lands on `d_0 (b^m - u_0^m)/(b - u_0)` and a run of `m` last digits on the same with `d_(k-1)` and `u_(k-1)`, which gives `J_a(l)` term by term. If some `A_a != 0`, `abs(J_a(l)) >= abs(A_a) abs(b)^(l-1) - abs(z_0) - abs(z_e)` grows without bound and a level breaks; if every `A_a = 0`, `J_a(l)` is periodic in `l` with period dividing 4 or 6. So unit steps at every level are a level-1 identity, the end of copy `a` meeting the start of copy `a+1`, together with a check of at most six levels, which is not implied: `Z[w]` base `1+w`, digits `1, w`, twists `1, w` meets the identity and has a step of norm 3 at level 2. This is the junction law above, read for radix words.

**Unless its attractor is one point, such a design is conjugate, map for map, to the walk of its twists, and `b = u_0 + ... + u_(k-1)`. Proved.** With every `A_a = 0`, `d_(a+1) - d_a = u_a z_e - u_(a+1) z_0`. If `z_0 = z_e`, this gives `d_a = (b - u_a) z_0`, so every map fixes `z_0` and the attractor is that point, as for `Z[i]` base `2+i`, digits `i, 1+i, 1+2i`, twists `1, i, -1`, whose twists sum to `i` (`junction`). Otherwise `h(x) = (x - z_0)/(z_e - z_0)` carries the digits to the partial sums `0, u_0, u_0 + u_1, ...` and the twists to themselves, and `z_e -> 1` forces the sum, so the attractor is the replacement curve of this note whose generator is the unit walk `u_0, ..., u_(k-1)` from `0` to `b`, every flag `F`. The word order runs that curve's vertices, up to a translate at each level, exactly when `d_0 = 0` or the twists are all equal: the level path is seeded at `h(0) = -z_0/(z_e - z_0)`, the walk's start exactly when `d_0 = 0`, and a seed `s != 0` moves the word `x` by `U(x) s`, the same for every word of a level exactly when the twists are equal. All 21 canonical pairs with unit steps have `d_0 = 0` or no twist, and the census asserts the sum on each (`census`).

**At the seven bases a pair without unit steps at every level keeps them through level 3 at most, so unit steps at level 4 give unit steps at every level there, and level 3 does not. Verified**, by the step law over all pairs (`census`, `junction`). The deepest is `Z[i]` base `1+i`, digits `0, 1`, twists `1, -1`, whose levels 1 to 3 are unit paths and whose level 4 is not, since `A_0 = 1 - 1/(2+i) != 0`; the same depth recurs on `Z[w]` at `2` and `1+w`.

### Glue and area

**Whether two words of one level land on one point, or share a place map, or sit in touching copies, is decided at every level by a finite automaton. Proved.** Two words that first differ at digits `a != a'` are compared through `delta`, their difference over the turn of the first, and `r`, the turn of the second over the first: they start at `delta = (d_a - d_(a'))/u_a`, `r = u_(a')/u_a`, and appending digits `e, e'` sends `delta -> (b delta + d_e - r d_(e'))/u_e` and `r -> r u_(e')/u_e`. The words land on one point at level `L` when `delta = 0` after `L - 1` steps, share a place map when also `r = 1`, and sit in touching copies when `abs(delta) <= 1`. Since `abs(b delta + d_e - r d_(e')) >= abs(b) abs(delta) - 2M`, with `M` the largest digit modulus, a state with `abs(delta) >= 2M/(abs(b) - 1)` and `abs(delta) > 1` never returns, so finitely many states matter and a breadth-first search settles every level at once. The level set is connected at every level exactly when, at every level, the graph on the digits joining touching copies is connected, by induction on the level, since each copy is a turned and moved copy of the level below; the set of states that reach a touch in exactly `n` steps is eventually periodic in `n`, so one piece at every level is decided too.

**With `k = N(b)` a design has positive area if and only if no two words of one length share a place map. Proved.** If two do, the `N(b)^L` copies of level `L` cover the attractor `K` with two of them equal, so `abs(K) <= (1 - N(b)^(-L)) abs(K)` and `abs(K) = 0`. If none do, seed the words at `x_0 = 1/3` on `Z[i]` or `1/5` on `Z[w]`, where `(r - 1) x_0` is in the ring for no unit `r != 1`: distinct maps then give `N(b)^L` distinct points of the lattice `b^(-L) x_0 R`, each owning a cell of area `covol(R) abs(x_0)^2 N(b)^(-L)`, all within a vanishing distance of `K`, so `abs(K) >= covol(R) abs(x_0)^2 > 0`.

**At `2+i` and `2+w` no twist moves the level set of the full residue system. Proved.** There the canonical residues are `0` and the units, a set every unit fixes, so level `L` is `b S + C` for the level-`(L-1)` set `S` and the residues `C`, whatever the twists, and the untwisted fill law gives distinct points. So all `4^5 = 1024` and `6^7 = 279936` full twists are plane-filling with distinct points at every level (`census`).

**The depths these need are exact and not small: at the seven bases the first revisit of a design with unit steps comes by level 5, the first shared place map of a full residue twist by level 7 and the first split of a canonical pair by level 13, and each is reached. Verified** by the automaton over every pair and walk (`census`, `curves`), and checked against brute force over the words of each level wherever the level fits (`junction`). The walk `1, w` at base `1+w`, two steps at 60 degrees of dimension `log 4/log 3`, is a curve through level 4 and lands two words on one point at level 5; `Z[w]` base `2` with the full residues twisted by `w^2, w^4, w^5, w` shares no place map through level 6 and shares one at level 7; `Z[w]` base `3`, code `223`, twists `1, w^4, w^3, w^5, w^2, 1, w^4` is one piece through level 12 and splits at level 13, checked by brute force to level 7.

### The census

**How the group acts on a twist. Proved.** For `h(x) = v x` with `v` a unit, `h phi_(d,u) h^(-1) = phi_(v d, u)`: the twist rides with its digit and never changes. For `h(x) = v conj(x)` with `conj(b) = eps b`, `h phi_(d,u) h^(-1) = phi_(eps^(-1) v conj(d), eps^(-1) conj(u))`: the mirror conjugates every twist and turns it by `eps^(-1)`. The scaled points move by `v`, or by a unit times `conj`, so every property of the level set is kept whenever the image digits are again the canonical residues. The order of the digits is not kept, so unit steps and curves are counted pair by pair. The elements that fix the canonical residues form the stabiliser, of order `2, 2, 4, 2, 2, 2, 6` at the seven bases in the order of the table; it cuts the `40353552` pairs at base `3` to `20181403` orbits, just under a factor of 2, the orbit walk agreeing with Burnside at every base (`census`).

**The code action of the radix dial, twists carried along, keeps one piece. Refuted.** It moves a digit to the canonical representative of its image class, which is a different design; at `2+i` and `2+w` it is the stabiliser and keeps everything. Witness: `Z[w]` base `2`, digits `1, w`, no twist, one piece at every level; conjugation sends the class of `w` to the class of `-1+w`, and the code `1, -1+w` is apart at level 1, while the true image `1, 1-w` is one piece. At base `3` that action would cut to `3371677` orbits.

**The census of canonical pairs with `k >= 2`, every property at every level. Verified**, by exhaustive enumeration with the step law and the automaton (`census`); orbits are under the stabiliser.

| ring | base | pairs | orbits | unit steps | curves | one piece | plane-filling | plane-filling with distinct points |
|---|---|---|---|---|---|---|---|---|
| `Z[i]` | `2` | 608 | 324 | 4 | 4 | 245 (131) | 128 (72) | 9 (6) |
| `Z[i]` | `1+i` | 16 | 8 | 2 | 0 | 12 (6) | 16 (8) | 10 (5) |
| `Z[i]` | `2+i` | 3104 | 816 | 0 | 0 | 1760 (468) | 1024 (280) | 1024 (280) |
| `Z[w]` | `2` | 2376 | 1217 | 8 | 5 | 1295 (667) | 324 (172) | 10 (6) |
| `Z[w]` | `1+w` | 324 | 165 | 6 | 0 | 179 (91) | 144 (72) | 14 (7) |
| `Z[w]` | `3` | 40353552 | 20181403 | 1 | 1 | 31610251 (15808376) | 82968 (41580) | 12 (8) |
| `Z[w]` | `2+w` | 823500 | 137775 | 0 | 0 | 702828 (117600) | 279936 (46956) | 279936 (46956) |

A cell counts pairs, with orbits in brackets. The canonical order makes unit steps rare: the 21 pairs with them are 10 untwisted straight segments and 11 twisted walks, the terdragon among them, all listed by `census`; the curves are exactly the 10 segments, every twisted walk landing two words on one point by level 5. At base `3` the full residues admit `6^9 = 10077696` twists, every one in one piece, of which `82968` have positive area and `12` keep every point distinct.

### The walks

By the step law every design with unit steps and more than one point is conjugate, map for map, to a walk: `k` unit steps `u_0, ..., u_(k-1)` from `0` to `b`, digits the partial sums, twists the steps, whose word order runs the vertices of the replacement curve with every flag `F`. A walk is an arc when the `k^L + 1` vertices of level `L`, its end `b^L` included, are distinct at every level; it is a radix walk when its partial sums are pairwise incongruent modulo `b`. Classes are up to reversal, which reverses the steps, and, where `conj(b) = eps b`, the mirror `u -> eps^(-1) conj(u)`; a turn of the plane moves the chord off `b`, so turns are not quotiented.

**No replacement curve with every flag `F` and order `N(b)` to a chord among the seven bases is point-avoiding: each lands two vertices on one point by level 4. Verified**, by the step law, the automaton and every walk (`curves`).

| ring | chord | order | walks | radix walks | radix and plane-filling | arcs | deepest first revisit, every walk |
|---|---|---|---|---|---|---|---|
| `Z[i]` | `2` | 4 | 16 | 4 | 4 | 0 | 2 |
| `Z[i]` | `1+i` | 2 | 2 | 2 | 2 | 0 | 4 |
| `Z[i]` | `2+i` | 5 | 50 | 5 | 1 | 0 | 2 |
| `Z[w]` | `2` | 4 | 34 | 10 | 0 | 0 | 3 |
| `Z[w]` | `1+w` | 3 | 6 | 6 | 6 | 0 | 3 |
| `Z[w]` | `3` | 9 | 116214 | 252 | 14 | 0 | 2 |
| `Z[w]` | `2+w` | 7 | 4235 | 56 | 1 | 0 | 2 |

Every chord of these norms is a turn of one in the table or of its mirror image, so on the triangular lattice the census of this note is settled for flags `F` alone at orders 3, 4, 7 and 9: a point-avoiding curve there runs some copy backwards or mirrored, as the two order-7 curves of the printer's line do. Below dimension 2 the arcs among all walks to the seven chords are the straight segments, the walk `1, w, w^5, 1` to `3`, which is the radix dial's Koch design, and the three-step walks `1, i, 1` to `2+i` and `1, w, 1` to `2+w`, of dimensions `2 log 3/log 5 = 1.3652` and `2 log 3/log 7 = 1.1292`; every other walk revisits a point by level 5 (`curves`).

**Seeded off its curve, a radix word order can fill the plane without revisiting a point. Verified** (`junction`). `Z[w]` base `1+w`, digits `-1, 0, 1`, twists `1, w^2, 1` has unit steps at every level by the step law, and by the automaton it shares no place map, never splits and lands no two words on one point at any level; brute force agrees through level 8. Its twists are the terdragon's walk and its seed sits at `h(0) = 1/2`, so its level path runs the midpoints of the terdragon's edges, which the terdragon never runs twice. The census of codes misses it because `-1` is not a canonical residue at `1+w`.

### The named curves

**The census finds the radix dial's named designs in these classes. Verified** (`named`, which asserts each design of `mrly_math::radix` equal to the census entry it names). The Koch design is the walk `1, w, w^5, 1` to `3`, the one arc class of four steps at base `3`; its digit `2` is not canonical and no turn of its digits is the canonical residues in canonical order, so a census of codes cannot hold it. The terdragon is the canonical pair code `7` at `1+w` with twists `1, w^2, 1`, the walk `1, w^2, 1`: unit steps at every level, positive area, and two words on one point already at level 2. The twindragon and the flowsnake are the untwisted full residues at `1+i` and `2+w`: positive area and distinct points at every level, and no unit steps, which fail at levels 3 and 2.

**The Hilbert curve and the Gosper curve are no radix word order. Proved.** A radix level 2 is `k` copies of level 1, each turned by a unit and read forward, and a turn keeps the sense of every corner. Hilbert's level 1 `0, i, 1+i, 1` turns right twice and the first quarter of its level 2, `0, 1, 1+i, i`, turns left twice; the Gosper curve's flag `R` reads a copy backwards, and its copies 1, 2 and 6 are no turn of its generator (`named`). With every flag `F` their generators are walks of the census: Hilbert's level 1 with the step that joins its copies is the walk `i, 1, -i, 1` to `2`, with positive area and a revisit at level 2, and the flowsnake direction sequence `0, 1, 3, 2, 0, 0, 5` is a walk to `2+w` that shares a place map and revisits at level 2.

## What stays open

- A proof at every order: the census suggests every gentle point-avoiding curve crosses by level 3, which would make the claim a finite statement about two levels; with every flag `F` the radix census decides crossing at every level by an automaton and finds no point-avoiding curve at orders 3, 4, 7 and 9, so what stays open there is a curve with some flag `R`, `M` or `MR`.
- Order 16 with `F, R` needs only an uncut `plain 16`; order 12 and 16 with mirrors need a compiled census; order 19 is the next beyond.
- The sharp-turn density minimum over replacement curves beyond order 7.
- Self-similar curves built from several motifs.
