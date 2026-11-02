---
title: Slices
lead: The diagonal slice of the solid cube: the `6n` census, centered-hexagonal vertices, and the splitting-prime rule.
figure: research-slices
slug: slices
---

Cut the solid cube of odd side `2k - 1` through its centre, perpendicular to the main diagonal, and the section is a regular hexagon tiled by `6*side^2` unit equilateral triangles. This page is the census of that mesh - triangles, edges, vertices, an [Euler characteristic](https://en.wikipedia.org/wiki/Euler_characteristic) that never moves - and of what the [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) designs do to it: which fills partition it, which fall into many pieces, and which pierce it with holes. The mesh itself is classical lattice geometry and no novelty is claimed for it; the designs are where the specific content lives. The back half leaves the slice for the surface census of the same designs in 3D. Lay these hexagons over each other at every odd side and they interfere; that stack, its exact ink laws and the arithmetic of its lines are [the hexagon moire](../../hexagon/notes/hexagon.md).

The generators are `six::census`, `six::topology`, `counts::six` and `three::census` in `crates/math`, each test named for the claim it pins, except the profile identity, which `cuts/lab/py/slice-ladder-controls` checks, and the slice ink, whose generator `demos::walsh_spectrum` and its fixture went with the `demos` crate; `mrly_num::boolean::walsh_spectrum` still reads the spectrum. The [slices demo](/three/) cuts the middle plane of any design at any odd side, draws the hexagon triangle by triangle, and reports the mesh census beside its closed forms, the carpet-and-net partition, and the pieces and holes over `k = 1..16`.

One boundary first. [The cuts page](../../cuts/notes/cuts.md) also slices a solid along `x + y + z = s` and also finds a hexagon, but everything else differs: there the solid is `bang dim 3, code 126` at base 2, fractal from the start, and the result - every slice a Sierpinski gasket, scheduled by the binary digits of the height - comes from a digit argument, with no mesh and no census. Here the solid is the plain filled cube of odd side, the cut is the single middle plane, and the result is a triangular mesh counted directly. The two pages share a plane and a hexagon and not one number.

## The solid slice, counted

Take the plane `x + y + z = 3*side/2` through the centre of `[0,side]^3`. It crosses the cells of the three middle diagonal layers, meeting each crossed cell in a triangle or a hexagon; cutting the hexagons into six from their centres leaves pieces that are all equilateral triangles of one size, tiling the section - a regular hexagon whose area is `6*side^2` triangles. The rim contributes `6*side` triangle sides, and the piece vertices are exactly the triangular-lattice points of the closed hexagon, `3*side^2 + 3*side + 1` of them. Every interior side is shared by exactly two triangles, so counting sides two ways gives `3*F = 2*I + B` at `F` the triangles, `I` the interior edges and `B` the boundary edges, which fixes the edge count, and the whole census follows. (Proved; each step is also asserted numerically in the build, not waved through.)

| count | in `k` | factored | at side 3 |
|---|---|---|---:|
| triangles | `24*k^2 - 24*k + 6` | `6*side * side` | 54 |
| boundary edges | `12*k - 6` | `6*side` | 18 |
| edges | `36*k^2 - 30*k + 6` | `6*side * (3k-1)` | 90 |
| interior edges | `36*k^2 - 42*k + 12` | `6*side * (3k-2)` | 72 |
| vertices | `12*k^2 - 6*k + 1` | `3*side^2 + 3*side + 1` | 37 |

Every count in the table carries the factor `6*side` except the last: the vertex count is odd and carries none. The carry demo runs the same cut at any base and dimension, where the count remembers only a carry, so `ceil(D/2)` past terms decide every term and the growth exponent misses the generic value, above it at odd dimensions and below it at even ones.

**Prior art.** The vertex count `12*k^2 - 6*k + 1` is [OEIS A154105](https://oeis.org/A154105), printed there as `12n^2 + 18n + 7` under an index shift of one: substituting `k = n + 1` gives `12(n+1)^2 - 6(n+1) + 1 = 12n^2 + 18n + 7`, and both forms return `7, 37, 91, 169, 271` over the first five terms. This lane also cites A154105 in [sequences](../../shared/sequences.md)'s A395241 row. `https://oeis.org/search?q=id:A154105&fmt=json` returns the name `a(n) = 12*n^2 + 18*n + 7.`, offset 0, and data `7, 37, 91, 169, 271, 397, ...`, so the name, the offset and the five terms above are all Verified at source.

Two consequences are exact, not asymptotic. The slice's own surface-to-volume ratio - boundary edges over triangles - is `6*side / 6*side^2 = 1/side` at every size, and `V - E + F = 1` identically - the slice is a topological disk at every `k`. (Proved; Verified with the Euler characteristic counted directly, never from the algebra, at every size built.) The polynomials are confirmed by direct census at `k = 1..8` (`mrly_math::six::census` test `the_five_closed_forms_match_the_census_to_eight`, which reads the Euler characteristic as 1 at every size it builds), by a blind quadratic fit through `k = 1..3` that reproduces `k = 4..10` with zero residual (test `a_blind_quadratic_fit_reproduces_the_wider_slices`), and by fresh builds at `k = 12, 16, 20` - at side 39: 9126 triangles, 13806 edges, 4681 vertices, 234 boundary edges (test `fresh_builds_at_the_wide_sides_hold_the_forms`). (Verified.)

One lemma serves every fill census below. **Lemma (Proved).** For any set of triangles of the mesh, the adjacency graph on them - one node per triangle, one edge per shared side - has exactly `E' - B'` edges, where `E'` and `B'` are the edge and boundary-edge counts of the sub-mesh those triangles span: every interior edge of the sub-mesh joins exactly two of its triangles, every boundary edge one, and no lattice edge lies in three. Verified by direct enumeration on twenty meshes - the solid and four design fills at base-3 levels one to four, sides 3 to 81 - with the adjacency count and the interior-edge count collected from different maps; the carpet's fill at side 81 reads `E' = 28188`, `B' = 6642`, adjacencies `21546` (`mrly_math::six::census` test `the_fill_adjacency_counts_the_sub_mesh_interior_edges`). The sub-mesh is the load-bearing word: read against the full hexagon mesh the identity is false - a fill-void edge is interior to the hexagon but joins no two fill triangles - and at side 27 the carpet's adjacency count is 2880 against the full mesh's `E - B = 6480` (test `the_lemma_needs_the_sub_mesh_and_not_the_hexagon`).

## The profile identity

For a binary array `X` let `P_X(t) = Sum over filled cells of t^(r+c)` be the anti-diagonal slice profile - the digit transform `hat F` of [mobius](../../mobius/notes/mobius.md) taken over the filled cells of `X` and read on the diagonal, `P_X(e(theta)) = hat X(theta, theta)`, written `P_X` on this page. Then

```
P_{A (x) B}(t) = P_A(t^(side_B)) * P_B(t)
```

because `r + c = side_B*(r_A + c_A) + (r_B + c_B)` under [the Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product). **Proved.**

- Stationary corollary: `P_level(t) = prod_{j=0}^{level-1} P(t^(base^j))`, the substitution product at one base.
- Mixed corollary: for an ordered word, `P_w(t) = prod_i P_i(t^(side_{i+1}...side_level))`, the mixed-radix substitution product.
- One identity therefore covers fractal slices, mixed-product slices, and the dimension ladder in [cuts](../../cuts/notes/cuts.md) and the base ladder in [spectra](../../spectra/notes/spectra.md). It is the generating function of a digit-restricted set, so every tool for those sets applies to it.
- **Verified** on all words of length 2 and 3 at base 2, dim 2, over the 15 nonempty `2x2` tiles, zero mismatches, `cuts/lab/py/slice-ladder-controls`.
- Take the whole profile and never one coefficient: the total and the centre coefficient are order-blind, and the profile as a whole is not.

Open, and cheap:

- Is the profile a rational series in the word, as components, boundary, Euler characteristic and holes all are? Its Hankel rank is the first thing to measure.
- Which slice functionals are order-blind? Total is, the centre coefficient is, peak and support are not. Classify them.
- Non-diagonal lattice planes: `x + y = c` decomposes because the coefficient vector is all ones, and a general primitive `alpha.x = c` does not. That is where the carry automaton enters, and nothing is known there.

## The vertex count is centered hexagonal

The closed hexagon of side `side` holds `CH(side+1)` lattice points, where `CH(m) = 3*m^2 - 3*m + 1` is the centered hexagonal sequence `1, 7, 19, 37, 61, ...` (OEIS A003215). With side `2k - 1` that is `CH(2k)`: the slice vertex count is the centered hexagonal number at even index. (Proved.)

The arithmetic attached is short and mostly classical. `12*k^2` and `6*k` are both divisible by 3, so the vertex count is `1 mod 3` for every integer `k`. (Proved.) And `CH(m) = m^3 - (m-1)^3 = m^2 + m*(m-1) + (m-1)^2`, a difference of consecutive cubes, so the centered hexagonal numbers that are [prime](https://en.wikipedia.org/wiki/Prime_number) are exactly the cuban primes, OEIS A002407 - a classical family, not a discovery of this page. What the slice adds is only the specialisation: its vertex counts walk the even-index half of that sequence.

The same form places the counts in [the Eisenstein integers](https://en.wikipedia.org/wiki/Eisenstein_integer) `Z[omega]`, the ring behind the hexagonal lattice on [the bases page](../../bases/notes/bases.md). A rational prime splits in `Z[omega]` exactly when it is `1 mod 3`, stays inert when it is `2 mod 3`, and 3 alone ramifies - standard algebraic number theory, cited rather than re-derived. Since every vertex count is `1 mod 3`, a prime vertex count is always a splitting prime, never inert and never the ramified 3, and it arrives with its own witness: `CH(m)` is already a value of the norm form at `(m, m-1)`. One convention note: the bases page writes the norm as `N(a + b*omega) = a^2 - a*b + b^2`; the two forms trade under `b -> -b`, so in that convention `CH(m) = N(m + (1-m)*omega)`.

The list itself. For `k = 1..20` the prime vertex counts are `7, 37, 271, 397, 547, 919, 1657, 1951, 2269, 4219`, at `k = 1, 2, 5, 6, 7, 9, 12, 13, 14, 19`; the other ten values, `91, 169, 721, 1141, 1387, 2611, 2977, 3367, 3781, 4681`, are composite. Continuing to `k = 40` adds `5167, 6211, 7351, 9241, 12097, 13669`. All sixteen primes are `1 mod 3`. (Verified by `mrly_math::counts::six` test `the_prime_vertex_counts_are_cuban_with_norm_witnesses`, with a norm-form witness exhibited for each - `4219 = 37^2 + 37*38 + 38^2`; the `CH(2k)` form and the `1 mod 3` residue run to `k = 40` in test `the_slice_vertex_count_is_centered_hexagonal_at_even_index`.) The entry `k = 19`, value 4219, is easy to skip and is prime.

## Four fills on one hexagon

The rest of the slice story needs the four historical families, and their naming needs care. The slice statements below are about the level-1 fill at odd side: the design's parity rule applied to the `side^3` grid directly. A triangle of the mesh belongs to a design when the cell it was cut from is filled.

| name here | rule | design | class | slice fill at side 3 |
|---|---|---|---|---:|
| carpet | at most one odd coordinate | `bang dim 3, code 23` | `bang dim 3, code 23` | 42 |
| net | at least two odd coordinates | `bang dim 3, code 232` | `bang dim 3, code 23` | 12 |
| tree | `x` and `y` both even | `bang dim 3, code 3` | `bang dim 3, code 3` | 18 |
| void | all coordinates one parity | `bang dim 3, code 129` | `bang dim 3, code 24` | 12 |

Carpet and net are one self-complementary class - that is why [the core page](../../core/notes/core.md) aliases both names to `bang dim 3, code 23` - and `bang dim 3, code 232` is its complement member, a different truncation of the same class, not a second design; its fill polynomial `4*k^3 - 9*k^2 + 6*k - 1` is the worked corollary on [the method page](../../method/notes/method.md). This page does not call the fourth family *void*, because the core page's alias *void* is the canonical `bang dim 3, code 24`, which fills `2*k^3 - 3*k^2 + k`, while `bang dim 3, code 129` fills `2*k^3 - 3*k^2 + 3*k - 1 = k^3 + (k-1)^3`, the centered cube numbers (OEIS A005898) - same class, different truncation, different polynomial. Nor is it the checkerboard: fill where `i + j + l` is even is `bang dim 3, code 105`, a different design again. This page calls it the void, after its two corners `(0,0,0)` and `(1,1,1)`. (Verified: orbit walks over the 48 signed permutations and cell-by-cell counts, `mrly_math::six::topology` test `the_four_families_fill_the_slice_and_name_their_classes`.)

**Carpet and net partition the hexagon.** Their corner rules partition `{0,1}^3` outright - every parity vector has popcount at most 1 or at least 2 - so every cell of any grid is filled by exactly one of the two, and in particular every crossed cell hands its whole section to exactly one. The two slice fills therefore partition the solid hexagon cell for cell:

```
carpet(side) + net(side) = 6*side^2
```

(Proved by the corner partition; Verified at every odd side 1..31, with disjointness and covering tested triangle by triangle: `42 + 12`, `72 + 78`, `204 + 90`, `210 + 276`, `486 + 240` at side 3, 5, 7, 9, 11, up to `3696 + 2070 = 5766` at side 31, `mrly_math::six::topology` test `carpet_and_net_partition_the_hexagon_triangle_by_triangle`.) Two cautions. This does not follow from the 3D complement identity - a volume identity says nothing about how a section decomposes; the extra fact, tested directly, is that the section decomposes along the same cells. And the split is not by the parity of `i + j + l`: the crossed cells lie in three consecutive diagonal layers and every layer splits between the two families - at side 5 the three layers split `15/3`, `7/12`, `15/3` - so the partition is by popcount, not by parity. (Refuted; Verified.) Measured by section area instead of by triangle, the same identity extends to even side, where the split is exactly half and half. (Verified at side 1..16 by `mrly_math::six::topology` test `the_layer_weighted_area_is_a_second_route_to_the_fill`.)

## Components and holes

The carpet's slice fill runs two regimes in alternation, both governed by the same centered hexagonal numbers as the vertex count. At odd `k` it falls into `CH((k+1)/2)` disjoint pieces with no holes; at even `k` it is one connected piece pierced by `CH(k/2)` holes. Over `k = 1..14` the component counts run `1, 1, 7, 1, 19, 1, 37, 1, 61, 1, 91, 1, 127, 1` and the hole counts `0, 1, 0, 7, 0, 19, 0, 37, 0, 61, 0, 91, 0, 127`. (Verified at `k = 1..14` by `mrly_math::six::topology` test `the_carpet_slice_counts_its_pieces_and_holes_two_ways`, which counts pieces by triangle adjacency and takes the hole count twice by routes sharing nothing: the piece count less the Euler number of the filled sub-mesh, and again as void regions the rim never reaches. `CH(6) = 91` and `CH(7) = 127` at `k = 11..14` are predicted by the law and come out right in both regimes, the law itself asserted in the same test.)

The other families sort cleanly. The tree's and the void's slices are hole-free at every `k` checked (`1..10`), and the net's carries the same sequence shifted by one in `k` - holes `1, 7, 19, 37` at odd `k = 3, 5, 7, 9` - so carpet and net are the two families that puncture the hexagon, in opposite phase. (Verified, `mrly_math::six::topology` test `the_other_families_puncture_in_opposite_phase`, both hole routes agreeing at every size.) `CH` is a classical sequence; its double occurrence in this slice - once as the vertex count, once as the component-and-hole law - is the fact.

## Surface: the face-count recurrence

Now the 3D families, no slice. Substitute a tile of side `side` to level `i` and count unit faces: `V(i)` visible, `H(i)` hidden between two filled cubes. Cells are multiplicative, `cells(i) = fc^i` with `fc` the tile's fill, so `V(i) + H(i) = 6*fc^i` at every level.

**The recurrence (Proved).** Substitution creates no new contacts inside a copy; the only new hidden faces arise across an interface where two filled cells of the tile meet face to face. Opposite faces of a copy carry the same fill pattern - checked as arrays, not merely as counts - and face fills are multiplicative, so an interface along axis `a` contributes `p_a^i` contacts at level `i`, with `p_a` the tile's face fill on that axis, and each contact hides two faces. Writing `adj_a` for the tile's face-adjacent filled pairs along axis `a`:

```
V(i+1) = fc*V(i) - 2 * sum_a adj_a * p_a^i
H(i+1) = fc*H(i) + 2 * sum_a adj_a * p_a^i
```

When every adjacency-carrying axis shares one face fill `l2`, the source term collapses to `2*W*l2^i` with `W = sum_a adj_a`, and the system closes: `H(i) = 2*W*(fc^i - l2^i)/(fc - l2)`, `V(i) = 6*fc^i - H(i)`. Both counts then live in the span of `fc^i` and `l2^i` - exactly the statement that `(V, H)` evolves by one fixed `2x2` matrix with eigenvalues `(fc, l2)`, unique whenever `W > 0` and `l2 != fc`.

| family | `fc` at side 3 | `W` | `l2` | `V(i)` |
|---|---:|---:|---:|---|
| carpet | 20 | 24 | 8 | `2*20^i + 4*8^i` |
| net | 7 | 6 | 1 | `4*7^i + 2` |
| tree | 12 | 8 | 4 | `4*12^i + 2*4^i` |
| void | 9 | 0 | none | `6*9^i` |

**Prior art.** The carpet row is not new. `2*20^i + 4*8^i` is [OEIS A332705](https://oeis.org/A332705) verbatim, "Number of unit square faces (or surface area) of a stage-n Menger sponge", formula contributed by Allan Bickle, Nov 2022. Its printed terms `72, 1056, 18048, 336384` are A332705(1..4), recomputed here from the closed form and matching exactly. The identification is recorded in [REFS.md](../refs.md) and belongs here, at the point of derivation. A direct request to `https://oeis.org/search?q=id:A332705&fmt=json` returns the name verbatim as quoted above, offset 0, data `6, 72, 1056, 18048, 336384, ...`, and the formula block *"From Allan Bickle, Nov 28 2022: a(n) = 2*20^n + 4*8^n"*. So the closed form, the attribution and the term alignment `72, 1056, 18048, 336384 = A332705(1..4)` are all Verified at source.

(Verified against brute-force face counts at side 3 to level 4, side 5 to level 3, side 7 to level 2, `mrly_math::three::census` tests `the_face_ledger_prints_the_family_closed_forms` and `the_face_matrix_fits_its_eigenvalues_and_predicts`; the carpet at side 3 reads visible `72, 1056, 18048, 336384` against hidden `48, 1344, 29952, 623616`. Matrices fitted exactly over the rationals from levels 1..3 predict the next level, and trace and determinant give `(fc, l2)` in all six fitted cases - e.g. carpet at side 3: `[[12, 4], [8, 16]]`, trace `28 = 20 + 8`, determinant `160 = 20 * 8`.) Across bases the second eigenvalue is a face count of the family: `side^2 - floor(side/2)^2` for the carpet, `floor(side/2)^2` for the net, `ceil(side/2)^2` for the tree - `8, 1, 4` at side 3, `21, 4, 9` at side 5, `40, 9, 16` at side 7, nine of nine. (Verified.) The tree is the case that makes the rule bite: its copies touch only along one axis, and it is that axis's face fill, 4 rather than 6, that drives the recurrence.

The void row is the degenerate case, and a second eigenvalue of 0 for it is wrong. (Refuted.) With no face-adjacent pair, `W = 0`, the hidden channel is never fed, and the substitution step is literally `fc` times the identity - eigenvalues `(fc, fc)` - while the face data sits on a single ray, so no `2x2` matrix is determined by it at all. The law for that family is the next section's, with no second eigenvalue to quote.

This matrix shares nothing but a name with the digit-borrow [transfer matrices](https://en.wikipedia.org/wiki/Transfer-matrix_method_(statistical_mechanics)) this lane uses elsewhere: those are automata on the binary digits of a position, this one is a `2x2` face ledger. Two arguments, one word.

## The whole skeleton, counted

The surface recurrence counts the faces open to void. The same Kronecker structure counts every unit face of the design at once: the vertices, edges, squares and cubes of the polyhedral complex of the union of closed unit cubes, each counted once however many cubes share it. That is exactly what the ledger renders as vertices, edges, faces and euler: `counts::faces` agrees with `three::census` on the level-one and level-two sponge, `64, 144, 96, -4` and `896, 2304, 1728, -80`, and with `two::census` on the carpet (test `the_rendered_skeleton_is_the_census`).

**The patterns.** Fix a tile `T` of side `n` in dimension `d` and let `D_L` be its `L`-fold Kronecker power, the cells of `[0, n^L)^d` whose every base-`n` digit vector lies in `T`. A pattern is a nonempty set `P` of cells of the `2^d` window `{0,1}^d`, normalised so that some cell of `P` has coordinate `0` on every axis; `N_P(L)` counts the translates `v + P` inside `D_L`. The single cell gives `N(L) = fc^L`, the pair along axis `a` gives the adjacent pairs of `D_L`, the full window gives its filled `2 x ... x 2` blocks.

**Faces from patterns (Proved).** A unit `k`-face is a window of type `B`, the `d - k` axes on which it is a point, together with the `2^(d-k)` cells around it, and the face belongs to the union exactly when one of those cells is filled. Inclusion and exclusion over the filled cells of the window, every nonempty subset `U` of the window cells counted once per window containing it, turns the count of `k`-faces into a signed sum of pattern counts:

```
F_k(L) = sum over B with |B| = d - k, sum over nonempty U in {0,1}^B of (-1)^(|U|+1) N_U(L)
chi(L) = sum over nonempty P in {0,1}^d, one per translation class, of (-1)^(|P|+1) N_P(L)
```

The second line is the nerve of the cover by closed cubes, and it is also the alternating sum of the first over `k`. In the plane the vertex count reads `4 fc^L - 2 N_x - 2 N_y - N_diag - N_anti + (the four L-trominoes) - N_square`, and on the level-one carpet that is `32 - 20 + 4 - 0 = 16`.

**The carry recurrence (Proved).** Write a translate `v + P` in `D_(L+1)` by its lowest digit `w` in `[0, n)^d` and its higher digits `u` in `D_L`. Adding `p` in `P` to `v` moves the lowest digit to `(w + p) mod n` and carries one into `u` exactly on the axes `W` where `w` sits at `n - 1` and `p` has a `1`. So `v + P` lies in `D_(L+1)` exactly when the wrapped pattern `w + P mod n` lies in `T` and the carried pattern `P ^ W`, the cells of `P` with every coordinate outside `W` set to `0`, lies in `D_L` at `u`:

```
N_P(L+1) = sum over Q of M[P][Q] N_Q(L),  M[P][Q] = #{w in [0,n)^d : w + P mod n in T, P ^ W(w) = Q}
```

The matrix `M` is the only tile data, read off the tile once. Its entries for the pair along `a` are the exposure data, `M[pair][cell] = adj_a` and `M[pair][pair] = S_a`, so `N_pair(L+1) = adj_a fc^L + S_a N_pair(L)`, and the exposure recurrence `V(L+1) = fc V(L) - 2 sum adj_a S_a^L` is the case `P = pair` of this one. `P ^ W` never has a larger support than `P`, so `M` is triangular in the support order, with the diagonal entry `lambda_P = M[P][P]` the spanning count of `P`: the digits `w` at `n - 1` on the whole support of `P` that still fit it wrapped. The single cell has `lambda = fc`.

**The roots (Proved).** The generating function `G_P(x) = sum N_P(L) x^L` satisfies `(1 - lambda_P x) G_P = N_P(0) + x sum over Q below P of M[P][Q] G_Q`, so by induction down the support order `G_P` is a rational function whose denominator is `prod over lambda of (1 - lambda x)^(c_lambda(P))`, where `c_lambda(P)` is the largest number of patterns of spanning count `lambda` on one carry chain from `P` (a path of nonzero off-diagonal entries), and whose numerator has degree below the denominator's degree plus `z`, with `z` the most patterns of spanning count `0` on one chain (a zero-count pattern adds no denominator degree, so `z = 0` gives degree below the denominator's). Every face count and the Euler characteristic therefore obey, from `L = 0`, the linear recurrence with characteristic polynomial `prod (x - lambda)^(c_lambda)` over the patterns their signed sums reach, and each is a sum of `lambda^L` times a polynomial of degree below `c_lambda` from `L = z` (code 123456789 has `z = 2` for its vertex count and leaves the power form at `L = 0` and `1`, `z = 1` for its edge count); when every `c_lambda` is `1` it is a plain sum of powers with rational coefficients fixed by the first terms. (Proved; `counts::faces` test `the_terms_obey_the_recurrence` checks the recurrence on twelve levels for six codes, `the_recurrence_matches_the_rendered_skeleton` the counts against a brute-force face census of the rendered design for nine codes across dims 2 and 3 and bases 2 and 3, to level 4 in the plane and at base 2 in space, level 3 at base 3 in space.)

**The classics (Proved).** The carpet and the sponge have the roots `1, 3, 8` and `1, 3, 8, 20`, all simple, so every count is a sum of powers (`counts::faces` test `the_classics_close`):

| count | carpet, side 3 | sponge, side 3 |
|---|---|---|
| vertices | `(40 + 56*3^L + 44*8^L)/35` | `(16320 + 25536*3^L + 28424*8^L + 20160*20^L)/11305` |
| edges | `(8*3^L + 12*8^L)/5` | `(192*3^L + 408*8^L + 420*20^L)/85` |
| squares | `8^L` | `2*8^L + 4*20^L` |
| cubes | | `20^L` |
| euler | `(8 - 8^L)/7` | `(192 - 38*8^L - 21*20^L)/133` |

The carpet vertices `16, 96, 688, 5280` are [OEIS A293143](https://oeis.org/A293143) from its second term, whose recurrence `a(n+1) = 8 a(n) - 8 (3^n + 1)` has the same roots; the sponge vertices `64, 896, 15616, 295808` are [A293144](https://oeis.org/A293144), given there by the recurrence `a(n) = 20 a(n-1) - 24 A293143(n-1)`, and it already carries the closed form above, which equals the page's form after clearing denominators (`64*255 = 16320`). The sponge edges `144, 2304, 42048, 810432`, squares `96, 1728, 33024, 648192` and Euler characteristics `-4, -80, -1408, -26432`, and the carpet edges `24, 168, 1272, 9960`, are absent from the OEIS dump at four terms (Conjecture until a live read). The sponge squares are the surface plus half the buried faces, `(exposed + 6*20^L)/2`, so the `(d-1)`-faces close by the exposure alone in every dimension: `F_(d-1) = d fc^L + exposed/2` (checked on six codes in `the_terms_obey_the_recurrence`). The carpet Euler characteristic says `1 - chi = (8^L - 1)/7` holes, one per cell removed along the way.

**The multiplicity bites (Proved).** The bound `c_lambda` is attained: `bang dim 3, code 123456789` at base 3 has vertex roots `0, 0, 1, 1, 2, 3, 5, 16`, and its vertex count `8, 57, 699, 10463, 164778, 2626182` carries the term `(11/30) L` beside its powers (test `a_repeated_root_carries_a_linear_term`), so the vertex count of a design is not always a sum of powers; the edge count of the same design has the simple roots `0, 1, 2, 3, 5, 16`. The double root `0` is transient and moves the first two terms only.

**What this gives the ledger.** `counts::skeleton(code, number, dimension, base, level)` returns every `k`-face count and the Euler characteristic from the tile alone, `counts::skeleton_recurrence(code, number, dimension, base, k)` the recurrence coefficients, `Placements::roots` the roots with their multiplicities, in dimensions one to four, where the window patterns fit a machine word.

## No hidden faces

**(Proved.)** Two cells share a unit face exactly when they differ by 1 in one coordinate - which flips exactly one bit of the parity vector. So a design ever hides a face if and only if two of its filled corners sit at Hamming distance 1. The void's corners are at distance 3, so no two of its filled cubes ever meet face to face, at any side and any level: every face is exposed, and

```
surface(k) = 6 * cells(k), cells(k) = 2*k^3 - 3*k^2 + 3*k - 1 = k^3 + (k-1)^3
```

with `surface(level) = 6 * 9^level` across fractal levels at base 3 - no hidden-face correction term, ever. (Proved; Verified by direct face counts at `k = 1..12`, in 2D as `edges = perimeter`, and at levels 1..4, `mrly_math::three::census` test `the_void_buries_no_face`.) Hamming distance is preserved by cube symmetry, so total exposure is a class property. The test is also exhaustive: sweeping all 256 designs at side 3, 5, 7, `surface = 6 * cells` holds for exactly the 35 designs whose filled corners are pairwise at Hamming distance at least 2 - the independent sets of the cube graph - and they form whole classes: `bang dim 3, code 0`, `bang dim 3, code 1`, `bang dim 3, code 6`, `bang dim 3, code 22`, `bang dim 3, code 24`, `bang dim 3, code 105`. (Verified by `mrly_math::three::census` test `total_exposure_holds_for_the_independent_corner_sets` and `mrly_math::bang::universe` test `total_exposure_names_the_independent_corner_sets`.)

## Corners

A last piece of 3D background, because it explains an asymmetry the complement identity leaves open. **(Proved.)** At odd side `2k - 1` a grid-corner coordinate `a` sees the cell indices `a-1` and `a`, clipped to the grid: the boundary values 0 and `side` see only an even index, the `2k-2` interior values see one of each parity. So a design touches every one of the `(side+1)^3` grid corners if and only if its rule contains the all-even corner - the grid corner at the origin can be touched by no other cell. Carpet, tree and the void all contain it, so their solids touch the whole grid, `8*k^3` corners. The net does not - it needs two odd coordinates - so a grid corner is touched exactly when at least two of its coordinates are interior:

```
net vertices = m^3 + 6*m^2 = 8*k^3 - 24*k + 16 = 8*(k-1)^2*(k+2), m = 2k-2
```

short of the full grid by `24*k - 16`. (Proved; Verified by brute force at `k = 1..20` in 3D and `k = 1..24` in 2D, where the same argument gives `4*k^2 - 4`, and the all-even criterion swept over all 256 designs at `k = 1..3`: exactly the 128 rules containing the all-even corner touch every corner, `mrly_math::three::census` tests `the_net_falls_short_of_the_grid_corners` and `the_all_even_rule_touches_every_grid_corner`, the 128 rules named again by `mrly_math::bang::universe` test `half_the_rules_hold_the_all_even_corner`.) Two footnotes. Odd side is load-bearing: at even side the boundary coordinate sees an odd index and every family loses corners. And this is the geometric half of the complement story on [the method page](../../method/notes/method.md): the carpet and its complement trade cells exactly, but not vertices, because complementing the rule loses the all-even corner that every boundary grid corner depends on.

## The slice ink, proved order by order

The Walsh spectrometer - the exact quasipolynomial that reads a 3D parity design's slice ink off its Walsh spectrum, level by level, `ink_dim(side) = S0 - (1/2) S3 s + [(2/3) S1 - (1/3) S2 s]/side + [(2/3) S2 - ((1/3) S1 + (1/2) S3) s]/side^2` at `s = (-1)^((3*side-1)/2)`, with `S_j` the design's level-`j` Walsh sum, proved in the [walsh-spectrometer](../../hexagon/claims/diagonal-slice-stack.md) lane - has all three orders proved: the constant terms and the `1/n` and `1/n^2` terms alike. **(Proved.)**

Fix odd `side` and let

```
P_side = { (x,y,z) in {0..4*side-1}^3 : x + y + z = 6*side - 2, z even }
```

be the plane the slice lives on. For a macro-parity triple `e` in `{0,1}^3` let `N_e(side)` count the points of `P_side` whose macro coordinates `(floor(x/4), floor(y/4), floor(z/4))` are congruent to `e` mod 2. Then `N_e(side)` depends only on the weight `k = |e|`, and with side `2h + 1`,

```
N_k(2h+1) = [t^(3h+1)] (1 + 6t + t^2) E_h(t)^(3-k) O_h(t)^k,
E_h(t) = 1 + t^2 + ... + t^(2h), O_h(t) = t + t^3 + ... + t^(2h-1).
```

The three ingredients are all elementary. Writing `x = 4u + r_x` and likewise for `y` and `z`, the residue sum `r_x + r_y + r_z` can only be `0`, `4` or `8`, with multiplicities `1, 6, 1` read off `(1 + t + t^2 + t^3)^2 (1 + t^2)` - the same `1:6:1` law as the constant terms, now kept at the lower orders instead of averaged away. Each macro variable ranges over `0..2h`, so a parity triple of weight `k` contributes `E_h^(3-k) O_h^k`. And `6*side - 2 = 0 mod 4` because the side is odd. The generating function proves on sight that only the weight matters and never which entries are odd.

Extraction is a finite binomial calculation, not an asymptotic: substituting the two geometric series and using `[t^m](1-t^2)^(-3) = C(m/2 + 2, 2)` for even `m >= 0` and `0` otherwise gives the four exact quadratic quasipolynomials, `(9*side^2 + 18*side + 21)/8` at `k = 0` and side `1 mod 4` down to `(9*side^2 - 18*side + 21)/8` at `k = 3` and side `3 mod 4`. Since slice ink is linear in the eight `N_e(side)`, dividing by the slice normalization can produce only constant, `1/side` and `1/side^2` terms, with the two residue classes carried by `s = (-1)^((3*side-1)/2)`. No approximation enters anywhere. (`demos::walsh_spectrum`, gone with the `demos` crate, was the generator: it read a code's four level sums and printed the law's fill at every odd side, and the `demos` fixture held those fills against `slice_census`, the hexagon counted triangle by triangle, at codes 23 and 11 over side 1..11, zero mismatches. **(Verified.)** Neither test code carries a weight-3 corner, so `N_3(side)` is untouched by the crate; the eight-triple split, the popcount aggregate against `C(3,k) N_k(side)` and the 28 layers to side 55 are the `walsh-spectrometer` lane's own check.) The spectrometer demo points the slice at any of the 256 designs, bars its four level sums, and plots this law against the ink counted triangle by triangle at every odd side, with a mystery mode that hides the code and asks you to read the recipe back.

Scope guard. This page is the base-3 census at dim 3 and nothing on it speaks to the slice DIMENSION at other dimensions. The sign law `sgn(dim_slice - (log(fill) / log(base) - 1)) = (-1)^(dim+1)` at `dim_slice` the dimension of the slice, Conjecture S, is on no note - [spectra](../../spectra/notes/spectra.md) names the lane that carries it, under WHERE THE REST LIVES - and the transfer matrix both it and Conjecture R run on is defined once in [cuts](../../cuts/notes/cuts.md). Nothing here should be read as evidence for or against either.

## Coprimality on the slice

One pointer. The arithmetic of the slice - which of a design's slice points are visible from the origin - lives on [the coprimality spine](../../coprime/notes/coprime.md) under "Coprimality on the slices", with the generator `coprime/lab/py/slice-coprimality`. The one-line summary: a slice point's gcd divides its height, so slice coprimality is finite arithmetic per height, the base prime peels the central slice one step off-centre, each foreign prime costs `1/p^2` rather than the solid's `1/p^3`, and the central slice's bill is the factorization of the repunit in the base, so its density never converges.

## Where the numbers live

Everything above is computed by `crates/math`: `six::census` carries the mesh census, its extrapolations and the adjacency lemma, `six::topology` the fills, the partition, the pieces and the holes, `counts::six` the closed forms and the prime vertex counts, and `three::census` the face recurrence, total exposure and corners, with `bang::universe` naming the corner-set classes; `apps/slices` puts them in the browser behind [the slices demo](/three/), calling `six.cut`, `six.census`, `six.components` and `six.holes` of `@mrlyprod/math`. The profile identity is checked by `cuts/lab/py/slice-ladder-controls`, the slice pointer is owned by `coprime/lab/py/slice-coprimality`, and the slice ink was rebuilt by `demos::walsh_spectrum`, whose law fills the `demos` fixture held against `slice_census` triangle by triangle; both are gone with the `demos` crate. Sequences that come out of this construction are held to the standard in [the sequence ledger](../../shared/sequences.md).
