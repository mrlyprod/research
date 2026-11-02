---
title: Crop
lead: The circle on a design: the corner disc count, its log-periodic main term, and the crossing shell read as a rooted tree whose transfer operator is derived from the geometry.
figure: research-crop
slug: crop
---

A crop lays an exact shape - a ball, a box, a diamond, any rational polytope - over a design's grid and classifies every cell as In, Cut or Out, in exact integer arithmetic with no floats: ball tests clear denominators and compare squared distances, half-space tests evaluate the linear form at two extreme corners. Cropping then zeroes the filled cells the shape rejects, keeping Cut cells on request, and `Shape::Anti` swaps In with Out while Cut stays. This page is the census of what a circle and a diamond keep of the carpet, what a sphere and an octahedron keep of the sponge, and the one open lane the Cut column points at: the dimension of a curved slice.

The generators are `mrly_math::shape` in `../crates` - `classify`, `regions`, `crop`, `refine`, `census` - of which `classify`, `crop`, `refine` and `census` are each pinned by tests named for their claims, including an independent `2^dim`-corner oracle, while `regions` carries no test of its own name and is exercised only through the others. `lab/rs/crop-counts` is the one pass that prints every number below. The crop demo draws a named shape over a design and counts the in, cut and out regions before anything is rendered.

## Two identities, by construction

**Partition (Proved).** Every cell is exactly one of Out, Cut, In - `classify` returns one region, `census` tallies cells and filled cells per region, and the three filled tallies sum to the design's fill. `crop(types, shape, true)` keeps the filled cells of In and Cut, `crop(types, shape, false)` keeps In alone, so the two crops bracket the boundary from both sides.

**Anti-crop complement (Proved).** `Shape::Anti` flips In and Out and fixes Cut, so `crop(types, Anti(shape), false)` keeps exactly the filled Out cells and `crop(types, Anti(shape), true)` keeps Out and Cut: the crop and the anti-crop with the complementary cut rule partition the filled set exactly, whichever side gets the boundary. Both identities are read off the definition of `classify` and asserted, both ways, on all 118 configurations `lab/rs/crop-counts` prints, and independently by the `mrly_math::shape` partition tests. (Proved; Verified.)

## The inscribed ball and diamond, level by level

Code 7 at `dim = 2` is the carpet, code 23 at `dim = 3` the sponge; the shape is centered at `(1/2, ..., 1/2)` with the inscribed radius `1/2`, so the ball touches the four or six face midpoints and the diamond is the inscribed cross-polytope. `filled_in` and `filled_cut` are the filled cells fully inside and crossing the boundary; `exposed_after` counts the exposed unit faces of the keep-cut crop, a face being exposed when its neighbour is empty or off the grid. Level 0 is the single filled cell. (Verified, `lab/rs/crop-counts`, every row.)

Carpet, `level = 0..5`, side `3^level`:

| `level` | ball in | ball cut | ball exposed | diamond in | diamond cut | diamond exposed |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 4 | 0 | 1 | 4 |
| 1 | 0 | 8 | 16 | 0 | 8 | 16 |
| 2 | 32 | 28 | 80 | 12 | 32 | 64 |
| 3 | 332 | 76 | 400 | 168 | 104 | 304 |
| 4 | 2908 | 204 | 2688 | 1596 | 320 | 1792 |
| 5 | 23900 | 580 | 20160 | 13560 | 968 | 12400 |

Sponge, `level = 0..4`, side `3^level`:

| `level` | ball in | ball cut | ball exposed | diamond in | diamond cut | diamond exposed |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 6 | 0 | 1 | 6 |
| 1 | 0 | 20 | 72 | 0 | 20 | 72 |
| 2 | 44 | 216 | 792 | 0 | 116 | 504 |
| 3 | 2320 | 1224 | 8400 | 132 | 476 | 2016 |
| 4 | 54800 | 7968 | 134448 | 4320 | 2612 | 18456 |

Two structural readings, both exact. The inscribed sphere never enters the sponge at level 1: `census` reads cells `[0, 26, 1]` - the one In cell is the empty centre, all 20 filled cells are Cut - so the keep-cut crop keeps everything and the strict crop keeps nothing (Verified, `mrly_math::shape` test and the `level = 1` row). And the inscribed octahedron holds no filled sponge cell fully inside through level 2: its deep interior is exactly where the sponge is empty. (Verified.)

## The radius sweep

At carpet level 4 (side 81) and sponge level 3 (side 27), radii `r = 1/24 .. 24/24`. The summary; every number is a printed line of `lab/rs/crop-counts`. (Verified.)

- Dead zones: the carpet crop is empty - In and Cut both zero - for `r < 1/6` under ball and diamond alike, the central hole's inradius in both norms, and the sponge diamond crop for `r < 1/3`; both radii are exact. The sponge ball's exact contact radius is `sqrt(2)/6 = 0.2357`, off the sweep grid: the level-3 cell `[13/27, 14/27] x [8/27, 9/27] x [8/27, 9/27]` is filled - digit triples `(1,0,0), (1,2,2), (1,2,2)`, each with at most one middle digit - and its nearest point to the centre is `(1/2, 1/3, 1/3)` at distance `sqrt(2)/6`, so the sweep reads empty at `r = 5/24 = 0.208` and first cuts at `r = 6/24`, the first grid radius past contact. (Proved by the witness cell; the sweep rows Verified.)
- First contact on the sweep grid: carpet ball at `r = 1/6` reads `in = 0, cut = 4`; sponge ball at `r = 6/24` reads `in = 0, cut = 60`; sponge diamond at `r = 1/3` reads `in = 0, cut = 12`.
- At the inscribed `r = 1/2` the sweep reproduces the level tables: carpet ball `2908 / 204`, carpet diamond `1596 / 320`, sponge ball `2320 / 1224`, sponge diamond `132 / 476`.
- Saturation: the carpet ball crop holds all `4096` filled cells with zero Cut from `r = 17/24`, the first sweep radius past the circumradius `sqrt(2)/2 = 0.7071`; the sponge ball from `r = 7/8`, past `sqrt(3)/2 = 0.8660`. The sponge diamond never saturates in the sweep - at `r = 1` it reads `in = 5356, cut = 1332` of 8000, since the cube's corners sit at `L1` distance `3/2`.
- The Cut column is not monotone in `r`: the carpet ball's runs `204, 184, 204, 144` across `r = 5/12 .. 13/24`. How this count fluctuates along `log r` is the [dimensions](../../dimensions/notes/dimensions.md) question below. (Verified for the values; nothing is claimed about the fluctuation.)

## Chapters

- [The circle count](circle-count.md): The corner disc count as one sequence with a log-periodic main term, the circle theorem with its exact integer defect, the measured exponents, the digits of the crossing shell, the 1923 bound on them, and the shell as a tree.
- [The transfer operator](transfer-operator.md): The transfer operator on the crossing shell, derived from the geometry, the product over scales and its decay, the resonance-tracking bound, the centre hole, and the routes that stop.

## The open lane: curved slices

The Cut column at fixed `r = 1/2` is a box count: `filled_cut` at level `level` counts the `3^-level` cells that meet both the circle and the pre-fractal at level `level`, so its growth exponent `log_3` of the ratio measures the dimension of the circle's trace on the carpet - an upper bound for the trace on the limit set, since a surviving cell need not meet the intersection itself. The printed ratios give `log_3(28/8) = 1.140`, then `0.909, 0.899, 0.951` for the carpet circle, and `2.166, 1.579, 1.705` for the sponge sphere. (Verified for the ratios; any limit is **Conjecture** - five levels decide nothing.)

The yardstick is the straight case. Furstenberg's slice conjecture, proved independently by [Shmerkin 2019](https://arxiv.org/abs/1609.07802) and [Wu 2019](https://arxiv.org/abs/1609.08053), concerns the intersection of a `xp`-invariant and a `xq`-invariant closed subset of the line with `p` and `q` multiplicatively independent - equivalently, irrational-slope slices of the product `A x B` - and bounds its box dimension by `max(dimension - 1, 0)`. That theorem is cited, not claimed, and the carpet is a single `x3`-invariant planar set, not such a product, so the theorem does not literally cover the carpet's straight lines either, let alone circles: for the carpet and the sponge, both the straight `dimension - 1` bound and any curved analogue are yardsticks by analogy, not covered cases. The comparison values are `dimension - 1 = 0.8928` for the carpet and `1.7268` for the sponge, and the printed exponents hover near both - suggestive and unproved. Whether the slice bound, or the generic-slice value from the Marstrand side, extends to curved slices - a circle on the carpet, a sphere on the sponge - appears uncharted for self-similar carpets: the arithmetic digit structure that drives the straight-line theory has no obvious action on a curve of no rational slope anywhere. That is the lane. (**Conjecture**, all of it.)

## The neighbours

- [slices](../../slices/notes/slices.md) and [cuts](../../cuts/notes/cuts.md) section these solids with single planes, one offset at a time, and get exact meshes and digit-scheduled gaskets. A sphere crop of the sponge sweeps all plane sections at once: the Cut shell at radius `r` meets tangent planes of every orientation, and the sweep in `r` runs through every offset, so the sphere's cut census aggregates the whole two-parameter family of plane sections the plane pages take one by one. The price is exactness - the plane pages get closed forms, the shell so far only counts.
- [dimensions](../../dimensions/notes/dimensions.md) predicts that a lattice set's counts oscillate along `log` of the scale with period `2*pi/ln 3`. `filled_in(r)` detrended by `r^dimension` along `log r`, and `filled_cut` along the level, are the same observable for the crop; the 24-point sweep above is too coarse to fold, and nothing is measured yet. (Open.)
- The trivial end, stated honestly: a grid-aligned polytope crop is digit counting, not geometry. A box whose walls sit on multiples of `3^-k` keeps exactly the cells whose coordinates lie in integer intervals at level `k`, and counting a digit-product set over an integer box is the restricted-digit machinery of `mrly_math::press` - the count factors along digit positions and the crop adds nothing. (Proved, read off the definitions.) The lane above starts precisely where this reduction stops: the ball is the simplest shape with no digit structure.

## Where the numbers live

`mrly_math::shape` carries the machinery - `Frac`, `classify`, `regions`, `named`, `crop`, `refine` with its 20-million-cell guard, `census` - with the corner oracle, the partition identities, the diamond closed form `2m(m-1)` and the sponge shell `[0, 26, 1]` pinned by its tests. `lab/rs/crop-counts` prints every table and sweep line above and asserts both identities on each; `lab/rs/circle-crop` prints every line of the circle count, its exponents and its bands, its transfer pass prints the shell's tree, the pattern operator with its certified Perron brackets and minorisation floors, and the gain ladder that turns `Psi` into a product over scales, and its derivation pass prints the frozen-slope alphabet, the fault rates of the derived pattern law against the shell, the centre rates against `1/3`, and the straight-line gain ladder that refutes the frozen gap.
