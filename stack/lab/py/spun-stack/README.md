# spun-stack

- What the line stack does when every layer is rotated: which rotations keep it exact, what the exact spun stack lights, and how bright each node is.
- Layer `n` at angle `theta` draws the lattice `(1/n) R_theta Z^2`. Two layers share a node other than the origin iff `cos` and `sin` of their relative angle are both rational, so the whole-degree schedules are decided by Niven's theorem.
- The exact spun stack indexes layers by nonzero associate classes of `Z[i]`: layer `z` is `z^{-1} Z[i]`, scale `|z|` and rotation `-arg z` in one multiplication. Its lit set is the Gaussian Farey set, its brightness a Gauss circle count.
- Every rotation and every node is exact: `Fraction` coordinates, Gaussian Euclidean gcd, and the cyclotomic field `Q(zeta_360)` for the rationality of `cos` and `sin` at whole degrees.

## THE FIXED-INCREMENT COROLLARY

- Under the fixed increment, layer `k` sits at `k theta`, so layers `j` and `k` share an exact lattice iff `(j - k) theta` is a multiple of 90. **Proved**, three lines from the dead-spin theorem.
- Write `theta/90 = p/q` in lowest terms. Then `(j - k) p/q` is an integer iff `q` divides `j - k`, so the layers fall into exactly `q` angle classes modulo 90, sharing happens inside a class and never across, and the sharing-pair count is `Sum_classes C(size, 2)`. **Proved.**
- If `theta/90` is irrational, `(j - k) theta/90` is never a nonzero integer, so no two layers share and the stack is dead everywhere but the origin. **Proved.**
- So the schedule lights up exactly at the Farey fractions of a quarter turn, and `q` - the Farey denominator - is the number of eyes. **Verified** by the table below.

## WHAT IT PRINTS

- `rational_angle_degrees` reduces `2 cos(d deg) = zeta^d + zeta^-d` and `2 sin(d deg) = zeta^(90-d) - zeta^(90+d)` modulo the cyclotomic polynomial `Phi_360` and returns the degrees where both remainders are constant: `0, 90, 180, 270`, four of 360.
- `spot_check_degrees` re-decides twelve of those degrees through `sympy.minimal_polynomial` and agrees with the cyclotomic route on all twelve.
- `pythagorean_hits` checks that all 68 rational points of the unit circle with denominator at most 60 are of the form `w^2/N(w)` for a Gaussian `w` in the box of side 12.
- `dead_spin_pairs` runs the whole-degree schedule `layer k at angle 0, 2, 3, 5, 7, 11, ...` for `k = 1..30`: 5 of the 435 pairs share a node other than the origin, every one of them at a relative angle of exactly 90 degrees.
- `shared_witness` exhibits the shared lattice for each of the 5: `(5, 7)` with `(26, 97)`, `(6, 11)` with `(27, 101)`, `(7, 13)` with `(28, 103)`, `(8, 17)` with `(29, 107)`, `(9, 19)` with `(30, 109)`. The shared set is `R_alpha (1/g) Z^2` with `g = gcd(m, n)`, of density `g^2` per unit area.
- `unit_square_shares` counts that lattice in the open unit square with the origin excluded, the convention the script prints: `1, 9, 49, 1, 9` for the five, read at angles 7, 11, 13, 17 and 19 degrees. The count is not `g^2` in general - it depends on the angle, and `share_count_spread` gives the whole range over degrees 1 to 89: `2, 3, 4` at `g = 2`, `8, 9, 10` at `g = 3`, and `41, 47, 48, 49, 50` at `g = 7`.
- `increment_classes` runs the fixed-increment schedule on the 28 odd scales `1..55`, layer `k` at `k theta` degrees, and prints the angle classes and the sharing-pair count for each `theta`.

| `theta` | `theta/90` | classes `q` | class sizes | sharing pairs |
|---|---|---|---|---|
| 0 | `0/1` | 1 | 28 | 378 |
| 10 | `1/9` | 9 | 4 and eight 3 | 30 |
| 11.25 | `1/8` | 8 | four 4 and four 3 | 36 |
| 15 | `1/6` | 6 | four 5 and two 4 | 52 |
| 18 | `1/5` | 5 | three 6 and two 5 | 65 |
| 22.5 | `1/4` | 4 | four 7 | 84 |
| 30 | `1/3` | 3 | 10, 9, 9 | 117 |
| 36 | `2/5` | 5 | three 6 and two 5 | 65 |
| 45 | `1/2` | 2 | 14, 14 | 182 |
| 60 | `2/3` | 3 | 10, 9, 9 | 117 |
| 67.5 | `3/4` | 4 | four 7 | 84 |
| `90 (sqrt2 - 1)` | irrational | none | all 1 | 0 |

- Every row is confirmed by a second route in the same call: the pairwise exact test `(j - k) theta / 90 in Z` agrees with the class formula on all twelve, and for the eight whole-degree increments a third route, the file's own Niven-free `rational_cos_sin` test on the relative angle, agrees as well.
- `near_miss` sweeps the other 430 pairs over every nonzero lattice vector in the box of side 6 and reports the closest a rotated node comes to a coincidence: `0.003390`, at scales 17 and 25 and angles 53 and 89.
- Layers are the nonzero associate classes of `Z[i]`, the zero layer excluded and the four units counted once, so one layer per class; `literal_stack` builds the spun stack at norm bound 50 by exact stacking of all 40 layers and finds 672 nodes, equal to the Gaussian totient sum `sum Phi(d)` over the same 40 classes.
- `closed_brightness` reads each node from its reduced Gaussian denominator alone as `g(floor(N/N(d)))` with `g(t) = sum_j (floor(t/(4j+1)) - floor(t/(4j+3)))`: 672 comparisons against literal stacking, 0 mismatches, the origin at 40.
- `circle_classes_jacobi` against `circle_classes_direct` for `t = 0..400`: equal throughout; the first twelve values are `1, 2, 2, 3, 5, 5, 5, 6, 7, 9, 9, 9`.
- `totient_sum` gives 672, 10608 and 168088 nodes at norm bounds 50, 200 and 800, ratios to `N^2` of `0.268800`, `0.265200`, `0.262638` against `pi / (8 zeta(2) G) = 0.260635`.
- `base_depth_check` nests the layers `c^-k Z[i]` for `c = 1 + i` to depth 8 and `c = 2 + i` to depth 4, and confirms that brightness equals `depth + 1 - address` at every one of the 256 and 625 deepest nodes.
- `base_c_overlap` over the box of side 10: `1 + i` puts all 440 nonzero points in the coarser layer, `3/2 + i/2` puts 220, and `sqrt(2) e^i` puts none, its closest approach `0.0604355`.

## RUN

- `uv run python research/stack/lab/py/spun-stack/spun_stack.py`
- From the repo root. One core, about seven seconds.
- Domain is the full source domain: all 360 whole degrees, all 435 prime-schedule pairs, all 378 pairs of the 28-layer increment schedule at twelve increments, the spun stack at norm bound 50 node by node, totient sums to norm bound 800.
- Nothing is written to disk.

## WITNESSES

- The fixed-increment corollary: `increment_classes` gives the class count `q`, the class sizes and the sharing-pair count at all twelve increments, by three agreeing routes.
- The dead-spin theorem: `rational_angle_degrees` gives the four rational rotations, `dead_spin_pairs` and `near_miss` give the 5 sharing pairs and the `0.003390` floor under the other 430, `unit_square_shares` and `share_count_spread` the window convention and the angle dependence of the count.
- The exact spun stack: `literal_stack` against `closed_brightness` and `totient_sum`, 672 nodes and 0 brightness mismatches at norm bound 50.
- The Gaussian twin of `floor(N/b)`: `circle_classes_jacobi` against `circle_classes_direct` on `t = 0..400`.
- The complex-base corollary: `base_depth_check` and `base_c_overlap`.
