# eisenstein-stack

- The hexagonal twin of `lab/py/spun-stack`: the exact spun stack whose layers are the associate classes of the Eisenstein integers `Z[omega]`, `omega = e^(2 pi i/3)`, norm `N(a + b omega) = a^2 - a b + b^2`, six units, unique factorisation.
- Layer `z` is the lattice `z^-1 Z[omega]`, scale `|z|` and rotation `-arg z` in one multiplication, one layer per nonzero associate class with the zero class excluded.
- Coordinates are taken in the basis `1, omega`, so the fundamental domain of `C/Z[omega]` is the unit square of that chart and the picture is its image, the fundamental parallelogram. Every node is a pair of `Fraction`s in that chart; the arithmetic is exact throughout, with nearest-integer Euclidean gcd in `Z[omega]` and the cyclotomic field `Q(zeta_360)` for the whole-degree rotations.
- The constant the base-3 page carries, `zeta_K(2) = zeta(2) L(2, chi_-3)`, is the constant in this stack's node count.

## WHAT IT PRINTS

- `field_rotation_degrees` reduces `zeta^d + zeta^-d` and `(zeta^d - zeta^-d)(zeta^120 - zeta^240)` modulo `Phi_360` and returns the degrees where both remainders are constant, that is where `cos` is rational and `sin` is a rational multiple of `sqrt 3`: `0, 60, 120, 180, 240, 300`, six of 360.
- `spot_check_degrees` re-decides twelve of those degrees through `sympy.minimal_polynomial` on `cos(d deg)` and `sin(d deg)/sqrt 3` and agrees with the cyclotomic route on all twelve; 90 and 270 are the instructive failures, rational cosine and irrational `sin/sqrt 3`.
- `rotation_hits` checks that all 58 solutions of `p^2 + 3 q^2 = r^2` with `r` at most 60, read as the pair `(cos, sin/sqrt 3) = (p/r, q/r)`, are of the form `w/conj(w) = w^2/N(w)` for an Eisenstein `w` in the box of side 20.
- `hex_classes_closed` against `hex_classes_direct` for `t = 0..400`: equal throughout; the first twelve values of the hexagonal circle count are `1, 1, 2, 3, 3, 3, 5, 5, 6, 6, 6, 7`.
- `literal_stack` builds the stack at norm bound 50 by exact stacking of all 31 layers and finds 630 nodes, equal to the Eisenstein totient sum `sum Phi(d)` over the same 31 classes from `totient_sum`.
- `closed_brightness` reads each node from its reduced Eisenstein denominator alone, computed by `reduced_denominator`, as `h(floor(N/N(d)))` with `h(t) = sum_j (floor(t/(3j+1)) - floor(t/(3j+2)))`: 630 comparisons against literal stacking, 0 mismatches, the origin at 31, which is also the maximum.
- `totient_sum` gives 630, 9606, 151020, 337026, 945486 and 2419950 nodes at norm bounds 50, 200, 800, 1200, 2000 and 3200, ratios to `N^2` of `0.252000`, `0.240150`, `0.235969`, `0.234046`, `0.236372` and `0.236323`; `main` prints each deviation from the limit scaled by `N/log N`, `+0.214, +0.186, +0.090, -0.198, +0.304, +0.438`, of both signs; read as bounded they were the first guess `N log N` for the error, Refuted by `lab/rs/node-error`, where the scale climbs to a root mean square of `0.8018` and the error reads `N^(5/4)`.
- `main` prints the constant's factors at 30 working digits through the Hurwitz zeta form `L(2, chi_-3) = (zeta(2, 1/3) - zeta(2, 2/3))/9`: `zeta(2) = 1.644934066848226`, `L(2, chi_-3) = 0.781302412896486`, `zeta_K(2) = 1.285190955484149`, the residue `pi/(3 sqrt 3) = 0.604599788078073`, and the limit `pi/(6 sqrt 3 zeta(2) L(2, chi_-3)) = 0.235217881630015`.
- `csl_zeta_ratio` expands `(1 + 3^-s)^-1 zeta_K(s)/zeta(2s)` by Dirichlet division to bound 100 and `csl_euler_product` expands `prod_{p = 1 (3)} (1 + p^-s)/(1 - p^-s)` over the same bound: equal, with nonzero coefficients `1` at 1, `2` at `7, 13, 19, 31, 37, 43, 49, 61, 67, 73, 79, 97` and `4` at 91.
- `draw` stacks the 774 lit nodes at norm bound 60, brightest 35 at the origin, into `eisenstein-stack.png`, dots sized by the square root of brightness inside the fundamental parallelogram.

## RUN

- `uv run python research/stack/lab/py/eisenstein-stack/eisenstein_stack.py`
- From the repo root. One core, about three seconds.
- Domain is the full source domain: all 360 whole degrees, the circle count on `t = 0..400`, the stack at norm bound 50 node by node, totient sums to norm bound 3200, the coincidence series to bound 100.
- Writes `eisenstein-stack.png` beside itself and nothing else.

## WITNESSES

- The Eisenstein address: `literal_stack` against `closed_brightness` and `totient_sum`, 630 nodes and 0 brightness mismatches at norm bound 50.
- The hexagonal twin of `floor(N/b)`: `hex_classes_closed` against `hex_classes_direct` on `t = 0..400`.
- The hexagonal rational rotations: `field_rotation_degrees` and `spot_check_degrees` give the six whole degrees, `rotation_hits` the 58 rational rotations as `w/conj(w)`.
- The constant: `totient_sum` at six norm bounds against the printed `pi/(6 sqrt 3 zeta(2) L(2, chi_-3))`, with the scaled deviation.
- The coincidence census: `csl_zeta_ratio` against `csl_euler_product`.
