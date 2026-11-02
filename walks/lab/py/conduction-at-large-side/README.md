# conduction-at-large-side

- Conduction on `bang dim 2, code 7` at odd side `N` and level `L`: the side-`N` tile voids a cell iff both digits are odd, the render is its `L`-th Kronecker power in the unit square, unit conductivity on the filled cells, potential 0 on the left edge and 1 on the right, every other boundary insulated.
- The network: one node per cell of the render refined `k` times per axis, conductance `2ab/(a+b)` between face neighbours of conductivities `a, b`, conductance `2a` to an electrode, so a solid square reads exactly 1; at `k = 1` it is the graph of the walk census. A symmetric render is solved on one quarter: the potential is odd about `x = 1/2` and even about `y = 1/2`. Sparse LU (SuperLU, COLAMD) in dim 2, Jacobi-preconditioned CG to `1e-11` in dim 3.
- `cell` computes the level-1 conductivity at infinite side on the unit square minus its top right quarter, between the left edge and the lower half of the right edge, at mesh `2M` per side, `M = 64..512`; Richardson with order `4/3`; then inclusion conductivities `z = 0.01, 1/9, 1/3, 3, 9, 100` by three-mesh Aitken against `sqrt((1 + 3z)/(3 + z))`, the products `sigma(z) sigma(1/z)`, and `z = 1e8` against `sqrt(3)`. It asserts the insulating value within twice its bar of `1/sqrt(3)`, every contrast within `6e-8` of the formula and every Keller product within `7e-8` of 1.
- `dual` multiplies the conductance with insulating voids by the one with voids of conductivity `1e8` at `(N, L) = (5, 1), (3, 2)`: Aitken on `k = 16, 32, 64`, the bar its shift from `k = 8, 16, 32`, and it asserts the product within its bar of 1.
- `side` prints `N (sigma_k(N, 1) - sigma_k(infinity, 1))` at `k = 4, 8, 16` for `N = 11, 21, 41, 81`, the infinite-side value at the matched mesh coming from `cell`, and `sigma(N, 1)/sigma(N, 2)` at `k = 1, 2, 4` for `N = 5..41`.
- `rate` prints `sigma(N, 1)/sigma(N, 2)` at `k = 2`, `N = 9..41`, and `N (sqrt(3) - rho(N, 1))`, and fits it by `r - a/N + b/N^(4/3)` with `r` free and by `sqrt(3) - a/N + b/N^q` at `q = 5/4, 4/3, 3/2`.
- `scale` prints the per-level ratio `rho(N, L) = sigma(N, L)/sigma(N, L + 1)` at `k = 1, 2` for `N = 3, 5, 7, 9, 11` and every level whose render has side at most 3200 at `k = 1` and 2700 at `k = 2`. It reads `rho(N)` as the deepest value at the mesh whose last level step is smaller, `k = 1` at `N = 3, 5, 7` and `k = 2` at `N = 9, 11`, with bar the larger of that step and the mesh gap at the deepest common level; it also prints the cut and short bounds `[beta_N, alpha_N]`, `d_w = log(fill rho)/log N`, the law `2 + log(3 sqrt(3)/4)/log N`, `(d_w - 2) log N` with the ends of its trap, `log(rho)/log(N^2/fill)` and its gap to `log 3/log(16/9)`, the drift from the first ratio to the deepest at `k = 1, 2`, and the table rows of the note. It asserts every ratio inside the bounds, which it prints as exact fractions.
- `dim3` computes the infinite-side cell of `bang dim 3, code 23`, the unit cube minus the four octants with at least two upper coordinates, at `M = 8, 16, 32, 64`, Aitken on two triples and Richardson on the last pair, and prints `log(1/sigma)/log 2`.

## RUN

- `uv run python research/walks/lab/py/conduction-at-large-side/conduction.py` runs every verb in about 32 s: `cell` 13 s, `dual` 0.5 s, `side` 3 s, `rate` 5 s, `scale` 6 s, `dim3` 5 s.
- One verb by name: `uv run python research/walks/lab/py/conduction-at-large-side/conduction.py scale`.
- Needs numpy and scipy; prints only, writes nothing.

## WITNESSES

- conduction.md, section "Conduction at large side": the cell value and its bar, the orders, the contrast and Keller readings, the duality products and bars, the level-1 boundary readings, the ratios `sigma(N, 1)/sigma(N, 2)` and the fits, the table, the drifts, the trap ends, the Archie readings and `-0.1412`, and the dim 3 value and exponent.
