# layer-spectra

- The spectrum view of [planar-block](../../../notes/planar-block.md), section Four-term progressions, bullet Graphs and spectra: a function `R: F_7^m -> F_7` is good when its graph in `F_7^(m + 1)` holds no four-term progression, that is when no progression `w, w + c, w + 2c, w + 3c`, `c != 0`, is sent to four values with both second differences `D_0 = R_0 - 2R_1 + R_2`, `D_1 = R_1 - 2R_2 + R_3` zero; the spectrum is the set of pairs `(D_0, D_1)` that occur. A sum of layers on complementary subspaces is good iff no choice of a spectrum point or `0` from each layer, not all `0`, sums to `0` in `Z_7^2`.
- No dependency: one file, `std` only, one thread.

## THE VERBS

- `univariate`: every function `F_7 -> F_7`, its spectrum, its degree by interpolation (`c_k = -sum h(x) x^(6 - k)`), whether it is good, whether its third difference `D_1 - D_0` never vanishes; the quintics `a (x + t)^5 + q(x)`, `q` quadratic, are built and asserted to be exactly the degree-`5` functions of the last kind; the distinct spectra are listed with their sizes.
- `layers`: depth-first search over the distinct good univariate spectra for the deepest zero-sum-free stack, repetition allowed, the reachable sums kept as a `49`-bit mask; a complete search within its node budget.
- `lift m`: the rank over `F_7` of the constraints `sum_x h(w + xc) x^j = 0`, `j = 0, 1, 2`, over the lines of `F_7^m`, which say that `h` has degree at most `3` on every line; the solution space is compared with the cubic polynomials, dimension `binom(m + 3, 3)`.
- `monomial`: the functions `Tr(a v^k)` on `F_49 = F_7(sqrt 3)` and `F_343 = F_7(cbrt 2)`, `a` over coset representatives of the `k`-th powers (the others are reparametrisations `v -> lambda v`), their spectra when good, and the zero-sum-free stacks `F_49 + F_343`, `F_343 + F_343`, `F_49 + F_49`, `F_49 + F_49 + F_49`.
- `plane`: the univariate functions with nonzero third differences, up to affine maps of `F_7`, scaling and the addition of a quadratic, each placed on the axis of `F_7^2` and completed, if it can be, to a function with nonzero third differences on every progression of the plane, by a depth-first search with forward checking on the `1176` four-point windows, `150000000` nodes a seed.

## RUN

- `cargo run --release --manifest-path research/progressions/lab/rs/layer-spectra/Cargo.toml -- all` from the repo root; `-- univariate`, `-- layers`, `-- lift`, `-- monomial`, `-- plane` run one part; all but `plane` finish in `2` seconds.
- Prints only, writes nothing; every check asserts.

## WITNESSES

- `568302` of the `823543` functions are good, `294, 2058, 8232, 69972, 487746` at degrees `2..6`; `45276` have a nonzero third difference on every progression, `2058, 14406, 28812` at degrees `3, 5, 6`, none at `4`; `323` distinct good spectra, from `3` points (the quadratics, `{(t, t): t in 2a QR}`) to `31`.
- `layers`: the deepest zero-sum-free stack of univariate spectra has `4` layers, two quadratic spectra and two of `14` points; stacks by depth `323, 341, 74, 6, 0, 0`.
- `lift`: dimensions `10` at `m = 2` (`168` constraints, rank `39`) and `20` at `m = 3` (`8379` constraints, rank `323`), the dimensions of the cubic polynomials.
- `monomial`: `F_49`, `240` functions, `28` good, spectra of `6, 42, 48` points for `9, 7, 12`; `F_343`, `2331` functions, `98` good, spectra of `42` and `48` points for `49` each, the norm `k = 57` among the `42`; stacks `441` for `F_49 + F_343`, for instance `k = 2` with `k = 57`, and `0` for the other three shapes.
- `plane`: `3` seeds, one cubic `(0, 0, 0, 1, 4, 3, 6)`, one quintic `(0, 0, 0, 1, 0, 2, 4)`, one sextic `(0, 0, 0, 1, 0, 2, 5)`; the cubic extends in `719308` nodes, `0.9` seconds; the quintic and the sextic are undecided at the cap of `150000000` nodes, about `195` seconds each; `389` seconds in all.
