# Burnol Residue

- The residue of `K(s) = sum n^(-s)` over the positive integers whose base-3 digits lie in `{0,1}` at the off-real poles `s_(0,k) = log_3 2 + 2 pi i k/log 3`, `k = 0..10`, enclosed in a complex rectangle by interval arithmetic (`mpmath.iv`, 128 bits) with every truncation bounded by a proved tail; the same engine on base 3 `{0,2}` (the residues scale by `2^(-s_(0,k))`) and on the full digit set `{0,1,2}`, where `K` is the Riemann zeta function and the off-real residues are exactly zero.
- The engine is the level recursion `Lev_(j+1)(w) = base^(-w) sum_(l >= 0) (-1)^l ((w)_l/l!) base^(-l) gamma_l Lev_j(w+l)` on the vector of shifts `w = s_(0,k) + i`, `i <= I`, run to `level`; the `l`-truncation is bounded by a hypergeometric tail and the levels past `level` by a geometric one; `R_k = lim_j Lev_j(s_(0,k)) = lambda_(0,k) log base` by Burnol's Proposition 5.1.
- A second certified enclosure by the functional-equation route: `R_k = 1 + sum_(m >= 1) (-1)^m ((s)_m/m!) 2^(-1) 3^(-m) K(s+m)` with each `K(s+m)` summed directly in intervals to a second `level` plus its tail, the `m`-series cut at `M` plus its tail; wider, and it must meet the first.
- Two controls in floats: the `k`-th Fourier coefficient of the log-periodic profile `Phi(u) = lim_j A(3^(u+j))/2^(u+j)` from exact counts `A(x)` on a grid, through `Res_(s_(0,k)) K = s_(0,k) c_k`, and Burnol's level-sum limit by direct enumeration.
- The column: `lambda_(m,k)` for `m = 1..3` by Burnol's Proposition 7.1 recurrence, checked at `m = 1, 2` against the closed forms `(s-1)/4` and `(s-1)(s-2)/64` times `lambda_(0,k)`.
- Safe rounding: every endpoint is converted to an exact fraction before printing, lower endpoints floored, upper endpoints ceiled.

## RUN

- `uv run python research/dimensions/lab/py/burnol-residue/burnol_residue.py`
- `uv run python research/dimensions/lab/py/burnol-residue/burnol_residue.py --band` adds `k = 2..10`.
- Prints only, writes nothing; the default run takes about 8 seconds, the band about 55.

## WITNESSES

- string-to-composition.md the arithmetic pole, certified, "The string above lives in [0,1]": Burnol's Proposition 5.1 limit formula and the statement that the paper does not study the off-real residues further.
- string-to-composition.md the arithmetic pole, certified, "Write G(u) = A(3^u)/2^u": the profile lemma: `Phi(u) = 2^(1-u)` on `[1 - s_0, 1]`, `Phi(1 - s_0) = 2^(s_0) = 1.548562`.
- string-to-composition.md the arithmetic pole, certified, "For Re s > s_0": the dictionary `Res_(s_(0,k)) K = s_(0,k) c_k`.
- string-to-composition.md the arithmetic pole, certified, "Let Lev_j(w) = sum n^(-w)": `lambda_(0,1)` in `[0.231891517689918, 0.231891517689919] + i [-0.501067414481069, -0.501067414481068]`, width `3.4e-18`, level tail `2.4e-19`, distance from zero at least `0.552125193`, `I = 66`, `level 40`.
- string-to-composition.md the arithmetic pole, certified, "The same certificate at k = 0..10" and its table: the band `k = 0..10`, `I` from 52 to 164, widths from `2e-18` to `4e-13`.
- string-to-composition.md the arithmetic pole, certified, "Two independent floating-point computations": the level-sum control `0.231891518 - 0.501067414 i` at distance `1.9e-10` against its bound `7.5e-10`, the Fourier control `c_1 = -0.082138842 - 0.049607499 i` at distance `7.2e-08`, the second enclosure `[0.23189098, 0.23189280] + i [-0.50106814, -0.50106632]`, the `{0,2}` enclosure `[0.135292875345483, 0.135292875345484] + i [0.329874091244466, 0.329874091244467]`, the zeta boxes around `1` and `0`.
- string-to-composition.md the arithmetic pole, certified, "Burnol's Proposition 7.1 recurrence carried down": the column `lambda_(1,1)`, `lambda_(2,1)`, `lambda_(3,1)` and the closed forms with `[t^1] = -1/4`, `[t^2] = 1/64`.
