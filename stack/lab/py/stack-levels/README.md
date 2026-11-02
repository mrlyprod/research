# stack-levels

- The correlation law of two stacked layers when the layers are different designs, or the same design at different levels, computed exactly.
- The general law it runs on: for 1-periodic `f, g` with Fourier coefficients `hat f(k), hat g(k)`, the mean of `f(mx) g(nx)` over the unit interval is `sum_j hat f(j n') hat g(-j m')` with `g = gcd(m,n)`, `m' = m/g`, `n' = n/g`, and the covariance is the same sum over `j != 0`.
- `check_parity` recovers the landed parity law from that sum by exact rational integration on the lcm grid: the master integral `gcd(m,n)^2/(mn)` at all 820 pairs to 40, reading `1/15, 1/3, 1/3` at `(3,5), (3,9), (5,15)`, and the layer covariance `(g^2-1)/(4mn)`, zero at all 159 coprime odd pairs.
- The same block prints the price of the tree's half-period convention: the 1-periodic odd strip `p(x) = 1` iff `floor(2x)` is odd carries `g^2/(4mn)`, nonzero at 159 of the 490 coprime pairs to 40, so coprime independence belongs to the odd sub-stack and not to the design.
- `check_carpet` computes the covariance of `f_level(mx)` and `f_level(nx)` for the base-3 Sierpinski shadow `f_level(x) = 1` iff none of the first `level` base-3 digits is 1, exactly, at level 1, 2, 3 and all 820 pairs to 40, by cell counting on the `lcm(m,n) 3^level` grid.
- It checks three laws there: the covariance depends only on `m/g, n/g`; the product `Cov * m' n'` depends only on `(m' mod 3^level, n' mod 3^level)`; and the covariance is zero exactly when `3^level` divides `m'` or `n'`.
- At level 1 it checks the closed form `Cov = (2/9) chi(m') chi(n')/(m'n')` with `chi` the nontrivial character mod 3, reading `-1/9` at `(1,2)`, `1/45` at `(2,5)` and `0` at `(1,3)`, on all 1600 ordered pairs.
- It prints the kernels `G_2` on the units mod 9 and the row `G_3(1,b)` on the units mod 27, and the non-separability witness `G(1,1)G(4,4) - G(1,4)^2 = 76/729`.
- `carpet_2d` cell-counts the 2D field `f_level(mu) f_level(mv)` against `f_level(nu) f_level(nv)` on the squared grid and confirms `Cov_2D = Cov_1D (M + (2/3)^(2 level))`, so the two zero sets coincide.
- `check_levels` verifies the level split `f_level(x) = prod_(i < level) f_1(3^i x)` as masks at level 1..4, the nesting `f_level(nx) <= f_(level-1)(3nx)`, and the cross-level zero law `Cov(f_level(mx), f_level'(nx)) = 0` exactly when `3^level | n'` or `3^level' | m'`, on all 14400 ordered triples with `level, level' <= 3` and `m, n <= 40`.
- `check_designs` verifies two closed forms at level 1: between two base-3 one-digit-removed shadows, `Cov = c/(27 m'n')` with `c` in `-6, -3, 3, 6` and zero exactly when 3 divides `m'n'`; between the base-2 strip and a base-3 shadow, `Cov = c/(18 m'n')` with `c` in `-3, 0, 3`. The strip against the Sierpinski shadow is exactly zero at all 1600 pairs.
- `check_gram` computes `det[gcd(m,n)^2/(mn)]` over odd `m, n <= 2K+1` by exact elimination and matches `prod J_2(k)/k^2` over the same odds at `K = 1..12`, `J_2` the Jordan totient; the value at `K = 12` is `11399736556781568/21994507608198125`.

## RUN

`uv run python research/stack/lab/py/stack-levels/stack_levels.py`

Two seconds, prints only, writes nothing.

## WITNESSES

- The general law and its parity specialisation, against the landed `gcd(m,n)^2/(mn)`: `check_parity`.
- The gcd law does not survive to base 3: the carpet covariance is a character-twisted gcd law at level 1 and a mod-`3^level` kernel above it, and coprime layers are correlated: `check_carpet`, `carpet_law_level1`, `carpet_kernel`.
- Levels are not new layers: `check_levels`, `level_product_mask`.
- Which design pairs are coprime-independent: `check_designs`, `design_law_33`, `design_law_23`.
- The layers of the parity stack are linearly independent: `check_gram`, `jordan2`, `exact_det`.
