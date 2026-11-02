# stack-dilations

- The parity stack's layers read as a system of dilated functions in `L^2(0,1)`, and the Dirichlet series that is its symbol.
- The square wave `s(x) = (-1)^floor(x)` is the odd 2-periodic extension of the constant `1` on `(0,1)`, so the layers `s(nx)` are the dilates `phi(nx)` of one function and the Hedenmalm-Lindqvist-Seip theory applies verbatim.
- Its sine coefficients against `e_n(x) = sqrt 2 sin(n pi x)` are `a_n = 2 sqrt 2/(pi n)` at odd `n` and `0` at even `n`, so the symbol is `S(s) = sum_n a_n n^-s = (2 sqrt 2/pi) (1 - 2^(-1-s)) zeta(1 + s)`.
- `check_dilation_shift` verifies the one identity the whole reading rests on, `<s(nx), e_k> = a_(k/n)` when `n | k` and `0` otherwise, by exact piecewise integration, max deviation `1.16e-15` over `n <= 8` and `k <= 60`: the dilate carries the Dirichlet series `n^-s S(s)`.
- `check_gram_two_ways` prints the Gram matrix of the layers two ways at all 78 pairs `m <= n <= 12`: `grid_integral` integrates `s(mx) s(nx)` exactly in rationals on the `lcm(m,n)` grid, `symbol_entry` reads the symbol sum `sum_j a_(j n') a_(j m')`, and they agree entry for entry, reading `1/15` at `(3,5)`, `1/3` at `(3,9)` and `0` at `(2,3)`.
- `symbol_series` sums the symbol numerically to 20000 terms and lands within its own tail bound of the closed form at every pair, largest gap `2.03e-05`.
- `check_blocks` shows the Gram entry vanishes whenever `v_2(m) != v_2(n)`, and that the block at `v_2 = a` is the odd Gram `gcd(m,n)^2/(mn)` itself, at all 2730 ordered pairs to 64: the full dilation system's Gram is a countable direct sum of copies of the odd one, so the odd stack carries the whole spectrum.
- `check_determinant` re-derives the Smith determinant `prod over odd k of prod over p | k of (1 - p^-2)` over the first `K` odd scales at `K = 1..13`, the value over the 13 odd scales `n <= 25` being `11399736556781568/21994507608198125`, the number `lab/py/stack-levels` prints.
- `check_inverse` verifies that the Dirichlet inverse of `a_n/a_1` is `mu(n)/n` at odd `n` and `0` at even `n`, at every `n <= 400`, so `1/S(s) = 1/(a_1 (1 - 2^(-1-s)) zeta(1 + s))` and the Mobius square wave, coefficients `mu(n) a_n`, has symbol `a_1^2/S(s)`.
- `a_n` itself is not totally multiplicative, `a_1 = 2 sqrt 2/pi`; `a_n/a_1` and its inverse are, which is what HLS Corollaries 5.3 and 5.8 need, and both verdicts are invariant under scaling `phi` by `a_1`. The prime sum those corollaries test already reads `1.9854` at `p <= 100000` and diverges by Mertens.
- Not-a-Riesz-basis needs no multiplicativity at all: HLS Theorem 5.2 with Theorem 3.1 demands `S` bounded on `Re s > 0`, and `S(s) = a_1 (1 - 2^(-1-s)) zeta(1 + s)` blows up as `s -> 0+`.
- `spectrum` computes the extreme eigenvalues of the odd Gram over the first `K` odd scales with a symmetric eigensolver at `K = 25, 50, 100, 150, 200`, printing `lambda_max`, `lambda_min`, the condition number and two normalisations.
- The window reads `lambda_max` from `2.01467` to `2.47224` and `lambda_min` from `4.393e-01` to `3.570e-01`, condition number `4.586` to `6.926`, with `lambda_max/(log N)^2` falling from `0.1330` to `0.0689` and `lambda_max/(log log N)^2` falling from `1.0910` to `0.7717`.
- The ceiling those numbers sit under is not printed here and belongs to the literature: the spectral norm of a gcd matrix at exponent one over any `k` distinct integers is of order `(log log k)^2`, Lewko-Radziwill 2014 Theorem 2, which settles the exponent-one spectral norm that Gal's 1949 bound on the gcd sum itself had left open.

## RUN

`uv run python research/stack/lab/py/stack-dilations/stack_dilations.py`

Under a second, prints only, writes nothing.

## WITNESSES

- The dilation reading of `levels-dilations-nodes.md`, "The layers as a dilation system": the symbol `S(s)`, the Gram identity, the block splitting, the failed Riesz condition and the completeness.
- The eigenvalue window it prints is the finite half of the layer Gram spectrum named as the open object in "What a breakthrough would look like".
- It does not re-derive the covariance law of `lab/py/stack-levels`; it matches it, and extends it to the even scales the odd convention drops.
