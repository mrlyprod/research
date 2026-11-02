# Digit Uniform Bound

- One bound on the digit transform of every one-missing-digit set in a base at once: it names the base and forgets the digit, so what it proves is a statement about `base` and not a list of bases.
- The target is the `l^1` exponent `alpha_1` of a one-missing-digit set against the threshold `1/4`, where the pair route turns a Mobius asymptotic on.

## THE DOMINATION

- The set is `F = {0..base-1}`, the digits of the code, less `{a0}`, `fill = base - 1`, and its normalised transform is `hat F(t) = (1/(base-1)) Sum_(a != a0) e(a t)`.
- The full residue sum is the Dirichlet kernel up to a phase: `Sum_(a<base) e(a t) = e(((base-1)/2) t) D_base(t)` with `D_base(t) = sin(base pi t)/sin(pi t)`.
- Subtracting the excluded term and taking absolute values, `|hat F(t)| = |e(((base-1)/2) t) D_base(t) - e(a0 t)|/(base-1) <= (|D_base(t)| + 1)/(base - 1) =: u_base(t)`, for every `t` and every `a0`.
- The right side names no digit, so `G(w) = sup_(cell w) |hat F| <= sup_(cell w) u_base` dominates all `base` excluded digits at once, and `|hat F| <= 1` gives the free sharpening `u_base -> min(1, u_base)`.
- The cost is exact: the triangle inequality throws away the phase `e((a0 - (base-1)/2) t)`, so `u_base` is attained only where that phase is antipodal to the kernel and the bound is strict everywhere else. The uniform bound is weaker than the per-digit one at every base, by construction.

## THE RUN DECOMPOSITION

- The shift sandwich reads `alpha_1 <= (1/N) log_base max_x Sigma_N(x)` with `Sigma_N(x) = Sum_(i < base^N) Prod_(j<N) |hat F(base^j (x + i/base^N))|`, and `Sigma_N <= Sigma_N^u`, the same sum with `u_base` in place of `|hat F|`.
- Expanding `Prod_j (|D_base(base^j t)| + 1)` over subsets `E` of `{0..N-1}` and telescoping each maximal run, `D_base(t) D_base(base t) ... D_base(base^(l-1) t) = D_(base^l)(t)`, gives the exact identity `(base-1)^N Sigma_N^u(x) = Sum_E Sum_(i<base^N) Prod_(runs (s,l) of E) |D_(base^l)(base^s (x + i/base^N))|`.
- Every term is a Lebesgue sum of a Dirichlet kernel over an arithmetic progression of points. Peeling the run of largest `s` first, whose factor depends on the fewest low digits of `i`, bounds each term by `base^N Prod_r lambda_(l_r)` with `lambda_l = L_(base^l)/base^l` and `L_M = max_x Sum_(j<M) |D_M(x + j/M)|`.
- The peel is checked term by term: every one of the `2^N - 1` run terms sits under `base^N Prod_r lambda_(l_r)` and `max_x Sigma_N^u` sits under `a_N (base/(base-1))^N`, at `base 3, 4` to `N = 4` and `base 5, 7, 10` to `N = 3`, slack about `12%`.
- Summing over `E` by first letter gives `a_N = a_(N-1) + Sum_(l<N) lambda_l a_(N-1-l) + lambda_N`, whose growth root solves `z = 1 + Sum_(l>=1) lambda_l z^(-l)`, and `alpha_1 <= log_base(z base/(base-1))`.

## THE LEBESGUE CONSTANT

- Writing the `M` points as distances `d_j` to the nearest integer and `|sin(M pi x)| = |sin(M pi d_j)|`, the sum is `s(theta) Sum_j 1/sin(pi d_j)` at offset `theta` in `[0, 1/2]`, with pairs at `(m + theta)/M` and `(m + 1 - theta)/M`.
- With `mu = pi(m + 1/2)/M` and `delta = pi(1/2 - theta)/M`, the pair identity `s = cos(M delta)` and `csc(mu - delta) + csc(mu + delta) = 2 sin mu cos delta/(sin^2 mu - sin^2 delta)` turn the pair into `2 csc mu` times `X = cos(M delta) cos delta sin^2 mu/(sin^2 mu - sin^2 delta)`.
- `X <= 1` is `sin^2 mu (1 - cos(M delta) cos delta) >= sin^2 delta`, whose left side grows with `mu`, so it is enough at `m = 0`, `u = mu_0 = pi/(2M) <= pi/4`, `delta = ru` with `0 < r < 1`, `M delta = pi r/2`; the half-angle identities `1 - cos(pi r/2) cos(ru) = 2 sin^2(ru/2) + 2 cos(ru) sin^2(pi r/4)` and `sin^2(ru) = 4 sin^2(ru/2) cos^2(ru/2)` turn it into `sin^2 u cos(ru) sin^2(pi r/4) >= sin^2(ru/2) (cos^2 u + cos(ru))`.
- `x cot x` decreases on `(0, pi)`, so `sin(pi r/4)/sin(ru/2)` decreases in `r` and is at least its value `1/(sqrt 2 sin(u/2))` at `r = 1`; hence the left side is at least `sin^2(ru/2) sin^2 u cos(ru)/(2 sin^2(u/2)) = sin^2(ru/2)(1 + cos u) cos(ru)`, which exceeds the right side by `sin^2(ru/2) cos u (cos(ru) - cos u) >= 0`. So every pair is largest at `theta = 1/2`, every `m` at once.
- The odd leftover at `pi/2 - delta` contributes `cos(M delta)/cos delta <= 1`, again largest at `theta = 1/2`.
- At `theta = 1/2` the sum is `2 Sum_(j < floor(M/2)) csc(x_j) + (M mod 2)`, `x_j = (2j+1) pi/(2M)`; split `csc = 1/x + g` with `g` nonnegative, increasing and convex on `(0, pi)`, `Int_0^(pi/2) g = log(4/pi)`, and use the harmonic bounds `H_N <= log N + gamma + 1/(2N)` and `H_N >= log(N + 1/2) + gamma`, both from the convexity of `1/x`.
- Even `M = 2n`: the `x_j`, `j < n`, are the midpoints of the `n` intervals of length `pi/M` partitioning `(0, pi/2]`, so `Sum g(x_j) <= (M/pi) log(4/pi)`, while `Sum 1/x_j = (2M/pi)(H_(2n) - H_n/2) <= (2M/pi)((1/2)(log M + log 2 + gamma) + 1/(2M))`; doubling, `L_M <= (2M/pi)(log M + gamma + log(8/pi)) + 2/pi`.
- Odd `M = 2n + 1`, `n >= 1`: `L_M = 2 Sum_(j < n) csc(x_j) + 1`, the `x_j` now the midpoints of `n` intervals of length `pi/M` inside `(0, pi/2)`, so `Sum g(x_j) <= (M/pi) log(4/pi)` as before. With `2n = M - 1` and `n + 1/2 = M/2` the harmonic bounds give `H_(2n) - H_n/2 <= log(M-1) + gamma + 1/(2(M-1)) - (1/2)(log(M/2) + gamma) <= (1/2)(log M + log 2 + gamma) - 1/M + 1/(2(M-1))`, using `log(M-1) <= log M - 1/M`. Multiplying by `4M/pi`, the two small terms give `-4/pi + 2M/(pi(M-1)) = -2/pi + 2/(pi(M-1))`, so with the `g` part and the leftover `1`, `L_M <= (2M/pi)(log M + gamma + log(8/pi)) + 1 - 2/pi + 2/(pi(M-1))`. For `M >= 5` the constant is at most `1 - 3/(2 pi) < 0.53 < 2/pi`; at `M = 3`, `L_3 = 2 csc(pi/6) + 1 = 5` against `3((2/pi) log 3 + gamma') + 2/pi > 5.62`.
- Hence `L_M <= M((2/pi) log M + gamma' ) + 2/pi` with `gamma' = (2/pi)(gamma + log(8/pi)) = 0.96252282676`, printed and used rounded up at `0.9625229`, so `lambda_l <= (2/pi) l log base + gamma' + (2/pi) base^(-l)` at every `l`, and the coarser `lambda_l <= (2/pi) l log base + c_0` with `c_0 = 0.97` needs `base^l >= 86`, since `0.97 - gamma' = 0.0074771` has to cover `(2/pi)/base^l`. The scan confirms the shape: `L_M/M - (2/pi) log M` falls from `0.965217` at `M = 4` to `0.962523`, and the maximising offset is `1/2` at every `M` tested.

## THE THRESHOLD

- With the exact hypothesis `lambda_l <= c_1 l log base + gamma' + c_1 base^(-l)`, `c_1 = 2/pi`, the root equation is `(z-1)^3 = c_1 (log base) z + gamma'(z - 1) + c_1 (z-1)^2/(base z - 1)`, and `alpha_1 < 1/4` follows from `z < base^(1/4)(1 - 1/base)`. The coarser `c_0 = 0.97` form drops the last term and is the same statement wherever `base^l >= 86`.
- Certified at 120 bits at `w = base^(1/4)(1 - 1/base)` for every `base` in `[125, 3000)`, and the chain fails at `base 124`; the certificate first checks that the constant it uses sits above the true `gamma'`, which the earlier printed `0.9625153` did not.
- Above `base 211` the closed-form cap `z <= 1 + sqrt(2 c_1 log base + c_0)`, which follows from `z >= 2` and `z <= 2(z-1)`, is below `base^(1/4)(1 - 1/base)` and both sides are monotone, so the tail needs no evaluation at all.
- The threshold is `base_u = 125`: `alpha_1 < 0.249980` at `base 125`, and it falls from there at every excluded digit, with no per-base spectral computation anywhere in the chain. The `c_0 = 0.97` form of the same chain gives `126`, one base worse, so the sharper constant strengthens that statement and does not contradict it.

## ANY NUMBER OF EXCLUDED DIGITS

- Nothing in the chain used `m = 1`. With `E` any set of `m` excluded digits, `|Sum_(e in E) e(e t)| <= m` gives `|hat F(t)| <= min(1, (|D_base(t)| + m)/(base - m)) =: u_(base,m)(t)`, again naming no digit, and the free cap `|hat F| <= 1` still applies.
- The run decomposition gains one weight. Expanding `Prod_j (|D_base(base^j t)| + m)` over subsets `E` of `{0..N-1}` carries `m^(N - |E|)` on the positions outside `E`, and each maximal run of `E` telescopes to the same Dirichlet kernel, so `(base-m)^N Sigma_N^u(x) = Sum_E m^(N - |E|) Sum_(i<base^N) Prod_(runs (s,l) of E) |D_(base^l)(base^s(x + i/base^N))|`, the identity checked term for term against the direct product.
- The peel is unchanged, since it bounds one run term at a time; summing over `E` by first letter gives `a_N = m a_(N-1) + m Sum_(l<N) lambda_l a_(N-1-l) + lambda_N` and `max_x Sigma_N^u <= a_N (base/(base-m))^N`.
- The growth root solves `z = m + m Sum_(l>=1) lambda_l z^(-l)`, so with `lambda_l <= c_1 l log base + gamma' + c_1 base^(-l)` the certificate is `(z - m)(z - 1)^2 <= m(c_1 (log base) z + gamma'(z - 1) + c_1 (z-1)^2/(base z - 1))` and `alpha_1 < e` follows from `z < base^e (1 - m/base)`. At `m = 1` this is the cubic of the section above, so the general form is a check on the special one and not a rewrite of it.
- The thresholds, certified at 120 bits on `[base_u, 3000)` with the chain failing at `base_u - 1`: against `1/4` it gives `base_u(1) = 125`, `base_u(2) = 649`, `base_u(3) = 1873`; against the weaker bar `1/3` it gives `32`, `105`, `230`. So every `base >= 649` clears `alpha_1 < 1/4` at every excluded pair, with no per-pair spectral computation anywhere. The `1/3` rung at `m = 1` lands at `base 32`, where `base^l = 32` is below `86`, so it is exactly the rung the coarser `c_0` form may not be used on.
- The cost of uniformity grows with `m` because the phase thrown away is `m` characters rather than one: the majorant loses `m/(base-m)` where the truth loses the interference between them, and the two-digit numeric floor is far below `650` (`lab/py/digit-transform-norms`, verb `pairfail`).

## THE WINDOW CEILING

- The same majorant run through the window machine measures how far the uniform idea can reach at all: `U(w) = min(1, (min(base, S_1(w) S_2(w)) + 1)/(base-1))` on the cell `[w/base^n, (w+1)/base^n)`, the transfer matrix `(M y)(v) = Sum_(c<base) U(v base + c) y((v base + c) mod base^(n-1))`, and `alpha_1 <= log_base rho`.
- The supremum is enclosed in closed form, never sampled. `|sin(pi t)|` is concave on `[0,1]`, so `sup_cell 1/|sin(pi t)| = 1/sin(pi delta(w))` exactly with `delta(w) = min(w, base^n - 1 - w)/base^n`; and `sup_cell |sin(base pi t)|` is `1` when the image interval `[w, w+1]/base^(n-1)` contains a half-integer and `max` of the two endpoint values otherwise. Both are monotone closed forms, so there is no sub-scan and no Lipschitz slack.
- The ceiling: the uniform window bound clears `1/4` first at `base 75` and at every base above it, and the bound moves by less than `0.001` from three window digits to four (`0.340029 -> 0.339085` at `base 21`, `0.297650 -> 0.297296` at `base 34`), so no refinement of the window rescues a base near `74`.
- So `base 75` is the floor of the digit-uniform route as a numeric fact, `base_u = 126` is the floor of what the route proves uniformly in `base`, and the difference is slack in the run decomposition, not in the majorant.

## THE CONSISTENCY CHECK

- Where both apply the uniform bound must be weaker than the per-digit ladder, and it is at every base tested, never stronger: `0.461532` against `0.323432` at base 9 missing `0`, `0.440960` against `0.350684` at base 10 missing `5`, `0.344129` against `0.283414` at base 20 missing `6`, `0.339085` against `0.250088` at base 21 missing `0`, `0.299529` against `0.253000` at base 33 missing `15`, `0.297296` against `0.249371` at base 34 missing `16`.
- The gap closes as the base grows, `0.1381` at `base 9` down to `0.0479` at `base 34`, which is the discarded phase costing `1/(base-1)`.

## THE CERTIFICATE

- The threshold inequality is `mpmath.iv` at 120 bits: `base^(1/4)`, `log base` and every product and difference are intervals, and a row passes only when the lower endpoint of the margin is positive.
- The window rows are floats with every operation pushed outward by `nextafter` and every sine inflated by `2^-45` in the unsafe direction, the Perron root bounded above by Collatz-Wielandt, `M y <= mu y` componentwise giving `rho <= mu`, with the float Perron vector as the test vector and a positive floor on it.
- The power iteration accepts a root only after fifty consecutive relative moves under `1e-13` and at least three hundred steps; a shorter test admits the transient cluster of near-zero cells, whose row sums are `base`, and returns a root too large by a factor of sixteen.
- Every printed exponent is `log_base mu` rounded up and asserted against `mu < base^e` before it prints, and the certified `L_M` rows are the exact maximum over the offset, not a fit.

## RUN

- `uv run python research/coprime/lab/py/digit-uniform-bound/ubound.py check` in about 25 seconds: the domination on nine bases and every distinct digit, the run identity against the direct product at five base-level pairs, the Lebesgue scan, the threshold and its certificate, and the consistency rows.
- `domination`, `identity`, `peel`, `lebesgue`, `threshold`, `window`, `least`, `compare` run the pieces; `least` is 6 seconds and `window` is 5.
- Prints only, writes nothing; every row that has a target asserts it, and the run ends by raising if any row is off.

## WITNESSES

- The domination `|hat F| <= (|D_base| + 1)/(base-1)` for every excluded digit, with the digit-slack column `1.000000` at `base 3` falling to `0.010045` at `base 200`.
- The exact run decomposition of `Sigma_N^u` into `2^N` Lebesgue sums of Dirichlet kernels at moduli `base^l`, and the peel bounding each of them.
- `L_M <= M((2/pi) log M + 0.9625229) + 2/pi`, with the maximum at the half-offset and the scanned constant falling from `0.965217` to `0.962523`.
- The threshold `base_u = 125`, certified at 120 bits on `[125, 3000)`, the chain failing at `124`; the `c_0 = 0.97` form of the same chain gives `126`.
- The window ceiling `base 75` at three digits, and the consistency rows against the per-digit ladder.
