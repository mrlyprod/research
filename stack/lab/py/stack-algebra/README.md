# stack-algebra

- Stacking a stack is Dirichlet convolution on the layer weights, and this study checks that by literal stacking, then reads off the closed form of every selected stack.
- The `w`-weighted line stack over scales `1..N` gives node `a/b`, denominator `b`, the brightness `sum_{k <= N/b} w(kb)`; the `u`-stack of the `v`-stack draws, at scale `kn`, the weight `u(k) v(n)`.
- Under the hyperbolic cut `kn <= N` the composite scale weight is exactly `(u * v)(m)` for every `m <= N`, so the two pictures agree node for node.
- Checked at `N = 60` on `u = v = 1`, `u = 1, v = mu`, `u = mu, v = mu` and `u = 1, v = n^-1`: 1102, 1102, 974 and 1102 nodes, zero mismatches by three routes that share no inner loop.
- The three routes are explicit copy placement, the line at `j/m` inside copy `i` of `k` landing at `((i-1)m + j)/(km)`; the composite scale weight built from the hyperbolic cut; and the node formula `sum_{k <= N/b} (u * v)(kb)`.
- The rectangular cut `24 x 40` is exact at every scale `m <= 24` and wrong above it, 16 of the 16 scales from 25 to 40 differing for `1 * 1`, `1 * mu` and `1 * n^-1` and 10 of 16 for `mu * mu`.
- `1 * mu = e`: the Mobius stack of the plain stack is one layer, weight 1 at scale 1 and nothing else.
- Closed-form node brightness of the five selections, checked against literal stacking at every `b <= N` for `N = 30, 61, 200, 501`, zero mismatches:
  - evens `floor(N / lcm(2, b))`
  - odds `0` at even `b`, else `ceil(floor(N/b) / 2)`
  - primes `pi(N)` at `b = 1`, `1` at prime `b <= N`, else `0`
  - squarefree `0` unless `b` is squarefree, else `sum_{d^2 <= N/b, gcd(d,b) = 1} mu(d) sum_{e | b} mu(e) floor(N/(b d^2 e))`
  - prime powers `sum_p floor(log_p N)` at `b = 1`, `floor(log_p N) - i + 1` at `b = p^i <= N`, else `0`
- The primes-only line stack lights 96 denominators at `N = 501`: `b = 1` and the 95 primes, every prime node at brightness exactly 1.
- The `s`-harmonic stack, weights `n^-s`: node `a/b` reads `b^-s H_s(floor(N/b))`, tending to `zeta(s)/b^s`; total node mass `sum_b phi(b) zeta(s) b^-s = zeta(s-1)` for `s > 2`, read at `N = 16000` as `1.644872`, `1.202057` and `1.082323` for `s = 3, 4, 5` against `zeta(2), zeta(3), zeta(4)`.
- The carpet layer covariance is rebuilt from scratch by exact rational integration on the lcm grid and matched to the gcd closed form at `(3,5), (5,7), (3,9), (15,21), (5,15)`.
- The primes-only carpet stack has zero covariance at every layer pair, so `L * Var` of the `L`-layer mean, with `L` the layer count, is the mean of the per-layer variances exactly, ratio to the independent value exactly `1` at every `L`; the values run `0.1429334753, 0.1595579958, 0.1833424270, 0.1869968711` at `L = 5, 10, 100, 1000`, rising to `3/16`.
- Per-layer exactness: `16 p^4 Var_p = 3p^4 - 4p^3 - 2p^2 + 4p - 1`, no breach over the first 1000 odd primes, so the limit constant is exactly `3/16` and `c = sqrt(3)/4 = 0.4330127`.
- The tree's full odd stack over the same estimator reads `L * Var = 0.2708541` at `L = 4000`, ratio `1.446579` and `c` factor `1.202738`, still climbing.
- The squarefree-odd stack is not uncorrelated: `Cov(15, 21) = 284/99225 = 0.0028621819`, Pearson `0.0165610084`, and 64087 of its first 1000 layer pairs share a factor, giving ratio `1.308596`.
- Davenport's expansion at `a = mu`, `sum_{n >= 1} mu(n)/n ({nx} - 1/2) = -sin(2 pi x)/pi`, truncated at `n <= 10^5` over five rational `x`: max error `5.49e-03, 1.37e-03, 2.08e-04` at cuts `10^3, 10^4, 10^5`.

## RUN

- `uv run python research/stack/lab/py/stack-algebra/stack_algebra.py`
- From the repo root. One core, under two seconds.
- Domain is the full source domain: `N = 60` for the convolution check, `N = 30, 61, 200, 501` for the selections, `N = 1000, 4000, 16000` for the harmonic stack, `L` up to 4000 layers for the variance, `n <= 10^5` for Davenport.
- Nothing is written to disk.

## WITNESSES

- The convolution theorem for weighted stacks and its hyperbolic cut, with the four checked pairs and their node counts.
- `1 * mu = e`, the Mobius stack of the plain stack collapsing to a single layer.
- The five selection closed forms and the primes-only stack's 96 lit denominators at `N = 501`.
- The `s`-harmonic node form and the total node mass `zeta(s-1)`.
- The primes-only carpet stack fading at exactly the independent rate, ratio `1` at every `L`, constant `3/16`, `c = sqrt(3)/4`.
- The squarefree-odd stack failing that property, `Cov(15, 21) = 284/99225`.
- Davenport's `mu` expansion checked numerically at truncation `10^5`.

## NOTE

- Every layer pair here is a carpet pair, so the gcd covariance law is used as a landed theorem and re-derived only at the five spot-check pairs.
- The exact rational `L * Var` at `L = 5` is `36324973301545304/254139019762753125`; beyond that the numerators run to thousands of digits and only the float is printed.
- The `{nx}` form of Davenport's identity differs from the sawtooth form by `(1/2) sum mu(n)/n`, which the script reads as `-0.00048723` at `n <= 10^5`.
