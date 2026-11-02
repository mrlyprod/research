# totient-constant

- The node count of the spun stacks is a totient sum over ideals, and this study proves its constant instead of fitting it.
- It serves `lab/py/spun-stack` (the Gaussian stack, node count 672, 10608, 168088 at norm bounds 50, 200, 800) and `lab/py/eisenstein-stack` (the hexagonal stack, 630, 9606, 151020, 337026, 945486, 2419950 at 50, 200, 800, 1200, 2000, 3200), whose limit constants were Conjecture and are now Proved.
- Everything is recounted from scratch here, twice: once by enumerating associate classes and factoring norms, once by the norm-indexed Dirichlet series, which is the route that also carries class number above one.

## THE THEOREM

- Let `K` be an imaginary quadratic field with ring of integers `O_K`, class number `h`, unit group of order `w`, and discriminant `D_K`. Write `N` for the ideal norm and `Phi(a) = |(O_K/a)^*| = N(a) prod_{p | a} (1 - 1/N(p))`.
- Proved: `sum_{N(a) <= N} Phi(a) = (rho_K/(2 zeta_K(2))) N^2 + O(N^(3/2))` with `rho_K = 2 pi h/(w sqrt |D_K|)` the ideal density, the sum running over nonzero ideals.
- At `h = 1` nonzero ideals are nonzero associate classes, which is the spun stacks' index set.
- Proved: for `K = Q(i)`, `h = 1`, `w = 4`, `|D_K| = 4`, the constant is `pi/(8 zeta(2) G) = 0.260634696495` with `G` Catalan's constant.
- Proved: for `K = Q(sqrt -3)`, `h = 1`, `w = 6`, `|D_K| = 3`, the constant is `pi/(6 sqrt 3 zeta(2) L(2, chi_-3)) = 0.235217881630`.
- The error `O(N^(3/2))` is the elementary one derived below; THE CIRCLE EXPONENT sharpens it to `O(N^(4/3))`, and the true size reads `N^(5/4)`, not the `N log N` that the deviations the two studies print were first read as (Refuted, `lab/rs/node-error`); THE MEAN SQUARE proves the mean square of the error is `O(N^(5/4))`.

## CONVENTIONS

- The sum is over nonzero integral ideals throughout; at class number one the unit group acts freely on nonzero elements, each class holds exactly `w` elements, and ideals and associate classes are the same monoid.
- `Phi` is the ideal totient, the order of the unit group of `O_K/a`, and it is multiplicative because the Chinese remainder theorem splits `O_K/a`.
- `mu_K` is the ideal Mobius function, `(-1)^k` on a squarefree product of `k` prime ideals and `0` otherwise.
- The Gaussian window is `a + b i` with `a >= 1` and `b >= 0`, one representative per class; the Eisenstein window is `a + b omega` minimal in its orbit under the six units, `omega` acting as `(a, b) -> (-b, a - b)` and `-1` as `(a, b) -> (-a, -b)`.
- Norm bound `N` throughout, never a radius: the height `|q| <= T` of a complex Farey fraction is the norm bound `N = T^2`.

## THE IDEAL COUNT

- Lemma, Proved: `A(t) = #{a : N(a) <= t} = rho_K t + E(t)` with `|E(t)| = O(sqrt t)` and `rho_K = 2 pi h/(w sqrt |D_K|)`.
- First the principal ideals. `O_K` is a rank two lattice in `C` of covolume `sqrt |D_K|/2`, so counting its points in the disc of radius `R = sqrt t` by attaching to each point its fundamental cell, of diameter `delta`, gives cells contained in the disc of radius `R + delta` and covering the disc of radius `R - delta`; dividing the two areas by the covolume gives `2 pi t/sqrt |D_K| + O(sqrt t)` lattice points, and removing the origin and dividing by `w` gives `2 pi t/(w sqrt |D_K|) + O(sqrt t)` principal ideals, with an effective implied constant, `4 pi delta/(w sqrt |D_K|)` at leading order.
- Then one class at a time. Fix an ideal class `C` and an integral ideal `a` in `C^-1`; the map `b -> a b = (x)` is a bijection from the integral ideals of `C` of norm at most `t` to the nonzero elements of `a` of norm at most `t N(a)`, taken up to units, since `(x) a^-1` is integral exactly when `x` lies in `a`.
- `a` is a lattice of covolume `N(a) sqrt |D_K|/2`, so the same cell argument gives `2 pi t/(w sqrt |D_K|) + O(sqrt t)` ideals in every class alike, and summing the `h` classes gives the lemma.
- This is the Gauss circle count for `Z[i]` and the hexagonal count for `Z[omega]`, whose closed forms `sum_j (floor(t/(4j+1)) - floor(t/(4j+3)))` and `sum_j (floor(t/(3j+1)) - floor(t/(3j+2)))` are Verified in `lab/py/spun-stack` and `lab/py/eisenstein-stack`; the lemma above is derived here from the lattice and cites neither.
- `rho_K` is also the residue of `zeta_K` at `s = 1` given by the class number formula, which is the cross-check: `pi/4` for `Q(i)` and `pi/(3 sqrt 3)` for `Q(sqrt -3)`.

## THE CONVOLUTION STEP

- The method is Dirichlet convolution followed by Abel summation, not a hyperbola split: the whole range of the ideal `b` is summed and the inner sum is estimated uniformly.
- Proved: `Phi = N * mu_K` as a Dirichlet convolution on ideals, since both sides are multiplicative and on a prime power `(N * mu_K)(p^k) = N(p)^k - N(p)^(k-1) = Phi(p^k)`; equivalently `sum_a Phi(a) N(a)^-s = zeta_K(s - 1)/zeta_K(s)`.
- Hence `S(N) = sum_{N(a) <= N} Phi(a) = sum_{N(b) <= N} mu_K(b) T(N/N(b))` with `T(X) = sum_{N(c) <= X} N(c)`.
- Abel summation on the lemma gives `T(X) = X A(X) - int_0^X A(u) du = (rho_K/2) X^2 + X E(X) - int_0^X E(u) du`, so `|T(X) - (rho_K/2) X^2| <= (5/3) C X^(3/2)` with `C` the constant of the lemma.
- The main term is `(rho_K/2) N^2 sum_{N(b) <= N} mu_K(b) N(b)^-2`, and the tail `sum_{N(b) > N} N(b)^-2 = O(1/N)` by Abel summation on `A(u) = O(u)`, so the main term is `(rho_K/(2 zeta_K(2))) N^2 + O(N)`.
- The error term is at most `(5/3) C N^(3/2) sum_{N(b) <= N} N(b)^(-3/2) <= (5/3) C zeta_K(3/2) N^(3/2)`, finite because `zeta_K(s)` converges for `s > 1`.
- Adding the two gives the theorem, with an effective implied constant `(5/3) C zeta_K(3/2) + O(1)`.
- Proved, conditional refinement: if the lemma holds with `E(t) = O(t^theta)` for some `0 < theta < 1`, the same three steps give error `O(N^(1 + theta))`; the elementary `theta = 1/2` is the one derived here, and no sharper circle exponent is claimed. The lower constraint is not decoration: at `theta = 0` the error layer needs `zeta_K(1)`, which diverges, and the conclusion drops to `O(N log N)`.

## THE CIRCLE EXPONENT

- Lemma, Proved from van der Corput's ellipse bound: `A(t) = rho_K t + P(t)` with `|P(t)| <= C_K t^(1/3)` for all `t > 0`, `C_K` depending on `K` alone.
- Class by class as in THE IDEAL COUNT: with `omega_1, omega_2` a basis of the ideal `a`, `N(x omega_1 + y omega_2) = N(a) f(x, y)` for a positive definite integral binary quadratic form `f = A x^2 + B x y + C y^2` with `4 A C - B^2 = |D_K|`, and the ideals of the class with norm at most `t` number `(L(t) - 1)/w` with `L(t) = #{(x, y) in Z^2 : f(x, y) <= t}`.
- The form `Q = (2/sqrt |D_K|) f` has `a c - b^2 = 1` in the shape `a t_1^2 + 2 b t_1 t_2 + c t_2^2`, and `f <= t` reads `Q <= x^2` with `x^2 = 2 t/sqrt |D_K|`, whose area `pi x^2 = 2 pi t/sqrt |D_K|` is the main term of THE IDEAL COUNT.
- Read in the arXiv source of arXiv:math/0410522, Section 3.1, the ellipse bound is `|P(x; E_2)| < 38 x^(2/3) + (3 a^(3/4) + 700 c^(3/4)) sqrt x + 11` for `x >= 1`, with `x` the linear scale, so `|L(t) - 2 pi t/sqrt |D_K|| <= C' t^(1/3)` once `t >= sqrt |D_K|/2`; below that, and for `t < 1`, the error is bounded by a constant of `K`, absorbed in `C_K`. The primary papers of van der Corput and Chaix are unread.
- Summing the `h` classes, dividing by `w` and removing the origin gives the lemma; for `t < 1` there is no ideal and `P(t) = -rho_K t`, of size at most `rho_K t^(1/3)`.
- Theorem, Proved: `sum_{N(a) <= N} Phi(a) = (rho_K/(2 zeta_K(2))) N^2 + O(N^(4/3))` for every imaginary quadratic field.
- Let `C = max(C_K, rho_K)`. Abel summation gives `T(X) = X A(X) - int_0^X A = (rho_K/2) X^2 + X P(X) - int_0^X P`, so `|T(X) - (rho_K/2) X^2| <= C X^(4/3) + (3/4) C X^(4/3) = (7/4) C X^(4/3)`.
- In `S(N) = sum_{N(b) <= N} mu_K(b) T(N/N(b))` replace each `T` by its main term: the replacements cost at most `(7/4) C N^(4/3) sum_b N(b)^(-4/3) <= (7/4) C zeta_K(4/3) N^(4/3)`, finite because `4/3 > 1`.
- What remains is `(rho_K/2) N^2 sum_{N(b) <= N} mu_K(b) N(b)^-2`, and `|sum_{N(b) > N} mu_K(b) N(b)^-2| <= sum_{N(b) > N} N(b)^-2 = int_N^infty u^-2 dA(u) <= 2 (rho_K + C)/N` by `A(u) <= (rho_K + C) u`, so it equals `(rho_K/(2 zeta_K(2))) N^2 + O(N)`.
- Adding, the error is at most `(7/4) C zeta_K(4/3) N^(4/3) + O(N)`.
- Over `Q(i)`, `A(t) = (sum_{n <= t} r_2(n) - 1)/4`, so Huxley's circle bound (equation (1.7) of the same source), `P(x) = O(x^(131/416) (log x)^(18637/8320))`, gives `O(N^(1 + theta))` for every `theta > 131/416` by the same steps; this is stated for `Q(i)` only, the other fields staying at `1/3`.
- The route bounds each Mobius layer in absolute value, and `P(t) = Omega(t^(1/4))` by Hardy (same source), so it cannot give an exponent below `5/4`; `lab/rs/node-error` measures `N^(5/4)` and `max |E|/N^(4/3)` at most `0.5425` Gaussian and `0.4395` Eisenstein over the windows `2^6 .. 2^22`.

## THE MEAN SQUARE

- Theorem, Proved: with `E(x) = sum_{N(a) <= x} Phi(a) - (rho_K/(2 zeta_K(2))) x^2` for real `x > 0`, `int_N^(2N) E(x)^2 dx <= C_K N^(7/2)` for all `N >= 1` and every imaginary quadratic field, so `(1/N) int_N^(2N) E^2 = O(N^(5/2))` and the root mean square of `E` over `[N, 2N]` is `O(N^(5/4))`.
- The same holds over the integers: `E(x)` differs from `E(n)` by at most `(rho_K/(2 zeta_K(2))) (2n + 1)` on `[n, n + 1)`, so `sum_{N <= n < 2N} E(n)^2 <= 2 int_N^(2N) E^2 + O(N^3) = O(N^(7/2))`.
- Layers, exact: write `T(y) = sum_{N(c) <= y} N(c)` for real `y > 0`, zero below 1, and `E_1(y) = T(y) - (rho_K/2) y^2`. The convolution `Phi = N * mu_K` gives `S(x) = sum_b mu_K(b) T(x/N(b))` over all ideals `b`, the terms with `N(b) > x` being zero, and `c_K = (rho_K/2) sum_b mu_K(b) N(b)^-2`, so `E(x) = sum_b mu_K(b) E_1(x/N(b))`, with no tail term. The sum converges absolutely because `abs(E_1(y)) <= B y^2` for `0 < y < 2`, `B = 1 + 2 rho_K`, the only ideal of norm below 2 being `O_K`.
- Lattice input, class by class: the ideals of a class with norm at most `t` number `(L(t) - 1)/w` with `L(t) = #{(x, y) in Z^2 : f(x, y) <= t}` for a positive definite integral form `f` of discriminant `D_K`, so `P(t) = A(t) - rho_K t = (1/w) sum_classes (P_f(t) - 1)` with `P_f(t) = L(t) - (2 pi/sqrt |D_K|) t`.
- The ellipse `Omega = {f <= 1}` is a convex domain with the origin inside, `C^infinity` boundary and non-vanishing curvature, and `L(t^2) = card(t Omega ∩ Z^2)`. Iosevich, Sawyer and Seeger (arXiv:math/0205127, Theorem 1.2.1, read at source) at `h = R` give `int_R^(2R) (L(t^2) - area(Omega) t^2)^2 dt <= C R^2` for `R >= 2`; the substitution `tau = t^2` has `dtau = 2t dt <= 4R dt` and gives `int_Y^(4Y) P_f(tau)^2 dtau <= 4 C Y^(3/2)` for `Y >= 4`, and `P_f` is bounded on `[1, 8]`.
- Summing over the classes by Minkowski, with the constant `h/w`: `int_Y^(2Y) P(t)^2 dt <= C_P Y^(3/2)` for `Y >= 1`. For `Q(i)` this is the circle's mean square, display (1.13) of arXiv:math/0410522, `int_0^X P^2 = C X^(3/2) + Q(X)` with `Q(X) = O(X^(5/4 + eps))` already by Cramer; no `eps` is lost on the main order.
- The integral term needs no bound of its own: splitting `[0, y]` dyadically down to `[0, 2)`, where `P` is bounded, gives `int_0^y P^2 <= C_3 y^(3/2)` for `y >= 1`, and Cauchy-Schwarz gives `abs(int_0^y P) <= y^(1/2) (int_0^y P^2)^(1/2) <= C_3^(1/2) y^(5/4)`.
- Hence `E_1(y) = y P(y) - int_0^y P` has `int_Y^(2Y) E_1^2 <= C_1 Y^(7/2)` for `Y >= 1`: the two terms have squared norms at most `4 C_P Y^(7/2)` and `2^(5/2) C_3 Y^(7/2)`, so `C_1 = 2 (4 C_P + 2^(5/2) C_3)`.
- Scaling: for an ideal `b` of norm `m`, the `L^2([N, 2N])` norm of `x -> E_1(x/m)` is `m^(1/2)` times the `L^2` norm of `E_1` on `[N/m, 2N/m]`; for `m <= N` the window starts at `N/m >= 1` and the norm is at most `C_1^(1/2) m^(-5/4) N^(7/4)`.
- For `m > N` the window lies in `(0, 2)` and the norm is at most `4 B m^(1/2) (N/m)^(5/2) = 4 B N^(5/2) m^-2`, and `sum_{N(b) > N} N(b)^-2 <= 2 (rho_K + C)/N` by THE CIRCLE EXPONENT.
- Minkowski over the layers: `norm(E) <= C_1^(1/2) zeta_K(5/4) N^(7/4) + 8 B (rho_K + C) N^(3/2)`, finite because `5/4 > 1`. Squaring gives the theorem.
- Why `5/4`: a layer of norm `m` is the `b = 1` term rescaled to the window `N/m`, of rms size `(N/m)^(5/4)`, so the layers sum exactly when the exponent exceeds 1, as `zeta_K(5/4)`; the route bounds the layers and says nothing about cancellation between them.
- Not proved here: the lower half, that `rms(E)` over `[N, 2N]` is at least a positive multiple of `N^(5/4)`, which needs the layers not to cancel the `b = 1` term (`rms(E - E_1)` is `0.2059` and `0.1239` of `N^(5/4)` against `rms(E_1)` at `0.3984` and `0.3218`, correlation `0.8775` and `0.9359`, `lab/rs/node-error`); and the pointwise bound `E(N) = O(N^(5/4 + eps))`, which contains the circle problem at `b = 1`.
- The exact layer identity is tested in `lab/rs/node-error`: `E(n)` against `sum_m layers(m) E_1(n/m)` within `0.01` at `n = 100, 250, 400` in both fields; the measured `rms(E)/N^(5/4)` over `2^12 .. 2^22` is at most `0.2750` and `0.2535`.

## THE CONSTANTS

- `zeta_K(s) = zeta(s) L(s, chi_D)` for a quadratic field, `chi_D` the Kronecker symbol of the discriminant, so `zeta_K(2) = zeta(2) G` for `Q(i)` and `zeta(2) L(2, chi_-3)` for `Q(sqrt -3)`.
- `c = rho_K/(2 zeta_K(2)) = pi h/(w sqrt |D_K| zeta(2) L(2, chi_D))`, which is `pi/(8 zeta(2) G)` and `pi/(6 sqrt 3 zeta(2) L(2, chi_-3))`.
- `main` prints, at 50 working digits through the Hurwitz form `L(2, chi_D) = |D_K|^-2 sum_a chi_D(a) zeta(2, a/|D_K|)` with Euler-Maclaurin and `pi` by Machin: `zeta(2) = 1.644934066848`, `G = 0.915965594177`, `L(2, chi_-3) = 0.781302412896`, `zeta_K(2) = 1.506703009923` and `1.285190955484`, `rho_K = 0.785398163397` and `0.604599788078`, and the constants `0.260634696495` and `0.235217881630`.
- Those two agree with the values `lab/py/spun-stack` and `lab/py/eisenstein-stack` print, `0.260635` and `0.235217881630015`, which are now derived rather than fitted.
- `class_number` derives `h` from the finite character sum `h = w (-sum_a chi_D(a) a)/(2 |D_K|)`, returning 1, 1, 2 and 3 at `D_K = -4, -3, -20, -23`, so the constant is assembled from `w`, `|D_K|`, `h` and `L(2, chi_D)` alone.

## THE GENERAL FIELD

- Proved: the theorem holds for every imaginary quadratic field, of any class number. The convolution step never sees the class number, and the lemma is proved above one ideal class at a time.
- Verified: at `D_K = -20`, `h = 2`, the constant is `0.378582556790` and the count at norm bound 102400 is 3969730264, ratio `1.000001116781`; at `D_K = -23`, `h = 3`, the constant is `0.425700255391` and the count is 4463281160, ratio `0.999885848149` (`norm_totient_sum`).
- Proved as an implication only: for a general number field, with ideals throughout, the convolution step is unchanged and gives `(rho_K/(2 zeta_K(2))) N^2 + O(N^(1 + theta))` for any ideal count `A(t) = rho_K t + O(t^theta)` with `0 < theta < 1`; that ideal count is not derived here, so nothing about degree above two is claimed outright.

## THE FAREY CONVENTION

- The lit set of the spun stack is the complex Farey set: Sayous arXiv:2407.04380 Section 2 defines `F_t = {pr_{O_K}(p/q) : p, q in O_K, p O_K + q O_K = O_K, 0 < |q| <= e^(t/2)}` inside `C/O_K` and states `card F_t ~ c_K e^(2t)` with `c_K = pi/(sqrt |D_K| zeta_K(2))`, attributed there to Cosentino, Ergodic Theory Dynam. Systems 19 (1999) 1437-1484, Theorem 4.
- Since `|q| <= e^(t/2)` is the norm bound `N = e^t`, the two conventions are: counting denominators as ideals gives `pi/(w sqrt |D_K| zeta_K(2)) = rho_K/(2 zeta_K(2))`, and counting them as elements gives `pi/(sqrt |D_K| zeta_K(2))`, exactly `w` times larger, `w = 4` for `Z[i]` and `6` for `Z[omega]`.
- Verified: `farey_set` builds the set `F_t` itself, layer by layer as the union of the lattices `q^-1 O_K` modulo `O_K` over all classes of norm at most `N`, deduplicating points in exact integer coordinates, and returns 672, 10608, 168088 for `Z[i]` and 630, 9606, 151020 for `Z[omega]` at norm bounds 50, 200, 800, equal in all six cases to the totient sum.
- `F_t` is written as a set of points of the torus, so its cardinality is the ideal one; the factor `w` is inherited from the cited Cosentino statement rather than chosen, and `c_K e^(2t)` is the size of the family indexed by element denominators, in which each point repeats once per associate of its denominator.

## WHAT IT PRINTS

- `classes` enumerates one representative per nonzero associate class in each field, 80414 Gaussian and 61914 Eisenstein classes of norm at most 102400.
- `totient` computes `Phi` from the factorisation of the norm alone, splitting rational primes into ramified, inert and split, a split `p` contributing `(1 - 1/p)^2` exactly when `p` divides both coordinates.
- `totient_sum` reads the prefix array built by `tables` and returns 672, 10608, 168088, 372872, 1045088, 2663864, 42663808, 683283504, 2732257376 for `Q(i)` and 630, 9606, 151020, 337026, 945486, 2419950, 38560362, 616981926, 2466558234 for `Q(sqrt -3)` at norm bounds 50, 200, 800, 1200, 2000, 3200, 12800, 51200, 102400, matching the published counts of `lab/py/spun-stack` and `lab/py/eisenstein-stack` at all nine bounds where they overlap.
- `main` prints `count/(c N^2)` at those bounds: `1.031329, 1.017516, 1.007684, 0.993494, 1.002445, 0.998113, 0.999097, 1.000066, 0.999746` for `Q(i)` and `1.071347, 1.020968, 1.003192, 0.995017, 1.004904, 1.004699, 1.000578, 1.000604, 1.000049` for `Q(sqrt -3)`, oscillating about 1 rather than approaching it from one side.
- `main` prints the deviation scaled by the derived error, `(count - c N^2)/N^(3/2)`, never larger than `0.119` in absolute value across both fields and all nine bounds, and the same deviation scaled by `N log N`, which reproduces the `+0.214, +0.186, +0.090, -0.198, +0.304, +0.438` of `lab/py/eisenstein-stack`. That scale is Refuted as the size of the error: `lab/rs/node-error` recounts to norm bound 8388607, where it climbs and `N^(5/4)` is flat.
- `convolution_sum` re-evaluates the same totient sum as `sum_{N(b) N(c) <= N} mu_K(b) N(c)` from an independently computed `mu_K`, agreeing at all nine bounds in both fields, which is the numerical form of `Phi = N * mu_K`.
- `farey_set` counts the complex Farey set directly, 672, 10608, 168088 and 630, 9606, 151020 at norm bounds 50, 200, 800.
- `norm_arrays` and `norm_totient_sum` recount by norm alone, with no lattice: the ideals of norm `n` number `sum_{d | n} chi_D(d)`, the ideal Mobius function by norm is `mu * (mu chi_D)` from `1/zeta_K = (1/zeta)(1/L)`, and `S(N) = sum_k m(k) T(N/k)`. This agrees with the class route at all eighteen bounds of the two class-number-one fields and supplies the tables at `D_K = -20` and `-23`.

## RUN

- `uv run python research/stack/lab/py/totient-constant/totient_constant.py`
- From the repo root. One core, about a second, standard library only.
- Domain is norm bounds up to 102400 for the counts, up to 800 for the direct Farey set, four discriminants for the general statement, and the constants at 50 working digits.
- Writes nothing.

## WITNESSES

- The theorem's constant: `main` against the counts of `totient_sum` at nine norm bounds in both fields, with the ratio and both scaled deviations.
- The recount: `totient_sum` against the published counts of `lab/py/spun-stack` and `lab/py/eisenstein-stack`, nine matches, no mismatch.
- The convolution `Phi = N * mu_K`: `convolution_sum` against `totient_sum`, eighteen agreements.
- Class number above one: `norm_totient_sum` at `D_K = -20` and `-23`, ratios `1.000001116781` and `0.999885848149` at norm bound 102400, with `class_number` deriving `h = 2` and `3`.
- The Farey convention: `farey_set` against `totient_sum`, six agreements, settling that the set convention carries `rho_K/(2 zeta_K(2))` and the element convention `w` times that.
- The circle exponent: the lemma and the `O(N^(4/3))` theorem above, with `lab/rs/node-error` meeting the bound with room at every window to norm bound 8388607.
- The mean square: the `O(N^(5/4))` theorem above, with the layer identity tested and the windows read in `lab/rs/node-error`.
