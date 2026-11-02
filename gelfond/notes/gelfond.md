---
title: The digit-sum twist on a design
lead: Under the per-digit certificate (W), which holds at every one-missing-digit base from `92317`, twisting a missing-digit design with two consecutive digits by its digit sum costs its transform a fixed share of its mass at every frequency once `norm((b - 1) alpha) > 0`, so the major arcs of the unconditional dissection need no L-function, the Mobius and von Mangoldt sums over the design save a power of `x`, and the primes of the design equidistribute in the classes of their digit sum mod `q` coprime to `b - 1`; with the design read in progressions mod every `N` coprime to `b`, the digit sum of its primes is uniform mod every `q` on the classes prime to `gcd(q, b - 1)`.
figure: research-gelfond
slug: gelfond
---

Fix a base `b >= 3`, a set `E` of `m >= 1` excluded digits, the rest `F` with `fill = b - m`, and the missing-digit design `S_F` of [mobius](../../mobius/notes/mobius.md): the positive integers whose base-`b` digits all lie in `F`. That page bounds `M_F(x) = sum_(n in S_F, n <= x) mu(n)` unconditionally only at the shape `A_F(x) exp(-c sqrt(log x))`, and [coprime](../../coprime/notes/coprime.md) carries the same shape for the primes of the design; in both, the ceiling is the major arcs, where `mu` and `Lambda` must be read in progressions mod powers of `b` through a zero-free region. This page twists the design by its digit sum and shows that the twisted transform loses a fixed share of its mass per pair of digits at every frequency, so the major arcs need no L-function and the saving becomes a power of `x`. Every number is printed by [lab/py/digit-sum-twist](../lab/py/digit-sum-twist/).

## The objects

- `s_b(n)` is the sum of the base-`b` digits of `n`, `D_k` the digit strings of length `k` over `F` read as integers in `[0, b^k)`, `g_F(t) = sum_(d in F) e(d t)` the digit symbol and `A_F(x)` the count of `S_F` up to `x`, all as on [mobius](../../mobius/notes/mobius.md).
- The twisted weight is `w(n) = 1_(S_F)(n) e(alpha s_b(n))` for a real `alpha`, and its transform at depth `k` is `hat w_k(t) = sum_(n in D_k) e(alpha s_b(n) + n t)`. With `n = sum_j d_j b^j` the phase is `sum_j d_j (b^j t + alpha)`, and the digits are independent, so `hat w_k(t) = prod_(j < k) g_F(b^j t + alpha)`. At `alpha = 0` this is the transform `hat F_k` of [mobius](../../mobius/notes/mobius.md).
- The twisted meters are `M_F(x; alpha) = sum_(n in S_F, n <= x) mu(n) e(alpha s_b(n))` and `psi_F(x; alpha) = sum_(n in S_F, n <= x) Lambda(n) e(alpha s_b(n))`.
- `beta = norm((b - 1) alpha)` is the distance from `(b - 1) alpha` to the nearest integer, and the lattice is the set `(b - 1)^(-1) Z` of `alpha` with `beta = 0`.
- `B_b(F) = sup_t sum_(r mod b) abs(g_F((t + r)/b))` is the one-step constant of [mobius](../../mobius/notes/mobius.md), a supremum over every shift `t`. The certificate (W) of that page's unconditional dissection is `B_b(F) <= fill b^(alpha_1)` with `alpha_1 < 1/5`, read from the per-digit constants `PB_b(m)` or `PB'_b(1, e_0)`; `alpha_b = log(fill)/log(b)`.
- `V_F = fill^(-1) sum_(d in F) (d - mean(F))^2` is the variance of the digit set, `K_b(x) = max(sin(pi b x)/sin(pi x), 1/sin(pi/b))` and `x_f = 1/(2(b - 1))`.

## Chapters

- [The two-position lemma to the power saving](two-position-lemma-to-regions.md): The two-position lemma and the twisted peel give the transform of the twisted design a fixed loss per pair of positions, which sums to a power saving, and the regions of the dissection say which hypothesis each needs.
- [Primes of a design by digit sum](primes-by-digit-sum.md): The lattice and the ceiling, the equidistribution of the primes of a design in the classes of their digit sum, the falsification of the claim past its range, and the comparison with the literature.
- [Progressions and every modulus](progressions-every-modulus.md): Read in progressions mod every N coprime to b, the design has its digit sum uniform mod every q on the classes prime to gcd(q, b - 1).

## Open

- **The pairing, which sets the exponent past about `10^5`.** The lemma books one deficit per two positions; the readings show close to one per position, `0.718676` per digit at `alpha = 1/2` against the proved `0.921947`. Since `sigma` binds from `100084` at `beta = 1/2` and from `93409` at `beta = 0.2` (the exponent bullet), the two-position constant, through regions B and C, limits `delta_0` at almost every admissible base, and a deficit per position would raise it there.
- **A lower base.** The wall is the per-digit (W). A twisted chain would need, in place of the Lebesgue sums of `D_(b^l)`, the grid `l^1` norms of the full-digit twisted kernel `prod_(j < l) D_b(b^j t + alpha)`, the discrete transform Mauduit and Rivat estimate; a twisted window would need a supremum over shift vectors. Either could bring the theorem toward the untwisted `584` or `115`.
- **The lattice at power strength.** The section Progressions and every modulus carries every lattice `alpha` and every `q` at `exp(-c sqrt(log x))`, through (M_N) at the moduli `N b^j`; on the lattice the major arcs keep full mass, so a power saving there asks a zero-free strip for the `L(s, chi)` mod `N b^j`, and the error term of the law at every `q` stays root-log while the off-lattice part saves a power. The off-lattice rate at small base is set by the frequency `l/(b + 1)`, where the locating step of the two-position lemma is attained.
- **The untwisted C1 at a power,** and with it whether the root-log ceiling of [mobius](../../mobius/notes/mobius.md) has any source besides the zero-free region at C2.
- **The exponent.** Region A's `1/5 - alpha_1` binds only on a window just above `92317`, and past it `sigma`, below `gamma/4`, binds; either way `delta_0` stays below `alpha_b/20`, and nothing here touches the sparse columns or the exponent `theta(F)` of [mobius](../../mobius/notes/mobius.md).
