---
title: The Mobius meter across digit designs
lead: The Mobius meter across digit designs: exact transfer between scaled columns, the 47-column cancellation census, a power saving on the dense columns under GRH, the pair route with its major-arc lemma and its `l^1` threshold `1/4`, and square-root cancellation as the open exponent.
figure: research-mobius
slug: mobius
date: 2026-09-01
---

# The Mobius meter across digit designs

The Mobius meter across digit designs: exact transfer between scaled columns, the 47-column cancellation census, a power saving on the dense columns under GRH, the pair route with its major-arc lemma and its `l^1` threshold `1/4`, and square-root cancellation as the open exponent.

- [notes/mobius.md](notes/mobius.md): The Mobius meter across digit designs
- [notes/meter-transfer-census.md](notes/meter-transfer-census.md): The meter, the transfer and the census
- [notes/digit-strings-across-divisors.md](notes/digit-strings-across-divisors.md): Digit strings across divisors
- [notes/grh-power-saving.md](notes/grh-power-saving.md): A power saving under GRH at large base
- [notes/pair-route.md](notes/pair-route.md): The pair route
- [notes/unconditional-dissection.md](notes/unconditional-dissection.md): The unconditional dissection
- [notes/half-interval-even-digits.md](notes/half-interval-even-digits.md): The half interval and the even digits
- [notes/digit-progressions.md](notes/digit-progressions.md): The digit progressions
- [claims/design-mobius-meter.md](claims/design-mobius-meter.md): The design Mobius meter
- [claims/digit-progressions.md](claims/digit-progressions.md): The digit progressions
- [claims/digit-restricted-mobius-exponent.md](claims/digit-restricted-mobius-exponent.md): The digit-restricted Mobius exponent
- [claims/digit-restricted-mobius-meter.md](claims/digit-restricted-mobius-meter.md): The digit-restricted Mobius meter
- [claims/digit-strings-across-divisors.md](claims/digit-strings-across-divisors.md): Digit strings across divisors
- [claims/even-digits.md](claims/even-digits.md): The even digits
- [claims/half-interval.md](claims/half-interval.md): Half interval
- [claims/moments-of-the-digit-transform.md](claims/moments-of-the-digit-transform.md): The moments of the digit transform
- [papers/first-base-below-a-quarter.md](papers/first-base-below-a-quarter.md): The First Base Below a Quarter
- [papers/half-interval-mobius.md](papers/half-interval-mobius.md): An Unconditional Mertens Bound on the Half Interval
- [papers/sparse-mertens-under-grh.md](papers/sparse-mertens-under-grh.md): A Power Saving for Mobius Sums on Missing-Digit Integers, under GRH
- [refs.md](refs.md): the references this topic cites
- [lab/py/base3-exponent/](lab/py/base3-exponent/): The `l^1` exponent of the digit transform of `{0,1}` at base 3, enclosed in interval arithmetic with its floor, its one-step and two-step bounds.
- [lab/py/design-meter/](lab/py/design-meter/): Reads the oscillation spectrum of the digit-design Mobius meter `M_F(x) = sum of mu(n)` over `n in S_F`, `n <= x`, with `F` the kept digits, normalised by `x^(alpha/2)` with `alpha = log_base fill`, and asks which frequencies it carries.
- [lab/py/digit-progressions/](lab/py/digit-progressions/): The numbers of the digit progressions on mobius
- [lab/py/even-digits/](lab/py/even-digits/): The numbers of the even digits on mobius
- [lab/py/interval-digits/](lab/py/interval-digits/): The numbers of the half interval on mobius
- [lab/py/mobius-dissection/](lab/py/mobius-dissection/): The numbers of the unconditional dissection on mobius
- [lab/py/mobius-region/](lab/py/mobius-region/): The three exponents a digit design hands the pair route, and the region of `(alpha, alpha_1, beta)` in which that route's inequalities all hold.
- [lab/py/mrly-pairing/](lab/py/mrly-pairing/): The identity that replaces `zeta(s) M(s) = 1` on a digit design `S_F`, with `F` its kept digits, and what the position pairing `M_F(base^level) = int_0^1 G_level(t) S_level(t) dt` costs when it is split by `base`-power denominator.
- [lab/rs/carry-free-mobius/](lab/rs/carry-free-mobius/): The design Mobius `nu_F` of the two-digit design `F = {0,1}` is base-free
- [lab/rs/mertens-numerology/](lab/rs/mertens-numerology/): Computes the kernel-sum constants attached to a digit set in base `base` that omits `m` of the `base` digits, and the exponent bookkeeping they force.
- [lab/rs/mobius-designs/](lab/rs/mobius-designs/): Regenerates every number on mobius
- [lab/rs/rho-decoupling/](lab/rs/rho-decoupling/): Computes `N_F(level
- [figures/paper-first-base-below-a-quarter.ts](figures/paper-first-base-below-a-quarter.ts): paper-first-base-below-a-quarter
- [figures/paper-first-base-below-a-quarter.json](figures/paper-first-base-below-a-quarter.json): paper-first-base-below-a-quarter census
- [figures/paper-half-interval-mobius.ts](figures/paper-half-interval-mobius.ts): paper-half-interval-mobius
- [figures/paper-half-interval-mobius.json](figures/paper-half-interval-mobius.json): paper-half-interval-mobius census
- [figures/paper-sparse-mertens-under-grh.ts](figures/paper-sparse-mertens-under-grh.ts): paper-sparse-mertens-under-grh
- [figures/paper-sparse-mertens-under-grh.json](figures/paper-sparse-mertens-under-grh.json): paper-sparse-mertens-under-grh census
- [figures/research-mobius.ts](figures/research-mobius.ts): research-mobius

## Lost

- Coefficient-free routes to the Type II sum: two box witnesses floor them at `x^alpha` on the dense cells, and the statement for every digit set stays a Conjecture (The meter and its yardstick).
- Weighting the frequency grid by denominator: it buys no exponent, since the `l^1` mass of the digit transform sits at the top denominator; the Refuted row in [claims/digit-restricted-mobius-meter.md](claims/digit-restricted-mobius-meter.md).
- The pair route at base 10: certified `beta > 1/4` at all ten one-missing-digit sets kills both branches of the criterion (The pair route, "Base 10 is refuted at every excluded digit").
- The pair route below `fill = base^(3/4)`: the Parseval floor `m_t >= 1 - alpha` shuts the gate (The pair route).
- A certificate of any shape `Sigma_i <= C lambda^i` at `L = base/3`: the spectral radius of `T_L` reads `alpha_1 = 0.207` at base `1000003` (`0.273` at `10007`, falling), so nothing of that shape clears `1/5` below about `2 * 10^6`; a reading on a 401-point grid, not a theorem (digit-progressions, verb `headroom`).
