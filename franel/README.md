---
title: The restricted Franel identity
lead: Franel's identity on the Farey fractions whose denominators lie in a digit design: the exponential sums are dilated Mertens sums, the square-root threshold is equivalent to their square-root cancellation in mean square over the dilates, the dilates are regular languages read by long multiplication, their masses need no Mobius input, and the design's repunits refute the uniform law that would give a rate in `d`.
figure: research-franel
slug: franel
date: 2026-09-23
---

# The restricted Franel identity

Franel's identity on the Farey fractions whose denominators lie in a digit design: the exponential sums are dilated Mertens sums, the square-root threshold is equivalent to their square-root cancellation in mean square over the dilates, the dilates are regular languages read by long multiplication, their masses need no Mobius input, and the design's repunits refute the uniform law that would give a rate in `d`.

- [notes/franel.md](notes/franel.md): The restricted Franel identity
- [notes/identity-and-recursion.md](notes/identity-and-recursion.md): The identity, the dilates and the recursion
- [notes/threshold-surrogate-repunits.md](notes/threshold-surrogate-repunits.md): The threshold, the surrogate and the repunits
- [notes/strict-set-and-average.md](notes/strict-set-and-average.md): The strict set and the average over designs
- [claims/franel-on-a-digit-design.md](claims/franel-on-a-digit-design.md): Franel on a digit design
- [papers/franel-converse.md](papers/franel-converse.md): The Franel Threshold Is a Mean Square
- [refs.md](refs.md): the references this topic cites
- [lab/py/design-averages/](lab/py/design-averages/): Averages `abs(M_S(x))^2`, `M_S(x) = sum_(n <= x, n in S) mu(n)`, over the three flat families of the Franel note's section The average over designs
- [lab/py/restricted-franel/](lab/py/restricted-franel/): Franel's identity on a digit design
- [lab/py/shifted-sums/](lab/py/shifted-sums/): Walks the shifted sums `T_s = sum_{w in B_s} mu(6w + 1)` of the Franel note's repunit section, `B_s` the `2^s` integers below `3^s` whose base 3 digits lie in `{0,1}`, to `s = 26`.
- [lab/rs/design-recursion/](lab/rs/design-recursion/): The design meter `M_F` at base 3 `{0,1}` to `3^27` beside `A`, its sum over the elements coprime to `3`, with the recursion `M_F(x) = A(x) - A(x/3)` of the Franel note's section The recursion in the base asserted at every rung; verbs `level`, `fake` and `pretend` for the mean square over a level, the multiplicative class and its pretentious members.
- [figures/paper-franel-converse.ts](figures/paper-franel-converse.ts): paper-franel-converse
- [figures/paper-franel-converse.json](figures/paper-franel-converse.json): paper-franel-converse census
- [figures/research-franel.ts](figures/research-franel.ts): research-franel
- [figures/research-franel.json](figures/research-franel.json): research-franel census

## Lost

- Square-root cancellation on average over a flat family of designs: not cheap, it costs the zero-free half plane `Re s > 1 - alpha/2` for `zeta` (The average over designs).
- The accepting-set law uniform in coprime `d <= Q^theta`: the design's own repunits refute it (The repunits; the Refuted row in [claims/franel-on-a-digit-design.md](claims/franel-on-a-digit-design.md)).
- The mean square as a route to `M_F` unbounded: for every `N` a bounded non-multiplicative fake agrees with `mu` on `[1, N]`, has its support and its recursion, so no Parseval identity, no use of the recursion and no finite segment alone decides it; multiplicativity must enter (The mean square over a level).
