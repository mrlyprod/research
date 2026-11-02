---
title: Cuts
lead: The slice sign law: the central diagonal cut of the `D`-dimensional middle-digit solid obeys a recurrence of order `ceil(D/2)` and grows strictly slower than a generic cut in every even dimension at bases 3 and 5, and at least as fast in every odd dimension at base 3.
figure: research-cuts
slug: cuts
date: 2026-08-28
---

# Cuts

The slice sign law: the central diagonal cut of the `D`-dimensional middle-digit solid obeys a recurrence of order `ceil(D/2)` and grows strictly slower than a generic cut in every even dimension at bases 3 and 5, and at least as fast in every odd dimension at base 3.

- [notes/cuts.md](notes/cuts.md): Cuts
- [papers/slice-sign-law.md](papers/slice-sign-law.md): The Slice Sign Law
- [claims/conjecture-s-even-half.md](claims/conjecture-s-even-half.md): Conjecture S: even half
- [claims/conjecture-s-odd-half.md](claims/conjecture-s-odd-half.md): Conjecture S: odd half
- [claims/conjecture-s-the-even-half-transient.md](claims/conjecture-s-the-even-half-transient.md): Conjecture S: the even-half transient
- [claims/diagonal-cuts-of-the-parity-solid.md](claims/diagonal-cuts-of-the-parity-solid.md): Diagonal cuts of the parity solid
- [claims/diagonal-slice-ladder.md](claims/diagonal-slice-ladder.md): Diagonal slice ladder
- [claims/slice-ladder-rate-and-dead-routes.md](claims/slice-ladder-rate-and-dead-routes.md): The slice ladder: rate and dead routes
- [claims/slice-sign-law-in-every-dimension.md](claims/slice-sign-law-in-every-dimension.md): Slice sign law in every dimension
- [refs.md](refs.md): the references this topic cites
- [lab/py/base3-transient-exhaustion/](lab/py/base3-transient-exhaustion/): Settles `level*(D)`, the last level with `V(level) < 0`, at every even carry parameter `D = 6..120` in base 3, middle-digit design, with the exhaustion proved rather than assumed.
- [lab/py/even-block-valuations/](lab/py/even-block-valuations/): Prints `det M_even` exactly and its `v_p` for the base-3 middle-digit diagonal slice at every `dim = 2..40`, and tests the route "Lemma H by one prime": a prime `p` with `v_p(det M_even) = 0` at every `dim`.
- [lab/py/pascal-shear/](lab/py/pascal-shear/): The level-`level` cell set of the 2D design with code `7`, tile `[[1, 1], [1, 0]]`, is `{(i, j) in [0, 2^level)^2 
- [lab/py/slice-sign-certificates/](lab/py/slice-sign-certificates/): Certifies every numerical step of the even-half theorems at bases 3 and 5 of the paper The Slice Sign Law, Theorems 6.1 and 7.1, and regenerates the exact certificates, determinants and anchors the topic's claim rows print.
- [lab/py/slice-ladder-controls/](lab/py/slice-ladder-controls/): Checks the anti-diagonal profile identity `P_{A (x) B}(t) = P_A(t^side_B) P_B(t)` on every base-2 word of length 2 and 3, counting mismatches.
- [lab/py/slice-ladder-rate/](lab/py/slice-ladder-rate/): Settles the decay rate `r_inf` of the slice excess `rho_dim - fill/3` of the base-3 central diagonal slice: the 3-adic cosine product `1/prod_(k>=2) cos(2 pi/3^k)` against the retired Richardson fit, from exact Perron roots to `dim = 120`.
- [figures/paper-slice-sign-law.ts](figures/paper-slice-sign-law.ts): paper-slice-sign-law
- [figures/research-cuts.ts](figures/research-cuts.ts): research-cuts

## Lost

- A fixed-epsilon spectral separation for the slice transfer matrix: none is expected, `lambda_1/abs(lambda_2)` reading `(dim+2)/(dim-2) + O(dim^-3)` and tending to 1 on the computed range, a conjecture; the row in [claims/slice-sign-law-in-every-dimension.md](claims/slice-sign-law-in-every-dimension.md).
- Lemma H by a single prime (`H det M_even != 0` from one nonzero residue): dead for every prime `p <= 59` on `dim = 2..40`, with `61, 67, 73, 97` dividing no determinant there and untested beyond; `v_p(det M_even) > 0` at `dim = 2..40` for `p = 2` at 39 of 39 dims, `p = 3` at 38, `p = 5` at 19 and `p = 7` at 13, and every prime `p <= 59` divides it at some `dim <= 30` (`dim = 5` for `p = 5`, `dim = 7` for `p = 7`, `dim = 30` for `p = 41`), `v_3` reaching 181 at `dim = 40`; the row in [claims/diagonal-slice-ladder.md](claims/diagonal-slice-ladder.md).
