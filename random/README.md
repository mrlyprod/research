---
title: Random
lead: Randomness from designs: the rounds, where every round that writes two cells or more is even and a census reads which splits reach the alternating group; the design ladder, three no-go rungs for generators built from designs alone; the mess, the nonlinear step, whose Sources all sit below the nonlinearity ceiling past five neighbours; four chapters on top: the moire limit, the schedules, the mixing and the machine; and four on reading a design: along the squares, along polynomials (the readings), along the primes and the major arcs.
figure: research-random
slug: random
date: 2026-10-26
---

# Random

Randomness from designs: the rounds, where every round that writes two cells or more is even and a census reads which splits reach the alternating group; the design ladder, three no-go rungs for generators built from designs alone; the mess, the nonlinear step, whose Sources all sit below the nonlinearity ceiling past five neighbours; four chapters on top: the moire limit, the schedules, the mixing and the machine; and four on reading a design: along the squares, along polynomials (the readings), along the primes and the major arcs.

- [notes/random.md](notes/random.md): Random
- [notes/design-ladder.md](notes/design-ladder.md): The design ladder
- [notes/machine.md](notes/machine.md): The machine
- [notes/major-arcs.md](notes/major-arcs.md): The major arcs
- [notes/mess.md](notes/mess.md): The mess
- [notes/mixing.md](notes/mixing.md): The mixing
- [notes/pascal-images.md](notes/pascal-images.md): Pascal images
- [notes/primes.md](notes/primes.md): Along the primes
- [notes/readings.md](notes/readings.md): The readings
- [notes/rounds.md](notes/rounds.md): The rounds
- [notes/schedules.md](notes/schedules.md): The schedules
- [notes/squares.md](notes/squares.md): Along the squares
- [notes/two-bases.md](notes/two-bases.md): The moire limit
- [claims/design-ladder.md](claims/design-ladder.md): The design ladder
- [claims/machine.md](claims/machine.md): The machine
- [claims/major-arcs.md](claims/major-arcs.md): The major arcs
- [claims/mess.md](claims/mess.md): The mess
- [claims/mixing.md](claims/mixing.md): The mixing
- [claims/primes.md](claims/primes.md): Along the primes
- [claims/readings.md](claims/readings.md): The readings
- [claims/rounds.md](claims/rounds.md): The rounds
- [claims/schedules.md](claims/schedules.md): The schedules
- [claims/squares.md](claims/squares.md): Along the squares
- [claims/two-bases.md](claims/two-bases.md): The moire limit
- [papers/design-generators.md](papers/design-generators.md): Necessary Exits from the Automatic Class: A Second Base or an Exercised Nonlinearity
- [papers/design-rounds.md](papers/design-rounds.md): Design Rounds on the Torus
- [papers/polynomial-readings.md](papers/polynomial-readings.md): Normal Is Not Random: Digit Designs Read Along Polynomials
- [refs.md](refs.md): the references this topic cites
- [lab/py/mess/](lab/py/mess/): The independent half of the mess study
- [lab/rs/design-ladder/](lab/rs/design-ladder/): The generator of the design ladder
- [lab/rs/machine/](lab/rs/machine/): The counter-mode machine against Rule 30
- [lab/rs/major-arcs/](lab/rs/major-arcs/): The major arcs summed by orbit trees against the census along the squares and the cubes
- [lab/rs/mess/](lab/rs/mess/): The mess as a Boolean function
- [lab/rs/mixing/](lab/rs/mixing/): The pairwise mixing of design rounds
- [lab/rs/primes/](lab/rs/primes/): The generator of the chapter Along the primes
- [lab/rs/readings/](lab/rs/readings/): The generator of the readings
- [lab/rs/rounds/](lab/rs/rounds/): The round census
- [lab/rs/schedules/](lab/rs/schedules/): The schedule study: block closures, local conditions and long prime cycles
- [lab/rs/squares/](lab/rs/squares/): The generator of Along the squares: designs along n^2 and the primes, certified Fourier exponents
- [lab/rs/two-bases/](lab/rs/two-bases/): The generator of the moire limit
- [figures/paper-design-generators.ts](figures/paper-design-generators.ts): paper-design-generators
- [figures/paper-design-rounds.ts](figures/paper-design-rounds.ts): paper-design-rounds
- [figures/paper-polynomial-readings.ts](figures/paper-polynomial-readings.ts): paper-polynomial-readings
- [figures/research-random-design-ladder.ts](figures/research-random-design-ladder.ts): research-random-design-ladder
- [figures/research-random-machine.ts](figures/research-random-machine.ts): research-random-machine
- [figures/research-random-major-arcs.ts](figures/research-random-major-arcs.ts): research-random-major-arcs
- [figures/research-random-mess.ts](figures/research-random-mess.ts): research-random-mess
- [figures/research-random-mixing.ts](figures/research-random-mixing.ts): research-random-mixing
- [figures/research-random-primes.ts](figures/research-random-primes.ts): research-random-primes
- [figures/research-random-readings.ts](figures/research-random-readings.ts): research-random-readings
- [figures/research-random-rounds.ts](figures/research-random-rounds.ts): research-random-rounds
- [figures/research-random-schedules.ts](figures/research-random-schedules.ts): research-random-schedules
- [figures/research-random-squares.ts](figures/research-random-squares.ts): research-random-squares
- [figures/research-random-two-bases.ts](figures/research-random-two-bases.ts): research-random-two-bases
- [figures/research-random.ts](figures/research-random.ts): research-random

## Lost

- One round's spectral gap as the mixing rate: every round that writes two cells or more has `sigma_2(P) = 1`, a difference on the fill passing as a point mass, and the guess that this marks an imprimitive group is refuted by all 147 alternating `3 x 3` cells; the rate is `abs lambda_2` of the chain; the rows in [claims/mixing.md](claims/mixing.md).
- A local spread for fast mixing: the light cone bounds the distance below until the cone covers the torus, so rounds grow with the side; the row in [claims/mixing.md](claims/mixing.md).
- Local conditions alone for a primitive group: a strongly connected reads digraph, spanning second differences and no affine direction, refuted by the `4 x 4` single rounds `0x0055`; the row in [claims/schedules.md](claims/schedules.md).
- An exercised nonlinearity at unbounded depth as an exit from the automatic class: refuted by rule 182 from one cell; the row in [claims/design-ladder.md](claims/design-ladder.md).
- The machine's 2.77 and 5.69 over untransported Rule 30 (5.69 at an earlier `R*`, 4.675 at the current one) as a mixing advantage: it is the grid's light cone against the ring's, and with the same slide Rule 30 comes within 1.06 and 1.18 at the antipodal load; the rows in [claims/machine.md](claims/machine.md).
- One spectral invariant as the ranking of designs by how far from random they read along the squares: the mirror pair `6:5`, `6:12345` has one spectral measure, yet along `n^2`, `n < 2^30`, it reads `X2 = 13152` against `714`; the row in [claims/squares.md](claims/squares.md).
- A polynomial reading with no near-head among the battery pairs as one that passes the battery: at `N = 2^26` Thue-Morse along `n^2 + n` fails the dilation test at `z = 14.42` against the control's `2.81`, and Rudin-Shapiro along `n^3 + n`, though white, fails its kernel test at `z = 17.10`; the row in [claims/readings.md](claims/readings.md).
- The `q = 3` term alone as the main term of `sum_(n<x) t(n) Lambda(n)`: it misses by `31.68 sigma` at `2^31`, where `T_6` cancels a third of it, and the integers prime to 6 still miss by `4.08 sigma` at `2^32`; the row in [claims/primes.md](claims/primes.md).
