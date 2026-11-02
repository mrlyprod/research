---
title: Random
lead: Randomness from designs: the rounds, all even once they write two cells or more, and a census of which splits reach the alternating group; the design ladder, three no-go rungs for generators built from designs alone; the mess, every Source below the nonlinearity ceiling past five neighbours; and nine chapters: pascal images, moire limit, schedules, mixing, machine, squares, readings, primes, major arcs.
figure: research-random
slug: random
---

A generator on a grid of `N` bits permutes the `2^N` grids, so two questions precede any statistical test, asked of every design split of the `3 x 3` and `4 x 4` torus: which permutations its rounds can make, and whether one round looks random.

A **round** acts on the torus `Z_n x Z_n`, `N = n^2` cells and `M = 2^N` grids, through a **split**: a cell set `A`, the fill or written side, with complement `B`, the void or read side. The **mess** writes `A` from `B`, each cell of `A` counting the live cells of `B` in a design neighbourhood and writing `[count in S]` for a Source `S`; the round is `R = t_c L F`, the Feistel step `F(x) = x + mess(x)`, which changes only `A` and reads only `B`, a spread `L` that moves every bit by a torus map, and a stamp `t_c(x) = x + c` by a design. [The rounds](rounds.md) proves every round a bijection and, once it writes two cells or more, an even permutation, and reads off a census of the `3 x 3` and `4 x 4` catalogs which splits reach the alternating group with every stamp and why the others miss. The figure runs the level-2 carpet split `0x135f` (carpet, `primes`, `slide`; an alternating cell of 258 cycles, `lab/rs/rounds` verb `word 7 7`) for 64 rounds from one live cell.

Each chapter names its generator, a `lab/rs` study that prints every number of its headline below.

## Headline results of the chapters

- **The rounds.** Proved: every design round is a bijection and, once it writes two cells or more, an even permutation, so rounds and stamps generate a subgroup of `A_M` unless a round writes one cell with a mess of odd weight; five design laws, four block systems and a product action, each force `G != A_M`. Verified: the group of a round and every stamp is `A_512` in 147 of 3240 nonlinear `3 x 3` cells and `A_65536` in 198 of 2142 at `4 x 4`, a law explains every other nonlinear cell, and a single round averages 136.43 cycles at `3 x 3` against `H_512 = 6.8165`. Conjecture: on every torus a nonlinear round generates `A_M` with every stamp exactly when no law holds.
- **The moire limit.** Proved: Thue-Morse correlations decay in binary runs, `abs C_t(h) <= 1.0093 . 0.7177^(r_2(h))`, so the moire XOR has `C_z(h) -> 0` and a Rajchman spectral measure though neither factor is, at the effective rate `abs C_z(h) <= 1.0093 . 3^(-G(h)/4)` for `h > 25`, `G(h)` of order `log log h/log log log h`, its constant not computed. Verified: the exact `max abs C_z` on every `[2^k, 2^(k+1))`, `k <= 24`. Conjecture: `sum_h C_z(h)^2` is finite.
- **The schedules.** Proved: a trap inside the written side that no mess breaks; a strongly connected second-difference graph with isolated live-read cells gives a primitive group on any torus, and a prime cycle in `(M/2, M-3]` makes a transitive group alternating or symmetric. Verified: `A_(2^25)` for one `5 x 5` pair. Refuted: strongly connected reads, spanning second differences and no affine direction do not force primitivity. Conjecture: they do for pairs and onto masks.
- **The mixing.** Proved: with uniform stamps the input difference is a doubly stochastic chain, irreducible exactly when the group is 2-transitive; a round that writes two cells or more has `sigma_2 = 1`; a local schedule needs `Omega(n)` rounds. Verified: a hyperbolic spread on the board covers the torus 4 to 10 rounds after the first at `L = 2..6`, where `shift` takes 5 to 65. Refuted: `sigma_2 = 1` exactly when imprimitive. Conjecture: `O(L + log 1/eps)` for that spread.
- **The machine.** Verified: at the antipodal load, given the same slide, Rule 30 comes within 1.06 and 1.18 of the slide machine on neighbour reads at sides 16 and 32, the torus-over-ring cone advantage is 2.77 and 4.675, and a three-offset XOR spread beats Rule 30 + slide by 1.70 and 2.40, winning at every load tested; ChaCha8 stays faster per byte. Proved: `L^n = I` over `F_2`, so every such round is a bijection. Refuted: the ring ratios 2.77 and 5.69 as a mixing advantage.
- **Along the squares.** Proved: the Fourier exponent `eta` of all 57 designs, bracketed within `7.80e-5`. Verified: the 30 statistics that light most fit `N^(-2 eta)` within `0.05`. Refuted: a spectral invariant ranking the designs, since the mirror pair `6:5`, `6:12345` shares one measure and reads `X2 = 13152` against `714`.
- **The readings.** Proved: a monomial reading keeps its design's dilation identity, so normal `t(n^2)` returns exactly `1`. Refuted: a reading with no near-head passing the battery, `t` along `n^2 + n` failing at `N = 2^26` with `z = 14.42` against `2.81`. Conjecture: a free Thue-Morse pair with linear remainder of slope divisible by 3 decays as `N^(-0.415)`.
- **Along the primes.** Proved: for Thue-Morse at `2^k` the `q | 6` part of the major-arc sum is `-3^floor(k/2)`. Verified: the sum cut at `q <= 1024` meets the census of the primes below `2^32` within `2.35 sigma` for all 55 distinct designs. Refuted: the `q = 3` term alone as the main term, `31.68 sigma` off at `2^31`.
- **The major arcs.** Proved: the major arcs cut to the denominators dividing `Q` are the local model `rho_Q(m)`, the tree of an orbit point `c/q` sums to the Gauss sum `G_k(-c, q)/q` times the design twisted by `c/q` against the `B`-adic density, a zero Gauss sum kills its whole tree, and for Thue-Morse at `1/3` the coefficient is `(4/3) e(-1/3)`. Verified: along `n^2`, `n < 2^32`, the observed bias has the model's sign in all 180 cells where the model exceeds 2 in `z`, from 21 designs; along `n^3`, `n < 2^33`, no lag passes `3.91` in `z`, and `3:12` holds a bias the model misses. Conjecture: the base-36 coefficients of the mirror pair `6:5`, `6:12345` are `3.6656` and `0.4377`, so the mirror split is a factor and not a zero, and the period-2 model is the leading term of the bias along `n^2` for every design with a nonzero coefficient on an arc below `1/4`.

## Chapters

- [The rounds](rounds.md): Every round that writes two cells or more is even; a census of the `3 x 3` and `4 x 4` catalogs reads which splits reach the alternating group, and five laws explain every nonlinear miss.
- [The design ladder](design-ladder.md): Three no-go rungs for generators built from designs alone: one base, two bases, and rounds. A polynomial reading of degree 2 or more leaves the automatic class but not the battery: centred monomials keep exact identities, and only Rudin-Shapiro along `n^2 + n` remains a conjectural fourth rung.
- [Pascal images](pascal-images.md): Rules 22, 126, 129, 161 and 182 from one cell as Pascal images, with the paper's Conjecture 9.4 closed.
- [The moire limit](two-bases.md): Thue-Morse times the ternary parity: correlations that tend to zero.
- [The mess](mess.md): The nonlinear step as a symmetric Boolean function: past five neighbours every Source sits below the nonlinearity ceiling.
- [The schedules](schedules.md): Which schedules generate a primitive group, for every `n`.
- [The mixing](mixing.md): Rounds to pairwise independence.
- [The machine](machine.md): The counter-mode machine against Rule 30 given the same slide.
- [Along the squares](squares.md): A digit-parity design read along the squares turns random at a rate its spectrum fixes, with a strength it does not.
- [The readings](readings.md): A normal reading of a design can still fail the battery.
- [Along the primes](primes.md): A design read along the primes is biased by its major-arc sum.
- [The major arcs](major-arcs.md): The circle method's major arcs, summed along orbit trees, predict a design's bias along the squares.

## Papers

- [Necessary Exits from the Automatic Class: A Second Base or an Exercised Nonlinearity](../papers/design-generators.md)
- [Design Rounds on the Torus](../papers/design-rounds.md)
- [Normal Is Not Random: Digit Designs Read Along Polynomials](../papers/polynomial-readings.md)

## Open

- **Rounds to random.** [Mixing](mixing.md) measures the rounds to pairwise independence; the rounds until a composite's cycle type matches `H_M`, `0.6243 M` and one fixed point stay open.
- **A complete law list.** Law IV was needed only at `4 x 4`; whether five laws close the nonlinear cells at `5 x 5`, where [schedules](schedules.md) proves 7 of 81 line pairs primitive, is open.
- **The converse as a theorem.** A proof that no law means `A_M` for nonlinear cells in both catalogs would run through O'Nan-Scott, like Caranti, Dalla Volta and Sala 2009b for AES-like ciphers.
