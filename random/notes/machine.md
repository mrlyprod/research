---
title: The machine
lead: The counter-mode design machine against Rule 30 given the same slide, in one harness under a nine-test battery; the slide machine beats it on neighbour reads by only 1.06 at side 16 and 1.18 at side 32, a three-offset XOR spread, invertible because L^n = I over F_2, beats it by 1.70 and 2.40 and wins at every load tested, and ChaCha8 stays the fastest per byte.
figure: research-random-machine
---

The rounds and the mess price one step of the design machine. This page runs the machine as a generator. A seed fills a torus, a 32-bit counter is written into it, rounds run, and a raw block of 64 cells is read out, so each counter gives one 64-bit output. The question is the one every counter-mode generator is asked: after how many rounds does the map from counter to output pass a fixed battery of tests, and what does that cost?

An earlier version of this page measured the machine against Rule 30 on a ring and reported ratios of 2.77 and 5.69 in the machine's favour. The 5.69 is `561 * 3/(8 * 37)`: the ring's 561 steps at side 32 against the slide machine's earlier six-test `R*` of 37, which is 45 under the nine tests and gives 4.675. Those ratios measure dimension, not mixing. A ring cell reads one neighbour each side per step, so a counter half a ring away needs at least `n^2/2` steps to reach the output, while the machine's spread carries the whole grid. The fair baseline gives Rule 30 the machine's transport. Here that is **Rule 30 + slide**: the Rule 30 step on the row-major line of the torus, then the machine's slide, with the machine's seed, counter cells and output block. The battery is also corrected: the full Walsh spectrum gets a bar calibrated for its `2^26` coefficients, and rank runs at 512 and 1024 so that a linear complexity of 256 shows up.

The answer has four parts.

- The slide machine beats Rule 30 + slide on neighbour reads by 1.06 at side 16 and 1.18 at side 32 under the antipodal load. That is a tie at 16 and a narrow win at 32.
- The reason is diagnosed. The slide moves one active cell to one active cell, so all growth of a difference comes from the Feistel step. That step gives about one new active cell per active cell placed uniformly on the grid, and growth is held by the two-dimensional cone.
- The change is a spread that XORs three translates of the grid. It is invertible for every triple of distinct offsets, because `L^n = I` over `F_2` when `n` is a power of 2. Over an enumerated family of 84 triples, the first four, in enumeration order, of eight triples tied at the top of the second screen all beat Rule 30 + slide at both sides. The one taken for the loads, the branch-3 machine, wins by 1.70 at side 16 and 2.40 at side 32. It wins at every load tested, also when the spread's own reads are counted.
- ChaCha8 in the same harness is faster per byte than every design machine.

One generator, `lab/rs/machine`, prints every number below. Each verb runs single-threaded, from under a second to under four minutes.

## The harness

The harness is the same for every generator.

- **The grid.** The torus `Z_n x Z_n` with `n = 16` (256 bits) or `n = 32` (1024 bits). Rule 30 untransported runs on a ring of `n^2` cells, so every generator carries the same state.
- **The seed.** splitmix64 from a fixed seed fills every cell. This expander sits outside every machine and is identical for all of them.
- **The counter.** A 32-bit counter is XORed into a 4 x 8 block. There are three loads. The **antipode** puts the block within two cells of the output centre plus `(n/2, n/2)`, as near 16 of its cells in the split as the window allows. **Upstream** puts its corner at `(n - 4, n - 8)`, just behind the output along the slide. **Below** puts it at `(n/2 + 2, 0)`, half a torus under the output. Rule 30 + slide uses the same cells. Rule 30 on the ring takes cells `n^2/2 + 16` to `n^2/2 + 47` (antipodal) or the last 32 cells (near).
- **The output.** The raw 8 x 8 block at the origin of the torus, or the first 64 ring cells: 64 bits per counter.
- **ChaCha8.** The RFC 8439 block function at 8 rounds, with key and nonce from splitmix64 on the seed. The counter goes in word 12, and the output is words 0 and 1: 32 bits in, 64 bits out, as for every other generator.
- **Bit-slicing.** One `u64` holds one cell for 64 consecutive counters, so a round is a fixed program of word operations, run once per 64 counters. ChaCha8 runs scalar, one counter at a time.

A machine round is `R = t_c L F`. The Feistel step `F` writes each cell `a` of the split with `x_a + g(count)`, the count over the live cells of the read side in its neighbourhood. The spread `L` moves bits. The stamp `t_c` flips a design. The stamp of round `t` is the split translated by `(t mod n, floor(t/n) mod n)`, so no two rounds below `n^2` share a stamp. Rule 30 writes `x_(p-1) + (x_p OR x_(p+1))` at every cell.

**The ChaCha in the crate is the RFC 8439 block function. Verified.** The crate's tests reproduce the quarter-round vector of section 2.1.1 and the ChaCha20 block of section 2.3.2. A second implementation, written in Python, reproduces the same block. It also agrees with the crate's ChaCha8 on all 8 counters the `chacha` verb prints, from 0 to `2^32 - 1`.

## The battery

Nine tests read the map from counter to output at a given round. Each has its sample size and bar stated in advance. Where the false-alarm rate of an ideal function can be computed, the bar sits near `4.5` standard deviations of the whole test.

- **Monobit.** All `2^26` output bits of counters `0..2^20`, `|z| < 4.5`.
- **Avalanche.** For `2^14` base counters and each of the 32 counter bits, the fraction of base counters whose output bit `j` flips when that counter bit flips. The bar is `4/2^7` on the largest deviation from `1/2` over the 2048 pairs.
- **Weight-one Walsh.** The correlation of each output bit with the constant and with each of the low 20 counter bits, over `2^20` counters, under `4.5/2^10`.
- **Full Walsh.** All `2^20` Walsh coefficients of each output bit over the low 20 counter bits, `M = 2^26` in all. The bar is `sigma (sqrt(2 ln 2M) + 1) = 7.118/2^10` with `sigma = 2^-10`.
- **Cubes.** For `d = 2..20`, the XOR of each output bit over the `2^d` counters that vary the low `d` counter bits, with the next `20 - d` at each value. Each dimension's count of ones has `|z| < 4.5`. A map of degree below `d` has every such sum zero.
- **Linear complexity.** Berlekamp-Massey on the first `2^16` stream bits, the 64-bit outputs in counter order, `|L - 2^15| <= 10`.
- **Rank 64.** The 64 x 64 matrices whose rows are 64 consecutive outputs, `2^14` of them. A chi-square on the classes 64, 63, 62 and at most 61 against the exact law must stay under 27.
- **Rank 512 and rank 1024.** Matrices of consecutive stream bits, 256 of side 512 and 64 of side 1024. The chi-square on the classes full rank, one short and at least two short must stay under 24, and no matrix may be more than 10 short.

**`R*` is the first round at which all nine pass and keep passing four more rounds.** Rule 30's step count is `T*`, with the same definition. The six-test `R*`, without the last three tests, is printed beside it. The headline is `C = R* n^2 k/64`, neighbour reads per output bit, with `k = 8` for the carpet neighbourhood and 3 for Rule 30. The strict count adds the spread's reads: one per cell for a map, three for a three-offset sum.

**The battery is calibrated. Verified.** splitmix64 and ChaCha8 in counter mode pass all nine tests at round 1. Monobit reads 1.07 and 0.29, avalanche 0.0199 and 0.0194, weight-one Walsh 3.79 and 3.54, full Walsh 5.84 and 5.76 in units of `2^-10`. Linear complexity is `2^15` and `2^15 + 1`. The rank chi-squares are 0.23 and 2.72 at side 64, 0.02 and 1.76 at 512, and 3.41 and 0.48 at 1024.

**The full Walsh bar. Verified.** The largest of `2^26` independent `|N(0, 1)|` has median 5.73 and exceeds 7.118 with probability `7.35e-5`. The earlier page set `4.5/2^10` for one coefficient and then had to drop the full spectrum from the gate, since splitmix64's maximum is 5.84. With the bar calibrated for the count, the full spectrum is back in the gate. It matters: it moved `R*` in 12 of the 21 machine runs on this page, by 1 to 8 rounds.

**The low bit of xoshiro256+ fails rank at 512 and 1024. Proved.** Bit 0 of `s0 + s3` is the XOR of the two low bits, and the state transition is linear over `F_2`. So the low-bit stream has linear complexity at most 256. Every row of 512 or 1024 consecutive bits is then a linear image of a 256-bit state, which caps every such matrix at rank 256. Measured: Berlekamp-Massey reads exactly 256, and all 256 matrices of side 512 and all 64 of side 1024 have rank exactly 256. The chi-squares are 1659.66 and 414.91. **That it fails the 64 x 64 test is Refuted:** its classes there are `4721, 9453, 2111, 99`, chi-square 1.84, because a recurrence of order 256 is longer than a 64-bit row.

## The catalog

The design cells are enumerated, never picked.

- **Splits.** The level-`L` design of each of the 14 base-2 letters, the Kronecker power `k^(x)L` at side `2^L`, and its complement: 28 splits.
- **Neighbourhoods.** Carpet, net and void.
- **Spreads.** The six maps `id`, `shift`, `drop`, `slide`, `turn` and `flip`.
- **Rules.** At each neighbourhood size `k`, the distinct Source sets of the mess census, united with every rule whose live and dead cases both differ from the parity on about half the counts. At `k = 8` that gives 14 Source sets and 252 rules at parity distance 4 or 5, 259 in the union. At `k = 4` it gives 24.

That is 504 geometries and 51576 cells at side 16. **456 of the 504 geometries have a structural cone that misses a (counter bit, output bit) pair through round 64. Verified.** The cone propagates, for each counter bit, the set of cells that can depend on it through every read and every spread, so it bounds every rule from above. An avalanche screen over 64 counters, then a second screen over `2^14` counters reading monobit, weight-one Walsh and rank, leave 30 cells. 10 of the 30 were shortlisted and run through the nine tests; the other 20 stop at the screens. The best of the 10 is the level-4 power of letter 14, the carpet neighbourhood, the slide spread and the mess `g = {1, 3, 6, 8}`, with `R*` 24. This page calls it the slide machine.

## The headline

Three lines carry the result. Lines (b) and (c) are printed by the `headline` verb; line (a) is printed by the `speed` verb, which only times.

- **(a) Time.** At `R*`, in ns per output byte, median of 7 runs with min to max on a shared laptop: the slide machine 25.4 (21.2 to 54.2), the branch-3 machine 15.6 (15.5 to 15.6) and Rule 30 + slide 26.8 (26.7 to 26.9) at side 16; 128.1 (127.9 to 128.5), 77.0 (76.8 to 77.3) and 215.7 (214.9 to 227.4) at side 32; ChaCha8 scalar 4.5 (4.5 to 4.5).
- **(b) Neighbour reads.** `C(Rule 30 + slide)/C(machine)` at the antipodal load is 1.06 and 1.18 for the slide machine and 1.70 and 2.40 for the branch-3 machine, at sides 16 and 32.
- **(c) Cone advantage, dimensional.** `C(Rule 30 untransported)/C(slide machine)` is 2.77 and 4.675 (`26928/5760`, printed as 4.67 by the verb) under the antipodal ring load and 1.52 and 0.81 under the near ring load.

Line (c) is not the result. It is the size of the advantage the geometry gives away: a 2-D torus with a moving frame against a 1-D ring whose counter is half a ring away. The result is line (b), and the win that counts is the branch-3 machine's.

| side | generator | `R*` at antipode, upstream, below | `C` at the antipode | ns/byte at `R*` |
|---|---|---|---|---|
| 16 | slide machine | 24, 24, 24 | 768 | 25.4 (21.2 to 54.2) |
| 16 | branch-3 machine | 15, 16, 15 | 480 | 15.6 (15.5 to 15.6) |
| 16 | Rule 30 + slide | 68, 84, 83 | 816 | 26.8 (26.7 to 26.9) |
| 16 | Rule 30 on the ring | 177 (antipodal), 97 (near) | 2124 | |
| 32 | slide machine | 45, 37, 49 | 5760 | 128.1 (127.9 to 128.5) |
| 32 | branch-3 machine | 22, 17, 22 | 2816 | 77.0 (76.8 to 77.3) |
| 32 | Rule 30 + slide | 141, 133, 142 | 6768 | 215.7 (214.9 to 227.4) |
| 32 | Rule 30 on the ring | 561 (antipodal), 97 (near) | 26928 | |
| any | ChaCha8 | 1 | | 4.5 (4.5 to 4.5) |

**Rule 30 + slide in the identical harness has `T*` 68, 84, 83 at side 16 and 141, 133, 142 at side 32. Verified.** Its cones are 40, 39, 41 and 81, 71, 81. With the slide, Rule 30 is no longer limited by the ring distance, but its step still spreads a difference along one line, and the output block spans 8 rows of it.

**Against Rule 30 + slide the slide machine wins by 1.06, 1.31, 1.30 at side 16 and 1.18, 1.35, 1.09 at side 32, at the antipode, upstream and below loads. Verified.** Counting each spread's reads, one per cell for both, the ratios are 1.26, 1.56, 1.54 and 1.39, 1.60, 1.29. **The reading of 2.77 and 5.69 as a mixing advantage is Refuted.** Given the same slide, Rule 30 comes within 6 percent of the slide machine at side 16 and 18 percent at side 32.

**The branch-3 machine wins by 1.70, 1.97, 2.08 at side 16 and 2.40, 2.93, 2.42 at side 32. Verified.** With the spread's three reads per cell counted, the ratios are 1.65, 1.91, 2.01 and 2.33, 2.84, 2.35. It beats Rule 30 + slide at every load tested at both sides, so this headline survives the counter being moved.

On time, each number is the median of 7 runs printed by `speed` on a shared laptop, with its min to max in the table; the slide machine at side 16 ranges from 21.2 to 54.2 ns per byte. The branch-3 machine is the fastest design machine at both sides. It is faster than Rule 30 + slide by 1.7 at side 16 and 2.8 at side 32. It is slower than ChaCha8 by 3.5 at side 16 and by 17 at side 32. ChaCha8 here computes a 512-bit block and the harness keeps 64 bits of it, so the full block would be cheaper still per byte. No design machine on this page is a faster generator than ChaCha8.

## Why the slide machine was slow

**At side 32 one Feistel step gives an active cell placed uniformly on the grid (an active read-side cell gets 1.30) `beta = 0.991` new active cells on average. Verified.** The count is exact. A written cell that reads `m` live-side cells flips when one input flips with probability `P(g(c) != g(c + 1))`, with `c` binomial on the other `m - 1` inputs. Summed over the 243 written cells and their 1048 reads and divided by the 1024 cells, that gives 0.991. The slide is a permutation: it moves one active cell to one active cell and adds nothing. So in the sparse limit the active count multiplies by `1 + beta = 1.991` per round, with almost all of it from the Feistel step.

The measurement is the mean active-cell count of a single-bit counter difference, averaged over the 32 counter bits and 1024 counters, round by round at side 32:

- **Slide machine:** 1.0, 1.7, 2.6, 4.1, 6.3, 9.4, 15.1, 21.4, 27.6, 34.3, 40.8, 49.6, 59.1, 69.7, 81.5, 94.0, 108.7, 124.5, 142.2, 161.0, 182.2, 204.7, 228.5, 252.7, 280.1 at rounds 0 to 24.
- **Rule 30 + slide:** 1.0, 2.0, 2.9, 3.4, 4.2, 4.7, 5.4, 5.8, 7.1, 9.2, 11.2, 14.1, 16.6, 18.8, 21.8, 24.8, 26.8, 29.8, 33.3, 36.3, 40.1, 43.9, 46.0, 50.5, 54.0.

The slide machine's multiplier starts at 1.70, 1.54 and 1.58, near the prediction less the overlap of new cells, and falls to 1.11 by round 24. The cause is the cone. In the frame that moves with the slide, the Feistel step reads only the eight neighbours, so after `r` rounds a difference lies within Chebyshev distance `r` of where it started: at most `(2r + 1)^2` cells. Growth is quadratic, not exponential, and at round 24 the difference covers 280 of the 1024 cells. Rule 30 + slide is worse still. Its step spreads along one line, so its count grows by about 2 per round, linear in `r`.

So the slide machine pays twice. Its cone is narrow, one cell per round in each direction. And its fill inside the cone is thin, since only written cells next to an active read cell can become active. A spread that writes each bit to three places attacks both at once.

## The branch-3 spread

The change is the spread `L(x) = rot(x, v0) + rot(x, v1) + rot(x, v2)`, the XOR of three translates of the grid by distinct torus offsets. As a bit program, cell `p` becomes `x(p - v0) + x(p - v1) + x(p - v2)`. It maps one active cell to three.

**For distinct offsets on the torus of side `n = 2^L`, `L^n = I`, so the round is a bijection. Proved.** The translations commute, and over `F_2` the square of a sum is the sum of the squares, since every cross term appears twice. So `L^2 = T_(2 v0) + T_(2 v1) + T_(2 v2)`, and after `L` squarings `L^n = T_(n v0) + T_(n v1) + T_(n v2) = 3 I = I`. Hence `L^-1 = L^(n-1)`. The Feistel step is an involution and the stamp is its own inverse, so `R^-1 = F L^(n-1) t_c`, and the round is a bijection for every split, mess and stamp. The same holds for any odd number of distinct offsets. In the group ring of a 2-group over `F_2`, an element is a unit exactly when its number of terms is odd.

**For the winning spread the `n^2 x n^2` matrix has rank 256 at side 16 and 1024 at side 32. Verified.** On a random grid `L^n = I`, and 8 rounds undone by `F L^(n-1) t_c` return all 64 random states.

### The family

**The family is every triple of distinct offsets in the box `{-1, 0, 1}^2`: 84 triples at each side, none excluded. Verified.** A triple is excluded when the difference of two of its offsets is a translation symmetry of the split, since then two of its terms would carry the same design. Letter 14's gasket has no nontrivial translation symmetry at either side, so none is excluded. Split, neighbourhood and mess are those of the slide machine, so the family varies the spread and nothing else.

The same two screens run at both sides. At side 16 one triple is structurally blind and 64 pass the second screen. At side 32 none is blind and 64 pass. 57 pass the second screen at both sides. They are ranked by the worse side of `C(Rule 30 + slide)/C`, and the first four, in enumeration order, of eight triples tied at the top of the second screen go through the nine tests at both sides:

| spread | `R*` at 16 | `R*` at 32 | ratio at 16 | ratio at 32 |
|---|---|---|---|---|
| `(-1,-1)(-1,0)(1,1)` | 15 | 22 | 1.70 | 2.40 |
| `(-1,-1)(0,-1)(1,1)` | 15 | 23 | 1.70 | 2.30 |
| `(-1,0)(-1,1)(1,-1)` | 15 | 21 | 1.70 | 2.52 |
| `(-1,0)(-1,1)(1,1)` | 17 | 21 | 1.50 | 2.52 |

**All four beat Rule 30 + slide at both sides. Verified.** The first three tie on the worse side at 1.70, and enumeration order took the first as the branch-3 machine for the loads and timings. The third is better at side 32. No one triple is the result: the family is. Under the six-test battery the four read 13, 14, 13 and 15 at side 16. The full Walsh spectrum failed in the rounds between, which shows the corrected test doing work right at `R*`.

The symmetry check: of the `8 n^2` torus symmetries, 2 fix the split, the identity and the transpose, but only the identity commutes with the winning spread, at both sides.

### The growth law

**At side 32 the branch-3 machine's predicted sparse multiplier is `3 (1 + beta) = 5.973`. Verified.** Its measured active count at rounds 0 to 24 is 1.0, 4.3, 8.6, 20.4, 29.7, 53.7, 84.5, 122.8, 161.5, 206.4, 255.0, 304.9, 359.9, 412.5, 454.5, 478.9, 489.9, 495.9, 500.5, 504.0, 506.9, 509.2, 511.0, 511.7, 511.9. A uniformly random difference has mean 512. The first multiplier is 4.30, then 1.99, 2.37, 1.46, 1.81. It stands at 255.0 at round 10, about half the random mean 512, and at 511.0 at its `R*` of 22.

The oscillation is the linear part. **`L^(n/2)` is the translation by `(n/2, 0)` for the winning spread. Proved.** By the squaring rule, `L^(n/2) = T_((n/2) v0) + T_((n/2) v1) + T_((n/2) v2)`, and `(n/2)(-1,-1)` and `(n/2)(1,1)` are the same point, so those two terms cancel. More generally, `L^t` is the product of `L^(2^j)` over the binary digits of `t`, each a sum of three translates, so a single cell spreads to at most `3^(s(t))` cells, `s(t)` the binary digit sum. **The impulse `|L^t e|` is exactly `3^(s(t))` for every `1 <= t < n/2` at sides 16 and 32. Verified.** It reads 3, 3, 9, 3, 9, 9, 27, then 1 at `t = 8` on side 16. So the linear part alone collapses at every power of 2. That is why the active count jumps in odd rounds and stalls in even ones. The Feistel step's new cells, injected every round, are what carry the count through the collapses.

The law that explains the win has two parts. The sparse phase multiplies by about `terms x (1 + beta)`: 1.99 for the slide and 5.97 for the sum, measured 1.70 and 4.30 in the first round. The dense phase is held by the cone. The sum's offsets and the reads each move a difference up to one cell per round, so its cone widens by two cells per round against one for the slide. Measured at side 32, the branch-3 difference reaches the random mean, 511.0 of 1024 cells, by round 22, its `R*`. The slide machine's covers only 280 cells at round 24, and its `R*` is 45. For the branch-3 machine, `R*` sits where the difference saturates, four rounds after its cone of 18.

## Falsification

The plan was three loads that must agree, two bad loads, a fold measured apart, and a symmetry check.

- **Three loads.** At side 16 they agree within a round: 24, 24, 24 for the slide machine and 15, 16, 15 for the branch-3 machine. At side 32 they do not: 45, 37, 49 for the slide machine and 22, 17, 22 for the branch-3 machine. **That `R*` does not depend on the counter load is Refuted at side 32.** For the branch-3 machine the upstream load starts the difference closer to the output, its cone is 12 against 18, and `R*` drops by 5. So part of its `R*` at the antipode is transport. For the slide machine, the earlier six-test reading was 37 at every load. Under the full Walsh test the antipode fails through round 44 and the below load through round 48.
- **The headline survives.** Rule 30 + slide is moved with the counter, and the ratio at each load uses the same cells for both generators. The branch-3 machine wins at all six (side, load) pairs, by 1.70 to 2.93. The slide machine wins narrowly at all six, by 1.06 to 1.35.
- **Bad loads.** For the slide machine at side 16, the all-zero fill gives `R*` 25 and the counter inside the output block gives 28, both unchanged under the nine tests. A counter bit written into the output appears there linearly at round 0, and it takes four extra rounds to wash out.
- **The fold.** XOR each output cell with the cell `(8, 0)` away, a pairing of width 2 that is no symmetry of the split and touches no counter cell. `R*` falls from 24 to 23: a fold gain of 1, where the six-test battery read 2.
- **Symmetry.** For the slide machine, 2 of the `8 n^2` torus symmetries fix the split and both commute with the slide; only the identity fixes every stamp through round `4n`, at both sides.

## The cone advantage, dimensional

Rule 30 on the ring keeps its own section, because it explains the earlier ratios.

**Rule 30's `T*` on the ring is its cone plus 2. Verified, in all four runs, unchanged under the nine tests.** Under the antipodal load the cone is 175 at side 16 and 559 at side 32, exactly `n^2/2 + 47`, the rightward ring distance from counter bit 0 to output bit 63. Under the near load it is 95 at both sides, and `T*` is 97. The right edge of a difference moves exactly one cell per step, because the rule is linear in its left neighbour.

**A lower bound on Rule 30 at every even side `n >= 10`. Proved.** A ring cell reads one neighbour each side per step, so after `t` steps an output bit depends only on cells within ring distance `t`. Counter bit `b` sits at cell `n^2/2 + 16 + b`, and both ring distances from there to output cell `b + 16` equal `n^2/2`. Before step `n^2/2` that output bit cannot depend on that counter bit, and avalanche fails. Hence `T* >= n^2/2`, and against any carpet machine cell the ratio is at least `3 n^2/(16 R*)`. At the slide machine's `R*` that is 2.0 at side 16 and 4.26 at side 32, against the measured 2.77 and 4.675.

So the untransported ratio grows like `n^2/R*` for geometric reasons alone, and it says nothing about the step. Under the near load, which removes most of the distance, the ratio is 1.52 at side 16 and 0.81 at side 32: there Rule 30 already beats the slide machine at side 32.

## The mess that won

The mess census found every Source thin: its nonlinearity is held under the count of its ones, and at `k = 8` no Source exceeds 0.55 of the symmetric ceiling 120. **The winning mess `{1, 3, 6, 8}` has full degree 8, 93 ones and nonlinearity 93, 0.7750 of the ceiling, against the primes' degree 6 and 0.3000. Verified.** It is not a Source. It is one of the 252 rules at parity distance 4.

**It is the parity on counts 0 to 4 and the complement of the parity on counts 5 to 8. Proved from its values.** So it is linear on every written cell that reads at most 4 cells of the read side. A written cell of the gasket reads 4.3 of them on average, 346 reads over 81 cells at side 16 and 1048 over 243 at side 32. This is why `beta` is so close to 1: on counts 0 to 4 every single flip changes the output, and only the step from 4 to 5 does not.

The obstruction named in advance was a rule too close to linear. **At the slide geometry the Life rule with primes birth and odds survival, and the linear odds mess, never pass the avalanche screen through round 64. Verified.** **No Feistel mess can be that close to the parity. Proved:** a Feistel mess writes `x_a + g(count)`, so its live case is the complement of its dead case, and it agrees with the parity on exactly half of the (centre, neighbourhood) pairs for every `g`. The Life rule agrees on 476 of 512 and is not invertible on its cell.

## The figure

The figure is a 4 x 4 array of 16 x 16 grids. Each row follows one generator for one seed from counter bit 0 flipped, at four rounds ending at its `R*` at side 16. Row 1 is the slide machine at rounds 6, 12, 18, 24. Row 2 is the branch-3 machine at 4, 8, 11, 15. Row 3 is Rule 30 + slide at 17, 34, 51, 68. Row 4 is Rule 30 on its ring, laid out row by row, at 44, 88, 133, 177. In every grid the blue cells are where the two copies differ outside the output, the orange cells where they differ inside it, and the line cells the output cells that still agree. The branch-3 difference spreads over most of the grid by its second tile. The slide machine's grows as a patch carried along the diagonal. Rule 30 + slide fills its rows one line at a time, and Rule 30 on the ring is a band that reaches the output only at the end. The figure rebuilds the split from its definition and asserts its 81 written cells and 346 reads against the generator.

## What died and what is left

- Dead: the 2.77 and 5.69 ratios as the result. They are the cone advantage of a 2-D torus over a 1-D ring, bounded below by `3 n^2/(16 R*)`. The 5.69 used the six-test `R*` 37 at side 32; with the nine-test 45 it is 4.675.
- Dead: load-independence of `R*` at side 32, for both machines.
- Dead: the full-spectrum Walsh test's exclusion from the gate. With the bar `7.118/2^10` it gates, and it moved `R*` in 12 of the 21 machine runs.
- Dead: the rank failure of xoshiro256+ at 64 x 64. It fails at 512 and 1024, with every matrix at rank 256.
- Dead: the Life rule with primes birth and odds survival as a round.
- Alive: the branch-3 win, 1.70 to 2.93 over Rule 30 + slide at every load tested. The open question is its scaling. Both spreads are local, so at the antipodal load the cone forces `R* = Omega(n)` and `C` grows at least like `n^3`. Rule 30 + slide grew from 68 to 141 as the side doubled, the branch-3 machine from 15 to 22.
- Alive: a non-local spread. The mixing chapter's hyperbolic spread breaks the linear light cone. A hyperbolic map is a permutation, so the doors are a hyperbolic map XORed with translates, or three hyperbolic images.
- Alive: wider families. The box `{-2..2}^2` holds 2300 triples, and the mess, split and neighbourhood can vary with the spread. On this page only the spread varied.
- Dead: the shortlist's six-test readings as a stand-in for the nine-test ones. All 10 shortlisted cells run through the nine tests, and the other nine read `R*` 24, 26, 27, 34, 35, 36, 36, 43 and 51, each equal to its six-test reading. None beats the slide machine's 24; one ties it, the same split, neighbourhood and spread with the complement of its mess.
- Alive: ChaCha8. It stays 3.5 to 17 times faster per byte than the best design machine here. A design machine that beats it would need far fewer rounds or far cheaper ones.
