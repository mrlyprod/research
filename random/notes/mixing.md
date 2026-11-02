---
title: The mixing
lead: With a uniform stamp after every round, the difference of two inputs walks a doubly stochastic chain that is irreducible exactly when the round's group is 2-transitive; one round always has second singular value 1, every catalog spread obeys a light cone of linear depth, and a hyperbolic spread on the board split breaks it.
figure: research-random-mixing
---

A generator built from design rounds is fed a seed and a counter, so the question after the group is how fast two different inputs stop looking related. This page measures it for the rounds of the [random spine](random.md): the same torus grids, splits, neighbourhoods, Sources and spreads, with a fresh uniform stamp after every round. Then the difference of two inputs is a Markov chain on the `2^N - 1` nonzero grids, and the distance from pairwise independence after `r` rounds is the distance of that chain from uniform. The chain is reducible exactly when the group of the round and its stamps is not 2-transitive, so the census of the spine reads off which cells mix at all: the 147 alternating cells at `3 x 3` and the 198 at `4 x 4`, and no other nonlinear cell. Among those, the median cell needs 10 rounds to come within `1/4` of pairwise independence and 22 to come within `0.01` at `3 x 3`; at `4 x 4` the cell of the spine's figure needs 11 and 19.

Three results are new against the natural expectations. The second singular value of one round is always 1, so the textbook spectral bound needs powers of the round. A 2-transitive group does not mix: twelve affine cells are 2-transitive and keep their differences deterministic forever. And depth is geometry: every catalog spread is an isometry of the torus, so a difference spreads at most one cell per round and any local schedule needs a number of rounds linear in the side, while a hyperbolic spread on the board split covers the torus in a number of rounds that grows by one or two per doubling of the side.

One generator, `lab/rs/mixing`, prints every number below; its `dense.py` rechecks the `3 x 3` census by an independent dense route.

## The chain

The round is the spine's `R = t_c h` with `h = L F`: the Feistel step `F(x) = x + mess(x)`, which writes the fill `A` from the void `B`, then the spread `L`. Here the stamp `c` is not the split's design but a fresh uniform element of `V = F_2^N` at every round, a round key. A schedule of `r` rounds with stamps `c_1, ..., c_r` maps an input `x` to `E(x) = t_(c_r) h_r ... t_(c_1) h_1 (x)`.

**The difference of two inputs is a Markov chain. Proved.** Fix inputs `x != x'`, write `x_t, x'_t` after `t` rounds and `d_t = x_t + x'_t`. For a map `h` of `V` and a difference `a`, write `D_a h(y) = h(y) + h(y + a)` and

```
P_h(a, b) = #{y : D_a h(y) = b} / 2^N        (the difference table over 2^N)
```

Then `d_1 = D_(x + x') h_1 (x)` is a fixed nonzero difference, and for `t >= 1`, conditionally on `d_1, ..., d_t`, the next difference has law `P_(h_(t+1))(d_t, .)`. The proof is one induction. `x_t = h_t(x_(t-1)) + c_t` with `c_t` uniform and independent of everything before it, so `x_t` is uniform and independent of `(d_1, ..., d_t)`, because each `d_s` with `s <= t` is a function of `x_(s-1), d_(s-1)` and of no stamp after `c_(t-1)`. Then `d_(t+1) = D_(d_t) h_(t+1) (x_t)` with `x_t` uniform and independent of the history, which is the kernel. This is the Markov-cipher property of Lai, Massey and Murphy 1991, cited in [refs](../refs.md) as unread; nothing here depends on their proof.

**The first round only relabels the start. Proved.** `d_1` is a function of the inputs and is not random, and every nonzero `d_1` occurs: for any `x` take `x' = h_1^(-1)(h_1(x) + d_1)`. The stamps only start working in round 2.

**The distance from pairwise independence. Proved.** After round `r` the output pair is `(x_r, x_r + d_r)` with `x_r` uniform and independent of `d_r`. A uniform pair of distinct grids has the same shape, a uniform first point and an independent uniform nonzero difference, so the total-variation distance between the law of `(E(x), E(x'))` and the uniform pair is the distance of `d_r` from the uniform law `U` on the `2^N - 1` nonzero differences. The worst input pair is the worst start:

```
d(r) = max over a != 0 of TV(e_a P_2 P_3 ... P_r, U),     d(1) = 1 - 1/(2^N - 1).
```

**Every kernel is doubly stochastic, and it need not be symmetric. Proved.** A row sums to 1 by definition. A column `b` sums to `#{(a, y) : h(y) + h(y + a) = b} / 2^N`, and for each `y` exactly one `a` works since `h` is a bijection, so it is 1 too; `b = 0` is never hit from `a != 0`. So `U` is stationary for every round and every schedule, and the period product `Q = P_2 ... P_(k+1)` of a schedule of period `k` is doubly stochastic as well. `P` is the transpose of `P_(h^(-1))`, not of itself.

**The same holds for three inputs, on pairs of differences. Proved.** Three distinct inputs keep the differences `(d_t, e_t)` from the first; after round 1 they are a fixed ordered pair of distinct nonzero differences, every such pair occurs, and the next pair is `(D_(d_t) h(y), D_(e_t) h(y))` with `y` uniform. The 3-wise chain lives on `(2^N - 1)(2^N - 2)` states, 260610 at `3 x 3`; it is doubly stochastic by the same count, and its first coordinate is the 2-wise chain, a lumping: every eigenvalue of `P` is an eigenvalue of the 3-wise kernel `P_3`, and every singular value of `P` is at most the matching one of `P_3`.

## The spectrum

**The singular-value bound. Proved.** For a doubly stochastic `Q` on `n` states the sum-zero vectors are invariant under `Q` and `Q^T`, and `sigma_2(Q)`, the norm of `Q` on them, is the second singular value. A start `e_a` has `|e_a - U|_2 <= 1`, and total variation is half the `l_1` norm, at most `sqrt(n) / 2` times the `l_2` norm. So for a schedule of period `k` with period product `Q`,

```
d(1 + k s) <= (1/2) sqrt(2^N - 1) sigma_2(Q)^s.
```

The rate of the curve is `|lambda_2(Q)|`, the largest modulus of an eigenvalue on the sum-zero space; `Q` is not normal, so `|lambda_2| <= sigma_2` and the two can differ.

**Irreducible means 2-transitive. Proved.** With every stamp allowed the rounds generate `G = <T, h>`, `T` the translations. Since `G` contains `T`, a `G`-orbit of ordered pairs of distinct grids is a union of translation classes, and a class is named by its difference; call `D(a)` the differences in the orbit of `(0, a)`. A path `a -> a_1 -> ... -> b` in the support graph of `P` is a word: at each step translate the pair so that its first point is the `y` the edge needs, then apply `h`. Conversely every element of the finite group `G` is a word in `h` and translations, translations keep the difference, and `h` moves it along an edge. So `D(a)` is the set the chain reaches from `a`, and `G` is 2-transitive, one orbit on distinct pairs, exactly when the chain is irreducible. A block system of `G` is a subspace `U` whose nonzero vectors form a closed class (the spine's subspace criterion), and a primitive group that is not 2-transitive splits into several orbitals: both make the chain reducible.

**A reducible cell never mixes. Proved.** A closed class `C` traps a start in it, so `d(r) >= 1 - |C| / (2^N - 1)` for every `r`. Since `P` is doubly stochastic, every class is closed.

**One round always has second singular value 1. Proved.** A difference supported on the fill, `a_B = 0`, passes deterministically: the mess reads only `B`, so `D_a h(y) = L a` for every `y` and the row of `a` is a point mass. Two such starts `a != a'`, which exist once `|A| >= 2`, give the sum-zero vector `e_a - e_(a')` mapped to `e_(La) - e_(La')` of the same norm, so `sigma_2(P) = 1`. The same count with the cells `p` whose first `j` images `p, sigma(p), ..., sigma^(j-1)(p)` all lie in the fill gives `sigma_2(P^j) = 1` once two such cells exist, and `d(j + 1) = 1 - 1/(2^N - 1)` once one exists: its start `e_p` is the point mass at `sigma^j(p)` after `j` steps, and no law on the nonzero differences is farther from `U`, since it gives some difference at least the uniform weight `1/(2^N - 1)`. Write `k*` for the least `j` with `sigma_2(P^j) < 1`; the bound above is then used with `Q = P^(k*)`.

**So the guess that `sigma_2 = 1` marks the imprimitive cells is false. Refuted.** It is 1 for every design round with two written cells or more, mixing or not; only a power of the round separates them.

**A 2-transitive group need not mix. Proved; Verified.** If the mess is linear the round is affine, `D_a h(y) = M a` for a fixed linear map `M`, and `P` is the permutation matrix of `M` on the nonzero differences: the difference is deterministic forever and `d(r) = 1 - 1/(2^N - 1)` for every `r`. Such a `P` is irreducible exactly when `M` permutes the `2^N - 1` nonzero vectors in one cycle, a Singer cycle, and then `G = T M` is 2-transitive. The census finds 12 such cells at `3 x 3`, the 12 primitive affine cells of the spine: four `odds` cells with the `net` neighbourhood under `shift`, four under `drop` and four with `void` under `slide`. Group-theoretic randomness, 2-transitivity included, says nothing about how fast; the nonlinear step is what makes the difference random.

## The census at 3 x 3

The kernel is built without a dense table. The mess reads only `B`, so `D_a h(y) = L(a + D_(a_B) mess(y_B))` and a row depends on `a` through `a_B` alone; for each `a_B` the law of `mess(y) + mess(y + a_B)` over the `2^|B|` voids is exact from two Walsh-Hadamard transforms, the autocorrelations of the component functions `u . mess`. At `3 x 3` every row of every one of the 4320 cells is compared with the direct table `h(y) + h(y + a)`: 0 rows differ. The round itself is the spine's, reproduced and recounted: the 4320 and 2856 cells, the mean cycle counts 136.43 and 16759.65 of the nonlinear cells, and the figure cell `0x135f` with 258 cycles, longest 24865, 3 fixed points.

**Exactly the alternating cells mix. Verified.** Of the 3240 nonlinear cells, 147 have an irreducible chain and 3093 a reducible one; of the 1080 `odds` cells, 12. By the section above the 3093 are the spine's 3087 imprimitive cells and 6 primitive products, and the 147 are its alternating cells, a second certificate for their 2-transitivity that shares nothing with Jordan's theorem: one search of the support graph. Every reducible nonlinear cell has a closed class of at most 255 differences, so `d(r) >= 1 - 255/511 = 0.500978` at every `r`: half the distance is never removed. The independent rebuild in numpy, which builds every kernel from the direct table and counts strong components, finds the same 147 and 12.

**The distance curves of the 147 alternating cells. Verified.** Every start is run, all 511 rows of `P^(r-1)`, until `d < 10^-6`.

| quantity | min | median | max |
|---|---|---|---|
| rounds to `d < 1/4` | 7 | 10 | 29 |
| rounds to `d < 0.01` | 13 | 22 | 76 |
| `k*` | 3 | 3 | 6 |
| `sigma_2(P^(k*))` | 0.449163 | 0.632930 | 0.829949 |
| `sigma_2(P^(k*))^(1/k*)`, per round | 0.765834 | 0.858834 | 0.963408 |
| `abs lambda_2(P)` | 0.586924 | 0.748207 | 0.931477 |

`d(2) = 1 - 1/511` in every cell, as the fill-supported start predicts, and the curve stays above `0.95` for at least three rounds before it falls geometrically. `k*` is 3 for 135 cells, 5 for six and 6 for six. The twelve slow ones are `binary` cells on the splits of the twelve Singer cells with the neighbourhood swapped: `void` under `shift` and `drop` where the Singer cells have `net`, and `net` under `slide` where they have `void`. They are the twelve slowest `binary` cells, `abs lambda_2` from 0.852 to 0.891 and 38 to 51 rounds to `0.01`. The slowest cell of all is `0x0054` with `net`, `fibonacci` and `shift`, and its `drop` twin, `abs lambda_2 = 0.9315`, 29 rounds to `1/4` and 76 to `0.01`; the fastest, 7 and 13, are twelve `primes` and `binary` cells on the splits `0x0007`, `0x0049` and `0x0054`.

The worst start is a difference supported on the fill or on one cell in 9539 of the 9744 round maxima, the light differences that the Feistel step carries without mixing.

The dense route recomputes every curve by dense matrix powers of the 511 by 511 kernel: on all 147 cells the rounds to `1/4` and to `0.01` agree with the sparse route. On every seventh cell, 21 cells, it also takes the full eigenvalues and singular values: `k*` agrees, `sigma_2(P^(k*))` agrees within `5e-7` and `abs lambda_2` within `7e-4`, the gap being the power iteration's. The power iteration needs its 400 burn-in steps: 40 steps leave errors up to `2.6e-2`.

The spectral bound is loose by a factor near two. For the median cell `0x004a` with `net`, `primes`, `shift`, `sigma_2(P^3) = 0.632067` gives `d(1 + 3 s) <= 11.3 x 0.632067^s`, below `0.01` from `s = 16`, round 49; the curve is there at round 22. The rate is what the bound misses: `0.632067^(1/3) = 0.858` per round against `abs lambda_2 = 0.752`.

## The alternating cells at 4 x 4

At `4 x 4` the kernel has 65535 rows and from 0.56 to 4.0 million nonzero entries, median 2.0 million, against `4.3 x 10^9` for a dense table.

**The irreducible nonlinear cells are the 198 alternating cells. Verified.** Of the 2142 nonlinear cells, 1802 fall to the spine's law I (a read-closed set, which makes the chain reducible without building it), 142 more are reducible by search and 198 are irreducible: the spine's 113 module, 8 complement, 18 parity and 3 product cells against its 198 alternating cells, count for count.

**Sampled curves. Verified on 23 cells.** Every ninth alternating cell in catalog order and the figure cell, 23 cells, each with 32 seeded kernel rows checked against the direct table. The distance is the worst of 18 starts, the 16 single cells, the whole fill and the whole void, so it is a lower bound for `d(r)`; the single-round bound with `sigma_2(P^(k*))` is the upper side. `abs lambda_2` is the largest modulus among the Ritz values of restarted Arnoldi on the sum-zero space: 80 steps with full reorthogonalisation, restarted from the Ritz vector until the unit Ritz vector `x` has `||P x - lambda x|| <= 1e-10`; the largest residual is `5.0e-11` and a second seed agrees within `9e-13` on every cell.

| quantity | min | median | max |
|---|---|---|---|
| rounds to `1/4`, 18 starts | 8 | 11 | 61 |
| rounds to `0.01`, 18 starts | 11 | 19 | 106, one cell past 120 |
| `k*` | 3 | 4 | 4 |
| `sigma_2(P^(k*))` | 0.237400 | 0.471190 | 0.820184 |
| `abs lambda_2(P)` | 0.473767 | 0.725209 | 0.951110 |

The figure cell of the spine, the level-2 carpet `0x135f` with the carpet neighbourhood, `primes` and `slide`, has `k* = 4`, `sigma_2(P^4) = 0.369446`, `abs lambda_2 = 0.651795`, and needs 11 rounds to `1/4` and 19 to `0.01`. The slow end is the same story as at `3 x 3`: `0x1212` with `net`, `fibonacci`, `shift`, `abs lambda_2 = 0.951`, 61 rounds to `1/4` and not within `0.01` after 120. The median cell needs about as many rounds at 16 bits as at 9; the light cone below says why that cannot last.

## Three points

The 3-wise chain runs on the 260610 ordered pairs of distinct nonzero differences at `3 x 3`, applied on the fly from the same derivative table; its first-coordinate marginal reproduces the 2-wise curve within `4e-16` at every round, the lumping checked. Every thirtieth alternating cell, five cells, is run from the eight starts `(e_p, e_q)` with `p` the first fill cell.

| cell | `k*` | `sigma_2(P_3^(k*))` | `sigma_2(P^(k*))` | `abs lambda_2(P_3)` | `abs lambda_2(P)` | 3-wise rounds to `1/4`, `0.01` | 2-wise |
|---|---|---|---|---|---|---|---|
| `0x004a` net primes shift | 3 | 0.632067 | 0.632067 | 0.7523 | 0.7524 | 11, 22 | 10, 22 |
| `0x006d` void primes shift | 3 | 0.626550 | 0.626550 | 0.7105 | 0.7105 | 10, 19 | 10, 19 |
| `0x000e` net binary drop | 3 | 0.651366 | 0.651366 | 0.7456 | 0.7482 | 11, 23 | 10, 21 |
| `0x0055` void binary drop | 6 | 0.614733 | 0.614733 | 0.8906 | 0.8907 | 24, 53 | 21, 50 |
| `0x005e` net binary slide | 3 | 0.626550 | 0.626550 | 0.7096 | 0.7096 | 10, 19 | 10, 20 |

**Three inputs mix at the rate of two on these cells. Conjecture.** The lumping proves `sigma_2(P_3^j) >= sigma_2(P^j)` and that every eigenvalue of `P` is one of `P_3`; on all five cells Lanczos finds no larger singular value, to six digits, and power iteration no larger rate, within its `3e-3`. The 3-wise distance from eight starts trails the 2-wise one from all starts by at most three rounds. A Lanczos value is a lower estimate, so equality is measured, not proved.

## Design stamps

The chain needs the stamp to be uniform on all of `V`. The spine's own stamp is the split's design, one fixed grid, and a generator that draws its stamps from a few designs has a key space of `k^r` after `r` rounds; the output pair then takes at most `k^r` values, which the uniform law on the `2^N (2^N - 1)` pairs of distinct grids weighs at most `k^r / (2^N (2^N - 1))`, so the worst pair's distance from it is `d(r) >= 1 - k^r / (2^N (2^N - 1))` whatever the rounds.

**Design-shaped stamps break the Markov property. Verified.** The witness is the first alternating cell, `0x004a` with `net`, `primes`, `shift`, and stamps uniform on the three neighbourhood designs `{495, 186, 341}`, every path enumerated exactly. Under uniform stamps the law of `d_2` given `d_1` is the kernel row, distance 0. Under the three designs, from the void start `e_0` and a uniform input, the law of `d_2` given `d_1` is at distance `0.1667` from the kernel row. From a uniform input pair, where the first step is stationary and agrees with the kernel, the law of `d_3` given `(d_1, d_2) = (0x084, 0x101)` is at distance `0.5` from the law of `d_3` given `d_2 = 0x101` alone: the history matters, so no kernel of any kind describes the difference. The reason is the proof above read backwards: a design stamp leaves `x_t` correlated with the differences that produced it.

## The light cone

A difference spreads only where the mess reads it. The mess of a fill cell reads the cells within `l_infinity` distance 1 of it, for all three neighbourhoods, and every spread of the catalog is an isometry of the torus in that metric: the four translations are, and so are `turn` and `flip`, which are linear with an orthogonal integer matrix.

**The light cone. Proved.** Suppose every round of a schedule reads within distance `m` and its spread maps every ball of radius `t` into a ball of radius `t + s`, and let `c = m + s`. A start `e_p` reaches after `k` Markov steps only differences supported in the ball of radius `c k` around some cell, so

```
d(r) >= 1 - (2^(w^2) - 1) / (2^(n^2) - 1),     w = min(n, 2 c (r - 1) + 1).
```

The proof is the support: `D_a mess(y)` is zero outside the fill cells that read `supp a`, all within `m` of it, and the spread adds `s`. A law supported on `S` is at distance at least `1 - U(S)` from uniform. While `2c(r - 1) + 1 <= n - 1` the bound exceeds `1 - 2^(1 - 2n)`, so any local schedule needs at least `1 + (n - 1)/(2c)` rounds to bring `d` to `1/2`, which at side `n = 2^L` is `Omega(2^L)`. Every catalog spread has `s = 0` and `c = 1`.

**A translation at any scale cannot break the cone. Proved.** A translation by any vector is an isometry of the torus in `l_infinity`, so a schedule that slides by `2^j` in round `j` has `s = 0` and `c = 1` like the plain slide: moving the whole cone does not widen it. A spread breaks the cone only if it stretches distances.

**A sum of shifts is not an isometry. A door.** The layer `1 + x^(2^j) + y^(2^j)` adds to a cell function its shifts by `2^j` along each axis. It is not a permutation of the cells: it sends a cell to three cells `2^j` apart, so `s = 2^j`, the radius bound `c k` is lost, and the support triples per layer. Two facts hold by hand on the torus of side `n = 2^L`. It is invertible: in characteristic 2 put `x = 1 + X` and `y = 1 + Y`, so `X^n = Y^n = 0` and `x^(2^j) = 1 + X^(2^j)`, and the layer is `1 + X^(2^j) + Y^(2^j)`, one plus a nilpotent, a unit. It cannot cover by itself, in `L` rounds or in any number: the product of the layers `j = 0, ..., L - 1` has coefficient 1 at `x^p y^q` exactly when `p AND q = 0` (pick one of the three terms for each `j`), which is the level-`L` carpet, `3^L` of the `4^L` cells; and every product of layers is a unit, so its support has odd size (the coefficient sum is the value at `x = y = 1`, which is 1) and is never the even-sized torus. Whether the layers together with the mess, which adds the readers between layers, cover the torus in `O(L)` rounds and then mix is open: no lab verb builds a spread that is not a permutation of the cells, so neither fact is run and neither the cover time nor the distance is checked.

**The exact cone.** The bound is far from tight because the Feistel step only spreads from void cells into the fill cells that read them. The exact reachable support `R -> sigma(R + {a in A : a reads R})` gives `d(r) > 1/2` for every `r <= T`, `T` the number of rounds after the first, one application of the map each, that the last start needs to cover the torus, since a support of `N - 1` cells or fewer carries less than half of `U`; `d` first falls below `1/2` at round `1 + T` at the earliest. On the level-`L` carpet `{i AND j = 0}`, the spine's design, `slide` covers in `T = 8, 15, 30, 61, 123` rounds at `L = 2..6`, about `2n`, four times the radius bound. At `n = 4` the cone is tight: the exact chain of the spine's figure cell first falls below `1/2` at round 9 = `1 + T`, and the board's `shift` and `shear` below at rounds 6 and 7, again `1 + T`.

**Every schedule needs depth `Omega(L)`. Proved.** A void cell is read by at most 8 fill cells, so the reachable support grows at most ninefold per round and covering `n^2` cells takes `T >= log_9(n^2) = 0.63 L` rounds, whatever the spread, as long as it permutes the cells; a layer that sends a cell to three would triple the growth and leave `log_27(n^2) = 0.42 L`, still linear in `L`.

**A cone that reaches a fixed set short of the torus is law I. Proved.** For a single spread the size of the reachable support grows or stays; if it reaches a fixed set `R* != V` then `sigma(R* + W) = R*` with `W` the fill cells reading `R*`, and since `sigma` is a bijection `W` lies in `R*` and `sigma(R*) = R*`: a read-closed spread-invariant cell set, the spine's law I, and the group is imprimitive. A spread orbit that lies inside the fill is the smallest case. A cone can also stall and move, as in the product cell `0x0054` with `void`, `primes` and `shift`, which never covers and has no law-I set.

## The multiscale spread

The first stretching spreads tried are the natural design ones on the carpet, and they all fall to the law above. The digit rotation `(rot i + 1, rot j + 1)`, `rot` the cyclic shift of the `L` binary digits, moves neighbours at scale 1 to scale `2^(L-1)`, and every start from `L = 3` on lies in a proper closed set. The cat map `(2i + j + 1, i + j + 1)` fixes the cell `(n - 1, 0)`, which lies in the fill: an orbit inside the fill at every `L`.

| spread | carpet `L = 2..6` | board `L = 2..6` |
|---|---|---|
| `slide` `(i + 1, j + 1)` | 8, 15, 30, 61, 123 | never |
| `shift` `(i, j + 1)` | | 5, 9, 17, 33, 65 |
| `scales`, slide by `2^(t mod L)` | 9, 22, over 64, over 128, over 256 | |
| `row` `(rot i + 1, j + 1)` | 7, 14, 27, 53, 103 | |
| `shear` `(i + j + 1, j + 1)` | 7, 14, 28, 54, 105 | 6, 9, 13, 21, 37 |
| `rotate` `(rot i + 1, rot j + 1)` | never | never |
| `cat` `(2i + j + 1, i + j + 1)` | an orbit in the fill | never |
| `hyper` `(2i + j + 1, 3i + 2j)` | | 4, 6, 7, 8, 10 |

Rounds after the first for the last of every start to cover the torus. "Never" is certified: some start lies in a proper set closed under the spread and under reading, law I; "over" is the cap of `4n` rounds with no closed set found.

**The board and the hyperbolic spread.** The board `{i + j even}` is the design whose word is the full letter at every level but the last and the diagonal letter 9 at the last; it is the split of the two-sided Feistel network, every fill cell reading its four orthogonal void neighbours. `shift` maps the fill onto the void and back, the classic alternating Feistel, and its cone is linear, `T = n + 1`. The spread `hyper(i, j) = (2i + j + 1, 3i + 2j)` keeps the alternation and stretches. **It is a bijection that swaps fill and void at every `L`. Proved.** Its matrix `[[2, 1], [3, 2]]` has determinant 1, so it is invertible modulo every `n`, and `(2i + j + 1) + (3i + 2j) = 5i + 3j + 1` has the parity of `i + j + 1`. Its eigenvalues `2 +- sqrt 3` stretch an unstable direction by `3.73` per round. Its cone covers the torus in 4, 6, 7, 8, 10 rounds from `L = 2` to `6`, while the board's `shift` needs 5, 9, 17, 33, 65.

**The exact chains at `n = 4`. Verified.** The cone only bounds the support; the chain decides. The board with `primes` is reducible under `shift`, `shear` and `hyper` alike, its distance stuck at `0.9844` or `0.9375`; on four read cells `primes` is the quadratic pattern `[w mod 4 in {2, 3}]` of [the mess](mess.md). With `binary` all three are irreducible:

| schedule at `n = 4` | `abs lambda_2` per round | rounds to `1/2` | `1/4` | `0.05` | `0.01` |
|---|---|---|---|---|---|
| carpet primes slide, the figure cell | 0.6518 | 9 | 11 | 15 | 19 |
| carpet primes scales, period 2 | 0.6762 | 9 | 11 | 16 | 20 |
| board binary shift | 0.6921 | 6 | 7 | 10 | 14 |
| board binary shear | 0.6296 | 7 | 8 | 10 | 14 |
| board binary hyper | 0.6921 | 6 | 7 | 10 | 14 |
| board fibonacci hyper | 0.8884 | 19 | 28 | 44 | 59 |

The distance is the worst of the 18 starts; `abs lambda_2` is the same Arnoldi, every residual at most `1e-10`. At `n = 4` the board's cone is already full after 4 to 6 rounds and the three board spreads tie; the difference shows only at larger sides.

**Sampled runs at `n = 4, 8, 16`. Verified.** Past `n = 4` the chain has `2^64` and `2^256` states; collision counting needs about `sqrt(2^N)` runs and is blind there, the collision count of 16384 runs falling to the uniform level by round 8 on the board and 20 on the carpet, long before the difference is random. The meter is the largest bit bias, `max over cells p of abs(Pr(d_r(p) = 1) - 2^(N-1)/(2^N - 1))`, a lower bound on `d(r)` since each bit is a test, from 16384 runs started at the cell the cone reaches last, at one fixed seed per side (`0x6d756c746900 + L`, `L = 2, 3, 4`, the same for the three schedules of a side); its sampling floor is about `0.02`, the printed tails running `0.007` to `0.018`, so the rounds read off are one sample per side, estimates of a lower bound, not ranges over seeds.

| schedule | bias below `0.2`, `n = 4, 8, 16` | bias below `0.05` |
|---|---|---|
| carpet primes slide | 10, 18, 45 | 14, 25, 55 |
| board binary shift | 6, 11, 20 | 8, 13, 23 |
| board binary hyper | 6, 8, 9 | 8, 9, 11 |

The local schedules double their depth with the side; the hyperbolic one adds one or two rounds per level.

The figure draws the cone of one difference on the `32 x 32` board, started at cell `(0, 0)`, 1 to 8 rounds after the first: the top two rows under `shift`, a patch that grows by at most one cell a side per round, the bottom two under `hyper`, which covers the torus 8 rounds after the first. Blue is a fill cell the difference can reach, orange a void cell, the rest is faint. The figure recomputes both cones and checks the cover and the light-cone bound of `(2k + 1)^2` cells `k` rounds after the first on the shift.

**The depth of the hyperbolic board is `O(L + log 1/eps)`. Conjecture.** The lower side is proved, `Omega(L)` for every schedule whose spreads permute the cells. The upper side is measured on three sizes for a lower-bound meter, and the falsification test, a growth faster than linear in `L` at `n = 4, 8, 16`, did not trigger: the increments are 1 and 2 rounds, and the exact cone keeps adding one or two rounds per level up to `n = 64`.

## The obstruction

A uniform spectral gap in `n` is not the right target: with `sigma_2` bounded away from 1 the bound above still costs `log(2^(N/2)) / log(1/sigma_2)` rounds, linear in `N = 4^L`. What a depth `O(L)` needs is that the difference decorrelates at every cell exponentially fast once the cone has passed, uniformly in `n`, and then a union over the `n^2` cells.

For a local schedule the chain on the plane is a random local map, a probabilistic cellular automaton driven by the uniform voids, and that uniform decorrelation is its exponential ergodicity. The contraction criteria that descend from Dobrushin 1971 give it when the total influence of the neighbours on a cell is below 1; here every cell's own old difference passes to its image with influence 1, so no such criterion holds at the first step. And ergodicity can fail outright: the abstract of Gacs 2001 recalls a two-dimensional automaton with every local transition of positive probability that keeps one bit forever, and builds a one-dimensional one. No general theorem applies; each design would be proved alone. The light cone is the same obstruction seen from outside: in infinite volume a finite difference never fills the plane, and on the torus the filling takes a time proportional to `n`.

The hyperbolic spread leaves the class: `hyper` is not a cellular automaton, it moves neighbours far apart, and the chain becomes a random perturbation of a toral automorphism. Its cone stays within a factor two of the counting bound `log_5(n^2)` of the board, whose void cells are each read by four fill cells, but the decorrelation behind the cone has no theorem to borrow.

## Open

- **A proof of the conjecture.** A depth `O(L + log 1/eps)` for the hyperbolic board, or for any design schedule, with a constant read off the round.
- **A design-native stretching spread.** `hyper` is a matrix, not a design; whether a spread built from the tree's own letters, a level permutation composed with a design-dependent shift, can stretch without an orbit inside its fill is open; every one tried here falls to law I.
- **A sum-of-shifts layer.** `1 + x^(2^j) + y^(2^j)` is invertible but not an isometry, and by itself never covers the torus; whether it, with the mess between layers, covers in `O(L)` rounds and mixes needs a verb that builds a spread that is not a permutation of the cells.
- **Three points.** Whether `sigma_2(P_3^j) = sigma_2(P^j)` for every alternating cell, so that 3-wise independence costs nothing beyond 2-wise.
- **The worst start.** At `3 x 3` the worst start is a fill or one-cell difference in 98 percent of the round maxima; a proof that the light differences are the worst would turn the 18-start curves at `4 x 4` into exact ones.
- **Rounds to random.** The spine asks when a composite's cycle type looks random; pairwise mixing is a different clock, and how the two compare on one schedule is open.
