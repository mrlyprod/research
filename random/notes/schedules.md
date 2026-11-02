---
title: The schedules
lead: Which schedules of design rounds generate a primitive group, for every n: a trap inside the written side that no mess can break, a theorem that turns two local checks, a connected second-difference graph and isolated read cells, with one live read, into primitivity on any torus, two cells showing that the obvious three checks are not enough, and the first 5 x 5 group the tree has read, A_(2^25).
figure: research-random-schedules
---

The [rounds](random.md) census reads the group of one design round and every stamp on the `3 x 3` and `4 x 4` tori by brute force: for every nonzero vector it computes the smallest block through it. That route ends at `4 x 4`, since a `5 x 5` grid has `2^25` states and as many vectors. This page asks for primitivity from the shape of a schedule instead, for every `n`. It finds one source of blocks that every mess obeys, a trap inside the written side, and counts it in closed form for the six spreads at every `n`. It proves that a schedule whose second differences connect the torus and whose read cells are isolated and not all dead generates a primitive group, on any torus and for any number of rounds, and it shows by two `4 x 4` cells that the three checks one would try first (connected reads, spanning second differences, no affine direction) are not enough. The theorem then certifies the checkerboard pair at every even `n >= 6` and seven of the 81 line pairs at `5 x 5`, one of which generates the alternating group on `2^25` states by a single prime cycle.

One generator, `lab/rs/schedules`, prints every number below; it rebuilds the round of the rounds census and does not share its code.

## The schedule

The grid is the torus `Z_n x Z_n`, a state `x` in `V = F_2^N`, `N = n^2`, `M = 2^N`. A **round** writes a split `A` from its complement `B`: a written cell `a` reads its neighbourhood cells in `B`, the set `N(a)` of size `m_a`, and the mess `m(x)` sets bit `a` to `f_a(x) = g(|x ∩ N(a)|)`, with `g` the indicator of a Source. The round without its stamp is `h = L F`, `F(x) = x + m(x)`, `L` the spread, a permutation `sigma` of the cells. A **schedule** is a list of rounds `h_1, ..., h_k`, each with its own split and spread, and its group is `G = <T, h_1, ..., h_k>`, `T` the `M` stamps. Two schedules carry the questions of this page: a **pair**, a split and its complement written in turn with the identity spread, the two-sided Feistel network; and a **single**, one split with one spread, among them the masks whose spread carries `A` onto `B`.

`T` is transitive, so the question is primitivity, and the subspace fact of Caranti, Dalla Volta and Sala 2009 in [refs](../refs.md), as the rounds page reads it, answers it: `G` is imprimitive exactly when a proper nonzero subspace `U` has `D_u h_j(x) = h_j(x + u) + h_j(x)` in `U` for all `u` in `U`, all `x` and every `j`. Such a `U` is a **block**. Every computation below rests on one identity:

```
D_u h(x) = L (u + D_u m(x))
```

Two consequences are used throughout. If `u` is in a block `U` then `L(D_u m(x) + D_u m(x'))` is in `U` for all `x, x'`, so `L D_c D_u m(z)` is in `U` for all `c, z`, and so is every higher derivative of `D_u m` taken along at least one more direction. And if `u` touches no cell that the mess reads, `D_u m = 0` and `L u` is in `U`.

A written cell is **live** when `g` is not constant on `0..m_a`; a dead cell writes a constant and never matters to a derivative. A cell is **live-read** in a round when some live written cell reads it, and **inert** in that round otherwise. Five local conditions are computed for every cell of every schedule:

- **R, connected reads.** The reads digraph on the cells, with an edge `p -> sigma_j(p)` for every round `j` and an edge `p -> sigma_j(a)` whenever a live `a` reads `p` in round `j`, is strongly connected.
- **S, spanning second differences.** For each live-read `p`, the span `S_j(p)` of the second differences `D_c D_(e_p) m_j(z)`, which lies in the coordinates of the live readers of `p`, is all of them.
- **W, no affine direction.** Every live window rule, `g` read on `m_a` variables, is neither affine nor self-complementary.
- **I, isolation.** Every live-read cell `q` of every round has a live reader `a` and a nonempty set `Q` of cells of `N(a)` other than `q` such that `a` is the only live written cell of that round reading all of `Q`, and `Delta^|Q| g`, read on the remaining `m_a - |Q|` variables, is neither affine nor self-complementary. Here `Delta g(t) = g(t) + g(t + 1)`.
- **G, the second-difference graph.** The digraph `Gamma` with an edge `p -> sigma_j(a)` whenever `e_a` lies in `S_j(p)`, and an edge `p -> sigma_j(p)` whenever `p` is inert in round `j` or `D_(e_p) m_j(0)` lies in `S_j(p)`, is strongly connected. Under S every reads edge is a `Gamma` edge, so R and S give G.

A symmetric function of `m` variables is affine when `g(t) = g(t + 2)` for all `t <= m - 2`, the constants, the parity and its complement, and self-complementary when `g(t) = g(m - t)` for all `t`.

## The trap inside the written side

**A subspace of the written side is a block of a round exactly when the spread maps it onto itself, whatever the mess. Proved.** If `W` lies in `V_A`, the mess reads only `B`, so `D_w m = 0` and `D_w h(x) = L w` for every `w` in `W`: the block condition is `L W = W`. The invariant subspaces of a coordinate permutation inside `V_A` all lie in `V_O`, `O` the union of the `sigma`-orbits that `A` contains, and `V_O` is one; so `V_O` is the largest such block, and a trap exists exactly when `A` contains a whole orbit. For a schedule the same proof gives a block in every subspace that lies in every written side and is kept by every spread.

The trap is a trapdoor in the plain sense: the rounds keep the cosets of `V_O`, and the bits outside `O` after any word of stamped rounds are a function of the bits outside `O` before it, because the mess reads only `B`, which `O` avoids, and `sigma` keeps the complement of `O`. The quotient is the same machine on the cells outside `O`, with `A` minus `O` written. This is law I of the rounds page in the case where the closed set lies inside `A`, and its proof never looks at the mess, so no Source and no neighbourhood can repair it.

**The number of proper splits of the `n x n` torus with a trap is `2^N - 1 - prod over orbits o of (2^|o| - 1)`, and the trapped fraction of all masks is `1 - prod (1 - 2^-|o|)`. Proved; Verified.** A mask avoids a trap when it misses a cell of every orbit, `2^|o| - 1` choices per orbit, and the full mask is the one trapped mask that is not proper. The six spreads have these orbits on `Z_n x Z_n`:

| spread | orbits | trapped fraction of all masks |
|---|---|---|
| `id` | `N` fixed cells | `1 - 2^-N` |
| `shift`, `drop`, `slide` | `n` orbits of `n` cells: rows, columns, diagonals | `1 - (1 - 2^-n)^n`, about `n 2^-n` |
| `turn` | `n^2/4` orbits of 4 for even `n`; the centre and `(n^2 - 1)/4` of 4 for odd `n` | `1 - (15/16)^(n^2/4)`; `1 - (15/16)^((n^2-1)/4) / 2` |
| `flip` | `n` diagonal cells and `n(n - 1)/2` pairs | `1 - 2^-n (3/4)^(n(n-1)/2)` |

`trap` checks the count by brute force over every mask at `n = 3, 4, 5`: shift, drop and slide trap 168 of the 510 proper splits at `3 x 3`, 14910 of 65534 at `4 x 4` and 4925280 of 33554430 at `5 x 5`; turn traps 286, 14910 and 22163806; flip 484, 64806 and 33495382; the identity traps all of them. It also computes the closure of a unit vector in every orbit of every `3 x 3` split with every spread, neighbourhood and nonlinear Source, 123930 closures, and finds each orbit inside `A` to be exactly its own block and no other unit closure to stay inside `A`, 0 faults.

So the spreads split in two as `n` grows. The translations trap a vanishing fraction, `0.000244` of all masks at `n = 16`, because an orbit is a whole row or diagonal. The turn and the flip trap almost every split, `0.983925` and `1 - 2^(-16) (3/4)^120`, about `1 - 1.6e-20` (`trap` prints `1.000000`), at `n = 16`, because their orbits have four cells or fewer and a random mask contains one with probability tending to 1. A design round with the turn or the flip spread is imprimitive for almost every split however strong its mess, and at `n = 32` the fraction left is below `10^-6` for both. A pair never has a trap, since `V_A` and `V_B` meet only in 0; this is the first reason to write both sides.

## The linear structures of a symmetric rule

**A symmetric Boolean function of `m >= 2` variables that is not affine has no nonzero linear structure except possibly the all-ones vector, and the all-ones vector is one exactly when `g(t) + g(m - t)` does not depend on `t`. Proved; Verified.** A linear structure is a `v != 0` with `D_v f` constant. Suppose `v` is one, with constant `c`, and `v` is neither 0 nor all ones; take `i` with `v_i = 1` and `j` with `v_j = 0`, and `tau` the swap of coordinates `i` and `j`. Since `f` is symmetric, `D_(tau v) f(x) = D_v f(tau x) = c`, so `tau v` is a linear structure with the same constant, and `D_(v + tau v) f = c + c = 0`. Now `v + tau v = e_i + e_j`, and by symmetry every `e_k + e_l` has a zero derivative. These span the even-weight vectors, so `f` depends only on the parity of `|x|`: it is affine. For the all-ones vector, `D_1 f(x) = g(m - |x|) + g(|x|)`. `lemma` lists the linear structures of every symmetric function of `m = 2..10` variables by brute force and finds no other; at `m = 10` there are 2044 non-affine functions, 60 with the all-ones structure.

So for a window rule that is neither affine nor self-complementary, `D_v f` is not identically zero for any `v` that meets the window. That is condition W, and the same lemma applied to `Delta^|Q| g` is the second half of condition I.

## Unit propagation

**If the second-difference graph is strongly connected, the only block containing a unit vector is the whole space. Proved.** Let `e_p` lie in a block `U`. If `p` is inert in round `j`, `D_(e_p) h_j = L_j e_p`, so `e_(sigma_j(p))` is in `U`. Otherwise `L_j D_c D_(e_p) m_j(z)` is in `U` for all `c, z`, so `L_j S_j(p)` lies in `U` and `e_(sigma_j(a))` is in `U` whenever `e_a` lies in `S_j(p)`; and if `D_(e_p) m_j(0)` lies in `S_j(p)`, then `L_j e_p = L_j(e_p + D_(e_p) m_j(0)) + L_j D_(e_p) m_j(0)` is in `U` as well. These are the edges of `Gamma`, so `U` contains `e_r` for every `r` that `p` reaches, which is every cell.

**Condition R is necessary. Proved.** If the reads digraph is not strongly connected, some nonempty proper set `O` of cells has no edge leaving it. For `u` in `V_O`, `D_u m_j(x)` lives on the live readers of `O`, so `L_j(u + D_u m_j(x))` lives on `sigma_j(O)` and the `sigma_j`-images of those readers, all in `O`: `V_O` is a block. This is law I of the rounds page for any schedule.

Over the 420 pair cells of both catalogs and the 2856 `4 x 4` singles below, no cell has a strongly connected `Gamma` and a proper unit closure, as the lemma requires. Condition S is not necessary: the 12 alternating `3 x 3` pairs all fail it.

## The isolation theorem

**For any schedule of design rounds on the `n x n` torus, any `n`, any spreads and any Sources, if the second-difference graph is strongly connected, every live-read cell is isolated and some cell is live-read, then the group of the rounds and every stamp is primitive. Proved.** Let `U` be a nonzero block.

1. **A read vector.** If some `u` in `U` lives on cells inert in every round, then `D_u h_j = L_j u` and `L_j u` is in `U` for every `j`. Repeat on the images. Either some vector of `U` has a live-read cell of some round in its support, or the union `O` of the supports of all these images is a nonempty set kept by every `sigma_j` and inert in every round. Such an `O` has no `Gamma` edge leaving it, so `O` is every cell, and then nothing is live-read anywhere, which is excluded. So `U` contains a `w` with `w_q = 1` for a cell `q` live-read in some round `j`.
2. **A unit vector.** Take `a` and `Q` from condition I for `q`, `k = |Q| >= 1`. The vector `D_Q D_w m_j(z)`, the derivative of `D_w m_j` along every `e_r` with `r` in `Q`, has `L_j D_Q D_w m_j(z)` in `U` for every `z`. Its coordinate at a live `a' != a` vanishes, since `a'` misses a cell of `Q` and `f_(a')` does not depend on it; a dead cell contributes nothing. At `a`, `D_Q f_a(z) = Delta^k g(|z ∩ (N(a) - Q)|)`, a symmetric function `f'` of the `m_a - k` cells of `N(a)` outside `Q`, so `D_Q D_w f_a = D_(w') f'` with `w'` the part of `w` on those cells, nonzero because it contains `q`. By the lemma `D_(w') f'` is not identically zero, since `f'` is neither affine nor self-complementary. Choose `z` with value 1: then `D_Q D_w m_j(z) = e_a` and `e_(sigma_j(a))` is in `U`.
3. **Everything.** By unit propagation, `U = V`.

The two sources of blocks the target named are the two steps. A block with no read vector is the trap of the last section, generalised to cells inert in every round, and step 1 kills it with `Gamma`. A block with a read vector is killed in step 2, but only through isolation: a read vector of large weight gives second differences on several written cells at once, and without a set `Q` that singles one of them out, nothing forces a unit vector into `U`. The next section shows that this is a real gap and not a weakness of the proof.

**For a pair whose mess derivatives span each written side, every block lies in a split block `X + Y`, `X` in `V_A`, `Y` in `V_B`, with the derivatives of `m_A` along `Y` in `X` and those of `m_B` along `X` in `Y`. Proved.** For `u` in a block `U`, `D_u F_A(x) = u + D_(u_B) m_A(x)` and `u` are in `U`, so `D_(u_B) m_A(x)` is in `U` and in `V_A`. Hence `X + Y`, with `X` and `Y` the projections of `U`, is a block containing `U`; if it were all of `V`, `U` would contain every derivative of `m_A` and of `m_B`, which span `V_A` and `V_B`, and `U = V`. Read dually, a proper block of a pair is a set of parities of the written side that the mess maps to functions of a set of parities of the read side, and back.

## Every even n

**The checkerboard pair, `A` the cells with `i + j` even, with the `carpet` or `net` neighbourhood and the `binary` or `fibonacci` Source, generates a primitive group with every stamp on the `n x n` torus for every even `n >= 6`. Proved; Verified.** Both neighbourhoods read the four orthogonal neighbours, which are the other colour, so `m = 4` everywhere, every cell is live, and the reads digraph is the plus graph between the colours, connected for every `n`. For isolation take `Q = {a + e_1, a - e_1}`: a cell reading both is orthogonally adjacent to both, which for `n >= 6` is `a` alone; `q = a + e_2` is outside `Q`, and `Delta^2 g(t) = g(t) + g(t + 2)` on two variables is non-affine and not self-complementary exactly when `g(0) != g(4)`, true for `binary` (0 and 1) and `fibonacci` (1 and 0). For spanning, the cell `p + 2 e_1` is read by the reader `p + e_1` of `p` and by no other reader of `p`, so `D_(e_(p+2e_1)) D_(e_p) m(z)` is `Delta^2 g(t) e_(p+e_1)`, `t` the live count among the two other reads of `p + e_1`; `Delta^2 g` on `0..2` is `1, 1, 0` for `binary` and `0, 0, 1` for `fibonacci`, not identically zero, so `e_(p+e_1)` lies in `S(p)`, and by symmetry so does each of the four readers. Every check compares two cells at most 4 apart in each coordinate, a cell of `Q` against a neighbour of another, and such cells coincide on the torus only when `n <= 4`; so the conditions are the same at every even `n >= 6`, and `family` computes all five at `n = 6, 8, 10` and finds them met.

At `n = 4` the picture differs: `a + 2 e_1 = a - 2 e_1`, the pair `Q` has two common readers, and the census finds the `4 x 4` checkerboard alternating for `binary` and `fibonacci` and imprimitive for `primes`, which keeps the parity of the whole torus. `primes` fails isolation at every `n`: `Delta^2 g` is constant on two variables and no other `Q` isolates a cell. The figure draws the second-difference graph of the checkerboard `carpet binary` pair at `n = 10`, recomputed from the Source, as the depth at which a unit vector reaches each cell: `1, 4, 8, 12, 16, 18, 16, 12, 8, 4, 1` cells at depths 0 to 10, the torus distance.

By the same count the stripes pair, `A` the even columns, with `carpet` and `primes` or `binary`, meets every condition at every even `n >= 6`, spanning as computed by `family` at `n = 6, 8, 10`: a written cell reads six cells, the diagonal pair `{a + (-1, -1), a + (1, 1)}` isolates it, and `Delta^2 g` on four variables is `1, 1, 1, 0, 0` for `primes` and `1, 1, 0, 0, 1` for `binary`, neither affine nor self-complementary. For `fibonacci` it is `0, 0, 1, 0, 0`, self-complementary, and isolation fails.

## The long prime

**A transitive group of degree `M` with an element that has a cycle of prime length `p`, `M/2 < p <= M - 3`, is primitive and contains the alternating group. Proved.** Every other cycle is shorter than `M/2 < p`, so no other length is divisible by `p`, and the power of the element by the product of the other lengths is a `p`-cycle `c`. If the group kept blocks of size `s`, `2 <= s <= M/2`, then `c`, of prime order `p`, would either move the `M/s <= M/2 < p` blocks in `p`-cycles, impossible, or fix every block, and then its one orbit of `p` points would sit in a single block of size `s < p`. So the group is primitive, and Jordan's theorem in [refs](../refs.md), with `p <= M - 3`, gives the alternating or the symmetric group; since every round that writes two cells is even, it is the alternating group.

The certificate needs no block search, so it reproves the alternating verdicts of the rounds census by a route that shares nothing with it. Words of one to four stamped passes through the schedule, at most 512 per cell, find such a cycle in every one of the 12 alternating `3 x 3` pairs (primes from 257 to 499 of 512 states), the 10 alternating `4 x 4` pairs, and 195 of the 201 primitive `4 x 4` singles; by the rounds census the other 6 are its 3 primitive products, which can never have one, and 3 alternating cells it certifies by a shorter lone prime.

## The census

The pair census rebuilds the rounds census's 12 `3 x 3` splits up to symmetry and complement and its 23 `4 x 4` words, with three neighbourhoods and four Sources; the 108 and 207 nonlinear cells give `A_512` in 12 and `A_65536` in 10, as there. The conditions read:

| `3 x 3` pairs | cells | | `4 x 4` pairs | cells |
|---|---|---|---|---|
| not R, imprimitive | 62 | | not R, imprimitive | 179 |
| R, not S, imprimitive | 34 | | R, not S, imprimitive | 18 |
| R, not S, alternating | 12 | | R, not S, alternating | 1 |
| R and S | 0 | | R and S, alternating | 9 |

**The 96 `3 x 3` pair cells that miss `A_512` fail connected reads or spanning second differences. Verified.** 62 have a read-closed set. The other 34 are connected but do not span: 32 have a proper mess module, and 2, `0x004e` with `net` and `0x001b` with `void`, both `binary`, keep a one-dimensional block spanned by a vector that is neither a unit nor all ones, a linear structure of every window it meets; both fail W. Of the 197 `4 x 4` pairs that miss `A_65536`, 179 are read-closed and 18 connected without spanning: 13 by the mess module, 3 by the parity of the whole torus, and 2, `0x134c` and `0x125a` with `carpet` and `primes`, by a single kept parity, of the two diagonal `2 x 2` tiles `0xcc33` and of the checkerboard `0xa5a5`, a design parity that is not the parity of a partition. No pair in either catalog meets the theorem's hypotheses: the 26 and 45 isolated pair cells are all read-closed, and no connected pair cell is isolated. For pairs the theorem's reach starts at `5 x 5`.

**The single census and the theorem agree on every `4 x 4` cell. Verified.** Over the 2856 single cells of the rounds census, 2142 nonlinear, the laws recomputed here split the imprimitive cells exactly as there, read-closed 1802, mess 113, complement 8, parity 18, and leave 201 primitive. The theorem's hypotheses hold in 59 of them, all primitive and 58 certified alternating by a long prime. The closure of the 16 unit vectors and the all-ones vector alone decides primitivity in every nonlinear `4 x 4` cell, single or pair: 238 nonlinear singles have every unit closure full, and each imprimitive one among them has a proper all-ones closure. At `3 x 3` that shortcut fails in exactly the two pairs above. By spread, the cells with every unit closure full are 0 of 225 for `id`, 0 of 540 for `flip`, 53 of 351 for `shift` and for `drop`, 30 of 315 for `slide` and 102 of 360 for `turn`, of which 49, 49, 27 and 76 are primitive.

## Three checks are not enough

**Connected reads, spanning second differences and no affine direction do not force a primitive group. Refuted.** Two `4 x 4` singles meet all three and are imprimitive. The split `0x0055`, four cells of the first two rows, with `net`, `primes` and `slide` keeps the parities of the eight cosets of a torus subgroup of order 2, `{p, p + (0, 2)}`; the split `0x0505` with `net`, `primes` and `turn` keeps those of the four cosets of a subgroup of order 4. Each block is the set of grids even on every coset, of codimension 8 and 4, and contains no unit vector, as unit propagation requires: the obstruction is a kept parity system, exactly the case step 2 of the theorem handles by isolation, and both cells fail isolation. These are the three checks one would try first; isolation is the repair, and with it the statement is the theorem above.

The two counterexamples are singles with a spread that does not carry the split onto its complement. On pairs and on masks carried onto their complement, the three checks have no counterexample in the catalogs: the 6 `4 x 4` pairs that meet them are all alternating, and of the 63 nonlinear `4 x 4` singles whose spread carries `A` onto `B`, the 14 that meet them are all alternating, while the other 49 are read-closed (23), mess-bound (19) or parity-bound (7).

**For a pair, or a single whose spread carries the split onto its complement, connected reads, spanning second differences and no affine direction force a primitive group. Conjecture.** It holds in all 20 such cells of the catalogs, and by the split-block proposition a counterexample pair would have to keep a parity system between the two sides with no unit vector in either half.

## Five by five

The `5 x 5` torus has no design catalog in the tree, so the census takes the simplest family on it: a split made of parallel lines `{(i, j) : j - s i in S}` of the affine plane over `Z_5`, with `s = 0` (columns), `s = 1` (diagonals) or `s = 2` (knight lines) and `S` one of `{0}`, `{0, 1}`, `{0, 2}`: 9 pairs and 81 nonlinear cells. The theorem needs no block search, so it reads all 81 in a fraction of a second.

**Seven of the 81 nonlinear `5 x 5` line pairs generate a primitive group with every stamp. Proved.** They are columns `{0, 2}` with `carpet primes`, `carpet binary` and `void binary`; diagonals `{0, 2}` with `carpet binary` and `net binary`; and knight lines `{0, 1}` and `{0, 2}` with `carpet binary`. Another 8 meet the three checks without isolation; the conjecture claims them, and nothing here decides them.

**The pair writing columns 0 and 2 and then the other three, with `carpet` and `primes`, generates `A_(2^25)`. Verified.** Primitivity is the theorem. The first word tried, two stamped pair rounds, has exactly one cycle of length divisible by the prime 1229311, of that length, so a power of it is a 1229311-cycle, and with primitivity Jordan's theorem gives the alternating group; independently, the 13th word tried, also of two pair rounds, has a cycle of prime length 31896037, between `2^24` and `2^25 - 3`, which by the long-prime lemma alone gives a primitive alternating group without the theorem. Both primes are confirmed by PARI. The decomposition runs on the fly over the `2^25` states with a 4 MB visited set, about two seconds a word, 36 seconds in all. The 13 words of two or three stamped pair rounds have 14 to 42 cycles, against `H_M = 17.91` for a uniform permutation of `2^25` points. This is the first `5 x 5` group the tree has read.

## Open

- **The conjecture for pairs.** Prove that connected reads, spanning second differences and no affine direction suffice for a pair, perhaps through the split-block proposition, or find a pair that keeps a parity system; the eight `5 x 5` line cells that meet the three checks without isolation are the first place to look, by a long prime for the positive side.
- **Alternating for every n.** The theorem gives primitivity for the checkerboard at every even `n >= 6`; the alternating group needs a prime cycle at each `n`, and a uniform construction of one, or a route through O'Nan-Scott as the rounds page suggests, would close the family.
- **Isolation on small tori.** Isolation never holds for a connected `3 x 3` or `4 x 4` pair, yet 22 of them are alternating; a sharper extraction of a unit vector, from read vectors of weight two, may cover them.
- **The turn and the flip.** Their trap swallows almost every split as `n` grows; a schedule that alternates them with a translation removes every trap, and the rounds-to-random question of the spine should be asked of such schedules.
