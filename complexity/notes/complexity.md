---
title: Complexity
lead: Boolean complexity of the catalog, and the Laplacian spectra of the fractals it builds.
figure: research-complexity
slug: complexity
---

A design is a [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) rule on the corners of a cube. Read its fill vector in corner order and it is the truth table of a Boolean function - that is the theorem of [the bijection page](../../bijection/notes/bijection.md), and it means the whole apparatus of Boolean complexity applies to the catalog without translation. Sensitivity, certificate complexity, decision-tree depth: these are properties of a design, computable from its filled corners with nothing but parity-cube combinatorics.

The question this page answers is how much of that complexity the *geometry* can see. The census already carries genus, degree, popcount and the odd-side fill polynomial. Do they determine how hard the function is? At dim 3 yes, but only because the catalog is too small for the question to bite. At dim 4 no - and the failure is exhibited by two named designs that draw indistinguishable fractals and differ in six of seven complexity measures.

A second sense of complexity closes the page: the spectra of the fractals the designs build. [The graph Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix) of the Sierpinski triangle that `bang dim 2, code 7` draws is degenerate to a degree that is itself a law: a third of its spectrum sits at eigenvalue 1, proved at every level by a matching count, while the next multiplicities are fitted to level 8; the base-2 flake of `bang dim 3, code 23` carries an exact interior band gap; and every fractal tested clusters its eigenvalues where a random matrix would repel them.

The generators are `six::census` and `six::topology` in `crates/math`, with `mrly_math::spectrum` for the eigensolver and the exponent, together with `lab/py/boolean-measures`, `lab/py/base-q-anf`, `lab/py/laplacian-degeneracy`, `lab/py/tree-matching`, `lab/py/flake-band-gap` and `lab/rs/spectral-spacings`.

## Seven measures on a finite catalog

Fix `dim`. There are `2^(2^dim)` designs and, up to cube symmetry, `A000616(dim)` classes - 22 at dim 3, 402 at dim 4. Both catalogs are small enough to measure exhaustively and exactly: no sampling, no approximation.

- `s`, **sensitivity**: the largest number of single-bit flips at one input that change the value.
- `bs`, **block sensitivity**: the largest number of pairwise disjoint blocks of coordinates whose flip changes the value.
- `C`, **certificate complexity**, `max(C_0, C_1)`: the fewest bits that pin the value down.
- `dt`, **decision-tree depth**: the exact minimum depth of a deterministic query tree.
- `deg`, **real degree**: the degree of the unique multilinear real polynomial representing `f`. This is *not* the `GF(2)` algebraic degree of [the core page](../../core/notes/core.md), and the two genuinely differ - the carpet `bang dim 3, code 23` has `GF(2)` degree 2 and real degree 3.
- `dnf` and `cnf`: the minimum number of terms in a DNF, and of clauses in a CNF, by exact prime-implicant cover.

All seven are constant on a hyperoctahedral orbit, so they descend to the classes unambiguously. (Proved: cube symmetry is permutation and negation of the input variables, and every measure above is defined by a property of the truth table that survives relabelling and complementing inputs - the output is never touched, which is what would swap `dnf` and `cnf`. Verified orbit-wide over all 256 designs at dim 3, `lab/py/boolean-measures`.) Hand cases pin the code down - the dictator has every measure 1, the 3-bit AND has `s = bs = C = dt = deg = 3` but `dnf = 1`, parity has `s = bs = C = dt = deg = 3` with `GF(2)` degree 1. (Verified.)

## The 3D catalog, fully measured

Twenty-two classes, every measure exact. (Verified: `lab/py/boolean-measures` recomputes all seven measures from the truth-table definitions and diffs the result cell by cell against the catalog files `measures_d3.csv` and `measures_d4.csv` - 286 cells at dim 3 and 5226 at dim 4, zero mismatches.)

| name | genus | gf2deg | pop | `s` | `bs` | `C` | `dt` | `deg` | `dnf` | `cnf` |
|---|---|---|---|---|---|---|---|---|---|---|
| `bang dim 3, code 0` | iso | -1 | 0 | 0 | 0 | 0 | 0 | -1 | 0 | 1 |
| `bang dim 3, code 1` | iso | 3 | 1 | 3 | 3 | 3 | 3 | 3 | 1 | 3 |
| `bang dim 3, code 3` | axis | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 1 | 2 |
| `bang dim 3, code 6` | comp | 2 | 2 | 3 | 3 | 3 | 3 | 3 | 2 | 3 |
| `bang dim 3, code 7` | comp | 3 | 3 | 2 | 2 | 2 | 3 | 3 | 2 | 2 |
| `bang dim 3, code 15` | axis | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `bang dim 3, code 22` | iso | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 4 |
| `bang dim 3, code 23` | iso | 2 | 4 | 2 | 2 | 2 | 3 | 3 | 3 | 3 |
| `bang dim 3, code 24` | iso | 2 | 2 | 3 | 3 | 3 | 3 | 2 | 2 | 3 |
| `bang dim 3, code 25` | comp | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 2 | 3 |
| `bang dim 3, code 27` | comp | 2 | 4 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| `bang dim 3, code 30` | comp | 2 | 4 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| `bang dim 3, code 31` | comp | 3 | 5 | 2 | 2 | 2 | 3 | 3 | 2 | 2 |
| `bang dim 3, code 60` | comp | 1 | 4 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| `bang dim 3, code 61` | comp | 3 | 5 | 3 | 3 | 3 | 3 | 3 | 3 | 2 |
| `bang dim 3, code 63` | comp | 2 | 6 | 2 | 2 | 2 | 2 | 2 | 2 | 1 |
| `bang dim 3, code 105` | iso | 1 | 4 | 3 | 3 | 3 | 3 | 3 | 4 | 4 |
| `bang dim 3, code 107` | iso | 3 | 5 | 3 | 3 | 3 | 3 | 3 | 4 | 3 |
| `bang dim 3, code 111` | comp | 2 | 6 | 3 | 3 | 3 | 3 | 3 | 3 | 2 |
| `bang dim 3, code 126` | iso | 2 | 6 | 3 | 3 | 3 | 3 | 2 | 3 | 2 |
| `bang dim 3, code 127` | iso | 3 | 7 | 3 | 3 | 3 | 3 | 3 | 3 | 1 |
| `bang dim 3, code 255` | iso | 0 | 8 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |

The genus column reads `iso 10, axis 2, comp 10` and the `GF(2)` degree histogram reads `-1: 1, 0: 1, 1: 3, 2: 9, 3: 8`, both matching the core page exactly. (Verified by `lab/py/boolean-measures` - a crossref, recomputed rather than copied.)

One structural fact is already visible: `dt > bs` happens. Three classes - `bang dim 3, code 7`, `bang dim 3, code 31`, and the carpet `bang dim 3, code 23` - have `bs = 2` but `dt = 3`. The query tree is deeper than block sensitivity, which is the first place the query chain comes apart. (Verified, `lab/py/boolean-measures`.)

## The pin family: where geometry wins

The `dim + 1` axis designs are the ones that say "these `r` named axes must be even". Such a design *is* the AND of `r` literals, and every classical measure of an `r`-bit AND is known in closed form.

**Proved.** For the design pinning `r` of the `dim` axes, `s = bs = C = dt = deg = deg_GF2 = r`, in every dimension. The function is `x_1 AND ... AND x_r`. At the input where all `r` literals hold, flipping any one of them changes the value and flipping anything else does not, so `s = r`; querying those `r` literals decides the function, so `dt <= r`. With `s <= bs <= C <= dt` - proved below - the chain closes: `r = s <= bs <= C <= dt <= r`. The multilinear real form is the monomial `x_1 * ... * x_r`, of degree `r`, and so is the `GF(2)` form. (Verified at dim 3 and dim 4, all `r`, `lab/py/boolean-measures`.)

"Every measure equals `r`" is wrong on the last two (Refuted): the DNF size is 1 for every `r` - one term is the whole function - and only the CNF size is `r`. (Verified: `dnf = 1`, `cnf = r`, at dim 3 and dim 4, all `r`, `lab/py/boolean-measures`.) The statement covers six measures, not eight.

The isotropic designs are tame for a related reason: a level-set is symmetric under permuting coordinates, so its measures are forced by its level structure. The compounds are where the geometry runs out - and by the core page's own asymptotic count, almost every design is a compound.

## Where geometry stops

Make "determine" precise. Group the catalog by a **key**, a tuple of geometric data. The key *determines* a measure when every group carries a single value of it. A group carrying two values is a certificate that the geometry cannot see the measure.

The finest key available is `(genus, GF(2) degree, popcount, fill polynomial)`, the fill polynomial being the odd-side polynomial of [the method page](../../method/notes/method.md), read in `k` at odd side `2k - 1` - which, by that page's corollary, is equivalent to knowing how many filled corners sit at each Hamming weight. That is as much as the fractal geometry knows about a design's fill at every side and every level.

**Verified** (`lab/py/boolean-measures`). At dim 3 this key determines all seven measures. It does so vacuously: the fill polynomial alone separates 21 of the 22 classes, and the single collision (`bang dim 3, code 15` against `bang dim 3, code 27`, both with polynomial `4k^3 - 4k^2 + k`) is broken by genus and degree. There is nothing to test.

**Verified** (`lab/py/boolean-measures`). At dim 4 the same key determines *none* of the seven. Nor does any coarser key: not genus, not `GF(2)` degree, not popcount, not the fill polynomial alone. The fill polynomial sends 402 classes to 183 distinct values, 94 of them shared by two or more classes; under the full key, 92 groups hold two or more classes, 81 of those split at least one measure, and there are 279 measure-splits in total.

The smallest witness, by code:

| name | genus | gf2deg | pop | fill polynomial | `s` | `bs` | `C` | `dt` | `deg` | `dnf` | `cnf` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `bang dim 4, code 27` | comp | 3 | 4 | `4k^4 - 4k^3 + k^2` | 3 | 3 | 3 | 3 | 3 | 2 | 3 |
| `bang dim 4, code 281` | comp | 3 | 4 | `4k^4 - 4k^3 + k^2` | 4 | 4 | 4 | 4 | 3 | 3 | 5 |

**Verified** (`lab/py/boolean-measures`). Two designs with identical genus, identical `GF(2)` degree, identical popcount and identical fill polynomial, separated by six of the seven measures - only the real degree agrees. Their ANFs are `1 + x + z + xz + yw + zw + xyw + xzw` and `1 + w + z + xy + xyw + xyz + xzw + yzw`, and theirs is one of 14 groups of size two at dim 4 in which block sensitivity splits - the one with the smallest code.

Two honest qualifications. First, the fill polynomial is not a symmetry invariant: its leading coefficient is the popcount and *is* invariant, but the lower coefficients are not, because a parity flip translates the infinite tiling by one unit along each flipped axis but does not preserve its truncation to a finite side. The column belongs to the canonical representative of the class. That only strengthens the result - the key used here is finer than any genuine invariant, and it still determines nothing. Second, this is a statement about dim 4, verified exhaustively there, and not a theorem about all `dim`.

The reading: the fill polynomial is a real geometric coordinate on the catalog, and it *bounds* complexity - it fixes the leading term, hence the popcount, and at an even side the popcount fixes the fill `popcount * (side/2)^dim` and so [the fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension) `log(fill)/log(side)`. At side 3 it does not: the class of `bang dim 2, code 7` fills 8, 7, 7 and 5 cells there ([the core page](../../core/notes/core.md)). It does not *compute* complexity. Past its resolution there is no geometry-to-complexity map, and two designs are the proof.

## The catalog on the inequality web

The standard relations, and where the designs sit on them.

**Proved.** `s <= bs <= C <= dt`. Single-bit flips are disjoint blocks, which gives the first. For the second, a certificate for an input must intersect every sensitive block of that input - otherwise the block could be flipped without leaving the certified subcube, changing the value - and the blocks realising `bs` are pairwise disjoint, so the certificate has at least `bs` bits. For the third, the path a depth-`dt` tree takes on an input is itself a certificate for that input, of at most `dt` bits.

**Verified** (`lab/py/boolean-measures`). `C = bs` for every class at dim 3 and dim 4 - 424 classes, no certificate gap anywhere. A seeded spot check of 50000 uniformly random designs at dim 5 finds no exception either. This is reported as verified for dim at most 4 and open beyond: nothing here rules out a certificate gap in higher dimension, and no general theorem is claimed.

**Verified** (`lab/py/boolean-measures`). `s = bs` for every class at dim 3. At dim 4 there is exactly one exception, `bang dim 4, code 7128`, with `s = 2` and `bs = 3` - orbit size 24, ANF `x + y + xy + xz + yw + zw`. It is the family's first separation between sensitivity and block sensitivity, and it is unique. In the dim 5 spot check, 216 of the 50000 designs have `s != bs`, every one of them at `s = 3, bs = 4`; nothing wider appeared in a random sample.

**Verified** (`lab/py/boolean-measures`). `deg <= s^2` holds across both catalogs, and exactly two classes with `s >= 2` meet it with equality: `bang dim 4, code 855` and `bang dim 4, code 1911`, both with `deg = 4`, `s = 2`, both of orbit size 48.

**Proved.** `bang dim 4, code 1911` is the textbook extremal function. Its ANF is `1 + xy + zw + xyzw`, and expanding `(1 + xy)(1 + zw)` over `GF(2)` gives `1 + zw + xy + xyzw` - the same polynomial. So the design is an AND of two 2-bit NANDs, which over the reals is an AND of two ORs: the standard construction attaining `deg = s^2`, found here as a single orbit of the catalog with a geometric address.

On the literature. Huang's theorem - that a Boolean function's sensitivity and degree are polynomially related, resolving the Sensitivity Conjecture - is read for this page: [Induced subgraphs of hypercubes and a proof of the Sensitivity Conjecture](https://arxiv.org/abs/1907.00847), Hao Huang, 2019. Its Theorem 1.4 is `s(f) >= sqrt(deg(f))`, which is the `deg <= s^2` used above; the paper notes the bound is tight for AND-of-ORs, which is exactly the function found here. Its Theorem 1.5 is `bs(f) <= s(f)^4`, and it reports the Nisan-Szegedy bound as `bs(f) <= 2*deg(f)^2`, adding that this was later improved to `bs(f) <= deg(f)^2` by Tal. Hanging `bs(f) <= deg(f)^2` on Nisan-Szegedy is a wrong attribution (Refuted): by Huang's account that is the later, sharper bound of Tal, so the inequality stands and only the credit moves. The papers of Nisan-Szegedy and of Tal are not opened for this page, so nothing here rests on them beyond what Huang's own text says.

## Average sensitivity is exact

The uniform distribution on designs - each of the `2^(2^dim)` codes equally likely - is on classes the orbit-size-weighted distribution, which is precisely the measure [Burnside](https://en.wikipedia.org/wiki/Burnside's_lemma)'s count integrates. So an average-over-designs statement is a statement about the class catalog, reweighted by orbit size.

Under that measure the value at each corner is an independent fair bit.

**Proved.** For a fixed input `x`, the point sensitivity `s(f, x)` - the number of neighbours of `x` that disagree with it - is `Binomial(dim, 1/2)`. The `dim` neighbours are distinct corners, so their bits are independent of `x`'s bit and of each other, and each disagrees with probability `1/2`.

The maximum over `x` does not concentrate, because neighbouring inputs share edges. **Verified** by exhaustive enumeration (`lab/py/boolean-measures`): the probability that a design attains the full sensitivity `dim` is `1/2, 5/8, 69/128, 18253/32768` at dim 1, 2, 3, 4 - between `0.50` and `0.63`, with no sign of climbing to 1.

What does concentrate is the average, the total influence `I(f) = 2^-dim * sum over x of s(f, x)`.

**Proved.** Under uniform-over-designs, `E[I] = dim/2` and `Var[I] = dim/2^(dim+1)`. Let `U` count the bichromatic edges of the hypercube `Q_dim`, the edges whose endpoints get different values. Each such edge is counted once at each endpoint in the sum over `x`, so `I = 2U/2^dim = U/2^(dim-1)`. The cube has `m = dim*2^(dim-1)` edges. The indicator of an edge `{a,b}` is `X_a xor X_b`, a fair bit; and for two distinct edges the two indicators are independent, since together they involve three or four distinct corners and are two linearly independent parities of independent fair bits. Pairwise independence is all the first two moments need: `E[U] = m/2` and `Var[U] = sum of the variances = m/4`, giving `E[I] = dim*2^(dim-2)/2^(dim-1) = dim/2` and `Var[I] = dim*2^(dim-3)/2^(2 dim - 2) = dim/2^(dim+1)`.

Asserting `U ~ Binomial(m, 1/2)` would prove it faster, and that step is false (Refuted). The edge indicators are *not* mutually independent: around any 4-cycle of the hypercube they XOR to 0, so some values of `U` are unreachable. (Verified, `lab/py/boolean-measures`: at dim 3 no design at all has exactly one bichromatic edge, where `Binomial(12, 1/2)` would give `U = 1` probability `12/4096`; the realised values of `U` are `0, 3, 4, 5, 6, 7, 8, 9, 12` and nothing else.) The moments survive because they only need pairwise independence, which does hold - but the Binomial distribution does not, and the proof above uses only what is true.

**Verified** by full enumeration of all `2^(2^dim)` designs in exact rational arithmetic at dim 1, 2, 3, 4 (`lab/py/boolean-measures`): mean `1/2, 1, 3/2, 2` and variance `1/4, 1/4, 3/16, 1/8`, matching the closed forms with no rounding.

The relative spread is `sqrt(2/(dim*2^dim))`, which goes to 0, so a uniformly random design has total influence `dim/2` up to a vanishing fraction. That `E[I] = dim/2` for a random Boolean function is standard; the closed-form variance and its tie to the hypercube edge count are what this page adds, and they are modest.

## Degree in any base

The `GF(2)` degree and the real degree above are base-2 notions. The normal form behind them generalizes to any base, and the generalization round-trips exactly.

A design at base `base` in dimension `dim` is a subset of the residue cube `{0,...,base-1}^dim`. For prime `base`, applying the inverse of the Vandermonde matrix on the points `0..base-1`, axis by axis over `GF(base)`, recovers the unique polynomial that takes the fill indicator's value on every cell, with every variable's exponent capped at `base-1`; at base 2 the inverse is `[[1,0],[1,1]]` and the transform *is* the classical XOR Mobius ANF. For composite `base` the Vandermonde is never invertible mod `base` - its determinant is the product of all `j - i`, so some prime below `base` divides it - and an integer Mobius transform along the product of chains `0 < 1 < ... < base-1` steps in instead, giving a base-agnostic integer degree.

**Proved.** Both transforms round-trip on *every* value table at a given base and `dim`: transform and evaluation are both linear, so exactness is a single matrix identity, `T E = E T = I` over `GF(base)` and `B M = I` over the integers, checked by `lab/py/base-q-anf`. That settles base 3 at dim 3 for all `3^27` value tables and so all `2^27` designs at once, where sampling could not. The same study shows the Vandermonde is not invertible mod the base at base 4, 6, 8, 9, 10, 12 and is at every prime through 13.

**Verified** (`lab/py/base-q-anf`). Design by design anyway: all 16 base-2 designs at dim 2, all 256 at dim 3 and all 65536 at dim 4 round-trip and match the classical XOR ANF coefficient for coefficient - after reversing the variable order, which is exactly how the two labelings differ, recorded so nobody rediscovers it. All 512 base-3 designs at dim 2 round-trip under both transforms, with `GF(3)` degree histogram `-1: 1, 0: 1, 2: 24, 3: 144, 4: 342`. No design has `GF(3)` degree 1: a nonconstant affine function over `GF(3)` takes three values on some line, and an indicator takes two.

A base-3 degree table for the four historical names, generalized as rules - void fills where all coordinates are equal, tree pins every coordinate but one to 0, carpet fills where the coordinate sum is at most 1, net where it is at least `dim - 1` - rules that reduce at base 2 to exactly the classes [the core page](../../core/notes/core.md) names. All eight entries reproduce. (Verified, `lab/py/base-q-anf`.)

| rule | `GF(3)` degree, dim 2 | `GF(3)` degree, dim 3 |
|---|---|---|
| void | 2 | 4 |
| tree | 2 | 4 |
| carpet | 3 | 6 |
| net | 4 | 6 |

Fractal labels on that table need care, and the distinction matters. The degree-3 "carpet" row is the rule keeping 3 of 9 cells; [the Sierpinski carpet](https://en.wikipedia.org/wiki/Sierpi%C5%84ski_carpet) keeps 8 of 9 - every cell but the centre - and its `GF(3)` degree is 4, not 3. The dim 3 row keeps 4 of 27 cells, not [the Menger sponge](https://en.wikipedia.org/wiki/Menger_sponge)'s 20 of 27; the sponge's own `GF(3)` degree happens to be 6 as well, so that label survives only by coincidence - 6 is the ceiling `dim*(base-1)`. (Verified on the actual shapes, `lab/py/base-q-anf`.)

One distinction to keep: the `GF(q)` degree and the integer degree are different invariants and do not order designs the same way - at base 3, dim 2, the tree has `GF(3)` degree 2 and integer degree 1, while the void has `GF(3)` degree 2 and integer degree 4. Any use has to say which it means. (Verified, `lab/py/base-q-anf`.)

## The other complexity: the triangle's spectrum

The second sense of complexity is spectral, and it belongs to the fractal rather than the rule. Take `bang dim 2, code 7` in dimension 2 - the parity rule filling three of the four corners of `{0,1}^2`, tile

```
1 1
1 0
```

- and substitute it into itself `level` times. The array at level `level` is `2^level` on a side with `3^level` filled cells, and a cell `(i, j)` survives exactly when `i AND j = 0`: the Sierpinski triangle. (Verified at level 1..8 by `lab/py/laplacian-degeneracy`: the [Kronecker](https://en.wikipedia.org/wiki/Kronecker_product) construction is executed literally and the `i AND j = 0` characterisation is then checked against it, never assumed.)

Let `G` be the graph on the filled cells at that level, with an edge between axis-nearest neighbours and no diagonals. The operator is the **normalised** Laplacian `I - D^-1/2 A D^-1/2`, at `A` the adjacency matrix of `G` and `D` its diagonal of degrees, whose spectrum lies in `[0, 2]`. (The combinatorial Laplacian `D - A` is a different operator with a different spectrum; the laws below are about the normalised one.)

The spectra demo diagonalises that operator live. It draws the integrated density of states of the triangle at any level to level 6, shades the low window it fits, and prints the distinct count, the degenerate classes, the repeated fraction, the multiplicity of 1 and of the `1 -/+ sqrt(30)/6` pair beside the eight largest clusters, with the combinatorial Laplacian and any other flat code one control away.

The spectrum is extraordinarily degenerate, and the degeneracy is structured.

| `level` | nodes | distinct | degenerate classes | repeated fraction | mult of 1 | mult of `1 -/+ sqrt(30)/6` |
|---|---|---|---|---|---|---|
| 1 | 3 | 3 | 0 | 0.0000 | 1 | 0 |
| 2 | 9 | 7 | 1 | 0.3333 | 3 | 1 |
| 3 | 27 | 17 | 3 | 0.4815 | 9 | 2 |
| 4 | 81 | 43 | 9 | 0.5802 | 27 | 4 |
| 5 | 243 | 111 | 25 | 0.6461 | 81 | 10 |
| 6 | 729 | 289 | 67 | 0.6955 | 243 | 28 |
| 7 | 2187 | 755 | 177 | 0.7357 | 729 | 82 |
| 8 | 6561 | 1975 | 465 | 0.7699 | 2187 | 244 |

(Verified by `lab/py/laplacian-degeneracy`, one generator: the fractal built by Kronecker product, the nearest-neighbour graph checked to be one component, the normalised Laplacian diagonalised densely, and the sorted spectrum clustered by consecutive gaps above `1e-9`. Every row of the table is printed by that run.)

**Verified** (`lab/py/laplacian-degeneracy`). The second family's eigenvalues are `1 -/+ sqrt(30)/6`, matching to within `4 * 10^-15` at every level from 2 to 8, the worst deviation measuring `3.61e-15`. Equivalently, the corresponding eigenvalues of `D^-1/2 A D^-1/2` are the roots of `6x^2 - 5`. This is a numerical identification at machine precision, not a derivation.

**Proved.** The multiplicity of eigenvalue 1 is `3^(level-1)` at every level `L >= 1`: `1, 3, 9, 27, 81, 243, 729, 2187` at level 1..8, one third of the entire spectrum pinned at a single value. Write `T_L` for the graph, `n = 3^L` for its nodes and `nu` for its matching number, the size of a largest set of pairwise disjoint edges.

- **The graph is a tree.** Put `h = 2^(L-1)`. Since `i AND j = 0` forbids both coordinates reaching `h`, the cells fall in three blocks: `A` with `i, j < h`, `B` with `i < h <= j` and `C` with `j < h <= i`. Dropping the top bit shows `B` is `T_(L-1)` shifted by `(0, h)` and `C` is `T_(L-1)` shifted by `(h, 0)`. Blocks `A` and `B` touch only across `j = h-1 | h`, which needs `(i, h-1)` filled, so `i AND (h-1) = i = 0`: the one edge `(0, h-1) - (0, h)`. Likewise `A` and `C` touch only in `(h-1, 0) - (h, 0)`. `B` and `C` never touch, since adjacent cells share a coordinate and the two blocks differ in both. `T_1` is the path `(0,1) - (0,0) - (1,0)`, so by induction `T_L` is a tree with `3^L - 1` edges.
- **The nullity of `A` is `n - 2 nu`.** Expand `det(xI - A)` over permutations. A term is nonzero only if each moved point goes to a neighbour; a cycle of length 3 or more would be a cycle of the graph, so only transpositions along edges survive, each worth sign `-1` and entries `(-1)(-1)`. Hence `det(xI - A) = sum_k (-1)^k m_k x^(n-2k)`, with `m_k` the number of `k`-edge matchings. The powers are distinct, nothing cancels, and the lowest present is `k = nu`: `0` is a root of multiplicity `n - 2 nu`, which for symmetric `A` is its nullity.
- **Eigenvalue 1 of the normalised Laplacian is the nullity of `A`.** A tree on `n >= 2` nodes has every degree at least 1, so `D` is invertible. Then `(I - D^-1/2 A D^-1/2) v = v` exactly when `A (D^-1/2 v) = 0`, and `v -> D^-1/2 v` is a bijection of kernels. The operator is symmetric, so the multiplicity of 1 is that kernel's dimension.
- **`nu(T_L) = 3^(L-1)`, and the corner `O = (0, 0)` lies in every maximum matching.** Induction on `L`. At `L = 1` the path has `nu = 1` and both its maximum matchings use an edge at `O`. For `L >= 2` let `M` be a matching of `T_L` and `M_A, M_B, M_C` its edges inside the blocks; the two link edges join `A` to the corners `O_B = (0, h)` and `O_C = (h, 0)`. If a link is in `M`, its copy misses its corner, so by the induction that copy holds at most `3^(L-2) - 1` edges of `M`: each link gains one edge and costs one. Hence `|M| <= 3 * 3^(L-2)`, and the three copies' maximum matchings together attain it. If a maximum `M` missed `O`, then `M_A` would avoid the corner of `A` (the links sit at `(0, h-1)` and `(h-1, 0)`, never at `O`), so `|M_A| <= 3^(L-2) - 1` and `|M| <= 3^(L-1) - 1`, a contradiction.
- **Conclusion.** The multiplicity is `n - 2 nu = 3^L - 2 * 3^(L-1) = 3^(L-1)`.

**Verified** (`lab/py/tree-matching`). The proof is recomputed rather than trusted, to level 10 (59049 nodes) for the combinatorics: the cells are rebuilt from `i AND j = 0` with no Kronecker step; the graph has `3^L - 1` edges and one component; the block structure above holds (each block equals `T_(L-1)` and exactly the two link edges cross); `nu` agrees between leaf stripping and a Hopcroft-Karp bipartite matching, and removing the corner lowers it by exactly one; and `n - 2 nu = 3^(L-1)`, so the multiplicities at levels 9 and 10 are `6561` and `19683`. At level 1..8 a dense diagonalisation of `A` and of the normalised Laplacian returns the same nullity and the same multiplicity of 1, and at level 1..4 the exact integer characteristic polynomial equals the matching expansion coefficient by coefficient.

**Conjecture.** The multiplicity of each member of the `1 -/+ sqrt(30)/6` pair is `3^(level-3) + 1`, holding at level 3..8: `2, 4, 10, 28, 82, 244`.

**Conjecture.** The cascade continues one step: a third family at `1 -/+ 0.988332421566` carries `3^(level-4) + 1`, holding at level 4..8: `2, 4, 10, 28, 82`. This is weaker than it looks. Six distinct classes at level 8 carry multiplicity 82, so the multiplicity alone does not pick the family out - it is tracked by its eigenvalue across levels, and only then is the multiplicity read. No closed form for that eigenvalue is known.

**Conjecture.** The distinct-eigenvalue count is `2*Fibonacci(2 level) + 1` and the degenerate-class count is `2*Fibonacci(2 level - 3) - 1`, fitting all eight and all seven available levels respectively. Both are existing OEIS entries, read live: `A192908` matches the distinct count at index `level+1` (its name is not a Fibonacci formula - the Fibonacci expression is a formula listed on the entry, contributed by Bruno Berselli), and `A069403`, whose name *is* `a(n) = 2*Fibonacci(2*n+1) - 1`, matches the degenerate-class count at index `level-2`. No novelty is claimed for either sequence. The fit is a fit, not a proof.

Two cautions on the multiplicities that remain Conjecture. They are integers read off a floating-point spectrum: the widest degenerate class at level 8 spans `5.4 * 10^-14`, the class at eigenvalue 1 is isolated by 0.442 on both sides (neighbours `0.5578103942` and `1.4421896058`), and the counts are flat across seven decades of clustering tolerance, `1e-12` to `1e-5` - convincing, and not a proof by exact arithmetic (`lab/py/laplacian-degeneracy`). And level 9 is out of reach by that method: 19683 nodes means a 3.1 GB dense matrix before the eigensolver's workspace. The multiplicity of eigenvalue 1 no longer needs it, being a matching count; the pair families do.

The pair law stays a Conjecture: its eigenvalue is not 0, so no matching count reaches it. The natural reason to expect a mechanism is the graph's self-similarity, and the honest state of that argument is that no spectral-decimation polynomial for this graph is known, so there is nothing yet to run an induction on.

## The flake's band gap

The carpet `bang dim 3, code 23` again, but at base 2, where its tile fills 4 of 8 cells - the class [the core page](../../core/notes/core.md) lists under carpet and net, and the tetrahedron flake of [cuts](../../cuts/notes/cuts.md). The fractal at level `level` is a `2^level` cube with `4^level` filled cells, a cell `(x, y, z)` surviving exactly when `x AND y = y AND z = x AND z = 0`. Its face-adjacency graph is a tree - `4^level` nodes, `4^level - 1` edges, one component - asserted at every level, never assumed.

The operator here is the **unnormalised** combinatorial Laplacian `D - A`. That is deliberately not the operator of the triangle section above: the two have different spectra, and nothing transfers between the two sections in either direction.

**Verified** (`lab/py/flake-band-gap`). The spectrum has a persistent interior band gap, read by a matrix-free route that never forms a matrix, counting eigenvalues by Sylvester's law on an `O(N)` tree elimination and bisecting both band edges; a dense eigensolver on the assembled matrix agrees to `1.1e-14` (Conjecture: the lab carries one generator and the cross-check figure has none):

| `level` | nodes | lower edge `lo(level)` | upper edge | strictly inside |
|---|---|---|---|---|
| 1 | 4 | 1.000000 | 4 | 0 |
| 2 | 16 | 1.827520 | 4 | 0 |
| 3 | 64 | 1.975680 | 4 | 0 |
| 4 | 256 | 1.996862 | 4 | 0 |
| 5 | 1024 | 1.999605 | 4 | 0 |
| 6 | 4096 | 1.999950 | 4 | 0 |

**Verified** (`lab/py/flake-band-gap`). The upper edge is exactly 4, not 4 to float precision: over exact rational arithmetic the elimination of `Lap - 4I` has its only zero pivot at the root, an integer null vector satisfies `(Lap - 4I)v = 0` with residual exactly 0, and exact eigenvalue counts put `3*4^(level-1)` eigenvalues below 2, none anywhere in `[2, 4)`, and a jump of exactly one at 4 - so 4 is a simple eigenvalue and the gap converges to `[2, 4]`.

Two cautions. The gap is *interior*: 4 is not the top of the spectrum, which climbs to about 5.7090 (at level 2 the largest eigenvalue is `3 + sqrt(5)` to twelve decimals - a numerical identification, not a derivation), so the spectrum is a band below 2, a hole, and a band from 4 to about 5.71, split `3*4^(level-1)` below and `4^(level-1)` at or above. And the reading "width `2 = fill - 2`" with `fill = 4` the design's popcount is arithmetic on one design, not a tested law - no second design with another popcount was run, so nothing supports `fill - 2` over any other way of writing 2.

**Conjecture.** The lower edge closes at a rate near `8^(-level)`. The defect `2 - lo(level)` and its successive ratios, carried to level 10 by the matrix-free inertia count (`lab/py/flake-band-gap`, which reaches level 11 and prints `12.986773` there):

| `level` | `2 - lo(level)` | ratio | `(2 - lo) * 8^level` |
|---|---|---|---|
| 1 | 1.000000000000 | | 8.000000 |
| 2 | 0.172480093133 | 5.7978 | 11.038726 |
| 3 | 0.024320402324 | 7.0920 | 12.452046 |
| 4 | 0.003138378624 | 7.7494 | 12.854799 |
| 5 | 0.000395339783 | 7.9384 | 12.954494 |
| 6 | 0.000049510285 | 7.9850 | 12.978824 |
| 7 | 0.000006191639 | 7.9963 | 12.984807 |
| 8 | 0.000000774043 | 7.9991 | 12.986289 |
| 9 | 0.000000096758 | 7.9998 | 12.986658 |
| 10 | 0.000000012095 | 7.9999 | 12.986750 |

The ratios climb monotonically and cleanly toward 8 and the scaled defect flattens, so the fitted form is `2 - lo(level) = c * 8^(-level) + o(8^(-level))` with `c` near 12.9868. It is a fit. No mechanism, no decimation map and no closed form for `c` or for the exact 4 is derived anywhere; `8 = 2^3` is suggestive of a volume scaling and nothing more. Spectral decimation is established mathematics for closely related self-similar trees (Vicsek sets, post-critically finite fractal trees), so the rate and the exact edge may fall out of known machinery - that check is not done here, and the constant is not recognised as anything.

## Spacings: the spectra cluster

One sharp question can be asked of any family of operators: unfolded to unit mean spacing, do the eigenvalues repel, like a random matrix's, or cluster, like an integrable system's? For every mrly fractal tested the answer is cluster - a negative result that closes a thread, and nothing in it concerns the Riemann hypothesis.

The objects (`lab/rs/spectral-spacings`): the carpet cell graph at level 3, 4 (512 and 4096 nodes), the Menger sponge cell graph at level 2 (400), and the sponge's diagonal-slice graph at level 2, 3 (306 and 2250) - the sponge sectioned on the plane `x + y + z = const`, which lives on a triangular lattice and is a genuinely different graph, not a projection. The controls run first, because a pipeline that cannot separate the two laws on known objects cannot be trusted on unknown ones: an Erdos-Renyi random graph reads random-matrix (GOE - 100% of its eigenvalues distinct, the KS test unable to reject GOE at `p = 0.24` and `0.99` while rejecting Poisson at `p = 0.000`) and a `20 x 20` square lattice reads clustered, `P(spacing < 0.5)` 0.54 to 0.59.

**Verified** (`lab/rs/spectral-spacings`). Every fractal reads clustered - on two unfolders, a degree-12 Chebyshev fit of the counting staircase and a local window of 21 levels, and on both Laplacians. `P(spacing < 0.5)` runs 0.50 to 0.69 under the polynomial unfolder across the five fractal rows - carpet level 3 0.52/0.50, level 4 0.58/0.55, sponge 0.67/0.69, slice level 2 0.58/0.60, level 3 0.67/0.66, combinatorial/normalised - and 0.49 to 0.67 under the window unfolder, against a GOE prediction of 0.1783; the band itself is unfolder-dependent, so no numeric band is a property of the spectra. What is stable is the verdict: GOE, and a fortiori GUE, is excluded by a wide margin everywhere, KS `p = 0.000` in every fractal row.

One caveat that matters: the KS distance is nearer Poisson than GOE in every row, and that is weaker than it sounds - on the square-lattice control under the window unfolder the two distances tie to four decimals, 0.4987 against 0.4987, with both laws rejected, and on the most degenerate fractal, the sponge, they read 0.6112 against 0.6040 and 0.6293 against 0.6171. The finding is that the spectra are *more clustered than Poisson*, not that they are Poisson.

The mechanism is the degeneracy this page already documented for the triangle. The random control has no repeated eigenvalue at all; on the combinatorial Laplacian the slice has 62.09% and 63.16% of its eigenvalues distinct at clustering tolerance `1e-9` at levels 2 and 3 (58.50% and 61.29% on the normalised Laplacian; the percentage drifts with the tolerance and must carry it), and its largest multiplicity grows from 12 at level 2 to 48 at level 3 on the combinatorial Laplacian, 18 to 66 on the normalised. The spacing verdict rests on graphs up to 4096 nodes; the 16578-node level-4 slice is confirmed as an object - exactly 16578 nodes and 21546 edges - but its spectrum is not computed.

## The slice's spectral dimension

The same slice graph carries one more spectral question: the [random-walk](https://en.wikipedia.org/wiki/Random_walk) spectral dimension `d_s`, read as twice the low-window slope of the integrated density of states of the normalised Laplacian on the giant component. The reading is `mrly_math::spectrum::spectral_exponent` on the spectrum `mrly_math::spectrum::laplacian_spectrum` returns, wrapped for the section by `mrly_math::six::topology::spectral_exponent`; the eigensolver beneath it is the crate's own dense symmetric pair, Householder tridiagonalisation then implicit QL with Wilkinson shifts.

**Verified**, the structure (`mrly_math::six::topology` test `the_carpet_slice_percolates_at_base_three`). At side 3 the carpet slice percolates: one connected piece at every level, with 42, 306, 2250, 16578 triangles at levels 1 to 4.

**Prior art.** That census is not uncharacterised. With the level-0 value 6 in front, `6, 42, 306, 2250, 16578` is [OEIS A299916](https://oeis.org/A299916) shifted by one: the level-0 term lands on A299916(1), so `census(level) = A299916(level+1)`. The index shift is the whole point and must be stated, because dropping it conflates two counts - A299916 counts hexagram *holes* of the `n`-th descending size in the same cross-section, and this lane counts filled mesh triangles. The gain is a closed form the census lacked: A299916's recurrence `a(n) = 9*a(n-1) - 12*a(n-2)`, signature `(9, -12)`, reproduces every printed term here (`9*42 - 12*6 = 306`, `9*306 - 12*42 = 2250`, `9*2250 - 12*306 = 16578`), and it gives the slice its own dimension, `log((9 + sqrt(33))/2)/log(3) = 1.8184`. Cite the geometry carefully: A299916 is arithmetic in its own terms and its Menger reading is a single 2018 comment, never the entry's definition. The geometry is nevertheless proved by exhaustion in [CLAIMS](../README.md). `https://oeis.org/search?q=id:A299916&fmt=json` returns: name `a(n) = A299914(2n+1).`, offset 0, data `1, 6, 42, 306, 2250, 16578, 122202, ...`, recurrence `a(n) = 9*a(n-1) - 12*a(n-2)`, and the Menger reading present exactly as this page describes it, as a comment - *"a(n) is the number of holes shaped like six-pointed stars, in descending size, found in the cross-section, in the shape of a regular hexagon, of a Menger Sponge. - Albert Säfström, Jul 25 2018"*. Offset, terms, recurrence and the comment-not-definition status are Verified at source, so `census(level) = A299916(level+1)` reads off the entry's own data.

The tree and void slices shatter, largest components pinned at 8 and 6 while the component count climbs. The net slice does neither - it is exactly 12 triangles forming one piece at *every* level, so "largest component stuck at 6 to 8 nodes" is wrong for it (Refuted): it is not a shattering family but one whose slice never grows. And at side 5 the carpet slice itself fails to percolate (20 components, giant 192 of 1164 at level 2), so percolation is base-dependent, not family-dependent. (Verified, `mrly_math::six::topology` tests `the_other_slices_shatter_or_never_grow_at_base_three` and `the_carpet_slice_fails_to_percolate_at_base_five`.)

**Conjecture** - an open observation, not a measurement. A reading of `d_s` about 1.45 to 1.55, robust across levels, is not supported (Refuted). Measured by `mrly_math::six::topology` tests `the_slice_exponent_climbs_with_the_size_it_is_read_at` and `the_deep_slices_hold_their_spectral_exponents`, the low-window exponent is still rising at every reachable size: 0.91 at level 1, 1.25 at level 2, 1.44 at level 3 at the 10 percent window, while a solid-sheet control of the same size reads 1.79 where the true value is 2. The level-3 and solid side 19 rows are the `#[ignore]` release test, run by `cargo test --release -- --ignored spectral`; level 1, level 2, the solid sheet at side 9 reading 1.61 and the base-5 carpet's giant piece of 192 reading 1.05 are pinned by the normal test, every exponent to 0.005. The slice is clearly sub-two-dimensional for transport, but the method's finite-size bias is comparable to every difference being claimed, so neither the value 1.45 to 1.55 nor any distinction from the classical 2D Sierpinski carpet is established: measured like for like, the classical carpet at 4096 nodes reads about 1.6 by the same method and its own square-grid control 1.84 (Conjecture: no lab study rebuilds those two comparators), and the literature offers no single settled value to compare against. Level 4, at 16578 nodes, is a node count here, not a spectrum; the fits sweep five IDOS windows, and none of it converges the exponent. The verdict is the "the slice's spectral dimension" row of the table below - *unconverged at reachable sizes; an open observation*. No page and no ledger line carries `1.45 to 1.55` as measured.

## What is new here

| result | status |
|---|---|
| pin family: six measures equal `r`, `dnf = 1`, `cnf = r` | known family (AND); the geometric handle and the DNF correction are the additions |
| `deg = s^2` extremal pair, `(1 + xy)(1 + zw)` | known construction, located as an orbit of the catalog |
| `C = bs` across both catalogs | empirical for dim at most 4, not claimed in general |
| `E[I] = dim/2` | known, the random-function average sensitivity |
| `Var[I] = dim/2^(dim+1)`, and the orbit-weighting bridge | modest, exact, enumeration-checked; the Binomial derivation is Refuted and replaced |
| the geometry under-determines complexity at dim 4, with witness | new structural fact about this family |
| the triangle's eigenvalue-1 multiplicity `3^(level-1)` | proved at every level: the graph is a tree, the multiplicity is `n - 2 nu`, and `nu = 3^(level-1)` by induction |
| the triangle's other multiplicities | extend a three-point conjecture to eight points; still conjectures |
| the normal form at any base and its degree table | round-trip proved by a matrix identity; the fractal labels checked on the shapes |
| the flake's `[2, 4]` band gap, upper edge exactly 4 | verified, certified in exact integer arithmetic |
| the `8^(-level)` edge rate, `c ~ 12.9868` | a fit over ten levels with no mechanism; conjecture |
| spacings cluster, GOE and GUE excluded | negative result, established up to 4096 nodes |
| the slice's spectral dimension | unconverged at reachable sizes; an open observation |

No claim about `P` versus `NP` is made or implied.

## Where the numbers live

`lab/py/boolean-measures` holds the generator `measures.py` and the two catalog files `measures_d3.csv` and `measures_d4.csv` it diffs against; `lab/py/base-q-anf` holds `anf.py` for the normal form at any base. `lab/py/laplacian-degeneracy` holds `degeneracy.py`, the tolerance sweep, and the reason level 9 stops for diagonalisation; `lab/py/tree-matching` holds `matching.py`, the tree and matching check to level 10 behind the eigenvalue-1 proof. `lab/py/flake-band-gap` holds `flake.py` for the band gap, with the exact-arithmetic certificates and the float sweep to level 11. `lab/rs/spectral-spacings` is the spacing pipeline, both unfolders and both Laplacians. The slice's spectral exponent is a crate witness: `mrly_math::spectrum` carries the eigensolver, the clustering and the log-log fit, `six::topology` reads them on the section, and the slice census itself is `six::census` and `six::topology` in `crates/math`. The spectra demo puts both spectral readings in the browser. The complex-dimensions side of the same fractals is [the dimensions page](../../dimensions/notes/dimensions.md).

The definitions of design, genus, code and fill polynomial are [the core page](../../core/notes/core.md); the identity that licenses this whole page is [the bijection page](../../bijection/notes/bijection.md); the odd-side fill polynomial is worked through end to end in [the method page](../../method/notes/method.md). The standard every sequence in this project has to meet is [the sequence ledger](../../shared/sequences.md) - none of the spectral sequences above are claimed as new entries there, since both are existing OEIS sequences and the remaining multiplicity laws are conjectures rather than verified generators.
