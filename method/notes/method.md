---
title: Method
lead: How the results here are produced and checked, worked through on the odd-side fill polynomial.
figure: research-method
slug: method
---

The design space is finite. In dimension `dim` there are `2^(2^dim)` designs and nothing else, so a claim about designs is a claim about a finite list, and the honest way to settle it is to walk the list. That one fact sets the method used on every page here: enumerate rather than sample, produce every number twice, pin every formula to something literally drawn, publish the code, and label each claim with what was actually established rather than with how sure it feels.

The [universe demo](/designs/) is the enumeration itself, run in the browser: every orbit per dimension and base, counted by [Burnside](https://en.wikipedia.org/wiki/Burnside's_lemma).

## Exhaust, do not sample

A design is a subset of the `2^dim` corners of the [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) cube, so dimension `dim` holds `2^(2^dim)` of them: 4, 16, 256, 65536 at dim 1..4. Designs related by a symmetry of the cube draw the same shape, and quotienting by the hyperoctahedral group `B_dim` of signed permutations leaves 3, 6, 22, 402 classes. (Verified three independent ways in [the core](../../core/notes/core.md) by `coprime/lab/rs/design-census`, and again by the hyperoctahedral walk of `integers/lab/py/fill-polynomials`, which reproduces all four.)

That is small enough to be brutal with. Where another project would test a formula on the named examples, the sweep here runs every design in the dimension: all 256 in 3D, all 65536 in 4D. A statement quantified over designs is then a finite check rather than an induction, and a counterexample cannot hide in the part of the space nobody drew - a real risk in this family, where the fraction of classes that are nameable as a level-set or an axis pin tends to zero, so almost every design is a compound nobody has drawn. (Proved in [the core](../../core/notes/core.md).)

Enumeration does run out. Past dim 4 the design count outgrows any list, and the counting moves to the Burnside average, which needs only the cycle counts of the group and never builds an orbit. That is the one place the method changes shape: from *checking every object* to *proving a formula and evaluating it*. See [the bijection page](../../bijection/notes/bijection.md), where the class count is proved equal to the number of NP-equivalence classes of Boolean functions, uniformly in `dim`. That this count is the OEIS entry A000616 is Verified there, not proved, and cannot be otherwise.

## Every closed form is pinned to a render

The standing rule for a counting formula is that it is never checked only against another formula. It is checked against an array built cell by cell from the parity rule and summed, with no arithmetic in common.

This is not a slogan; it is the shape of the test suites. The fill-class census compares its closed form against an independently rendered array for every design it censuses, at every side 1..12, and raises on any disagreement. The second generator behind the worked example below builds the `side^dim` grid and counts. In the shipped code the same pattern holds: `crates/math/src/counts` carries the fill engine and its hexagonal projection formulas, and their tests compare each closed form against a rendered cell - `fill_matches_rendered_sum` against a built array, `pro_and_cut_match_census` against an actually projected one. (Verified: `cargo test -p mrly-math` passes, those two among its 515 unit tests.)

A formula proposes; the render disposes. Everything downstream - dimensions, polynomials, densities - inherits its credibility from that comparison.

## What a sequence has to survive

[The sequence ledger](../../shared/sequences.md) states the standard in full and records how each entry met it. In short, a sequence ships only with two independent generators that share no code and no method, a b-file diffed against both over the widest range they can reach, an independent re-verification by a second reader who reruns from the published directory and invents their own checks, and a novelty pass against a local OEIS dump at several windows, shifts and transforms plus the live entry where it matters.

The standard earns its cost by catching things. Re-verification of one entry catches a wrong cross-reference for `27^n` - A001024, which is `15^n` - and replaces it with A009971, powers of 27. (Verified: both entries read from the live OEIS, A001024 reading "Powers of 15" and A009971 "Powers of 27".) Another entry passes on its terms, b-file, closed form and novelty, and still ships marked defective, because its draft's own program block seeds the grid one level too high and so skips a term; the data was right and the generator printed was not. Neither error is visible to a reader who only checks that the numbers look plausible.

## One engine, then a census

The organizing move is to write one script that reads a design's invariants off its definition, then sweep it over the whole space, rather than to derive a formula per named family. The named designs stop being special cases and become rows.

The fill law is what makes this possible. With `E = ceil(side/2)` and `O = floor(side/2)` the number of even and odd residues available to one coordinate, a design `F` fills

```
fill(F, side) = sum over c in F of E^(dim - w(c)) * O^(w(c))
```

cells of the `side^dim` grid, `w(c)` being the number of odd coordinates of the corner `c`. (Proved: a cell is filled exactly when its parity vector is a filled corner, and for a fixed corner each coordinate independently has `E` or `O` admissible values, so the corner contributes that product; distinct corners contribute disjointly. The fill at a level is this raised to that level, proved in [the core](../../core/notes/core.md) from [the Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product).) The law takes a design as data, so one implementation covers the entire space - it is the engine in `crates/math/src/counts/counting.rs` and the engine behind every census in `lab/`.

Three censuses run on that principle, and their honesty is uneven in a way worth stating plainly. `coprime/lab/rs/design-census` sweeps 58 designs across bases 2 and 3 and dimensions 2 and 3, validating each cell by cell and checking its orbit counts against a Burnside average, and writes a 59-line csv. (Verified.) The same study measures 763 designs against a predicted coprime density and flags none. (Verified: 763 design lines, of which 522 are spanning and so inside the density claim and 241 are degenerate and excluded by construction; every per-case summary reports zero flagged.) But that `OK` verdict is numerical agreement inside a flat tolerance at a shallow level, not a proof - the deepest level is capped, and the solid half of that census is the exact finite-level identity, not the limit. `integers/lab/py/fill-polynomials` is the sweep behind the worked example below.

## The three tags

The tags are a claim about evidence, not about confidence.

- **Proved.** A proof is given or restated on the page, and the reader can follow it without running anything. A finite computation in exact or certified interval arithmetic, with its code in the tree and named on the page, counts as a proof, as it does for the four colour theorem. The design count, the fill law, the bijection onto Boolean functions up to NP-equivalence.
- **Verified.** Recomputed from scratch, usually more than once and usually including a route that shares no method with the first. Some claims can never be more than this: that a particular OEIS entry counts what we say it counts is a fact about the entry, checkable against its stated definition and its terms, not provable from our side.
- **Conjecture.** Neither - a fitted pattern or a limit with a proved mechanism and an unclosed step. The coprime densities are conjectures with proved Euler factors. The `3^(level-3) + 1` multiplicity of an eigenvalue pair of the triangle's [Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix) on [the complexity page](../../complexity/notes/complexity.md) fits six levels and is tagged as such even though it looks inevitable, while its sibling law for eigenvalue 1 stood on the same kind of fit until a matching count proved it.

A claim that cannot be re-established does not ship; what cannot be reproduced is dropped rather than softened.

## False friends

Each of these looks like evidence and is not.

- **A `1/zeta(2)` in a density is not a zeta connection.** Every coprimality density carries it, by [Mobius](https://en.wikipedia.org/wiki/M%C3%B6bius_function) inversion over squares. RH is about `zeta` in the critical strip, not about its special values, and `6/pi^2` is not a bridge - it appears in [farey](../../farey/notes/farey.md), [pi](../../pi/notes/pi.md) and [bases](../../bases/notes/bases.md) without any of them containing RH content. This is the likeliest trap on the tree.
- **Fill polynomials are not Jensen polynomials.** Griffin, Ono, Rolen and Zagier is about real-rootedness of polynomials built from Taylor coefficients of `xi`; this tree's come from parity counting. The resemblance is the whole of the connection.
- **Catalan, `phi`, `pi` and `L(2, chi_-3)` are eigenvalues, not omens.** A substitution rule with characteristic polynomial `x^2 - x - 1` produces `phi`, and the tree's own mass laws say so. Constants appearing is the expected behaviour of any sufficiently rich combinatorial system.
- **A Collatz in a Perron certificate is not the `3n + 1` map.** It is Collatz-Wielandt (after Collatz 1942 and Wielandt 1950), the min-max characterization of the Perron root, which `mrly_math::automaton::Automaton::perron` states as "For any nonnegative `T` and any `v > 0`, `min_u (T v)_u / v_u <= rho <= max_u (T v)_u / v_u`", read through the [Perron-Frobenius article](https://en.wikipedia.org/wiki/Perron%E2%80%93Frobenius_theorem). The `3n + 1` map has one topic page, [beneath](../../beneath/notes/beneath.md), section The carry of a Collatz step, which reads the carry of one step as a four-state transducer and never claims the conjecture; its ledger rows and `beneath/lab/rs/carry-skeleton` read that same section.
- **Exhaustive verification is not proof, and one story settles it.** The strong [Mertens](https://en.wikipedia.org/wiki/Mertens_function) conjecture `|M(n)| < sqrt(n)` was believed, verified far, and is FALSE - Odlyzko and te Riele, 1985. Put that in front of anyone proposing to promote a pattern on the strength of exhaustive verification, including this tree's own Conjecture tags.

## Before spending a pass

Five questions, answered in writing before the work starts.

1. Is the object already standard under a different name?
2. Is the proposed theorem stronger than a computation, and can its quantifiers be written in one sentence?
3. What known theorem would make it immediate, and which hypothesis fails?
4. What would falsify it at the smallest nontrivial level?
5. If it is true, who besides this project would cite it, and for what use?

If the answer to the last is only "it resembles RH", keep it as an exhibit. If it is "it gives a reusable digit-restricted sieve or automaton lemma", it belongs in the proof queue.

## What counts as an anomaly

The mixed-product universe is infinite, so brute enumeration alone will manufacture false stories. Every claimed anomaly carries all six of these or it is not one.

1. A formal schedule and canonical factor names.
2. A baseline that fixes every trivial multiplicative invariant.
3. Two independent implementations, one direct-tensor and one coordinate or digit based.
4. An invariance audit under harmless reorderings and symmetries.
5. A stated search universe and a negative-control family.
6. A proof target that would explain the observed difference.

A changed fill, a changed dimension or a familiar constant is not an anomaly. It is an anomaly only when a fixed-control comparison violates a justified structural expectation.

## What a negative result has to name

A negative result is final only when every candidate is counted or a proof covers the space. Testing exactly one candidate state for the connected-component count of a mixed Kronecker word - the four-corner partition of the running product - finds it exact at length 2 and wrong on 20 of 216 words at length 3, and "the component count is not a finite-state function of the code sequence" does not follow: a rank-4 linear representation exists ([connectivity](../../connectivity/notes/connectivity.md)), so that claim is **Refuted**. What the failed candidate actually bounds is the naive geometric state, which does grow like `2^(level-1)` - a true theorem, wearing a false hat. The rule this leaves: name the class of descriptions a negative result rules out, or it rules out nothing.

## Worked example: the odd-side fill polynomial

One theorem, taken through the whole procedure.

**Theorem (Proved).** Fix a dimension `dim` and a design `F` with popcount `p = |F|`. At odd side `2k - 1`, the fill is an integer polynomial in `k` of degree at most `dim` whose coefficient of `k^dim` is `p` - so the degree is exactly `dim` for every non-empty design, and only the empty one falls short. *Proof.* At side `2k - 1` there are `E = k` even residues and `O = k-1` odd ones, so the fill law reads `fill(F, 2k-1) = sum over c in F of k^(dim-w(c)) * (k-1)^w(c)`. Each summand is a product of `dim` linear integer factors, hence a monic integer polynomial of degree `dim` in `k`; a sum of `p` of them is an integer polynomial whose `k^dim` coefficient is `p`, and it is identically zero exactly when `p = 0`.

**Corollary (Proved).** Popcount is invariant under cube symmetry, so the leading coefficient is a class invariant ([the fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension) `log(fill)/log(side)` is one only at even side, where the fill is `popcount * (side/2)^dim`, and not at side 3, where the class of `bang dim 2, code 7` fills 8, 7, 7 and 5 cells, [the core](../../core/notes/core.md)). The lower coefficients are not. A parity flip shifts the infinite tiling by one unit (the tiling of `F` xor `t` is the tiling of `F` translated by `t`) but does not preserve the truncation to `side` cells, so members of one class can fill differently at a fixed side; the polynomial in the table below belongs to the canonical representative, the smallest code in the class.

**Corollary (Proved).** The caveat has an exact witness in the table's own Menger row. `bang dim 3, code 23` is self-complementary - complementing its corner set is a cube symmetry, which is why carpet and net name one class in [the core](../../core/notes/core.md) - but complementing does not commute with truncating. The complement member `bang dim 3, code 232` fills exactly the cells the canonical member leaves void, so its polynomial is the sponge's own void count,

```
(2*k - 1)^3 - (4*k^3 - 3*k^2) = 4*k^3 - 9*k^2 + 6*k - 1
```

as exact polynomials, and vice versa: the canonical polynomial is `bang dim 3, code 232`'s void count. Same leading coefficient 4, the class popcount; different tail. The orbit's eight members carry four distinct polynomials - `4*k^3 - 3*k^2`, `4*k^3 - 5*k^2 + 2*k`, `4*k^3 - 7*k^2 + 4*k - 1`, `4*k^3 - 9*k^2 + 6*k - 1` - one class filling four ways at a fixed odd side. (Proved by expansion; Verified by cell-by-cell counts at `k = 1..9`, each pair summing to `(2k-1)^3`, and by an orbit walk over the 48 signed permutations, `mrly_math::six::topology` test `the_four_families_fill_the_slice_and_name_their_classes`.) The complement's fill sequence `0, 7, 44, 135, 304, ...` is `(k-1)^2 * (4*k - 1)`, OEIS A395241 - one truncation of the sponge's own class, not a new design.

**Corollary (Proved).** The coefficient vector fixes the fill at every side and every level - the ledger's sense of two designs drawing the same fractal. The `dim + 1` polynomials `k^(dim-w) * (k-1)^w` are linearly independent - setting `k = 0` kills every term but `w = dim`, then dividing by `k` and repeating kills the rest in turn - so the polynomial determines how many filled corners carry each Hamming weight, which is the popcount profile that fixes the fill at every side and every level. That is the same lemma the fill-class identity rests on in [the ledger](../../shared/sequences.md), and it makes the number of distinct polynomials in dimension `dim` equal to `Prod_{w=0..dim} (1 + C(dim,w))`, which is the closed form of A129824.

**Verified** (`integers/lab/py/fill-polynomials`). Two generators sharing no code and no method: one sums the closed form over filled corners and interpolates with exact rational arithmetic, the other builds the `side^dim` grid cell by cell, tests each cell's parity vector, and fits by finite differences. Both run over every design - all 4, 16 and 256 at dim 1, 2, 3 - and agree term for term at `k = 1..8`, the grid count matching the closed form on 256 of 256 designs at dim 3 over the six odd sides the census renders. The closed-form sweep extends to all 65536 designs at dim 4; in every dimension the polynomial fitted on `k = 1..dim+1` predicts the true fill out to `k = 10`, the coefficients come out integral, and the leading coefficient equals the popcount with no exceptions. Distinct polynomials number 4, 12, 64, 700 at dim 1..4, which is A129824 at index `dim` (read live: offset 0, terms `2, 4, 12, 64, 700, ...`). The lower coefficients split 4 of 6 classes at dim 2, 20 of 22 at dim 3 and 400 of 402 at dim 4, independently reproducing a caveat the fill-class census records; the leading coefficient splits no class anywhere. A third route agrees: `coprime/lab/rs/design-census`, written separately, carries the level-1 fill at side 1..12 for each 3D class representative, and its six odd columns match both generators on all 22 rows.

The full 3D table, one row per class, ordered by popcount. (Verified as above.)

| design | popcount | fill at side `2k - 1` |
|---|---:|---|
| `bang dim 3, code 0` | 0 | `0` |
| `bang dim 3, code 1` | 1 | `k^3` |
| `bang dim 3, code 3` | 2 | `2*k^3 - k^2` |
| `bang dim 3, code 6` | 2 | `2*k^3 - 2*k^2` |
| `bang dim 3, code 24` | 2 | `2*k^3 - 3*k^2 + k` |
| `bang dim 3, code 7` | 3 | `3*k^3 - 2*k^2` |
| `bang dim 3, code 22` | 3 | `3*k^3 - 3*k^2` |
| `bang dim 3, code 25` | 3 | `3*k^3 - 3*k^2 + k` |
| `bang dim 3, code 15` | 4 | `4*k^3 - 4*k^2 + k` |
| `bang dim 3, code 23` | 4 | `4*k^3 - 3*k^2` |
| `bang dim 3, code 27` | 4 | `4*k^3 - 4*k^2 + k` |
| `bang dim 3, code 30` | 4 | `4*k^3 - 5*k^2 + k` |
| `bang dim 3, code 60` | 4 | `4*k^3 - 6*k^2 + 2*k` |
| `bang dim 3, code 105` | 4 | `4*k^3 - 6*k^2 + 3*k` |
| `bang dim 3, code 31` | 5 | `5*k^3 - 5*k^2 + k` |
| `bang dim 3, code 61` | 5 | `5*k^3 - 6*k^2 + 2*k` |
| `bang dim 3, code 107` | 5 | `5*k^3 - 7*k^2 + 3*k` |
| `bang dim 3, code 63` | 6 | `6*k^3 - 7*k^2 + 2*k` |
| `bang dim 3, code 111` | 6 | `6*k^3 - 8*k^2 + 3*k` |
| `bang dim 3, code 126` | 6 | `6*k^3 - 9*k^2 + 3*k` |
| `bang dim 3, code 127` | 7 | `7*k^3 - 9*k^2 + 3*k` |
| `bang dim 3, code 255` | 8 | `8*k^3 - 12*k^2 + 6*k - 1` |

The two ends are classical and forced. `bang dim 3, code 1` has one filled corner and fills `k^3`, the cubes; `bang dim 3, code 255` is the solid cube of side `2k - 1` and fills `(2k-1)^3`, the odd cubes. (Proved by the theorem, both endpoints; Verified against the live OEIS - A000578, `a(n) = n^3`, and A016755, `a(n) = (2n+1)^3`.) [The Menger sponge](https://en.wikipedia.org/wiki/Menger_sponge) is one interior row, `bang dim 3, code 23`, filling `4*k^3 - 3*k^2`; at side 3, which is `k = 2`, that reads 20, and the celebrated dimension `log(20)/log(3)` is one evaluation of an ordinary row. Two rows are identical: `bang dim 3, code 15` and `bang dim 3, code 27` are distinct symmetry classes - different shapes, related by no cube symmetry - that carry the same polynomial and so fill identically at every side and level. The two classifications cut across each other rather than refining one another - a polynomial can be shared by two classes, and fill is not constant within one - so the 22 classes carry 21 distinct polynomials while the 256 designs carry 64. (Verified, `coprime/lab/rs/design-census`.)

The proof is valid in every dimension and the sweep covers every design; the dim 3 statement checked to `k = 6` on the 22 class representatives is the special case.

**Not claimed.** The other twenty polynomials have not been through the novelty procedure - no dump grep, no shifted windows, no live search - so nothing is asserted about whether they already sit in the OEIS under other names. They are the output of a sweep, not a claim of new sequences, and they are printed here in that spirit.

## The coin: the fill as a biased expectation

The fill law reads as probability. A uniform digit at odd side `N` is a coin: it is even with probability `(1 + 1/N)/2`, so under `chi(b) = (-1)^b` its mean is `mu = 1/N`, and the fill ratio of a design is the chance that `dim` independent such coins land on a filled corner. At even side the coin is fair. Let `f` be the indicator of the filled corners `F`, `w = |F|`, `chi_S(c) = (-1)^(sum over i in S of c_i)` and `hat f(S) = 2^-dim sum over c in F of chi_S(c)`, the Fourier coefficient of [O'Donnell, Analysis of Boolean Functions](https://arxiv.org/abs/2105.10386), Definition 1.2 and Proposition 1.8. Write `W_j = sum over |S| = j of hat f(S)` for the level-`j` sum of the coefficients. It is a signed sum, not the book's Fourier weight `W^j[f] = sum over |S| = j of hat f(S)^2` (Definition 1.19), which returns below.

**Theorem (Proved).** At every odd side `N`, `fill(F, N) = sum_(j=0..dim) W_j N^(dim-j)`. So `fill(F, N)/N^dim = sum_S hat f(S) mu^|S| = T_mu f(0)` at `mu = 1/N`: the noise operator of Definition 2.46 and Proposition 2.47 applied at the all-even corner, which for `N >= 3` is the `p`-biased mean of section 8.4 at `p = (N-1)/(2N)`, written there `E[f^(p)] = T_mu f(1, ..., 1) = f(mu, ..., mu)` in the proof of the Margulis-Russo formula, with `+1` standing for an even digit. *Proof.* With `E = (N+1)/2` and `O = (N-1)/2` the fill law gives a filled corner `c` the count `2^-dim prod_i (N + chi(c_i))`. Expanding the product over the subsets `S` of the axes gives `2^-dim sum_S chi_S(c) N^(dim-|S|)`, and summing over `F` gives `sum_S hat f(S) N^(dim-|S|)`. The `2^dim W_j` are integers; `W_0 = w/2^dim`, the `W_j` sum to `f` at the origin, and `W_dim` is the polynomial at `N = 0`. Substituting `N = 2k - 1` returns the worked example's polynomial: `bang dim 3, code 23` has `W = (1/2, 3/4, 0, -1/4)` and fills `(2N^3 + 3N^2 - 1)/4`, which is 20 at side 3, and `bang dim 2, code 7` has `W = (3/4, 1/2, -1/4)` and fills `(3N - 1)(N + 1)/4`. (Verified, `lab/py/walsh-fill`: 339045 checks of the expansion at odd sides `1..9` and 271236 of `w N^dim` at even sides `2..8`, no mismatch, on literal grids for every code at `dim 1..3` and 2000 seeded codes at `dim 4`, and for all 65536 codes at `dim 4` through a literal count of the cells at each corner.) At even side `mu = 0` and the ratio is `W_0` exactly, which is Even sides are silent, Proved under The staircase on [magic](../../magic/notes/magic.md): every even side already sits at the `N -> infinity` limit of the odd ones.

**The level sums are the Krawtchouk transform of the profile (Proved).** `2^dim W_j = sum_w a_w K_j(w)`, with `a_w` the filled corners of `w` odd coordinates and `K_j(w) = sum_i (-1)^i C(w, i) C(dim-w, j-i)` the Krawtchouk polynomial of Exercise 5.28, because `sum over |S| = j of chi_S(c)` depends on `c` only through `|c|` and equals `K_j(|c|)`. The `W_j` see a design only through its profile, and by the worked example's last corollary, that the polynomial fixes the profile, they fix it back, so this matrix is invertible; it is the classical `K^2 = 2^dim I`, Verified at `dim 1..4`.

**The mirror is `N -> -N` (Proved).** Flipping every coordinate sends `chi_S(c)` to `(-1)^|S| chi_S(c)`, so the mirror `F'` has `fill_F'(N) = (-1)^dim fill_F(-N)`. With `n = (N+1)/2` this is the identity `P_F(1-n) = (-1)^dim P_F'(n)` of [pi](../../pi/notes/pi.md), Pi on the staircase, the roots pair exactly at a palindromic profile, and nothing more.

**`N -> -N` swaps fill and void exactly at a balanced profile (Proved).** With `void_F(N) = N^dim - fill(F, N)`, `fill_F(-N) = +- void_F(N)` as polynomials iff the sign is `(-1)^dim` and `a_j + a_(dim-j) = C(dim, j)` for every `j`, equivalently `W_0 = 1/2` and `W_j = 0` at every even `j >= 2`. *Proof.* `fill_F(-N)` has coefficients `(-1)^(dim-j) W_j` and the void has `[j = 0] - W_j`. An empty or full design never swaps; otherwise both leading coefficients are nonzero, of signs `(-1)^dim` and `+`, so the sign is `(-1)^dim`. Then the coefficients agree identically at odd `j`, force `W_0 = 1 - W_0` at `j = 0`, and force `W_j = -W_j` at even `j >= 2`. In profile terms the mirror fills like the complement, so by the profile corollary the two share a profile, `a_(dim-j) = C(dim, j) - a_j`. Choosing the two levels `j` and `dim - j` together, by Vandermonde, the swapping designs number `prod over j < dim/2 of C(2m_j, m_j)`, `m_j = C(dim, j)`, times `C(m, m/2)` at the middle level `m = C(dim, dim/2)` in even `dim`: 2, 4, 40, 2800 at `dim 1..4`, 46 of the 273 nonempty codes at `dim <= 3`. (Verified on every code at `dim 1..4`; the sign `-(-1)^dim` never occurs.) The self-dual designs, whose mirror is their complement, the 0/1 form of an odd function in Exercise 1.8, swap by construction; the majority rules of odd `dim` are among them, which is why `bang dim 1, code 1` and `bang dim 3, code 23` swap. They are only 2, 4, 16, 256 of the swapping designs, so the criterion is neither self-duality nor the parity of `dim`: `bang dim 3, code 27`, corners `000, 001, 011, 100`, swaps without being self-dual, `bang dim 2, code 3` swaps in even `dim`, `bang dim 3, code 1` fails in odd `dim`, and `bang dim 2, code 7` fails because its weight-1 level holds both corners.

**The roots are zeros of the biased mean (Proved).** `R(mu) = sum_j W_j mu^j = mu^dim fill(1/mu)` is, for `|mu| <= 1`, the mean of `f` under the `mu`-biased coin, and a root `r` of the staircase polynomial `P_F(n)`, the fill at side `2n - 1` as a polynomial in `n = (N+1)/2`, is a zero of `R` at `mu = 1/(2r - 1)`; a root `r = 1/2` is a drop in the degree of `R`, `W_dim = 0`. For `|mu| < 1` every corner has positive probability, so `R > 0` there on a nonempty design. That is the staircase law's real roots in `[0, 1]`, read as: the biased mean of a nonempty design vanishes nowhere on the open interval `|mu| < 1`, and at its ends only where a point mass misses `F`, `r = 0` being `mu = -1` on the all-odd corner and `r = 1` being `mu = 1` on the origin. At `bang dim 2, code 7`, `R = (3 - mu)(1 + mu)/4`, zeros `mu = 3` and `-1`, roots `r = 2/3` and `0`. (Verified on the 77 profiles at `dim 1..3`, worst residual `2.7e-15`.)

**The drift is half the slope of the log biased mean at the fair coin (Proved).** On a nonempty design the staircase law's `drift = dim/2 - mean`, `mean` the average odd count over the filled corners, equals `W_1/(2 W_0) = (1/2) d/dmu log R` at `mu = 0`, because `2^dim W_1 = sum over c in F of (dim - 2|c|)`. Over the bichromatic edges of the cube, `2^dim W_1` counts those filled at the even end minus those filled at the odd end, so `W_1 <= I/2`, with `I = 2^-dim sum_x s(f, x) = U/2^(dim-1)` the sensitivity total influence of [complexity](../../complexity/notes/complexity.md), `U` the bichromatic edges, which is four times `sum_S |S| hat f(S)^2` for the 0/1 indicator; equality holds exactly on the down-sets, the designs closed under making an odd coordinate even. There `drift = I/(4 W_0)`, the Margulis-Russo formula (8.9) at the fair coin; that equality holds nowhere else is the edge count, not the formula. (Verified: the drift identity on all 65808 nonempty codes at `dim 1..4`; equality on all nonempty down-sets, 2, 5, 19, 167 at `dim 1..4`, and on no other code.)

**Noise stability is the fill summed over the flips by filled corners (Proved).** Write `F + c` for the design whose corners are `x + c mod 2`, `x` in `F`. At odd `N`, `sum over c in F of fill_(F+c)(N) = 2^dim sum_j W^j[f] N^(dim-j)` with the book's Fourier weights, so `Stab_(1/N)[f] = 2^-dim sum over c in F of fill_(F+c)(N)/N^dim` (Definition 2.42, Fact 2.48, Theorem 2.49); over all `2^dim` flips the sum is `w N^dim` at every side, so the fill averaged over the flips is `w (N/2)^dim`, the even-side formula, at odd sides too. *Proof.* The flip by `c` multiplies `hat f(S)` by `chi_S(c)`, and `sum over c in F of chi_S(c) = 2^dim hat f(S)` while the sum over every corner is `2^dim [S = empty]`. The sponge's orbit shows it: its members carry the worked example's four polynomials, `(4N^3 + 6N^2 - 2)/8` and its flips `(4N^3 + 2N^2 + 2)/8`, `(4N^3 - 2N^2 - 2)/8`, `(4N^3 - 6N^2 + 2)/8` by corners of 0, 1, 2, 3 odd coordinates, with multiplicities `1, 3, 3, 1`, averaging `N^3/2`. (Verified by literal counts at odd sides `1..9`, every code at `dim 1..3`.)

**A level costs `log2(2^dim/w)` bits at infinite side (Proved).** On a nonempty design at odd side `N >= 3`, knowing that a uniform cell is filled is worth `-log2 R(1/N)` bits per level, and `-log2 R(1/N) = log2(2^dim/w) - (W_1/W_0)/(N ln 2) + O(N^-2)`. *Proof.* `R(1/N) = W_0 (1 + (W_1/W_0)/N + O(N^-2))` with `R > 0` at `|mu| < 1`, and `log2(1 + x) = x/ln 2 + O(x^2)`. At every even side the cost is `log2(2^dim/w)` exactly. The correction is `-2 drift/(N ln 2)`, a rebate when the drift is positive and a surcharge when it is negative. The limit costs are 1 bit at `bang dim 1, code 1`, `log2(4/3) = 0.415037499` at `bang dim 2, code 7` and 1 bit at `bang dim 3, code 23`, and at side 3 they pay `0.584962501`, `0.169925001` and `0.432959407` in that order; `bang dim 1, code 2`, drift `-1/2`, pays `1.584962501` at side 3 against 1 in the limit. (Verified to side 729.)

**Threshold rules (Proved).** The rule "at most `t` odd coordinates" fills `2^-dim sum_(w <= t) C(dim, w) (N+1)^(dim-w) (N-1)^w` at odd `N`, with `2^dim W_0 = S(dim, t) = sum_(w <= t) C(dim, w)` and `2^dim W_1 = (t+1) C(dim, t+1)` by the telescoping `(dim - 2w) C(dim, w) = (w+1) C(dim, w+1) - w C(dim, w)`, so its drift is `(t+1) C(dim, t+1)/(2 S(dim, t))`. The majority, `t = floor(dim/2)`, has ratio `1/2` in odd `dim` and `1/2 + C(dim, dim/2)/2^(dim+1)` in even `dim`, `3/4, 11/16, 21/32, 163/256` at `dim 2, 4, 6, 8`, since the levels pair `w <-> dim - w` and the middle level is left over. At `dim 4` "at most one odd", `bang dim 4, code 279`, fills `(5N^4 + 12N^3 + 6N^2 - 4N - 3)/16`, ratio `5/16`; "at least two odd" is its complement, ratio `11/16`, the majority's ratio because it is the majority's mirror. The mirror keeps `w`, so the rules "at least `s` odd" add no limit ratio, and `S(dim, t)/2^dim = 1 - S(dim, dim-1-t)/2^dim` pairs each ratio with its complement. (Verified at `dim 1..8`, every `t`, by the Walsh transform of each rule.)

**Shared limit ratios (Proved; Verified).** The ratio `S(dim, t)/2^dim` equals the single-corner ratio `2^-m` exactly when the Hamming ball of radius `t` in `{0,1}^dim` holds `2^(dim-m)` points, the sphere-packing equality: at `t = 1` that is `dim = 2^r - 1`, and `t = 2` reduces to `x^2 + 7 = 2^k` with `x = 2 dim + 1`, met at `dim = 90` by `181^2 + 7 = 2^15`. Over `0 <= t < dim <= 2048`, in exact integers, the value `1/2` is held by exactly the odd majorities and 26 further values are shared, each by exactly two rules: the eleven `2^-m` of radius 1 at `dim = 2^r - 1`, `r = 3..11`, of `(dim, t) = (23, 3)`, the binary Golay parameters, and of `(90, 2)`; their complements `1 - 2^-m`; and two pairs off the powers of two, `S(274, 52) = 8 S(271, 51)` and `S(1871, 357) = 8 S(1868, 356)`, with their complements. (Verified, `lab/py/walsh-fill`.)

**Lower-order-free designs (Proved).** `fill(F, N) = W_0 N^dim + W_dim`, with no term between, exactly when `F` holds a fraction `alpha` of every level of even weight and a fraction `beta` of every level of odd weight; then `fill = (alpha (N^dim + 1) + beta (N^dim - 1))/2`. *Proof.* Such a profile fills `alpha ((E+O)^dim + (E-O)^dim)/2 + beta ((E+O)^dim - (E-O)^dim)/2` with `E + O = N` and `E - O = 1`. Conversely a polynomial `A N^dim + B` is met by the profile with `alpha = A + B` and `beta = A - B`, and the profile corollary makes it the only one. Level 0 is one corner, so `alpha` is 0 or 1; in odd `dim` the level `dim` forces `beta` to 0 or 1 as well, which leaves the empty and full designs and the two parity designs `(N^dim +- 1)/2` of [pi](../../pi/notes/pi.md), Pi on the staircase. In even `dim`, `beta C(dim, w)` need only be an integer at every odd `w`, so there are more: for example `bang dim 2, code 11`, corners `00, 01, 11`, fills `(3N^2 + 1)/4`, the design of pi's flip-closed refutation, and `dim 4` holds 140 lower-order-free designs, 136 beyond the empty, full and parity ones. The condition is on level sums, not on the spectrum: none of the 896 bent designs at `dim 4` qualifies, `x_1 x_2 + x_3 x_4` having profile `(0, 0, 2, 4, 0)`. (Verified: 4, 8, 4, 140 lower-order-free designs at `dim 1..4`, each set equal to the criterion.)

**Not claimed.** The expansion, the noise operator, biased means, Krawtchouk polynomials, Margulis-Russo and the sphere-packing equality are classical, read at the definitions cited; what is stated here is their exact dictionary on designs, sides and fills. Whether the two pairs off the powers of two begin an infinite family is open.

## Where the scripts live

Each result on these pages has a study under `lab/` holding its generators and a README saying what it computes, how to run it, and which page lines it witnesses - including the checks that do not resolve in the project's favour. The rule: prose without the code is a claim, not a result.
