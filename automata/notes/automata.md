---
title: The automata
lead: An elementary cellular automaton is a three-dimensional design read as a rule and Life a nine-dimensional one: what the identity buys, and what it cannot.
figure: research-automata
slug: automata
---

An elementary [cellular automaton](https://en.wikipedia.org/wiki/Cellular_automaton) is a three-dimensional [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) design read as a rule, and Life is a nine-dimensional one. A rule table on `k` binary inputs is a subset of `{0,1}^k`, and a design is a subset of the corners of `{0,1}^dim`; put `k = dim` and they are one object. This page states what that identity buys - a class lattice, a Walsh profile, a genus and a dimension for every rule, the Moore neighbourhood as the level-1 carpet tile in every dimension, and a decoupling lemma that says which masks carry a rule space of their own - and what it cannot: geometry does not decide surjectivity, Life is not a product of two elementary rules, and no census here can reach universality.

The generators are `lab/rs/automata-census`, `lab/py/life-census` and `lab/py/sibling-census`. Two demos run the same objects through wasm: the [wolfram app](/wolfram/) draws the 256 rules from one seed beside their cube classes, and the [life app](/life/) runs Life and its siblings on any mask of the tree.

## 1 The identity

Wolfram's rule `N` sends the neighbourhood `(l, c, r)` to `1` exactly when bit `4l + 2c + r` of `N` is set ([Wolfram 1983](../refs.md)). The design `bang dim 3, code N` of [core](../../core/notes/core.md) fills corner `i` of `{0,1}^3` exactly when bit `i` of `N` is set, corner `i` carrying the digits `(x0, x1, x2)` with `i = 4 x0 + 2 x1 + x2`. Under `(x0, x1, x2) = (l, c, r)` the two are the same subset of `{0,1}^3`, bit for bit. (Proved by the two definitions; Verified on all 2048 rule cells over the 256 codes, `lab/rs/automata-census`.)

So every invariant of a design is an invariant of a rule with no translation step: popcount, `GF(2)` degree, genus, Walsh spectrum, fill polynomial and dimension, all from [core](../../core/notes/core.md), [bijection](../../bijection/notes/bijection.md) and [complexity](../../complexity/notes/complexity.md). The rest of this page asks which of them the dynamics respects.

## 2 The class lattice

Four groups act on the 256 truth tables. `R` is the left-right reflection, order 2. `H = <reflection, conjugation>` is Wolfram's equivalence, order 4, where conjugation complements every input and the output together. `B3` is the cube group of [bijection](../../bijection/notes/bijection.md), the signed permutations of the three axes, order 48. `B3 x Z2` adjoins the output complement, order 96. Each count is taken twice, by orbit walk and by [Burnside](https://en.wikipedia.org/wiki/Burnside's_lemma), and the two agree:

| group | order | classes | the count in print |
|---|---|---|---|
| `R` | 2 | 160 | - |
| `H` | 4 | 88 | [Wolfram 1983](../refs.md), derived in [Li and Packard 1990](../refs.md) |
| `B3` | 48 | 22 | [A000616](../refs.md) at `n = 3`, the NP classes of [Harrison 1963](../refs.md) |
| `B3 x Z2` | 96 | 14 | [A000370](../refs.md) at `n = 3`, the NPN classes of [Harrison 1963 IEEE](../refs.md) |

(Verified, `lab/rs/automata-census`; the 22 and the 14 are the dim 3 row of [bijection](../../bijection/notes/bijection.md), so they are Boolean-function counts and nothing about dynamics is claimed for them.)

`H` and `B3` are incomparable subgroups of `B3 x Z2` and meet exactly in `R`. (Proved.) Conjugation complements the output, so it lies outside `B3`; the axis swap `l <-> r` is the reflection, so `R` lies in both; and `B3` contains the flip of one input, which is not in `H`. Hence `H` meets `B3` in `R`, neither contains the other, and the two reductions `256 -> 88` and `256 -> 22` are two branches joined only at the 14, not a chain. The witness is rule 137: the `H` class of 110 is `110, 124, 137, 193`, and 137 lies outside the 24-rule `B3` orbit of 110, so a cube class is not a union of Wolfram classes; the popcounts on that `H` class are `5, 5, 3, 3`, so popcount is not an `H` invariant. (Verified, `lab/rs/automata-census`.)

The free 48-element group is not a symmetry of the lattice. A permutation of `(l, c, r)` that is not the reflection moves the centre off the centre, and the symmetries of the line are translations and reflections, which is why the automata literature classifies by `H` and its neighbourhood-scaling extensions and not by `B3` ([Schaller and Svozil 2025](../refs.md)). A classification of the 256 rules by `B3`, with the 22 representatives and the lattice above, is not found in print; the count itself is the classical NP count and is not new.

## 3 Geometry is not dynamics

Two exact dynamical properties are decided for every rule. Surjectivity of the global map is decided by determinising the labelled de Bruijn graph and asking whether the empty subset is reachable, which terminates on at most `2^4` subsets; injectivity by the pair graph on ordered pairs of de Bruijn nodes with equal output label. Both are decidable in one dimension ([Amoroso and Patt 1972](../refs.md)), and the de Bruijn route is [Sutner 1991](../refs.md); both are undecidable in two dimensions ([NKS note, surjectivity](../refs.md)), which is why nothing below is said about Life.

- Exactly 30 rules are surjective; the exact preimage-balance test on every word of length 1 to 12 shrinks `70, 34, 30, 30, ...` and equals the walk's set from length 3 on. The count is the known one ([NKS note, surjectivity](../refs.md)). (Verified, two generators, `lab/rs/automata-census`.)
- Exactly 6 rules are injective, `15, 51, 85, 170, 204, 240`, every one surjective, so 6 are reversible. The count is on record without the rule numbers ([NKS page 436](../refs.md)); the numbers are this study's. (Verified.)
- The reversible six are exactly the `B3` orbit of 204, the degree-1 single-axis designs `c, l, r` and their complements, so reversibility is constant on every one of the 22 cube classes. (Proved that the orbit of 204 is those six: a signed axis permutation sends the single-axis function `c` to one of `l, c, r` or its complement, and nothing else; Verified that the six reversible rules are that orbit and that reversibility is constant on all 22 classes, 0 mixed.)
- Surjectivity is not constant on the cube classes. Of the 22, 3 are wholly surjective, 18 wholly non-surjective and exactly 1 is mixed: the class of 30, of size 24, surjective on `30, 45, 75, 86, 89, 101, 106, 120, 135, 149, 154, 166, 169, 180, 210, 225` and not on `54, 57, 99, 108, 147, 156, 198, 201`. The witness pair is 30 surjective and 54 not, one orbit apart. (Proved by witness; Verified.)
- Surjectivity is constant on all 88 `H` classes, 0 mixed. (Verified.) So the lattice symmetry preserves both properties and the cube symmetry preserves only one; the geometry of the design fixes reversibility and does not fix surjectivity.

The Walsh transform on this page is taken on `(-1)^f`, `W(S) = sum_x (-1)^(f(x) + S.x)`, so the weight-zero level sum is `Sigma_0 = 8 - 2 popcount` on every rule. (Proved: `W(empty) = sum_x (-1)^f(x)` counts the zeros minus the ones.) Section 6 uses the `0/1` convention instead and says so. The multiset of `|W|` at each character weight is a `B3` invariant on all 256 rules, and the signed level sums `Sigma_j` are not: the `B3` orbit of 110 carries 6 distinct `Sigma` vectors. (Verified.)

The card of rule 110, the rule [Cook 2004](../refs.md) proves Turing universal. Its `B3` orbit is the 24 rules `61, 62, 91, 94, 103, 110, 118, 122, 124, 155, 157, 167, 173, 181, 185, 188, 199, 203, 211, 217, 218, 227, 229, 230`, none surjective, none reversible, all of degree 3, popcount 5 and genus compound, with one Walsh amplitude profile. No invariant of the design separates 110 from a class-mate; the single-seed diagrams of section 5 separate all 24, and the occurring neighbourhood set splits 122 and 218, which meet only the 5 neighbourhoods `00110111`, from the other 22, which meet all 8. (Verified, `lab/rs/automata-census`.)

## 4 The additive bridge

Where the rule is affine over `GF(2)` the space-time diagram from one seed is a design of the tree. Rule 60 read rightward from the seed is `bang dim 2, code 13`; rule 102 read leftward is `bang dim 2, code 14`; rule 90 in the sheared coordinate `j = (t + i)/2` is `bang dim 2, code 13`; each matches cell for cell at levels 1 to 8 and each is the unique match among the fill-3 codes `7, 11, 13, 14`. (Verified, two renderers, `lab/rs/automata-census`.) The Laurent polynomials are the ones of [Martin, Odlyzko and Wolfram 1984](../refs.md): rule 90 is `x + x^-1`, rule 150 is `x + 1 + x^-1`; that paper observes rule 90's dimension `log2 3` and does not name Pascal. Rule 60 is `x_i + x_(i-1)`, the binomial recurrence over `GF(2)`, so its rows from one seed are the rows of [Pascal's triangle](https://en.wikipedia.org/wiki/Pascal's_triangle) mod 2, [A047999](../refs.md), by induction on the row. (Proved.)

Rule 150 is not a design of the tree but an XOR substitution. Its rows are `(1 + x + x^2)^t` over `GF(2)`, and `(1 + x + x^2)^2 = 1 + x^2 + x^4` because the cross terms come in pairs, so `(1 + x + x^2)^(2t) = (1 + x^2 + x^4)^t`: row `2t` is row `t` spread by two, and row `2t + 1` is row `2t` xor its two unit shifts. (Proved; Verified on rows 0 to 128.) Its row populations are [A071053](../refs.md). (Verified, 65 terms.)

Let `P(k)` be the live cells in rows `0 .. 2^k - 1` and `B(k)` the adjacent live pairs summed over the same rows. Then `P(k+1) = 4 P(k) - 2 B(k)` and `B(k+1) = 2 P(k) - 2 B(k)`. (Proved.) Row `2t` has the cells of row `t` at even positions, `|r_t|` of them and no adjacent pair. Row `2t + 1` keeps those and adds the odd position `2i + 1` exactly when one of `i, i + 1` is live in row `t`, which happens `2(|r_t| - B_t)` times, so it has `3|r_t| - 2B_t` cells; its adjacent pairs are `(2i, 2i + 1)` with `i` live and `i + 1` dead and `(2i - 1, 2i)` with `i` live and `i - 1` dead, `2(|r_t| - B_t)` in all. Summing over `t < 2^k` gives the system. The matrix has trace 2 and determinant `-4`, so `P(k+1) = 2 P(k) + 4 P(k-1)` with `P(0) = 1`, `P(1) = 4`, hence `P(k) = c (1 + sqrt5)^k + (1 - c)(1 - sqrt5)^k` with `c = (5 + 3 sqrt5)/10`, equivalently `P(k) = 2^k F(k+2)` with `F` the Fibonacci numbers, and the growth exponent is exactly `log2(1 + sqrt5) = 1.6942419`. (Proved; Verified at `k = 0..12`, `P(12) = 1544192`, with `B(0..11) = 0, 2, 4, 16, 48, 160, 512, 1664, 5376, 17408, 56320, 182272`.) Nothing in this paragraph is new: `P(k)` is [A087206](../refs.md), whose entry already carries the recurrence and `2^n Fibonacci(n+2)`, and whose sum relation to A071053 is on the A071053 entry; the eigenvalue method is [Willson 1987](../refs.md) and the number `log2(1 + sqrt5)` is stated in [NKS note, rule 150](../refs.md), the theorem being Willson's and the number the note's. The study re-derives them.

The 16 affine rules `a0 + a1 l + a2 c + a3 r` split by the weight of the linear part: 2 constants, 6 of one variable, 6 of two, 2 of three. Under `B3` that is five classes, the constants 0 and 255 being fixed by every input flip; under `B3 x Z2` the constants merge and there are four. (Proved by the count; the 16 are Verified by the degree sweep of `lab/rs/automata-census`.)

## 5 The single seed

One live cell on a line of width `2T + 1 + 2 pad` with constant-0 boundary, evolved `T` steps and cropped to the central `2T + 1` window; every depth is run at `pad = T` and at `pad = 2T` and the two must agree cell for cell, which is the stability check.

- The 256 rules give 143 distinct diagrams, 89 up to left-right reflection, with class-size histogram `{1: 118, 2: 13, 4: 4, 8: 4, 16: 4}`, the same at `T = 64, 128, 256` and at both pads. (Verified, `lab/rs/automata-census`.)
- The occurring neighbourhood set, the set of the 8 neighbourhoods a rule ever meets, is measured on the whole evolved padded line and is the same at the three depths, with sizes `{4: 64, 5: 16, 6: 16, 7: 36, 8: 124}`. (Verified.)
- The key of a rule is its occurring set together with the rule restricted to it, measured on that same region. Equal key implies equal diagram, 0 failures over all 65536 ordered pairs: two rules that agree wherever either one is ever asked cannot differ on the seed. (Proved, and the check is the exhaustive one.) The converse fails: the key takes 152 values against 143 diagrams, the 16 rules `23, 31, 55, 63, 87, 95, 119, 127`, all the blinking diagram, and `151, 159, 183, 191, 215, 223, 247, 255`, all ones from `t = 1`, sitting in 2 diagram classes that carry 4 and 7 keys, 11 in all. The split is a boundary artefact. Inside the cropped window the eight rules of either group meet only `000, 001, 010, 100, 111` and agree on all five, so they share one window key; at the ends of the padded line they meet the boundary neighbourhoods `011` and `110` and what those seed, and differ there. Restricting the occurring set to the window only coarsens a key, and equal window keys still force equal window diagrams, because the one cell outside the window that a window edge reads is a background cell whose neighbourhood also occurs inside the window at the same step; so on the cropped window the key takes exactly 143 values, one per diagram. (Proved.)

## Chapters

- [Life, the siblings and the atlas](life-and-the-atlas.md): Life as a nine-dimensional design: the dictionary, the tower, why Life is not a product, three sibling automata, the waves and the retired atlas.
- [The encyclopedia](encyclopedia.md): every basic recipe of mrlylife, named so the name regenerates the run, on the infinite plane as an island or as wallpaper, read whatever its fate.

## Where the numbers live

`lab/rs/automata-census` is the 256-rule census: the class lattice, surjectivity and reversibility, the Walsh profiles, the single-seed diagrams and the additive bridge. `lab/py/life-census` is the `2^18` census and the two XOR rules on the Moore neighbourhood. `lab/py/sibling-census` is the dictionary, the tower, the decoupling census, the composite census and the three siblings. `lab/py/kernel-waves` is the Riesz product and the kill test of section 13, `lab/py/carpet-nodal` its nodal domains; `lab/rs/encyclopedia` is the engine of the encyclopedia, which replaced the atlas of section 14 and replays its mini preset as a control; the [chladni app](/chladni/) shows the Chladni figures of a design. The [wolfram app](/wolfram/) and the [life app](/life/) run the same objects in the browser. The Boolean-function side of the same designs is [bijection](../../bijection/notes/bijection.md) and [complexity](../../complexity/notes/complexity.md); the designs themselves are [core](../../core/notes/core.md).
