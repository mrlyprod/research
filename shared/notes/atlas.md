---
title: The sequence atlas
lead: Every integer sequence the tree writes, walked once and read against itself and the OEIS: 99245 rows, 12486 distinct sequences, a handful of laws that explain the collisions, and the absences ranked.
figure: research-sequences-atlas
slug: atlas
---

MrlyMath writes integer sequences everywhere. The fill law counts the cells of a design, the exposure recurrence counts its surface, the renders count its vertices, edges, faces and Euler characteristic, the hex slice counts triangles, holes and pieces, the carry ladder counts a diagonal cut, the press keeps the integers a digit rule allows, the arc tracer counts loops, and the life sources print their own counts. Each topic met its sequences one at a time. The atlas meets them all at once.

The atlas is one Rust crate, `research/shared/sequences`. It walks every generator, writes every sequence it can reach within a fixed budget, folds equal sequences together under a small shift, links classes that are a simple transform of each other, and searches each class against a stripped dump of the OEIS. This chapter reads what came out: what a sequence of the tree is, how many there are, which ones turn out to be the same sequence under different names and why, which ones are a transform of another, and which ones the OEIS does not hold.

The figure is the whole census on one square. Each cell is one distinct sequence, laid out row by row in order of growth rate (the log of the largest term over the term count), slowest first. One ink marks the sequences that carry two or more names, one the single-name sequences the OEIS holds, one the single-name sequences it does not; sequences too short to search stay blank. The slow corner is crowded with collisions and finds; the fast tail is almost all absences, which is the expected shape: the OEIS is dense where sequences are small and sparse where they explode.

## The object

A sequence of the tree is a list of exact integers with three labels: a name, the generator that printed it, and a measure. The name is a key such as `sequence_dim=3_code=23_measure=surface_axis=side`; the generator is a crate function such as `mrly_math::counts::Exposure::at`; the measure is a short slug such as `surface/side`. A row may also carry a closed form (when the generator proves one) and a stated OEIS id (when a topic already cited one).

### The designs and the two axes

Most rows come from the designs. A design is a code: a set of corners of the residue cube `(Z/base)^dim`, and a cell of the grid is filled when its coordinates reduced mod the base land on a chosen corner. The atlas walks every design of the nine spaces the ledger keeps:

- dim 1, 2, 3, 4 at base 2; dim 1, 2 at base 3; dim 1, 2 at base 4; dim 1 at base 5.
- 1282 designs up to symmetry: 3, 6, 22 and 402 at base 2 (A000616), 4 and 26 at base 3, 6 and 805 at base 4, 8 at base 5.

Each design is read along two axes:

- The level axis: the tile at side `max(base, 3)` is raised to its Kronecker power, level 1, 2, 3 and so on, and the measure is read at each level. At base 2 the carpet `bang dim 2, code 7` fills `8^L` cells at level `L`.
- The side axis: the design is drawn once at growing side, and the measure is read at each side. At base 2 the sides are the odd numbers from 3; at side `2k-1` the sponge fills `4k^3 - 3k^2` cells.

### Orbit mates

Symmetry folds codes into orbits, and the ledger keeps one code per orbit. That is safe for the infinite design but not for the sequences: the level axis at base 2 reads a 3 by 3 tile, and the tile breaks the symmetry. The carpet code 7 fills `8^L`, its orbit mate code 11 fills `7^L`, and code 9 fills `5^L`. The atlas therefore walks every code of each orbit and keeps a mate's rows whenever they differ from the representative's. A mate whose rows agree adds nothing and is dropped before the census.

### The measures

The design rows carry 36 measures, counting each axis apart. Grouped by the generator that writes them:

- Closed: fills and voids by the fill law (`mrly_math::counts::fill`), surface by the exposure recurrence (`mrly_math::counts::Exposure::at`); 54502 and 30708 rows.
- Convolved: the peak and the number of nonzero heights of the diagonal profile (`mrly_math::counts::profile_of_tile`); 5217 rows.
- Rendered: perimeter, vertices, edges, faces, Euler characteristic and hidden faces of the drawn design (`mrly_math::two::census`, `mrly_math::three::census`, `mrly_math::three::hidden`); 6742 rows in the plane, 614 and 132 in space.
- The hex slice of dim 3, base 2: triangles, gaps, the profile fills and voids (`mrly_math::counts::cut_fills` and its siblings, 160 rows each), holes, pieces and the slice Euler characteristic (`mrly_math::six::holes`, `components`, `euler`).

Beside the designs, the atlas walks every other sequence-shaped generator in the tree:

- The carry ladder at bases 3, 5, 7 and dims 2 to 6 (`mrly_math::counts::ladder::ladder`, 15 rows), its fill and its sign.
- The solid slice closed forms (`mrly_math::counts::six::solid_slice_*`, seven rows), the centered hexagonal numbers and the grid triangle counts.
- The design counts, bracelets and fill classes (`mrly_math::bang::baseq`).
- The press at small bases and dims (`mrly_math::press::members`, 160 rows).
- The arc loops drawn and the proved loop laws (`mrly_math::arcs::draw`, 41 rows; `mrly_math::arcs::law`, 15 rows).
- The 23 life sources, each read twice: as printed, and with its zeros and ones dropped.
- The rows of the ledger page `../sequences.md` as rendered, 93 of them, so that the page and the walk can be compared.

### The budget

Every rendered term costs cells, so the walk sets a budget: 16 terms a sequence, 600000 cells a rendered term, 2 seconds a sequence. A closed or convolved row reaches all 16 terms; a rendered level row stops at the cell budget, often at 4 to 6 terms. 4278 rows end under 8 terms. A short row joins a class only when exactly one class contains it, so the budget never invents a collision.

## The census

The crate prints its census on every run. The numbers below are that run's.

- Rows: 99245, of which 99152 are generated and 93 are read from the ledger page.
- Distinct sequences: 12486, after folding equal terms under a shift of up to 3.
- Classes with two or more names: 5764.
- Transform relations between classes: 1176.
- Classes long enough to search: 10831; with an OEIS record: 1391; with none: 9440.
- Absent classes with 8 or more terms that are not constant, linear or geometric: 9383. This is the ranked list of candidates for new OEIS entries.
- Ledger rows disagreeing with the walk: 0. Loop laws disagreeing with the drawn loops: 0.

The 5764 collisions split by kind:

- Explained, 3412: one measure shared by codes or mates, a source read twice, a loop law and its drawing, or plane surface and perimeter, which are the same count.
- Trivial, 85: constant, linear or geometric sequences, which collide by accident.
- One design under two measures, 177: a single code whose two different measures print the same sequence.
- Across topics, 16: a design measure equal to a press, an arc count, a source or a slice closed form.
- Across codes and measures, 2074: two different codes under two different measures.

The last kind is the largest because the side axis is polynomial. Every side-axis row at base 2 is a polynomial in `k` by the fill law, and low-degree polynomials with small coefficients are few, so unrelated codes meet. The atlas registers these meetings and leaves them; the interesting ones are the 177 and the 16.

## One sequence, many names

A collision is a question: why are these two sequences equal? Some answers are trivial bookkeeping, some are short proofs, and some are only verified. This section gives the answers the atlas found, the proved ones first.

### Two double counts

The commonest one-design collisions join a rendered count to the closed ones. They are double counts, and they hold for every design.

- In the plane: `edges = 2 fills + surface/2`.
- In space: `faces = (6 fills + surface)/2`.

Proof. A unit edge of a filled square lies in one or two filled squares. Count the pairs (square, edge of that square): each square has 4 edges, so the count is `4 fills`; each interior edge is counted twice and each boundary edge once. The boundary edges are exactly the exposed sides that the surface counts (in the plane the surface is the perimeter). So `4 fills = 2 interior + surface`, and `edges = interior + surface = 2 fills + surface/2`. In space the same count runs over faces: each cube has 6 faces, each face lies in one or two cubes, the faces in one cube are the exposed ones, so `6 fills = 2 interior + surface` and `faces = (6 fills + surface)/2`.

The census bin of the figure, `research/shared/census/src/bin/research-sequences-atlas.rs`, checks both laws against every row of the atlas: 843 plane designs and 3448 terms on the level axis, 842 and 13472 on the side axis, 78 solid designs and 312 terms on the level axis, 80 and 1280 on the side axis, every term agreeing. The consequence is that edges in the plane and faces in space are closed wherever fills and surface are, on both axes. For the carpet the edges are `2 x 8^L` plus half the perimeter, and the perimeter obeys the proved recurrence `a(L) = 11 a(L-1) - 24 a(L-2)`; the atlas finds the carpet's edges in one class with the surface of a base-4 code, `bang dim 2, base 4, code 45211`, whose exposure recurrence is the same.

### Edges against surface in space

In space the edges are not a double count of the faces, because a unit edge can lie in up to four cubes. The atlas found 26 solid designs on the side axis whose edges are exactly twice their surface, the sponge and the empty design among them, and 52 whose edges are not. The rule that separates them is a local count.

Proof. Around a unit edge of the grid sit four cell positions. Say `j` of them are filled. Count the exposed faces that contain the edge: 2 when `j = 1`; 2 when `j = 2` and the two cubes share a face; 4 when `j = 2` and the two cubes are diagonal, touching only along the edge; 2 when `j = 3`; 0 when `j = 4`. Each exposed face has 4 edges, so summing over the edges of the solid counts every exposed face 4 times:

- `4 surface = 2 edges + 2 diagonal - 2 full`,
- `edges - 2 surface = full - diagonal`,

where `full` counts the edges inside four cubes and `diagonal` the edges inside exactly two diagonal cubes.

Now take a base-2 design on the side axis at odd side `n`. A cell is filled when its parity vector is a corner of the code. An interior edge parallel to axis `a` at height `h` along `a` sees four cells whose other two coordinates run over two consecutive integers each, so their parities run over all of `{0,1}^2`, while their coordinate along `a` has parity `h mod 2`. The four positions therefore copy the code's slice at parity `h mod 2` across axis `a`. An edge on the boundary of the box sees at most two cells inside it, and those two share a face, so it is never full or diagonal. Each axis has `(n-1)^2` interior lines of edges; along each line `(n+1)/2` heights are even and `(n-1)/2` are odd. Write `s(a,p) = +1` when the slice at parity `p` across axis `a` is full, `-1` when it is a diagonal pair, `0` otherwise. Then

- `edges - 2 surface = (n-1)^2 sum_a ((n+1)/2 s(a,0) + (n-1)/2 s(a,1))`.

The right side vanishes at sides 3 and 5 only if `sum_a s(a,0) = sum_a s(a,1) = 0`, and then it vanishes at every side. So `edges = 2 surface` at every odd side exactly when the full and diagonal slices cancel at each parity.

For the sponge, `bang dim 3, code 23`, the corners are the parity vectors with at most one odd coordinate. Its slice at parity 0 holds three corners and its slice at parity 1 holds one, across every axis, so no slice is full or diagonal and `edges = 2 surface` exactly. The census bin checks the whole formula, not only its zero case, against all 78 solid designs the atlas renders on the side axis, at sides 3 to 33: 1248 terms, every one equal to the formula, and the formula picks out exactly the 26 designs the atlas found.

This is a closed form for the edges of a solid design on the side axis, proved here: the surface is closed by the exposure recurrence, and the edges are the surface twice plus a correction read off six slices of the code.

### The sponge's surface, and two neighbours

The side axis at base 2 is simple enough that a surface can be counted by hand. Exposure is six faces per cube minus two per adjacent pair.

- The sponge at side `n = 2k-1` fills `4k^3 - 3k^2` cells (the fill law). Two cells adjacent along `x` have `x` coordinates of opposite parity, so both are filled only when the other two coordinates are even: `k^2` choices, times `2k-2` positions along `x`, times 3 axes. The surface is `6(4k^3 - 3k^2) - 12(k-1)k^2 = 6k^2(2k-1)`, six times A015237 at `k`; the edges are twice that.
- `bang dim 3, code 126` keeps every corner but the all-even and the all-odd one, and fills `6k^3 - 9k^2 + 3k` cells. A pair along `x` is filled only when its other two coordinates have mixed parity: `2k(k-1)` choices, times `2k-2`, times 3. The surface is `6(6k^3 - 9k^2 + 3k) - 24k(k-1)^2 = 12k^3 - 6k^2 - 6k`, twelve times A160378 at `k`.
- `bang dim 3, code 136` keeps the two corners whose last two coordinates are odd, so it is a set of separate rods of length `n` along the first axis, `(k-1)^2` of them, each with surface `4n + 2`. Its surface is `2(k-1)^2(4k-1)`, and the sponge's voids are `(2k-1)^3 - (4k^3 - 3k^2) = (k-1)^2(4k-1)`, A395241 at `k-1`. The rods expose exactly twice as many faces as the sponge leaves empty.

The atlas reaches each of these three surfaces at sides 3 to 33 (`mrly_math::counts::Exposure::at`, 16 terms each), and each agrees with its closed form at every term.

### A Fibonacci peak

The diagonal profile of a design counts its filled cells on each plane `x_1 + ... + x_dim = s`. The atlas records the largest of these counts at each level, the peak. One peak is the Fibonacci sequence: `bang dim 2, base 4, code 1291` has peak `1, 2, 3, 5, 8, 13, 21, 34, 55` at levels 1 to 9, and the atlas puts it in one class with the Fibonacci source. That is a theorem.

Proof. The code holds the corners `(0,0), (0,1), (0,3), (2,0), (2,2)` of the 4 by 4 residue square, whose weights `x + y` are 0, 1, 3, 2, 4: one cell on each weight from 0 to 4. The profile at level `L` is the product of `W(t^(4^j))` over `j < L` with `W(t) = 1 + t + t^2 + t^3 + t^4`, so its coefficient `r_L(s)` counts the ways to write `s = d_0 + 4 d_1 + ... + 4^(L-1) d_(L-1)` with every digit between 0 and 4. Reduce `s` mod 4. When `s` is not a multiple of 4 the last digit is forced, and `r_L(s) = r_(L-1)((s - d_0)/4)`. When `s = 4m` the last digit is 0 or 4, and `r_L(4m) = r_(L-1)(m) + r_(L-1)(m-1)`. Of the two consecutive integers `m - 1` and `m`, at most one is a multiple of 4, so one of the two terms is itself a single term at level `L-2`. Writing `M_L` for the peak, `M_L <= M_(L-1) + M_(L-2)`, with `M_0 = M_1 = 1`, so `M_L <= F(L+1)`. For the other direction, suppose `{r_L(m), r_L(m-1)} = {F(L+1), F(L)}`. Then `r_(L+1)(4m) = F(L+2)`, and one of `4m + 1` or `4m - 1` carries `r_L(m)` or `r_L(m-1)`, whichever is `F(L+1)`, so the consecutive pair at `4m` carries `{F(L+2), F(L+1)}`. The chain starts at `L = 1`, `m = 1`, where both counts are 1. So the peak at level `L` is `F(L+1)` for every `L >= 1`.

The same argument works for any base `q >= 2` and any tile of side `q` with one cell on each of `q + 1` consecutive weights and no other cell: the only freedom is the one extra digit `q`, which can carry once, so the peak is Fibonacci whatever the base. Outside base 2 the atlas reads the peak of one code per orbit, and the diagonal is not symmetric under the orbit, so other tiles of this kind are not all walked; each has the Fibonacci peak by the same proof.

### One function under two names

The solid slice closed forms in `mrly_math::counts::six` carry two functions, `solid_slice_core_edges` and `solid_slice_interior`, whose bodies both return `36k^2 - 42k + 12`. The atlas puts their 16 terms in one class, and so this is a sequence of the tree under two names in the strictest sense: the same code. Their class matches A033570 after division by 6, and `solid_slice_edges`, `36k^2 - 30k + 6 = 6(2k-1)(3k-1)`, is six times A033568.

### The sponge's slice, read on the side axis

The hex slice of the sponge on the side axis switches with the side mod 4. These are verified, not proved:

- Gaps at sides `4m + 3` equal `solid_slice_edges(2m+1) = 6(2m+1)(3m+2)`, at sides 3 to 31, 8 terms.
- Pieces are the centered hexagonal number `3m(m+1) + 1` (A003215) at sides `4m + 1`, and 1 at sides `4m + 3`, at sides 3 to 33.
- Holes are A003215 at `m` at sides `4m + 3`, and 0 at sides `4m + 1`, at sides 3 to 33.

So the slice is one connected piece with a hexagonal number of holes at every fourth side, and a hexagonal number of pieces with no holes at the sides between. The switch is a parity of `(n+1)/2`, and none of the eight searched transforms undoes it, which is why the sponge's hole count is absent from the OEIS as a whole while each half is a classic.

### Across topics

The 16 classes that join two topics are mostly small:

- The hex slice holes of `bang dim 3`, codes 126, 62 and 30 at levels 1 to 4 are 1, 2 and 4 times `(6^L - 1)/5` (A003464, A005610, A125687), the base-6 repdigits that the press of codes 2, 4 and 16 keeps. The holes are rendered at 4 levels only, so this is verified at 4 terms; the repdigits are the press's own definition.
- The arc loops of base 3, code 85 at levels 1 to 6 are `4^n - 3^n` (A005061), the voids of `bang dim 1, base 4, code 7`.
- The Mersenne numbers `2^n - 1` (A000225) are the base-2 press of code 2 and the arc loops of base 2, code 9, whose law the arcs topic proves.
- The level-axis voids of one-dimensional designs reappear as heights of the diagonal profile in dimension 3: `3^n - 2^n` (A001047), `(3^n - 1)/2` (A003462), `(4^n - 1)/3` (A002450).
- Three base-5 side rows of dimension 1 are the press of base 6 under a digit rule, A039215, A037465 and A007706.
- The Fibonacci peak above.

## The transforms

Beyond equality the atlas looks for one class that is a simple transform of another: the negative, the first differences, the partial sums, `a - 1`, `a + 1`, `a/2`, the even-indexed terms and the odd-indexed terms. A transform relation needs 7 agreeing terms and joins classes with different labels. The run finds 1176:

- `a/2`: 432; the odd-indexed terms: 329; the even-indexed terms: 236.
- `a - 1`: 74; `a + 1`: 46.
- Partial sums: 29; differences: 19; the negative: 11.

In 1129 of the 1176 both classes hold a side-axis row, and in 1070 every name on both sides is a side-axis row. The side axis at base 2 is polynomial, and halving or bisecting a low-degree polynomial with small coefficients lands on another one often enough that most of these relations are coincidences between unrelated codes, the same effect that fills the 2074 cross-code collisions.

146 relations stay inside one design or cross a topic, and these carry the structure:

- `a/2`, 44 of them: 21 halve the edges of a solid design on the side axis to its surface, the `edges = 2 surface` case of the slice rule; others halve the edges of one code to the surface of another, such as the edges of code 22 to the surface of code 126.
- The even-indexed terms, 51, and the odd-indexed terms, 23: the hex slice counts that switch with the side mod 4 bisect to a classic, as the sponge's gaps, pieces and holes do above, and a few side rows of dimension 3 bisect to plane rows.
- `a - 1`, 14, and `a + 1`, 7: fills and voids that differ by one cell, and Euler characteristics one off a fill count.
- The negative, 2: the Euler characteristic of the star `bang dim 2, code 6` on the side axis is minus its fill count.

Three measures are now known to be transforms of the others on whole families: plane edges and solid faces by the double counts, solid edges by the slice rule. The double counts never show as relations, because they need two measures at once; the slice rule shows only in its zero case.

## The OEIS face

The search reads a stripped dump of the OEIS, one line per record. For each class it builds four 5-term windows at offsets 0 to 3 for eight transforms: `a`, `-a`, the differences, the partial sums, `a - 1`, `a + 1`, `a/2`, and `a` divided by the gcd of its terms. A record is a hit when it contains a window and agrees with the class on every term both carry.

### What checks

Every OEIS id that the tree already states checks by value. 116 rows of the atlas state 70 distinct ids, from the ledger page, the life sources and five crate generators (the carry ladder, the centered hexagonal numbers, the design counts, the bracelets and the fill classes), and the stated record contains each row's terms; the run prints 0 failures. That is the strongest statement the atlas makes about the existing pages: no topic cites a record that does not hold its sequence.

1391 classes meet a record. By the first transform that hits: 682 as printed, 275 by their differences, 244 after division by the gcd, 86 halved, 42 negated, 30 and 19 one off, 13 by partial sums. The gcd transform is where the solid side axis lives: the sponge's surface is 6 times A015237, code 126's is 12 times A160378, and the solid slice edges are 6 times A033568, and none of these is in the OEIS as printed.

### What is absent

9440 searched classes meet no record, and 9383 of them have at least 8 terms and are neither constant, linear nor geometric. A dump absence is a conjecture until a live read; it is also mostly expected. By the name each line of the list prints first, 4303 of the long absences are level-axis surfaces, and a surface is a linear recurrence with a code-specific characteristic polynomial, so every new polynomial is a new sequence. Another 2403 are side-axis surfaces and 863 side-axis voids, polynomials in `k` whose coefficients depend on the code.

The value of the list is in its head. These eight are the most beautiful absences by the run's judgment, each a natural object with a clean law:

- The carry ladder at base 3, dim 4: `1, 6, 132, 1848, 29040`, order 2 by the proved order law of the cuts topic. It already has an OEIS draft, so it heads the list as drafted.
- The carry ladder at base 5, dim 3: `1, 18, 414, 9702, 227646`, the base-5 sibling of A299916; its 16 terms fit `a(n) = 25 a(n-1) - 36 a(n-2)`.
- The carry ladder at base 3, dim 5: `1, 30, 1000, 35700, 1321600`, order 3 as the order law says; its 16 terms fit `a(n) = 60 a(n-1) - 914 a(n-2) + 3120 a(n-3)`.
- The sponge's hex slice holes on the side axis: `1, 0, 7, 0, 19, 0, 37`, the parity switch above.
- The perimeter of the base-3 carpet, `bang dim 2, base 3, code 495`, at odd sides: `16, 28, 44, 72, 92, 116, 160`; its 16 terms fit the characteristic polynomial `(x-1)(x^3-1)^2`, a quasi-polynomial of period 3.
- The level-axis surface of `bang dim 3, code 126`: `60, 840, 13920, 244560`, with the proved recurrence `a(L) = 23 a(L-1) - 90 a(L-2)`, roots 18 and 5.
- The integers whose base-6 digits are all 1 or 5, the units mod 6: `1, 5, 7, 11, 31, 35, 43, 47` (the press of dim 1, base 6, code 34).
- The hex slice holes of `bang dim 3, code 126` on the side axis: `1, 7, 19, 31, 55, 73, 109`; the 16 terms fit `(x-1)^3(x+1)^2`, two quadratics interleaved by parity.

Each fit above is a conjecture until a transfer matrix proves it; the code 126 surface recurrence is already proved by the exposure recurrence, and the ladder orders by the cuts order law. These, with twelve more from the head of the absence list, are the first candidates for new OEIS entries.

## What the atlas changes

- The ledger page and the walk agree on every row they share: the page holds no sequence the generators do not print.
- Three measures that were render-only are closed on whole families: plane edges and solid faces everywhere, solid edges on the side axis at base 2. The renders stay as checks.
- The sponge's surface and edges on the side axis, and the surfaces of codes 126 and 136, have proved closed forms with OEIS names after a scale.
- One crate redundancy is named: two solid slice functions are the same polynomial.
- The absence list is ranked, so one list replaces each topic searching on its own.

## Open questions

- The hex slice holes of codes 126, 62 and 30 match the base-6 repdigits at 4 levels. Is the hole count of code 126 at level `L` exactly `(6^L - 1)/5`, and what in the slice makes the base 6?
- The sponge's slice gaps at sides `4m + 3` equal the solid slice edges. Is there a map from the gaps of one to the edges of the other?
- The slice rule closes solid edges on the side axis at base 2. The level axis is a Kronecker power, and the same local count should close edges there through the tile's edge configurations; the vertices need the vertex configurations, eight cells around a point.
- The 2074 cross-code collisions are mostly coincidences of low-degree polynomials. Which of them survive past 16 terms, and is any one a theorem?

## Regenerate

- Run: `OEIS_STRIPPED=<stripped dump> OEIS_NAMES=<names file> cargo run --release --manifest-path research/shared/sequences/Cargo.toml` from the repository root, about 80 seconds.
- Output: `data/research/shared/sequences/atlas.json` (rows, classes, relations, failed ids), `atlas.md` (one table sorted by generator) and `absences.txt` (the long absences, one line each).
- Without the two environment variables the walk, the classes and the relations still run, and the search is skipped.
- The figure: `cargo run --release --manifest-path research/shared/census/Cargo.toml --bin research-sequences-atlas` reads `atlas.json`, checks the two double counts and the slice rule against every row, and writes `research/shared/figures/research-sequences-atlas.json`.
- Claims: `../claims/atlas.md`.
