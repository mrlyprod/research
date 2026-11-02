---
title: Integers
lead: The census the other way round: which integers the whole registry writes, which it never writes, and which it writes thousands of times.
figure: research-integers
slug: integers
---

Which integers this work writes, which it never writes, and which it writes many times - taken over the whole registry rather than over one favourite sequence. The [sequences](../../shared/sequences.md) ledger asks whether a given sequence is known; this page asks the opposite question, what the union of every sequence the registry holds covers, and answers it by census.

The generator is `lab/rs/integer-census`, one pass over `ledger::keys`, which prints its definition before any table and writes `rows.csv`, `multiset.csv` and a manifest into a directory given on the command line. Every number below is a line of that run. The registry it walks is the one [sequences](../../shared/sequences.md) is rendered from and the [sequences demo](/sequences/) searches live; the closed forms it replays are the fill law and the exposure recurrence of the sequence-census paper. The integers demo reads that union integer by integer: which of the first thousand the designs write, how many rows write each, and which the pinned window misses. The [plot demo](/sequences/) draws a row of the same ledger rather than listing it, with the smallest linear recurrence its terms satisfy, its characteristic polynomial and its growth beside it.

## The definition

The census is only as good as its window, so the window is pinned and printed, never assumed.

- A registry row is one `(design, measure, axis)` key of `ledger::keys` over the four cost tiers.
- A row's rendered window is its first `min(48, B)` terms, `B` the leading terms whose footprint fits `100000` cells, under the ledger's own budget of `100000` cells a term.
- A term's footprint is 1 cell for a closed measure, `number^dimension + level * span` for a convolved measure, `number^(dimension * level)` for a grid measure.
- A row whose rendered terms are strictly increasing stops at the first term above `100000`; the count of rows truncated this way is printed, never assumed to lose nothing.
- Row `R` writes `n` iff `n` is a term of `R` inside `R`'s rendered window and `1 <= n <= 100000`.
- Multiplicity counts rows, not `(row, index)` pairs: a row writing `n` at several indices counts once.
- An integer `n` appears iff some row writes it, and is missed iff no row writes it.

## The census

- **The registry is 18066 rows and every tier count is derived twice. Verified.** 7692 closed, 5044 convolved, 2665 side grid, 2665 level grid, each matched against an independent count over `SPACES`, `ledger::designs` and `Measure::applies`; no row goes unread. Truncation is declared and counted: 5529 rows stop at the ceiling, 6802 at the 48-term cap, 5735 at a cell budget, and 390 write no integer in the window at all.
- **The two readings of multiplicity differ, and the honest one is smaller. Verified.** There are 347308 `(row, integer)` incidences against 360703 `(row, index, integer)` incidences, so 13395 times a row writes the same integer twice and the row reading refuses to count it twice. 29144 rendered terms are at or below zero - [Euler characteristics](https://en.wikipedia.org/wiki/Euler_characteristic), the voids of a solid design - and are excluded and reported, never folded in.

| window | never | once | multiple | written | share written |
|---|---|---|---|---|---|
| `1..=1000` | 41 | 31 | 928 | 959 | 0.9590 |
| `1..=10000` | 3589 | 765 | 5646 | 6411 | 0.6411 |
| `1..=100000` | 88867 | 2897 | 8236 | 11133 | 0.1113 |

| decade | width | missed | miss density |
|---|---|---|---|
| `1..9` | 9 | 0 | 0.000000 |
| `10..99` | 90 | 0 | 0.000000 |
| `100..999` | 900 | 41 | 0.045556 |
| `1000..9999` | 9000 | 3548 | 0.394222 |
| `10000..99999` | 90000 | 85278 | 0.947533 |

- **The written set is finite, so the miss density tends to 1. Proved.** A row renders at most 48 terms, so whatever the ceiling, the registry writes at most `48 * 18066 = 867168` integers. The registry is a fixed finite object and the integers are not: past `867168` the census is almost all miss, and no growth of the ceiling changes that. This is the one statement on the page that survives any change of window.

## The miss set

- **Every integer to 268 is written and the first miss is a prime. Verified.** 269 is missed, and `1..268` is the longest written run in the window; the longest missed run is 447 wide, on `95265..95711`, with 95264 and 95712 both written. The ceiling itself, 100000, is written by 103 rows. The first thirty misses are `269, 362, 422, 443, 446, 487, 502, 538, 607, 611, 618, 626, 643, 653, 659, 668, 677, 691, 698, 701, 709, 723, 758, 773, 787, 797, 803, 835, 857, 878`.
- **The miss set is not a union of residue classes. Refuted.** All 2079 classes mod `2..64` hold a written integer on `10000..100000`, exhaustively; no modulus in that range separates written from missed.
- **The bias is divisibility and smoothness, not congruence. Verified.** On the tail `10000..100000` the written count by residue mod 12 runs `1175, 440, 145, 194, 715, 176, 420, 224, 531, 358, 229, 116`, a ratio of `10.13` between residue 0 and residue 11, and mod 6 it runs `1595, 664, 676, 552, 944, 292`, a ratio of `5.46`. Sorted by greatest prime factor the written share on the same band falls `0.5798`, `0.1406`, `0.0506`, `0.0313`, `0.0117` across the bands `1..10`, `10..100`, `100..1000`, `1000..10000`, `10000..100000`. Of the 9592 [primes](https://en.wikipedia.org/wiki/Prime_number), 750 are written and only 158 of the 8363 above `10000`; the first missed prime is 269, the first missed integer.
- **Every cube, fourth, fifth and sixth power is written; the squares are not. Verified.** Cubes `46/46`, fourth powers `17/17`, fifth `10/10`, sixth `6/6`, all of `1..100000`. Squares run `176/316`: every square to `98^2 = 9604` is written and `99^2 = 9801` is not, and the largest written square is `97969 = 313^2`, written by exactly one row, `sequence dim 4, code 28662, voids, side`, whose closed form `4k^4 - 8k^3 + 8k^2 - 4k + 1` is the square of the centered square numbers.
- **The square frontier is the cap, not arithmetic. Verified.** Row multiplicity at `96^2, 97^2, 98^2, 99^2, 100^2` is `321, 19, 480, 0, 123`: the dense square families are exhausted, not excluded. Deepening the window to 96 terms, in the section below, writes at least 228 of the 316 squares and moves the first missed square from `9801` to `38809 = 197^2`. Oddness excludes nothing: `97969` is odd and written, and `9801` is missed for want of depth.

## The depth of the window

The whole miss set is a statement about the rendered window, and the study measures how much of one rather than asserting it is harmless.

| rendered window | written | missed | first miss |
|---|---|---|---|
| 8 terms | 5263 | 94737 | 269 |
| 32 terms | 8749 | 91251 | 269 |
| 48 terms | 11133 | 88867 | 269 |

- **More than half the written set arrives past the head. Verified.** 5870 of the 11133 written integers appear only past term 8 and 2384 only past term 32, so a census read off the ledger's own 8-term heads sees less than half of what 48 terms see, and its miss set starts `269, 281, 302, 311` rather than `269, 362, 422, 443`.
- **A row's written column is rebuilt from its head and the stop rule alone. Verified.** 3608 rows have a head whose finite differences terminate at order 6 or less, by degree `207, 1104, 569, 518, 1202, 0, 8`. Extending the head by Newton forward differences and applying the pinned stop rule reproduces the row's `written` column exactly for 1306 of 1306 ceiling-stopped rows and 1325 of 1333 cap-stopped rows; the 8 failures are exactly the rows whose head degree reads 6, which eight terms cannot certify. The 969 budget-stopped rows carry no rendered length in their head and are not testable this way, which is said rather than hidden.
- **Deepening the cap moves every window-relative number except the longest missed run. Verified.** Extending only the 1325 cap-stopped rows the rebuild reproduces, out to 96 terms, gives a strict lower bound on the 96-term census: at least 11898 integers written, the first miss moved from 269 to 362, the longest written run at least 361, at least 228 of 316 squares. 269 becomes written; the run `95265..95711` does not move.
- **That a missed integer is written by no row at any depth is Conjecture.** 6802 rows are cut by the cap and their deeper terms are not rendered here; the 96-term reading is a lower bound, not a census, and the true frontier of the written set is not known at any depth.

## The champions

| rank | integer | rows | rank | integer | rows |
|---|---|---|---|---|---|
| 1 | 16 | 2858 | 6 | 64 | 2176 |
| 2 | 9 | 2811 | 7 | 3 | 1951 |
| 3 | 4 | 2559 | 8 | 6 | 1883 |
| 4 | 12 | 2303 | 9 | 8 | 1790 |
| 5 | 36 | 2270 | 10 | 33 | 1777 |

- **The whole top of the census is small and mostly a power. Verified.** All twenty champions lie below 65 - ascending, `2, 3, 4, 6, 7, 8, 9, 12, 14, 15, 16, 18, 20, 21, 24, 25, 33, 36, 49, 64` - and they carry 39007 of the 347308 incidences, a share of `0.1123`. The 366 perfect powers of `1..100000` carry 58906 incidences, a share of `0.1696` against a density of `0.003660`: `46.34` times their weight.
- **The champions are not the divisor-rich integers. Refuted.** On `1..1000` the mean row count is `193.42` over all integers, `995.26` over the squares and `920.58` over the perfect powers, but only `170.60` over the 413 integers with at least eight divisors - below the overall mean. Being a small perfect power is what a champion is; being highly divisible is not, and reads slightly against it.
- **Multiplicity is not driven by each row's first term. Refuted.** Dropping every row's first rendered term removes 17036 of the 347308 incidences, `4.9%`, and changes nothing that matters: the written set stays 11133, the never counts stay `41`, `3589`, `88867` in all three windows, and the leaders stay `36` at 2212, `64` at 2112, `16` at 2000, `9` at 1999 - the same integers in a different order.
- **The multiplicity spectrum is neither geometric nor a power law. Refuted.** With `S(m)` the count of integers written by at least `m` rows, `S(1) = 11133` and `S(2) = 8236` give a ratio `0.7398`, which predicts `S(64) = 6.312e-5` against the observed 977 - wrong by seven orders. The spectrum takes 410 distinct values with a maximum of 2858.
- **The effect is a property of a measure column, not of a design. Verified.** `euler.side` writes 1 in 695 of its 859 rows, `peak.side` writes 12 in 809 of its 1261, `heights.side` writes both 9 and 33 in 765 of its 1261. A champion is an integer that one reading of the geometry returns for most designs at once.
- **The one champion that is neither small-smooth nor a power is an offset. Verified.** 33 ranks tenth and is `3 * 11`. The eight integers below 100 that `heights.side` writes most often are `9, 17, 25, 33, 41, 49, 57, 65`, every one of them `1 mod 8`: the column runs arithmetic progressions whose common difference is a power of two, and `33 - 1 = 2^5`. The arithmetic of the champion is the arithmetic of the step, not of the integer.
- **The closed tier carries the census and the grid tiers carry its tail. Verified.** Of the 11133 written integers the closed tier covers 7628 with 3983 exclusive, the side grid 6203 with 2603, the level grid 1826 with 541, the convolved tier 792 with 130. Above `30000` there are 2174 written integers and the closed tier covers 1853 of them.
- **The tail is not a few dominant families. Refuted.** Restricted to `30000..100000` the rows' written sets collapse to 953 distinct families, and 875 of those families own a tail integer no other family writes - between them 2005 of the 2174 tail integers. So every cover of the written tail needs at least 875 families, and the tail is a wide superposition rather than a handful of dominant sequences.

## Against the OEIS

Every search below is exhaustive on both sides: every window of the census sequence is indexed and every record of a local copy of [the OEIS stripped dump](https://oeis.org/stripped.gz) is walked against that index, so no sampling of offsets is involved. The dump read holds 398817 records. Under the standing caveat of [sequences](../../shared/sequences.md) a dump is a snapshot, so every absence here is **Conjecture** and needs a live re-read before it is repeated.

- **The miss set is new to the OEIS only in its dense head. Conjecture.** No record carries any 4-term window of the miss set at offsets `0..416`; the first hit is at offset 417, in [A049537](https://oeis.org/A049537). Above that the miss set does hit, and the hits are near-interval records rather than identifications: 852 hits at window length `k = 4`, 130 at `k = 10`, 37 at `k = 15`, 15 at `k = 20`, the 20-term witnesses being [A112820](https://oeis.org/A112820) and [A118471](https://oeis.org/A118471), each a sequence that runs a block of consecutive integers through a region the census misses wholesale. The head is the informative part and is clean: the string `,269,362,422,443,` appears nowhere in the dump.
- **The write-once set is absent at every offset. Conjecture.** The 2897 integers written by exactly one row have no hit at any offset of any record at `k = 4, 10, 15, 20` - a cleaner absence than the miss set's, because the once set is thin where the miss set is an interval complement.
- **The champion set meets two records for exactly twelve terms. Verified.** Ascending, the twenty champions open `2, 3, 4, 6, 7, 8, 9, 12, 14, 15, 16, 18`, which is a window of [A100290](https://oeis.org/A100290) and of [A336231](https://oeis.org/A336231) and of no other record. All three part at the thirteenth: A100290 gives 21, A336231 gives 19, the census gives 20. Both records are binary-digit conditions, which is the right neighbourhood - the registry's designs are corner subsets of a [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) cube - and neither is the champion set.
- **The written-per-decade run meets one record and parts at the next term. Verified.** The written counts by decade are `9, 90, 859, 5452, 4722` with 100000 itself, summing to 11133. The prefix `9, 90, 859` sits inside [A209631](https://oeis.org/A209631) alone, an exponential-transform array, which continues 6689 where the census gives 5452.
- **No recognizable family is systematically missed. Verified.** 173 records hold at least ten distinct integers of `1..100000` and lie wholly inside the miss set, the longest being [A361796](https://oeis.org/A361796) at 41 terms. At a miss density of `0.88867` a 41-term run of misses has probability about `10^-2.1`, which 398817 records make ordinary: the census excludes nothing a catalogue would recognise, it just runs out of depth. The tour demo runs the other way in a dozen cards, each drawing a design beside the sequence it counts and the OEIS record that holds the terms.

## What is left

- Whether any integer of `1..100000` is written by no row at any depth. The 96-term reading is a lower bound and already moves at least 765 of the misses across, 269 among them; the honest frontier needs a cap the dimension-2 side grid can pay for, and that tier costs `cap^3`. **Conjecture.**
- Whether the miss set has any arithmetic characterisation at all. No modulus to 64 separates it, no growth order does, and it is closed under nothing; the only theorem on offer is the finiteness bound above. **Conjecture.**
- Whether the 953 tail families are 953 rules. The families are de-duplicated by written set and not by generating rule, so two rules with equal truncated value sets merge and 953 is a lower bound on the number of rules, never an upper one. **Conjecture** that no bounded union of named families reaches the written tail.
- Whether the multiplicity function, `a(n)` the number of rows writing `n`, is worth an entry. It is absent from the dump, but it is a function of the registry's own shape - the tier mix, the cap, the ceiling - and not of `n` alone, so it is a reading of this instrument and not a sequence of the integers. **Conjecture** that no reparametrisation of it is submittable.

## Chapters

- [The field ladder](field-ladder.md): The fill polynomial as a product of norm forms, one per factor of the weight enumerator, and the ladder of number fields the designs reach: every number field appears at some dimension.

## The rest of the tree

- [README](../../README.md) is the front door: the parity cube, [the Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product), and the index of every page.
- [sequences](../../shared/sequences.md) is the ledger this page is the complement of: which sequences are known, against which integers are reached.
- [method](../../method/notes/method.md) - how a claim here is produced and checked, worked through on the odd-side fill polynomial.
- [CLAIMS](../README.md) - where every line above is tagged with its witness and its refutation attempt.
- [REFS](../refs.md) - every sequence id above resolved to a canonical URL.
