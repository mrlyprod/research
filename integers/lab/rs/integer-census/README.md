# Integer Census

- Which integers MrlyMath writes, which it never writes, and which it writes many times, taken over the whole `ledger` registry rather than over one favourite sequence.
- The registry is 18066 rows, one row a `(design, measure, axis)` key of `ledger::keys` across the four cost tiers: 7692 closed, 5044 convolved, 2665 side grid, 2665 level grid. Every tier count is re-derived from `SPACES`, `ledger::designs` and `Measure::applies` and matched against the enumerator, so the registry size is checked, not trusted.
- The pinned definition, printed by the binary before any table: a row's rendered window is its first `min(48, B)` terms, `B` the leading terms whose footprint fits 100000 cells, under the ledger's own budget of 100000 cells a term; a term's footprint is 1 cell for a closed measure, `number^dimension + level * span` for a convolved measure, `number^(dimension * level)` for a grid measure; a row whose rendered terms are strictly increasing stops at the first term above 100000; row `R` writes `n` iff `n` is a term of `R` inside `R`'s rendered window and `1 <= n <= 100000`.
- Multiplicity counts rows, never `(row, index)` pairs. The two readings differ: 347308 `(row, integer)` incidences against 360703 `(row, index, integer)` incidences, so 13395 times a row writes the same integer twice and the honest reading refuses to count it twice.
- The window is `1..=100000`. Terms at or below zero - the Euler characteristics, the voids of a solid design - are counted (29144 of them) and reported, never folded into the census.
- Truncation is declared and counted, never assumed harmless: 5529 rows stop at the ceiling, 6802 at the 48-term cap, 5735 at a cell budget, and 390 rows write no integer in the window at all.
- The `1..=100000` miss set is a statement about the rendered window. That a missed integer is written by no row at any depth is **Conjecture**: 6802 capped rows have terms this study does not render.
- The window dependence is measured, not assumed: the census is re-read at 8 and 32 rendered terms as well as 48, and the cap is pushed to 96 on the rows a model can extend, for a strict lower bound on a deeper census.
- The model is a Newton forward extension of a row's 8-term head, licensed only where it reproduces that row's whole `written` column from the head and the pinned stop rule alone. Budget-stopped rows carry no rendered length in their head and are declared untestable rather than passed.
- Two self-checks that do not read the sweep they check: every closed-form row is replayed against its own closed form - `f^level`, `g^level - f^level`, the fill polynomial at the odd-side index `k`, the exposure recurrence - over 34604 rendered terms with zero mismatches; and the ledger's own documented terms (`8, 21, 40, 65` for the carpet side fill, `8, 64, 512` and `1, 17, 217` for the carpet level fill and void, `20, 81, 208, 425` for the sponge side fill, `72, 1056, 18048` for the sponge level surface) are found at the head of their named rows and inside the multiset.
- The histogram is folded once and read twice: the sum of the multiset against a second walk over the rows, and the leading champion's multiplicity against a linear scan for its writers.

## RUN

- `cargo run --release --manifest-path research/integers/lab/rs/integer-census/Cargo.toml -- <directory>`
- Three and a half minutes on 8 cores at 4 threads, peak memory under 100 MB.
- Prints the census tables; writes `rows.csv`, `multiset.csv` and `MANIFEST.md` into the directory given, or into a temporary directory it names on stdout. Nothing is written beside this file.
- `rows.csv` carries the key, the tier, the closed form or `none`, the first 8 rendered terms and the full list of integers the row writes, so it rebuilds `multiset.csv` on its own.

## WITNESSES

- The registry: 18066 rows, `closed 7692`, `convolved 5044`, `side 2665`, `level 2665`, each agreeing with its independent derivation, none unread.
- The windows: `1..=1000` never 41, once 31, multiple 928; `1..=10000` never 3589, once 765, multiple 5646; `1..=100000` never 88867, once 2897, multiple 8236. The share written falls `0.9590`, `0.6411`, `0.1113`.
- The champions, by rows written: `16` at 2858, `9` at 2811, `4` at 2559, `12` at 2303, `36` at 2270, `64` at 2176, `3` at 1951, `6` at 1883, `8` at 1790, `33` at 1777. The leader `16` is written by 2858 of the 18066 rows, and every champion of the top twenty is below 65.
- Every integer to 99 is written, and 100000 itself is written by 103 rows.
- The first missed integers: `269, 362, 422, 443, 446, 487, 502, 538, 607, 611`. The least is 269, a prime.
- The miss density by decade: `0` on `1..9`, `0` on `10..99`, `0.045556` on `100..999`, `0.394222` on `1000..9999`, `0.947533` on `10000..99999`. The registry writes almost every small integer and almost no large one.
- The closed-form replay: 34604 terms, zero mismatches, over every row the ledger hands a power, a difference of powers, a fill polynomial or the exposure recurrence.
- The depth: 5263 integers written at 8 rendered terms, 8749 at 32, 11133 at 48, so 5870 of the written set arrive only past term 8 and 2384 only past term 32.
- The model: 3608 rows have a head of finite-difference degree at most 6, by degree `207, 1104, 569, 518, 1202, 0, 8`; the rebuild passes 1306 of 1306 ceiling-stopped rows and 1325 of 1333 cap-stopped rows, the 8 failures being exactly the degree-6 detections, and the 969 budget-stopped rows are untestable.
- The deeper window: extending the 1325 rebuilt cap rows to 96 terms writes at least 11898 integers, moves the first miss from 269 to 362 and the first missed square from `9801` to `38809`.
- The arithmetic: squares `176/316`, cubes `46/46`, fourth powers `17/17`, fifth `10/10`, sixth `6/6`; primes `750/9592` with 158 of the 8363 above 10000; all 2079 residue classes mod `2..64` hold a written integer on `10000..100000`; the 366 perfect powers carry 58906 incidences against a density of `0.003660`.
- The spectrum: 410 distinct multiplicities to 2858; the top twenty carry 39007 of the 347308 incidences; the closed tier covers 7628 written integers with 3983 exclusive; above 30000 the written sets collapse to 953 families of which 875 own a tail integer alone.
- The page these witness is [integers.md](../../../notes/integers.md).
