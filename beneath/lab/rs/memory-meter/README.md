# memory-meter

- The Mobius meter of a memory rule: `M_W(x) = sum of mu(n)` over the integers `n <= x` whose binary digit word the rule accepts, against that set's own mass `A_W(x)`.
- Every width-`1`, width-`2` and width-`3` rule at dim `1`, base `2` is read: `4`, `16` and `256` codes, `276` rules in all, the width-`1` codes being the memoryless row.
- An integer `n >= 1` is its minimal base-`2` word, coarsest digit first, no leading zero; `0` is excluded from every sum. The rule reading is `mrly_math::memory::Rule`: a window of `k` digits is `w = sum_j c_j 2^(k - j)`, bit `w` of the code is set when that window is allowed, a word is accepted when every one of its `k`-windows is allowed, and a word shorter than `k` is accepted.
- No exponent is fitted. Printed per rule and phase are `A_W(x)`, `M_W(x)`, `max abs M_W(t)` over `t <= x`, the ratios `M_W/sqrt(A_W)` and `max abs M_W/sqrt(A_W)`, with `kappa`, `rho` and `card W` from the crate.

## THE METHOD

- One ascending pass over `n <= 2^30` carrying the **window profile** of `n`, the bitmask of the `2^3` windows its word contains, by `profile(n) = profile(n >> 1) | bit(n mod 8)` for `n >= 4`; `n` is accepted by `W` exactly when `profile(n)` is a subset of `W`.
- The width-`2` and width-`1` profiles are induced from the width-`3` one: for a word of three digits or more every `2`-window is a prefix or a suffix of some `3`-window and every digit lies in one, so one `u8` per `n` carries all three widths. The words of one and two digits are entered by hand.
- The mass per profile is a `256`-bucket count and `A_W` is its subset-sum transform at each phase; the meter is carried per rule so the running maximum is exact, and the subset-sum transform of the per-profile `mu` sums is asserted equal to it at every phase.
- Phases are `x = floor(2^(level + j/4))` for `level 8..30` and `j = 0..3`, `89` of them, `2^30` last.
- `mu` is `mrly_num::factor::mobius_sieve`, the crate's linear Mobius sieve; `mrly_num::sieve` carries the Sierpinski word and no Mobius, so nothing there is reused.
- `rho` and `kappa` come from `mrly_math::memory::perron` and `kappa`, which split the digraph into strongly connected components and make each component's Perron root exact against its integer characteristic polynomial, so every printed root is an algebraic integer of the right minimal polynomial.

## THE CONTROLS

- The full line, `k = 1` code `3`, is the Mertens function: its meter reads `-1, 1, 2, -23, -48, 212, 1037, 1928` at `10^1..10^8`, asserted, which is [A084237](https://oeis.org/A084237).
- The memoryless base-`3` designs are enumerated directly from the same sieve and asserted against [mobius-designs](../../../../mobius/lab/rs/mobius-designs/): digits `{0,1}` read `(M, max abs M) = (11, 105)` at `level 14`, `(149, 173)` at `level 16` and `(-30, 312)` at `level 18`, digits `{1,2}` read `(-1461, 1582)` at `level 18`. Digits `{0,2}` at `level 20` wants `3^20`, past `2^30`, so it is printed at `level 14, 16, 18` and not pinned. The four pinned pairs are read at source in [design-meter](../../../../mobius/lab/py/design-meter/), which computes them and cites [mobius-designs](../../../../mobius/lab/rs/mobius-designs/) as their census.
- The profile recurrence is asserted against a direct digit recount on all `276` rules below `2^20`, and profile containment against `Rule::accepts` on all `276` rules below `2^12`.
- Code `7` at `k = 2`, the golden rule forbidding `11`, opens `1, 2, 4, 5, 8, 9, 10, 16, 17, 18, 20, 21`, which is [A003714](https://oeis.org/A003714) without its zero, and its mass is the Fibonacci number. Code `11`, forbidding `10`, opens exactly the Mersenne numbers [A000225](https://oeis.org/A000225) without its zero, one per level.
- Leading zeros: a rule is **zero-closed** when prepending one zero to the word changes no membership below `2^20`. Both masses are printed per rule and the criterion is asserted code for code.

## RUN

- `cargo run --release --manifest-path research/beneath/lab/rs/memory-meter/Cargo.toml`
- About thirty seconds on one thread, `23s` of it the sieve; peak resident set about `2.15 GB`, one `i8` and one `u8` per integer. Prints only, writes nothing.
- An optional depth argument runs a shallower table at the same phases: `--manifest-path research/beneath/lab/rs/memory-meter/Cargo.toml -- 28` takes about seven seconds and `0.6 GB`, and its rows are the first `81` phases of the deep run verbatim.

## READS

- Every ratio below is a reading at a named phase, never a fit and never a sweep-wide claim unless it says so. `rho` and `kappa` are read off the transfer matrix of the **word** language while `A` and `M` are read off the **integer** set; on a rule that is not zero-closed those are different objects, and the columns sit side by side for that reason.
- The census is the `53` rules of all three widths holding at least `10^4` integers below `2^30`; the floor is part of the census and is fixed before any reading. The same rules under a floor of `64` number `106`.
- The full line at phase `30.00`: `A = 1073741824`, `M = -10374`, `max abs M = 11173`, `M/sqrt(A) = -0.316589`, `max abs M/sqrt(A) = 0.340973`. Three codes carry it, `k = 1` code `3`, `k = 2` code `15` and `k = 3` code `255`.
- The golden rule, `k = 2` code `7`, `kappa = 0.098239`, `rho = 1.618033989`: `A = 2178309`, `M = 551`, `max abs M = 716`, ratios `0.373329` and `0.485125` at phase `30.00`. Code `14`, forbidding `00`, reads `A = 3524576`, `M = -466`, ratios `-0.248218` and `0.454355`.
- The three named `k = 3` least codes at phase `30.00`: `23` supergolden, `kappa = 0.115204`, `A = 125492`, ratios `-0.448837` and `0.731125`; `54` plastic, `kappa = 0.260981`, `A = 13579`, ratios `-0.360425` and `0.677943`; `127` tribonacci, `kappa = 0.056639`, `A = 98950096`, ratios `-0.432878` and `1.239625`.
- The normalised peak over the census **at phase `30.00`** spans `[0.293624, 1.239625]`, least on `k = 3` code `125` and largest on `k = 3` code `127`. **Over all `89` phases** the same quantity spans `[0.500000, 2.169240]`, the top on `k = 3` code `232` at phase `16.75`; a last-phase span is not a sweep-wide span.
- Read against the full line at the same phase, the factor over the census is `[0.861136, 3.635552]` at phase `30.00` and reaches `6.375774` on `k = 3` code `190` at phase `12.75` over all phases.
- Grouped by `kappa` over the census the means read `0.340973` on `[0, 10^-9)` with `3` rules, then `0.622171` on `6`, `0.581446` on `14`, `0.644840` on `12` and `0.645966` on `18` for the bands `[10^-9, 0.1)`, `[0.1, 0.2)`, `[0.2, 0.3)` and `[0.3, 0.5)`. Every band of positive coupling sits above the `kappa = 0` band, and among them the means are not monotone in `kappa`, so the coupling does not order the spread. Restricted to `k = 3` the same bands read `0.340973, 0.698386, 0.581446, 0.644840, 0.645966` on `1, 4, 14, 12, 18`.
- The full line's own normalised peak runs `[0.272410, 0.500000]` over the grid, the ceiling at phase `8.00`, which is the grid's first point. That ceiling is a property of where the grid starts and not of the full line: below the grid the same ratio reads `1.000000` at `x = 1`, `0.894427` at `5`, `0.832050` at `13`, `0.718421` at `31` and `0.565685` at `200`. Every band statement is therefore read beside the same-phase factor, which needs no grid.
- The falsification fires. Five of the `53` census rules never enter that band at any phase where they hold `10^4` elements, all of them above it: `k = 3` codes `159`, `182`, `190`, `218` and `250`, holding `211116`, `13607`, `31535`, `59860` and `4126645` integers. Under the floor of `64` the count is `32`, the other `27` all holding under `500` elements.
- `16` of the `53` census rules attain their sweep-wide normalised peak in the last quarter of the phases, from `24.75` on. No rule at any width with at least `1000` elements has either ratio rise at every one of the last eight phases; that test asks `max abs M` to grow about `9%` per quarter-level across two levels, so an empty answer carries little, and the late-peak count is the informative statistic.
- Zero-closed under one prepended zero: `3` of `4`, `8` of `16` and `64` of `256`, exactly the codes allowing the digit `0` together with the empty code, the codes allowing `01`, and the codes allowing both `010` and `011`. Under any number of prepended zeros the word language agrees with the integer set on `2` of `4`, `4` of `16` and `16` of `256` codes, so at `k = 3` the two readings part company on `240` of `256`.
- Equal mass is not the same set. `k = 2` code `14` and `k = 3` code `126` both hold `28655` integers below `2^20` and both carry `rho = 1.618033989`, but they share only `1077`: the symmetric difference is `55156`, and `4` is the least integer in code `126` and not in code `14`.

## WITNESSES

- [beneath](../../../notes/memory-meter.md), The memory meter - the meter table, the control band, the `kappa` bands and the zero-closed criterion.
- research/claims/ the memory-meter rows.

## COLUMNS

- `control` the pinned lines; `classes` and `reps` the `88` orbit representatives at `(1,3)` under `G_(1,3)`; `rule` one line per code with `card W`, `rho`, `kappa`, the zero-closed flag and both masses, the last-phase reading and a track of the normalised peak at `level 8, 12, 16, 20, 24, 28, 30`; `row` one line per code and phase; `top`, `bottom`, `span`, `factor`, `latepeak`, `kappaband`, `band` and `rho` the printed bands at two mass floors, `gridstart` the full line below the grid, `pair` the two equal-mass rules, `reading` the word language against the integer set; `climbing` the monotone scan; `run` the depth and the runtimes.
