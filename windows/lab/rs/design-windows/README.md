# Design Windows

- Counts the `k x k` windows of the subshift of every plane code at bases 2 and 3, decides which subshifts are of finite type, and checks both against brute force, for [windows](../../../notes/windows.md).
- A code at base `b` is a bitmask over the `b x b` cells of level 1 in row-major order; the substitution sends a 1 to that tile and a 0 to the empty tile.
- `L_2`, the `2 x 2` render windows, is the closure of the level-1 windows under one substitution step, stopped at the first repeat.
- The census: for `k <= b^n + 1` the `k x k` windows are the blocks of `sigma^n(B)`, `B` in `L_2`; each block is named by the quadruple of the names of its two `(k-1)`-blocks on the diagonal and its two other corners, so equal names are equal blocks and the table size is `p(k)`.
- One code per orbit of the eight symmetries of the square; the census asserts every image of every orbit's code gives the same count for `k <= 12`.
- A boundary code, all filled cells in one boundary line of the box, has the subshift `{0}` and prints `p = 1`, with the render-window count after the word `edge`.

## VERBS

- `count <base> <top> [code]`: one line per orbit, `fill`, the size of `L_2`, the level it settles at, `periodic@k` at the first `k` with `p(k) <= k^2/2` or `above`, the counts `p 1..top`, the exact quadratic when one holds from `k <= 3` to `top`; the carpet's and the gasket's formula checks; the class tallies and the number of distinct sequences.
- `scan <base> <level> <top> [code]`: renders each orbit's code with `mrly_math::bang::factory::create`, asserts it equals the study's own `sigma^level(1)` and that code `7` at base 2 drawn at side 3 is base-3 code `495`, counts distinct `k x k` blocks of the render by their rows, and prints the gaps against `count`.
- `kind <base>`: `sft empty`, `sft full`, `sft boundary`, or `not-sft` with a witness of fast growth of `L_j(X_2)`: `blocks` (two blocks of side at most 2 freely arranged), `lines` (a branching window rule on pictures constant along an axis or a diagonal), `xor` (`L_2` is the eight blocks with one corner the sum of its neighbours); then `recognised m` for the least `m <= 3` at which every occurrence of `R_m` in a picture `sigma^m(B)`, `B` in `L_2`, is at one of the four aligned offsets, else `unrecognised`; it is run only off the empty, full and boundary codes.
- `pairs <base> <code>`: `L_2` as rows `ab/cd`.

## RUN

- `cargo run --release --manifest-path research/windows/lab/rs/design-windows/Cargo.toml -- count 3 82`: under 7 seconds.
- `count 3 244 495`: the carpet to `k = 244`, about 2 seconds; `count 2 129`: under a second.
- `scan 3 6 12`: about 46 seconds; `scan 2 7 17` under a second; `scan 3 7 16 495`: about 8 seconds.
- `kind 2` and `kind 3`: under a second each. Prints only, writes nothing.

## WITNESSES

- [windows](../../../notes/windows.md), "The objects": the symmetry check, the `102` and `6` orbits, the carpet as code `495`.
- "Counting through the substitution": the settle levels `1, 2, 3` on `6, 75, 21` orbits at base 3; the scan agreement, `0` gaps in `108` orbits and on the carpet at level 7.
- "The carpet": `p(k)` to `244`, `L_2` of size `10`, `formula holds for 2 <= k <= 244`.
- "Growth": the `recognised 1` certificate on `91` orbits at base 3 and the gasket, `unrecognised` on the line orbits; the exact quadratics; `formula holds for 1 <= k <= 129` for the gasket.
- "Products and lines": the counts of codes `27`, `45`, `325`, `63`, `365`.
- "Not of finite type": every verdict and witness; the counts `p(6) = 242`, `p(24) = 4570`, `p(72) = 42970`.
- "The census": the base-2 table, `83` distinct sequences, the tallies `26`, `10`, `476`.
