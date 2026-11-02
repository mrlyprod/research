# carry-census

- The numbers of [progressions](../../../notes/carry-criterion-census.md): the carry criterion for three-term progressions in a digit design, brute force as its control, and the census of the largest progression-free design per base to `40`.
- A digit set is a `u64` mask over the digits `0..b-1`; its signed set `S = F + F - 2F` is never built, only the six memberships `1, -1, b, -b, b + 1, -b - 1` are read off the `u128` mask of `F + F`.
- No dependency: one file, `std` only.

## THE METHOD

- `free(base, mask)`: `F` progression-free in the integers and none of the four closing patterns `b, -1`; `b, -b - 1, 1`; `-b, 1`; `-b, b + 1, -1` in `S`. `exits` reads whether `b` or `-b` lies in `S`, the two edges out of carry `0`.
- `families(base)`: every `F` with `0 in F` that is progression-free in the integers, by adding digits in increasing order and refusing a digit `d` that closes `2g - d in F`; `89442` sets at base `27`, `293820` at bases `2..27` together.
- `brute(base, mask, level)`: the sorted level, a bitset, and every pair `x < z` of equal parity tested for its midpoint; shares nothing with the criterion.
- `cyclic(n)`: the largest subset of `Z_n` holding `0` with no `f + g = 2h` mod `n` besides `f = g = h`, by the same increasing search with the congruence tested directly; the convention of A304884, where `f = g != h` counts at even `n`.

## THE CONTROLS

- `control`: at bases `3..12` every family member against brute force to level `5`; asserts `0` mismatches and prints the deepest first level, `3`.
- `verify`: at bases `2..27` every family member against brute force at level `3`, the decisive level; asserts `0` mismatches on `293820` sets.
- `census`: at bases `28..40` every maximiser against brute force at level `3`; every half-base `F` asserted free; the cyclic maximum asserted equal to the no-carry maximum at every base.
- `cargo test`: the asymmetry `-5 in S`, `5` not, at `F = {0, 1, 3}`; `(3, {0, 1})` free to level `8`; `(5, {0, 1, 3})` closing by a carry with the witness `0, 3, 6` at level `2`; `(4, {0, 2})` free with both exits; the base-`19` design free with no exit; the criterion against brute force at level `3` at bases `3..9`.

## RUN

- `cargo run --release --manifest-path research/progressions/lab/rs/carry-census/Cargo.toml` from the repo root, `23` seconds, most of it the families and the cyclic search at bases `33..40`. Prints only, writes nothing; every check asserts.
- `cargo test --release --manifest-path research/progressions/lab/rs/carry-census/Cargo.toml`, `9` tests, under a second after the build.

## READS

- The control table, the verify table, and the census table: per base the half-base count `r_3(ceil(b/2))`, the cyclic maximum, the no-carry maximum, the largest design with an edge out of `0`, the best, the number of maximisers, the least maximiser and `b^(log 2/log 3)`; then every maximiser at a base where the best passes the half-base count, with its exits; the bases where the best passes `b^(log 2/log 3)`; and the sequence of best fills at bases `2..40`.

## WITNESSES

- `0` mismatches in `control` and in `verify`; the deepest first level `3`.
- The half-base count is passed at bases `19`, `35`, `37`, `39` by `18`, `54`, `40`, `40` digit sets, every one with no exit.
- The dimension `log 2/log 3` is passed at base `37` only, by ten digits, the least set `{0, 1, 3, 7, 17, 24, 25, 28, 29, 35}`.
- Best fills at bases `2..40`: `1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 6, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 9, 8, 10, 8, 10, 9`, the cyclic maxima A304884 at `3..40`.
