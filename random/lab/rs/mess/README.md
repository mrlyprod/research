# Mess

- The mess as a Boolean function: a written cell XORs `g(count)`, `count` the live cells among `k` neighbours and `g` the indicator of a Source of `mrly_life::source`. The 22 fixed Sources fall into 14 distinct sets, because the fills of a design are the voids of its complement; `k = 1..64` covers the design fills `2..9, 16, 25, 36, 49, 64`.
- Per Source and `k`: the ANF vector by binary domination (Lucas), the degree, the ones count in exact `u128`, the Walsh spectrum by the Krawtchouk formula in `i128`, the nonlinearity `2^(k-1) - max|W|/2`, the autocorrelation by a direct double sum, the differential uniformity `(2^k + max|autocorrelation|)/2` and the correlation-immunity order.
- Asserted on the way: Parseval's identity to `k = 40`, a brute-force truth table with the fast Walsh-Hadamard and Moebius transforms at `k <= 14`, the bounds `NL <= min(ones, 2^k - ones)` and `NL <= ceiling`, and the laws odds = parity, evens = its complement, binary = Thue-Morse, primes at even `w`.
- The ceiling is `2^(k-1) - 2^(k/2-1)` for even `k` and `2^(k-1) - 2^((k-1)/2)` for odd `k`.
- The frontier walks all `2^(k+1)` symmetric functions for `k = 1..26` in Gray-code order, updating the spectrum and the ANF mask by one Krawtchouk row per step, and keeps per degree the count, the best `max|W|` and how many functions reach it, and the same for the balanced ones. It asserts that the ceiling is met at every `k`.
- The typical draws: 2000 symmetric functions per `k` in `8, 16, 24, 32, 48, 64`, once with each `g(w)` a fair coin from a printed seed and once from `Source::Random(0..2000)`; median, minimum and maximum `NL` over the ceiling, the share of full degree and of balanced.
- The dropped modes: the 22 Sources rebuilt through `counts()` with 0 and 1 kept, with 0 dropped and with both dropped; per mode the sets, the affine classes, the tight cells `NL = min(ones, 2^k - ones)`, the full-degree cells and the best class over the ceiling, with the brute-force table rerun; and which Sources contain 0 and which contain 1.
- Rule 30's local map `l XOR (c OR r)` and the 3-variable majority through the brute-force table.
- Writes `census.json` and `frontier.json` into the study's data folder, which the Python study reads.

## RUN

- `cargo run --release --manifest-path research/random/lab/rs/mess/Cargo.toml`
- About eight seconds on one thread; prints 598 lines and writes two JSON files of 96 KB and 13 KB.

## WITNESSES

- mess.md the object: 22 Sources in 14 classes, the eight twin pairs, the design fills, 10 Sources with 0 and 14 with 1, and the dropped modes: drop 0 gives 13 sets, 576 of 684 tight cells, 362 of 744 full-degree cells, best 0.7862; drop 0 and 1 gives 13 sets, 578 of 741, 368 of 806, best 0.7803.
- mess.md the degree by Lucas: the odds `c(w) = [w = 1]`, the evens `c(w) = [w <= 1]`, the binary `c(w) = popcount(w) mod 2` with degree `k` at 33 of 64 values, the primes `c(w) = [w = 2 mod 4]` at even `w`, the primes degree sequence to `k = 64`, full degree at 334 of 744 class cells.
- mess.md the spectrum: 980 brute-force cells at `k <= 14` agree.
- mess.md the ceiling and the frontier: the ceiling met by exactly 4 functions and `ceiling - 1` at full degree by exactly 8, for `k` up to 26; the table at `k = 4, 8, 12, 16, 20, 24, 26`; the best cubic and the best at degree `k - 1`; the balanced counts `2, 2, 4, 2, 8, 2, 16, 6, 32, 2, 64, 2, 144, 14, 256, 2, 512, 2, 1024, 6, 2048, 2, 4096, 50, 8192, 6`.
- mess.md the Sources against the frontier: the primes equal the quadratic pattern for `k < 5`, the best class over the ceiling `0.7867` at `k = 12`, the leaders by `k`, the table at `k = 8` and `k = 64`, `NL = min(ones, 2^k - ones)` on 632 of 684 cells, the share of ones, the lowest differential `0.7206`, no balanced non-affine cell past `k = 3`, the typical medians, Rule 30 and the majority.
