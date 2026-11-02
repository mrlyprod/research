# Mess check

- The independent half of the mess study: it shares no code with the Rust census and reads only its two JSON files.
- It rebuilds the 22 Sources from definitions, the primes by sieve, the Fibonacci and powers of two directly, the odd squares, and the eight design families from the fill law `fill(F, n) = sum over c in F of E^(2-w(c)) O^(w(c))` with `E = ceil(n/2)`, `O = floor(n/2)` at the odd sides, and compares them with the library's sets.
- Route one, a brute-force truth table at `k <= 20`: the table `g(popcount(x))`, the fast Walsh-Hadamard transform for the spectrum, the Moebius transform for the degree and the transform of `W^2` for the autocorrelation; seven quantities per class and `k` compared with the census.
- Route two, exact integers at `k <= 64`: the Krawtchouk matrix by polynomial products, the spectrum, Parseval's identity, the autocorrelation by the Wiener-Khinchin sum, and the ANF by binomial parity `C(w, j) mod 2`, which tests Lucas's theorem itself instead of the digit domination the Rust census uses.
- Route three, the frontier at `k <= 20`: all `2^(k+1)` symmetric functions as a bit matrix, the spectrum as a matrix product with the Krawtchouk matrix, the ANF as a product with the binomial-parity matrix, and per degree the count, the best `max|W|`, how many reach it, and the balanced counterparts, compared with the Rust frontier.
- Route four, the dropped modes: the Sources with 0 dropped and with 0 and 1 dropped, rebuilt from definitions, then the exact route at every `k <= 64` for the sets, the affine classes, the tight cells, the full-degree cells and the best class over the ceiling, with the undropped mode beside them as a control.
- It writes the ANF vector `c(w)` to `w = 255` and the degree sequence to `k = 64` of every class to its data folder, for a lookup against the OEIS dump.

## RUN

- `cargo run --release --manifest-path research/random/lab/rs/mess/Cargo.toml` first, then `uv run python research/random/lab/py/mess/mess.py`
- About seven seconds; prints and writes one text file; holds one `2^20` table at a time.

## WITNESSES

- mess.md the spectrum: 1960 brute-force cells over 14 classes at `k <= 20` with 0 disagreements, 6272 exact cells at `k <= 64` with 0 disagreements and no Parseval or Lucas failure.
- mess.md the ceiling and the frontier: 1380 frontier cells at `k <= 20` agree with the Rust enumeration.
- mess.md the object: the dropped modes recomputed from definitions, 13 sets, 576 of 684 tight cells, 362 of 744 full-degree cells and 0.7862 for drop 0, 13 sets, 578 of 741, 368 of 806 and 0.7803 for drop 0 and 1.
- mess.md the degree by Lucas: the first 105 terms of the binary `c(w)` equal A010060; the primes' ones counts equal A121497 at `k = 1..33`.
