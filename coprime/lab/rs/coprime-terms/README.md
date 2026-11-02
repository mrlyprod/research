# coprime-terms

- Exact `A(level)`, the points of a base-3 mask design at level `level` with coprime coordinates, without walking the `fill^level` points: the Menger sponge to `level = 18`, the Sierpinski carpet and the Vicsek plus to `level = 20`.
- Mask form: with `mask(x)` the set of base-3 positions where `x` has digit 1 (a digit away from 1 for the plus), a point lies in the level set iff its coordinate masks are pairwise disjoint.
- Mobius over the common divisor with the base prime peeled: `A(level) = W(level) - W(level-1)` where `W(level) = Sum_{m < 3^level, gcd(m,3) = 1} mu(m) (N_level(m) - 1)`, since `N_level(3m) = N_{level-1}(m)` whenever the zero digit vector is filled; the plus has `A(level) = W(level)` and no `-1`.
- `N_level(m)` is a pairwise-disjoint count over the masks of the `Y = ceil(3^level / m)` multiples of `m`, on at most `2^level` masks whatever the box; the sponge splits its moduli by `Y` over five counters, each cheapest in its band by measurement.
- `tail`, `Y <= 3`: closed forms `N - 1 = 3 + 4 [mask(m) = 0]` for `Y = 2` and `6 + 7 [mask(m) = 0] + 7 [mask(2m) = 0] + 6 [mask(m), mask(2m) disjoint]` for `Y = 3`, so the top two bands of `W(level)` are Mertens sums over automatic sets.
- `bitset`, `Y <= (512 level 2^level)^(1/3)`: one `Y`-bit row per multiple marks its disjoint partners and `N = 6 T + 3 z (Y - 1) + z`, `T` the popcount sum of row pairs above the diagonal and `z` the zero-mask multiples; about `Y^3 / 1000` cycles against `Y^3 / 100` for the nested loops it replaces.
- `rows`, `U <= 6 sqrt(level 2^level)` distinct masks: one `u16` zeta transform `g` of the mask histogram `c`, then `2 Sum_{a < b disjoint} c_a c_b g(complement(a | b))` plus the zero-mask diagonal, the disjoint pairs read off bitset rows built four masks per NEON lane inside buckets keyed by the top six bits.
- `cube`, the rest: `u32` ranked zeta transforms truncated at rank `K ~ level/2` with the rare heavier masks folded back by subset enumeration, then `N = Sum_T <c_hat(complement T), Delta[sq(T)]>` per set `T`, `sq` the square of the rank polynomial and `Delta` its binomial difference table, in wrapping `u64` where `Y^3 < 2^64` and `u128` otherwise.
- `residue`: the automaton on `(Z/m)^3` for `m <= 16`; the earlier nested-loop, zeta and ranked-convolution counters stay as references inside the gate.
- Cost: `3^level` moduli a level; the `bitset`/`rows` crossover costs `3^level (level 2^level)^(2/3) / 8` and the `rows`/`cube` balance `3^level level 2^(level/2)` cycles, so the level ratio in `terms/menger.txt` is `4.32` from `level = 17` (`28.317 s`) to `18` (`122.294 s`) on eight threads, `3.7x` and `3.8x` the `104 s` and `463 s` of visibility-parity-and-slices.md:116, and tends to `4.76` only past `level = 22`.
- `check` gates the engine: direct enumeration to `level = 6` for the sponge and `8` for the plane designs, every counter agreeing with every other on every modulus to level 8, each counter pinned alone agreeing with auto to level 8 (`bitset`), 10 or 11, 45 probed moduli at the sample level across counters, two pinned level-19 probes for the `u64` cube gate (`m = 440`) and the `u32` rows branch (`m = 3^8 + 1`), the `3m` peel identity, and the fill `N_level(1) = fill^level`.
- `terms/` holds the three ladders with seconds per level and their b-files.

## RUN

- `cargo run --release --manifest-path research/coprime/lab/rs/coprime-terms/Cargo.toml -- terms menger 18 8`
- Arguments: design (`menger`, `carpet`, `vicsek`), top level, threads.
- `-- check menger 6 14` runs the gate battery; `-- profile menger 16 8` prints seconds per `log2 Y` band and counter for one level; `cargo test --release --manifest-path research/coprime/lab/rs/coprime-terms/Cargo.toml` pins the stored terms.
- The sponge ladder to `level = 14` takes 1 s, to `level = 18` about 2.5 min (`terms/menger.txt`).

## WITNESSES

- visibility-parity-and-slices.md:116: the sponge ladder `A(10) .. A(18)`, ending `215134797774716879278017`, with the 104 s and 463 s timings of the previous engine.
- research/claims/:161: the same ladder and the gate battery.
- A399364: DATA `level = 0..17` and the b-file `terms/menger_bfile.txt` to `level = 18`.
