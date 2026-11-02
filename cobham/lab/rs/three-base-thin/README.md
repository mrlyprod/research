# three-base-thin

- Computes the three-base thin set `E = {n : every base-3 digit <= 1, every base-5 digit <= 2, every base-7 digit <= 2}` exactly, and the wider set with the base-7 bound at `3`, which is `{k : binomial(2k, k) is prime to 105}`, [A030979](https://oeis.org/A030979).
- Exact integer arithmetic throughout, `num-bigint` above `128` bits; no sampling, no fit.
- The walk enumerates the base-3 side, whose members are the subset sums of distinct powers of `3`, from the top power down. Once the powers `3^k` and above are chosen the remaining addition is at most `(3^k - 1)/2`, so with `j` least such that `5^j > (3^k - 1)/2` the high part `floor(n / 5^j)` of every `n` still reachable is one of two consecutive integers, and the branch dies when neither of them has all its base-5 digits inside the bound; base `7` cuts the same way. Membership in the base-3 set holds by construction and is never tested.
- Also prints the three dimensions of `E`, `log_3 2`, `log_5 3` and `log_7 3`, the base-7 dimension `log_7 4` of the wider set, the three pair budgets `dim_i + dim_j - 1`, the two triple budgets `sum dim - 2`, and the criterion `A/(p-1) + B/(q-1)` of Erdos, Graham, Ruzsa and Straus for each pair. A dimension prints to the nearest six places, a budget rounds up, an effective exponent truncates down.

## VERBS

- `budget` prints the four dimensions, the three pair budgets, the triple budget for the thin set and for the wider one, and the five values of the infinitude criterion. Runtime under `0.01` s.
- `control` rebuilds the members below `10^8` twice, once by the pruned walk and once by scanning every integer against all three digit rules, at both base-7 bounds, and rebuilds the `23` terms [A030979](https://oeis.org/A030979) publishes. Runtime `1.9` s.
- `seven M` prints the members of `E` below `7^M` with the node count. Runtime under `0.01` s at `M = 17`.
- `reach level` prints the members of `E` below `3^level` with the node count. Runtime `3.02` s at level 20000, `13.60` s at level 40000 and `61.48` s at level 80000, the last on about `3` GB, which is the wall: the stored powers cost `Theta(level^2)` bits and the node count grows near `21 level`.
- `wide level` and `ten X` do the same at base-7 bound `3`, below `3^level` and below `10^X`, printing the count, the largest member and the effective exponent `log(count)/log(height)`. Runtime `0.09` s at `X = 70` and `39.32` s at `X = 140`.

## RUN

```
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- budget
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- control
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- seven 17
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- reach 80000
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- ten 70
cargo run --release --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml -- ten 140
```

- `cargo test --manifest-path research/cobham/lab/rs/three-base-thin/Cargo.toml`, `7` tests, `1.18` s after the build.

## READS

- `E` below `7^17 = 232630513987207`: `0, 1, 3186, 3187, 20007`, `333` nodes.
- `E` below `3^80000`, a height of `38170` decimal digits: the same five, `1710789` nodes.
- The wider set below `10^70`: `1374` members, the length of the table [A030979](https://oeis.org/A030979) calls complete to `10^70`. Below `10^140`: `216020` members.
- Dimensions `0.630930`, `0.682606` and `0.564575` for `E`, and `0.712414` for the base-7 side of the wider set; pair budgets `0.313536`, `0.195505`, `0.247182`; triple budgets `-0.121889` thin and `0.025951` wide; criterion `1.000000` at `(3, 5)` and `0.833333` at `(3, 7)` and at `(5, 7)`.
- Effective exponents of the wider set `0.044828` at `10^70` and `0.038103` at `10^140`, both above `0.025951` and falling.

## WITNESSES

- `gasket-three-base-and-two-frequencies.md:24` the three dimensions of `E`
- `gasket-three-base-and-two-frequencies.md:25` and `gasket-three-base-and-two-frequencies.md:26` the pair budgets and the triple budget
- `gasket-three-base-and-two-frequencies.md:28` and `gasket-three-base-and-two-frequencies.md:29` the infinitude criterion at each pair
- `gasket-three-base-and-two-frequencies.md:31` the members below `7^17`
- `gasket-three-base-and-two-frequencies.md:32` the members below `3^80000`
- `gasket-three-base-and-two-frequencies.md:33` the scan control below `10^8`
- `gasket-three-base-and-two-frequencies.md:34` the fourth dimension `log_7 4 = 0.712414` and the wider budget `0.025951`
- `gasket-three-base-and-two-frequencies.md:36` the `23` published terms rebuilt, the `1374` below `10^70`, the `216020` below `10^140` and the two effective exponents
