# sumset-density

- Computes the sumset `S = A + B` of the base-3 design `A = {0, 1}` and the base-4 design `B = {0, 1}`, the object of Erdos problem 125, exactly as a bit array up to `3^K`, and reads its density `D(x) = card(S meet [1, x])/x`.
- The bit array holds one bit per integer in `[0, 3^K]`; the members of `A` are set directly, and each power `4^j <= 3^K` is folded in by one shift-or pass, `S |= S << 4^j`, so `floor(K log_4 3) + 1` passes over the array build the whole sumset, `14` at `K = 17` and `18` at `K = 22`. The other order, `B` direct and the powers of `3` shifted, is the control.
- The scan reads `card(S meet [1, x])` at `x = 3^k`, `4^m`, `floor(3^k/2)`, `floor(4^m/3)` and at the centres `d = (3^k - 1)/2 + (4^m - 1)/3`, and the maximum and minimum of `D(x)` over each window `[3^k, 3^(k+1))`, each window `[2^j, 2^(j+1))` and over all of `[1, 3^K]`, each with its location. Maxima and minima are compared exactly by cross multiplication; a density prints truncated to six places; an exponent `log(card)/log(x)` prints truncated.
- The energy verb computes the additive energy `E(k, m) = card{(a, b, a', b') in A_k^2 x B_m^2 : a + b = a' + b'}` exactly, as `sum over t of 2^(z_3(t) + z_4(t))` with `t` ranging over the `3^m` integers with base-4 digits in `{-1, 0, 1}`, `z_4(t)` the number of zero digits there and `z_3(t)` the number of zero digits in the balanced ternary expansion of `t` on `k` digits, or nothing when `abs(t) > (3^k - 1)/2`; a test checks it against the histogram of the representation function at three pairs. It prints, per pair `(k, m)` with `4^m` within a factor `3` of `3^k`, the fill `card(S meet [0, d])/(d + 1)`, the energy, the flat energy `4^(k+m)/(d + 1)`, their ratio `Q` rounded up, and the Cauchy-Schwarz bound `1/Q` truncated down.
- The ladder verb computes the same energy without the bit array: `t = u 4^12 + v` with the `3^12` low strings `v` held in a table and the high strings `u` shared over threads, and `z_3(t)` read eight balanced ternary digits at a time from a `3^8` table of zero counts and carries. A test checks it against the digit-string energy at every `k <= 11`, `m <= 14` and against the representation histogram at `(9, 7)`. It prints `Q` rounded up and `1/Q` truncated down for every pair `6 <= k <= K` with `1/3 < 4^m/3^k < 4`, which holds every pair within a factor `3` and the whole chain `1 <= 4^m/3^k < 4`, and least-squares fits of `Q` against `k` over the distinct pairs, the distinct chain pairs and the clean pairs, from `k = 6, 11, 16`: the power model `3^(eta k)` fitted in `log_3 Q` with its standard error, the same model fitted in `Q` by a golden-section search on `eta`, and the saturating model `a + b 3^(-s k)` with `s = log_3 2 - 1/2`, all scored by rms in units of `Q`. A pair with `2 4^m < 3^k + 5` or `3^(k+1) < 4^m + 5` is flagged `copy`: `A_k + B_m` is then two disjoint translates of `A_(k-1) + B_m` or of `A_k + B_(m-1)`, and `E(k, m)` is exactly twice that pair's energy, so it is printed but kept out of the extremes and the fits; the chain statistics are printed once more with the copies in. The extremes print the lower end truncated and the upper rounded up. A cache file, named as the second argument, holds `k m E` lines; the verb reads it, computes only the pairs it lacks and appends them.
- `constants` prints the constants of the average bound on the note: `alpha = log_3 2`, `gamma = 1/(2 - 2 alpha)`, the bound `6^(1/(2 gamma)) gamma/(gamma - 1)` on `E(abs(b - b')^(alpha - 1))`, `30^alpha` and their total, each rounded up.
- The phases verb reads, for each `k` on the chain `1 <= 4^m/3^k < 4`, the energy of the real atoms `a + sigma b` at `G` scalings `tau = sigma 4^m/3^k` spread evenly in `log tau` over `[1, 4)`, with a nearest-integer window and with the tent `(1 - abs(Delta))_+`; both equal `E(k, m)` at `sigma = 1`, which a test checks. It prints the rank of the lattice value among the `G` readings, the same after multiplying by the support length `1/2 + tau/3`, the rank and the ratio to the mean within `1/50` of `tau_k` in `log_4 tau` under both kernels, the integral over `[1, 4]` and the maximum with its `tau`.

## VERBS

- `density K` builds `S` to `3^K` and prints the readings above. Runtime `0.09` s at `K = 17`, `1.3` s at `K = 20`, `13` s at `K = 22` on a `4` GB bit array; `K = 23` wants `12` GB.
- `energy K` builds `S` to `3^K` and prints the energy table for every pair with `d <= 3^K`. Runtime `0.3` s at `K = 17`, `39` s at `K = 22`.
- `ladder K [cache]` prints the energy table and the fits for every pair with `6 <= k <= K`, `K <= 29`, since `k = 30` overflows the `u128` ratio. Runtime `0.4` s at `K = 22`; at `K = 29` the energies take `139` to `178` s on `8` threads, most of it the pair `(29, 23)`, and a run from a full cache takes under a second.
- `constants` runs instantly.
- `phases K G` prints the phase table for `8 <= k <= K`. Runtime `47` s at `K = 17`, `G = 4096` on `8` threads; `K = 18` adds about a minute.
- `control` compares the shift-or array against a double loop over `A x B` at `3^13`, the two shift orders against each other at `3^17`, the first `58` non-members against the terms of [A367090](https://oeis.org/A367090), and tests the reflection `x -> d - x` on `S meet [0, d]` at every centre below `3^17`. Runtime `0.08` s.

## RUN

```
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- density 17
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- density 22
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- energy 22
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- ladder 29 energies.tsv
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- constants
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- phases 17 4096
cargo run --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml -- control
```

- `cargo test --release --manifest-path research/cobham/lab/rs/sumset-density/Cargo.toml`, `11` tests, `1.8` s after the build.

## READS

- `card(S meet [1, 3^22]) = 26666749554`, `D(3^22) = 0.849772`.
- `D(3^k)` at `k = 4..22`: `0.975308, 0.835390, 0.858710, 0.887517, 0.908855, 0.864959, 0.778472, 0.837186, 0.858264, 0.874244, 0.814704, 0.763392, 0.831183, 0.858962, 0.881342, 0.792352, 0.767893, 0.831191, 0.849772`.
- `D(4^m)` at `m = 3..17`: `0.968750, 0.843750, 0.860351, 0.897460, 0.859313, 0.791305, 0.837238, 0.868845, 0.806823, 0.783585, 0.838184, 0.875988, 0.785523, 0.793552, 0.845272`.
- Maximum of `D` over `[3^k, 3^(k+1))` at `k = 5..21`: `0.913419, 0.903768, 0.912038, 0.931596, 0.913781, 0.875566, 0.875469, 0.881621, 0.908274, 0.885045, 0.865671, 0.882855, 0.886340, 0.910650, 0.874408, 0.865858, 0.872186`; minimum: `0.835390, 0.852729, 0.887517, 0.858945, 0.778468, 0.778472, 0.822506, 0.858264, 0.806430, 0.763391, 0.763392, 0.815887, 0.858962, 0.785230, 0.767893, 0.767875, 0.818358`.
- Minimum of `D` over `[1, 3^22]`: `0.763391` at `x = 3^15 - 1`; the maximum is `1` at `x <= 61`.
- The reflection `x -> d - x` fixes `S meet [0, d]` at every clean centre and at no mixed centre with `d >= 449`.
- `Q(k, m)` over the `27` pairs with `6 <= k <= 22` lies in `[1.467705, 2.060586]`, the maximum at `(16, 12)`; the Cauchy-Schwarz bound `1/Q` is at least `0.485298` on every one of them, against fills between `0.834213` and `0.928391`. Fitted as `3^(eta k)` on two endpoints, `Q` grows at `eta = 0.015852` from `k = 6` to `k = 22` and at `eta = 0.001987` from `k = 11` to `k = 22`, both rounded up.
- `Q` over the `43` pairs with `6 <= k <= 29` and `1/3 < 4^m/3^k < 4`: twelve are gap copies, `(6, 4), (7, 5), (11, 8), (12, 9), (16, 12), (21, 16), (26, 20)` and the chain levels `(9, 8), (14, 12), (19, 16), (24, 20), (28, 23)`; over the `31` others `Q` lies in `[1.638124, 2.004783]`, the lower end truncated, the maximum at `(26, 21)`; the whole chain, `24` levels with copies, lies in the same interval with mean `1.861840`. Least squares of `log_3 Q` on `k`: `eta = 0.003915`, standard error `0.001243`, over the `31`; `eta = -0.001109`, standard error `0.001254`, over the `25` with `k >= 11`, where the power model fitted in `Q` has `eta = -0.001093` and rms `0.065794` and the saturating model rms `0.065127` with `b = 0.282604`, all in units of `Q`.
- `constants`: `E(abs(b - b')^(alpha - 1)) <= 7.398167`, `30^alpha <= 8.549875`, total `189.335292`, each rounded up.
- Phases at `k = 8..17`, `G = 4096`: the local rank of `h_k(tau_k)` averages `0.8110` under both kernels, the local ratio lies in `[1.0090, 1.0525]`, the integral over `[1, 4]` rises from `3.937434` at `k = 8` to `4.572550` at `k = 17`, and the maximum sits next to `4^4/3^5` or `4^5/3^6`.

## WITNESSES

- `object-s-base-3-plus-base-4.md` section "Object S: the base-3 design plus the base-4 design", every number there.
