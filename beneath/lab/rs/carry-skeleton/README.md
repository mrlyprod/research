# carry-skeleton

- One step of `m n + 1` in base 2 split into its `GF(2)` skeleton and its carry, and the carry read as an object of the memory dial of [beneath](../../../notes/memory-dial.md).
- Digits run least significant first, so `2n` is a shift towards the higher digits and the ripple carry runs in the same direction as the skeleton's dependency.
- The skeleton of `3n + 1` is `n xor 2n xor 1`; the carry word is `defect(n) = (3n + 1) xor (n xor 2n xor 1)` and the local carry count is `d_loc(n) = popcount(defect(n))`.
- Every rule reading here is `mrly_math::memory::Rule`: `Rule::new(1, k, code)` at window width `k`, with the window `(c_1, ..., c_k)` read as `w = sum_j c_j 2^(k - j)`, first digit most significant, and `Rule::allowed(w)` true when bit `w` of the code is set.
- `rho` and `kappa` are `mrly_math::memory::perron` and `mrly_math::memory::kappa`, which split the digraph into strongly connected components and make each component's Perron root exact against its integer characteristic polynomial.

## THE METHOD

- `carries(n)` runs the four-state Mealy transducer directly: state `(n_(i-1), q_i)`, `q_0 = 1`, `q_(i+1) = MAJ(n_i, n_(i-1), q_i)`, and its output word is asserted equal to `defect(n)` digit for digit.
- The substitution `M = 2n + 1` is checked as the exact identity `M xor 2M = 2 (n xor 2n xor 1) + 1`, and `M xor 2M` is checked against rule 60 rebuilt digit by digit as `M_i xor M_(i-1)`.
- The zero-carry code of an odd multiplier `m` is recomputed from `supp m`: the difference set `{|j - j'| : j, j' in supp m}` is built, the window width is `deg m + 1`, and a window is allowed when no two of its `1` digits sit at a distance in that set. No code is copied.
- A word shorter than the width holds no window and `Rule::accepts` takes it, so every integer is padded with leading zeros past the width before the set equality is asserted; leading zeros never close a forbidden pair.
- The carry density is read two ways: the exact integer balance of the stationary vector of the four-state chain, and the mean of `d_loc` over every `n < 2^level`.
- The depth is the longest run of consecutive carry-on positions, sampled on uniform digit strings with a `splitmix64` stream seeded in the source, so every depth row reruns identically.

## THE CONTROLS

- `carries` against `defect`, the `M` substitution and the rule-60 rebuild: `0` mismatches each over every `n < 2^18`.
- The zero-carry set against the rule: `0` mismatches over every `n < 2^16` at `m = 3, 5, 7, 9, 11, 15`.
- `kappa` on codes `7` and `23` asserted against `0.098239` and `0.115204`, the two values [beneath](../../../notes/memory-dial.md) already prints, reached here from the multiplier and not from the rule.
- The carry-free map asserted non-growing and asserted to reach `1` from every `n < 2^20`.
- The six refutation witnesses pinned as explicit integer pairs, and the worst-case depth `level + 1` pinned at `level 1, 2, 4, 8, 16, 32, 40`.

## RUN

- `cargo run --release --manifest-path research/beneath/lab/rs/carry-skeleton/Cargo.toml` from the repo root. `7.53s` on one thread, `7.43s` of it the depth sampling. Prints only, writes nothing.
- `cargo test --manifest-path research/beneath/lab/rs/carry-skeleton/Cargo.toml`, `12` tests, under a second after the build.

## READS

- The identity: `0` mismatches on all three checks below `2^18`, and the radius witness holds at every length `2..40`.
- The zero-carry rules, one row per odd multiplier, every code recomputed from the difference set.

| `m` | `supp m` | difference set | `k` | code | `card W` | `rho` | `kappa` |
|---|---|---|---|---|---|---|---|
| 3 | `0, 1` | `1` | 2 | `7` | 3 | `1.618034` | `0.098239` |
| 5 | `0, 2` | `2` | 3 | `95` | 6 | `1.618034` | `0.167412` |
| 7 | `0, 1, 2` | `1, 2` | 3 | `23` | 4 | `1.465571` | `0.115204` |
| 9 | `0, 3` | `3` | 4 | `22015` | 12 | `1.618034` | `0.201999` |
| 11 | `0, 1, 3` | `1, 2, 3` | 4 | `279` | 5 | `1.380278` | `0.115524` |
| 15 | `0, 1, 2, 3` | `1, 2, 3` | 4 | `279` | 5 | `1.380278` | `0.115524` |

- `11` and `15` share a difference set, so they share a rule: the zero-carry set depends on the difference set of `supp q` and not on `m`.
- The density: the stationary vector is `(2, 1, 1, 2)/6` on the states `(0,0), (0,1), (1,0), (1,1)`, all four balance residuals exactly `0`, carry-on mass exactly `3/6`. The mean of `d_loc` over every `n < 2^level` reads `4.332031, 5.333008, 6.333252, 7.333313, 8.333328, 9.333332, 10.333333, 11.333333` at `level 8, 10, 12, 14, 16, 18, 20, 22`, which is `level/2 + 1/3 - (-1)^level/(3 * 2^level)` exactly: the integer form `3 sum = 3 level 2^(level-1) + 2^level - (-1)^level` has residual `0` at every `level 8..22`, and `mean/level` falls `0.541504, 0.533301, 0.527771, 0.523808, 0.520833, 0.518518, 0.516667, 0.515152` towards `1/2`.
- The carry-free map: `0` digit-count increases and `0` values failing to reach `1` over every `n < 2^20`, with `T_free(1) = 1`.
- The six refutations, each the least clashing pair below `2^16` and the number of unordered pairs of integers below `2^16` that share the statistic and disagree on `d_loc`. A statistic is refuted as soon as one such pair exists; the count is how many there are.

| statistic | least clashing pair | `d_loc` | disagreeing pairs below `2^16` |
|---|---|---|---|
| popcount | `1, 2` | `2, 0` | `256217518` |
| longest run of `1` digits | `1, 2` | `2, 0` | `417177932` |
| `v_2` | `1, 3` | `2, 3` | `659301399` |
| digit count | `2, 3` | `0, 3` | `662864636` |
| `(popcount, v_2)` | `3, 5` | `3, 4` | `88157572` |
| all four at once | `19, 25` | `3, 4` | `9331881` |

- The last row subsumes every pair: two integers agreeing on all four statistics already disagree on `d_loc`, so no pair of them is a summary either.
- The depth: the carry-on block of the transfer matrix is `[[0, 1], [1, 1]]`, characteristic polynomial `x^2 - x - 1`, the same as `transfer(Rule::new(1, 2, 7)) = [[1, 1], [1, 0]]`, so the survival rate per digit is `phi/2` and the prediction is `log_(2/phi) level` at `2/phi = 1.236067977`.

| `level` | samples | mean depth | standard error | `log_(2/phi) level` | offset | increment per quadrupling |
|---|---|---|---|---|---|---|
| 16 | 200000 | `6.338` | `0.008` | `13.082` | `-6.744` | - |
| 64 | 200000 | `12.116` | `0.011` | `19.623` | `-7.507` | `5.778` |
| 256 | 200000 | `18.462` | `0.013` | `26.164` | `-7.702` | `6.346` |
| 1024 | 200000 | `24.967` | `0.013` | `32.706` | `-7.739` | `6.504` |
| 4096 | 200000 | `31.481` | `0.014` | `39.247` | `-7.766` | `6.514` |
| 16384 | 200000 | `38.029` | `0.014` | `45.788` | `-7.759` | `6.548` |

- The predicted increment per quadrupling is `6.541119`; the read increments rise towards it from below, the last reading `6.548` at a standard error of `0.014` per mean, and the offset settles near `-7.76`, so the base is supported and not established. The mean is a sample mean over `200000` strings and no exponent is fitted.

## WITNESSES

- [beneath](../../../notes/memory-dial.md), The carry of a Collatz step - the identity, the unbounded radius, the zero-carry table, the two couplings, the density, the carry-free cycle theorem, the refutation line and the depth conjecture.
