# primes

- The generator of the chapter Along the primes: one binary, five verbs, no dependencies, up to 8 threads.
- A design is `f_{b,D}(n)`, the parity of the number of base-`b` digits of `n` in `D`, `2 <= b <= 6`, `D` a nonempty subset of `{1..b-1}`: 57 designs, named `b:D`, so `2:1` is Thue-Morse. Its sign is `x(n) = prod_i v(n_i)`, `v(d) = (-1)^[d in D]`. Beside them run Rudin-Shapiro and a coin, the top bit of splitmix64 at `n`; one `u64` word carries all 59, read through the digit tables of the squares study.
- The bias along the primes is `Psi(x) = sum_(n<x) x(n) Lambda(n)`, in units of `sigma(x) = sqrt(sum_(n<x) Lambda(n)^2)`, and unweighted `sum_(p<x) x(p)`, in units of `sqrt(pi(x))`.
- The major-arc sum is `M_Q(x) = sum_(q<=Q) T_q(x)`, `q` squarefree, `T_q = (mu(q)/phi(q)) sum_((a,q)=1) S(a/q; x)`, `S(alpha; x) = sum_(n<x) x(n) e(n alpha)` summed exactly by a digit walk: the product of the tile polynomials at `b^j alpha` gives `S(alpha; b^j)`, and the digits of `x` from the top give the rest. Rudin-Shapiro walks its two-state automaton the same way. `Q = 1024` unless a verb says otherwise.
- The unweighted model is the partial sum `M(x)/log x + int M(t)/(t log^2 t) dt`, product-integrated between the quarter-octave checkpoints with `M(t)/t` linear on each step.

## VERBS

- `census K`: a segmented sieve of Eratosthenes over the odd `n < 2^K`, 8 threads, spans of `2^20`. At every checkpoint `round(2^(c/16))` and at `1339484198` (the first `2^26` primes) it records `pi`, `theta`, `sum Lambda^2`, and for each stream `sum x(p)` and `Psi`; beside them the sieve model, `n` prime to 210 binned by its least prime factor up to 1021, with `sum x(n)` and `sum x(n)/log n` per bin, so that `M_z = (P(z)/phi(P(z))) sum_(n<x, (n,P(z))=1) x(n)` reads off at `z = 7, 11, 13, 17, 31, 47, 73, 127, 257, 1021`. It writes the table to the data folder; the other verbs read it and run the census when it is missing. `census 32` takes 4 to 8 s.
- `thue K`: Thue-Morse at every `x = 2^k`, `k = 12..K`: `Psi`, the terms `T_3`, `T_6`, `T_5`, `T_10` and the rest of `M_Q`, the residual against `T_3` alone, the `q | 6` part, the `q | 30` part and `M_Q`, the normalised coefficient `G = -Psi/x^(log_4 3)` against the model's, the unweighted residual, the closed forms `-3^floor(k/2)` and `-5^floor(k/4) (1, 1/2, 0, 5/2)[k mod 4]`, the largest `T_15 + T_30` and the largest sum of the terms `q | 210` with `7 | q`, the largest residual where `Q <= sqrt(x)` and where not, the range of `A(4x) - 3 A(x)` for the `q | 6` part `A` over `x <= 2^22`, the growth of the best orbits on `q = 255` and `65535` with `T_65535` at `2^K`, the range of the coefficient of the `q | 6` part over one period, its minimum refined on 4096 points, and the squares chapter's `z_B` along the first `2^26` primes against the model. About 9 s.
- `law K`: every stream at `x = 2^K` and at the 49 quarter-octave checkpoints from `2^20`: `z`, the model, the residual `r`, its largest absolute value and rms, the largest `T_q`, the `q | rad(b)` part, the unweighted reading; the counts of designs under each cut; the iid-sign sd of `r`, `sqrt(1 - theta H/sum Lambda^2)` with `H = sum_(q<=Q) mu^2(q)/phi(q)`, at `2^20`, at `2^K` and pooled over the checkpoints, against the pooled rms and the rms at `2^K`; the sign-twin pairs along the odd primes and the independent streams; the coin; and the sieve model against the expansion, with its density at 1021. About 7 s.
- `ladder K Q`: `M_Q` for `Q` from 6 to the given `Q` at `2^(K-4)`, `2^(K-2)` and `2^K`, the rms and largest residual over the 54 distinct non-folded designs and Rudin-Shapiro, the iid-sign floor for each `Q`, the ratio of the rms to it and the excess variance, and the largest terms of eight designs with the growth exponent `g` of the best orbit through `a/q`, the mean of `log_b abs(P(e(c/q)))` over the cycle of `c = a b^j mod q`. `ladder 32 4096` takes about 5 s.
- `exact`: at `x = 6^k`, `k <= 12`, for the 31 designs of base 6, the `q | 6` part against `3 (v(1) + v(5)) P(1)^(k-1)` and the `q in (5, 10, 15, 30)` part against `-(3/4) sum_(a=1..4) (v(1) zeta^a + v(5)) P(zeta^a)^(k-1)`, `P(z) = sum_(d<6) v(d) z^d`, `zeta = e(1/5)`; for `6:23` and `6:14` the latter is `-(3/4) L_(3k-1)` and `-(15/4) F_(3k-2)`; and Thue-Morse's `T_7` at `8^m`, `m <= 6`, turned by the phase of `-i sqrt 7`.

## RUN

```
cargo run --release --manifest-path research/random/lab/rs/primes/Cargo.toml -- census 32
cargo run --release --manifest-path research/random/lab/rs/primes/Cargo.toml -- thue 32
cargo run --release --manifest-path research/random/lab/rs/primes/Cargo.toml -- law 32
cargo run --release --manifest-path research/random/lab/rs/primes/Cargo.toml -- ladder 32 4096
cargo run --release --manifest-path research/random/lab/rs/primes/Cargo.toml -- exact
cargo test --release --manifest-path research/random/lab/rs/primes/Cargo.toml
```
