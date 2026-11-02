# design-recursion

- The design meter `M_F(x) = sum of mu(n) over n in S_F, n <= x` at base 3 `{0,1}` to `3^27`, beside `A(x)`, the same sum over the elements coprime to `3`, and the recursion of [the Franel note](../../../notes/identity-and-recursion.md)'s section The recursion in the base, `M_F(x) = A(x) - A(x/3)`, asserted exactly at every `3^n`.
- Verb `level`: the mean square over a level of the Franel note's section The mean square over a level, to `B_27`.
- Verb `fake`: the multiplicative class `mu^2 psi` of the same section, to `3^24`: `mu`, a greedy member and one random member.
- Verb `pretend`: the pretentious members of the Franel note's section The pretentious members, to `3^24`: the modified principal member and eight real-character members, with the plain and the squarefree-weighted character sums over every level.
- The default run prints, at every `n`, `M_F(3^n)`, `A(3^n)`, the running maximum `max abs M_F(y)` over `y <= 3^n`, the reading `log max / log 3^n` and the same over `log 2^n`, the least-squares slope of `log max` on `log x` over the five rungs `n-4..n` with the bracket `[min, max]` of the reading over them, and `abs A(3^n)/(n + 1)`, the floor the recursion puts under `max_(j <= n) abs M_F(3^j)`.

## METHOD

- The `k`-th element of `S_F` is the base 3 reading of the binary digits of `k`, so the walk over `k = 1..2^n` is ascending and the unit digit is the parity of `k`: `A` sums the odd `k`.
- `mu` by trial division over the primes up to the least `b` with `b^3 > 3^level` (`6562` at level `24`, `19684` at `27`): the cofactor has no prime factor below `b` and is at most `3^level < b^3`, so it is `1`, a prime, a prime square or a product of two distinct primes; below `b^2` it is prime, above it an integer square root tells `p^2` and Miller-Rabin on the twelve witnesses `2..37` tells `p` from `pq`. No Pollard rho is needed.
- Four threads over chunks of `2^16` consecutive `k`; a chunk returns its sum, its sum over odd `k` and the minimum and maximum of its local prefix, and the sequential merge recovers the exact running maximum since `abs(offset + prefix)` over a chunk peaks at an endpoint of its prefix range.
- `level` stores `mu` on the `2^level` strings of `B_level` (index `k`, `mu(0) = 0`) on all the machine's threads, then walks it once for `Q_n`, `W_n`, the walk of `U`, and the law of cosines `sum P^2 = sum O^2 + sum O(k/2)^2 - 2 sum O(k) O(k/2)` with its cosine `rho_n`, asserting `P(k) = O(k) - O(floor(k/2))` at every `k`; at every `n` it folds the two ladders (residues mod `3^j` on the low index bits, translates on the top bits) and runs a fast Walsh transform in `i32`, asserting Plancherel and both ladders exactly, and prints the rank of the trivial twist and the largest twist.
- `fake` factors every element coprime to `3` by trial division to `3^(level/2)` on all threads; `psi` is greedy along the ascending coprime squarefree elements (`+1` at every new prime but the largest, the largest stepping `A` toward `0`; an element with no new prime is forced) or the sign bit of `splitmix64(p)`; the meter of each member is walked with `f(3m) = -f(m)`.
- `pretend` reuses the factoring of `fake`; a member's value on a coprime squarefree element is the product of `psi` over its primes, `psi(p) = chi(p)` off the modulus `q` and `+1` at the primes dividing `q` (`-1` at `2` for `chi4-`), the meter is walked with `f(3m) = -f(m)`, and the plain sums `sum_(m in B_n) chi(m)` and their squarefree-weighted versions are accumulated over the ascending walk; it asserts the principal member against the per-level squarefree count, the `chi_4` closed form `2^((n-1)/2)` at odd `n` and `0` at even `n`, and the bound `(q - 1) 2^(o - 1)`, `o` the order of `3` mod `q`, at `q = 5, 7, 11`.
- Self-checks at every run: `mu` against a linear Mobius sieve on all `2^15 - 1` elements below `3^15`, and the recursion `M_F(3^n) = A(3^n) - A(3^(n-1))` with `A(3^n) = sum_(j <= n) M_F(3^j)` at every `n`, both asserted.

## RUN

- `cargo run --release --manifest-path research/franel/lab/rs/design-recursion/Cargo.toml` from the repository root walks to `3^24` in `5.5 s` on four threads; `-- 27` walks to `3^27` in `73 s`. Prints only, writes nothing.
- `-- level 27` runs the level verb to `B_27` in `51.5 s`, `-- level 24` in `4.0 s`; `-- fake 24` runs the fake verb in `5.5 s`; `-- pretend 24` runs the pretend verb in `6.4 s`.
- `cargo test --release --manifest-path research/franel/lab/rs/design-recursion/Cargo.toml`: thirteen tests, `mu` against the sieve to `200000`, the two-large-prime cofactors, four strong pseudoprimes rejected, the integer square root, the trial bound, the ascending enumeration, the parallel walk against a sequential one at `3^18`, the recursion to `3^16`, `M_F(3^n)` with its running maximum against the `mobius` census ([mobius-designs](../../../../mobius/lab/rs/mobius-designs/)) at its rungs `4` to `16`, the level walk against the default walk to `3^16`, the fast Walsh transform against the literal twisted sums at `B_10` with every ladder, the fake verb's `mu` column against the level walk at `3^14`, and the plain character sums at `B_12`: the `chi_4` closed form, the bounds at `q = 5, 7, 11` and the squarefree-weighted sum mod `5` against the literal sum.

## WITNESSES

- The Franel note, The recursion in the base: every number of its Verified paragraph, `M_F(3^n)`, `A(3^n)` and the running maximum at `n = 12..27`, the readings `0.311823, 0.307223, 0.301266` at `n = 12, 24, 27`, the bracket `[0.2995, 0.3148]` at `n = 27`.
- The `mobius` census agrees at all thirteen of its rungs `3^4..3^24`, `M_(0,1)(3^24) = -1886` with running maximum `3296` among them; PARI (`moebius` on the literal elements, and `A` as the literal shifted sum over `3w + 1`, `w` a padded string) agrees at every `n <= 14`.
- The Franel note, The mean square over a level: every number of its Verified paragraph, from `level 27` and `fake 24`.
- The Franel note, The pretentious members: every number of its Verified paragraph, from `pretend 24`.
- The three claims rows on the recursion in the base in [franel-on-a-digit-design](../../../claims/franel-on-a-digit-design.md).
