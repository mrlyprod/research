# oeis-terms

- Regenerates the two b-files this tree submitted to the OEIS, each by more than one method.
- A396934: pairs `0 <= i, j < 2^n` with `i AND j = 0` and `gcd(i, j) = 1`, for `n = 0..20`.
- The row walk indexes a pair by `m = i + j = i OR j`, keeps the odd submasks of odd rows with `gcd(i, m) = 1`, and doubles; eight threads split the rows.
- Cross-checks: the complement submask walk with `gcd(i, j)` to `n = 14`, and Pascal's triangle mod 2 by the additive recurrence to `n = 12`.
- Prints `a(20)/3^20`, the limit `16/(3*Pi^2)`, and their gap.
- A398348: `n x n x n` binary arrays up to `D_n^3` semidirect `S_3`, order `48*n^3`, by Burnside for `n = 1..14`.
- The cycle walk maps every cell under `(x_0, x_1, x_2) -> (eps_t * x_{p(t)} + s_t mod n)` and walks the cycles.
- Cross-checks: fixed points of the powers of the affine map `x -> M x + t` to `n = 8`, and direct orbit enumeration of all `2^(n^3)` colourings to `n = 2`.
- Uses `mrly_num::factor::coprime` and `num-bigint` for the Burnside sums.

## RUN

- `cargo run --release --manifest-path research/bijection/lab/rs/oeis-terms/Cargo.toml`
- An optional argument lowers the top A396934 level, e.g. `-- 16`.
- The full domain runs in under half a minute on eight cores.

## WITNESSES

- shared/sequences.md:65: `0, 2, 4, 12, 34, 122, 362, 1130, 3406, 10506, 31550, 95260`, and :68 and :71 the b-file to `n = 20`, ending `1884174908`.
- shared/sequences.md:67: `16/(3*Pi^2) = 0.5403796` and `a(20)/3^20 = 0.5403761`.
- coprime.md:113: `2, 4, 12, 34, 122, 362` from `n = 1`, `a(20)/3^20 = 0.5403760862` against `0.5403796461`, gap `-3.6e-06`.
- shared/sequences.md:80: `2, 22, 111618, 6005363762644688, 7089215977519836239803174210135872, 10157435539019790383692007859901914095646506996125324171134976`, `a(7)` with 100 digits and `a(8)` with 150, and :83 the b-file to `n = 14`.
- shared/sequences.md:86: the affine-power second route to `n = 8` and orbit enumeration at `n = 1, 2`.
- bijection.md:184 and README.md:47: `2, 22, 111618, 6005363762644688, 7089215977519836239803174210135872`, with bijection.md:188 the b-file to `n = 14`.
- Not regenerated here: the `2^27` flood fill of shared/sequences.md:86, and the A255016 convention check on the same line.
