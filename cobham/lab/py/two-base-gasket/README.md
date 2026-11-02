# two-base-gasket

- Computes `C(N) = card(A cap B cap [0, N)^2)` exactly, where `A` is the base-2 gasket `{(x, y) : x AND y = 0}` and `B` is the base-3 gasket, the pairs whose base-3 digit pairs all lie in `{(0,0), (0,1), (1,0)}`.
- Exact integer arithmetic throughout; no sampling, no fit.
- The fast count walks the base-3 exponents from the top down, choosing for each whether `3^i` joins `x`, joins `y` or neither, and prunes a branch when the high halves of `x` and `y` cannot be made disjoint: with `B_k` the bit length of `3^k` and `k` exponents still unread, the remaining addition to each coordinate is below `3^k`, so each final high part `x >> B_k` is its current value or one more, and a branch dies when all four combinations of the two high parts share a bit.
- Also prints the naive planar budget `dim A + dim B - 2` and the axis lower bound `2^(m+1) - 1`, with the bound truncated down and the budget rounded up.
- The recurrence test is exact over `Q`: Bareiss integer Hankel determinants, and one rational Gauss-Jordan per order on every equation the terms supply, so a reported absence is a proof that no such recurrence fits the terms given, and nothing else.

## VERBS

- `budget` prints `dim A = log_2 3`, `dim B = log_3 3`, the budget `dim A + dim B - 2` rounded up, the axis exponent `log_3 2` truncated down, and the excess. Runtime under `0.05` s.
- `control M` recomputes `C(3^m)` for `m = 0..M` twice: once by the pruned walk and once by testing every pair `(x, y)` in `[0, 3^m)^2` against both digit rules with one shared routine, and prints whether the two agree and whether the axis bound holds. Runtime `0.2` s at `M = 6`.
- `terms LO HI` prints `C(3^m)`, the axis bound `2^(m+1) - 1`, the ratio to the level below and `log_3 C(3^m) / m`, asserting the axis bound at every level. Runtime `3` min `13` s for `LO HI` equal to `0 24` on one core, the cost per level rising by a factor near `2.1`; `m = 25` costs about `3.5` min more and is not run.
- `hankel HI [RMAX]` recomputes `C(3^m)` for `m = 0..HI` and prints them, prints the Hankel determinants `det H_k` of the matrix with entries `a[i+j]` for `k = 1..RMAX+1`, and for each order `r = 1..RMAX` solves the system `a(n) = c_1 a(n-1) + ... + c_r a(n-r)` over every `n` with `r <= n <= HI`, printing the number of equations, the equations spare after the `r` unknowns, and either `none` for an inconsistent system or the coefficients. A consistent order also gets its characteristic polynomial and the modulus of its dominant root from `gp -q`. `RMAX` defaults to `12`, which is the largest order `HI = 24` leaves a spare equation for, order `r` wanting `2r+1` terms and so order `13` wanting `HI = 26`. Runtime `3` min `15` s at `HI = 24`, the census being the whole cost and the algebra under `0.01` s.
- `selftest` runs the algebra of `hankel` against answers known in advance: `bareiss` on five determinants computed by hand, reading `3`, `-2`, `-6`, `0` and `4`, the third needing a pivot swap and the fourth singular; the Fibonacci numbers, first fitted at order `2` with characteristic polynomial `x^2 - x - 1`; `2^n + 3^n + 1`, first fitted at order `3` with dominant root `3.0000000` from `gp -q`; and `n^3 + 2^n`, fitted at no order below `5` and at `5` with `(x - 1)^4 (x - 2)`. Eleven checks, each printing its expected value beside what it got. A solver that answered `none` to everything would fail the three positive fits, and one that fitted anything would fail `n^3 + 2^n` below order `5`. Runtime `0.05` s.

## KUMMER

- `kummer.py` is the variant in which the two gaskets are the Kummer conditions: `K(3^m) = card{(i, j) : 0 <= i, j < 3^m, gcd(binomial(i + j, i), 6) = 1}`, the pairs with no carry in base 2 (`i AND j = 0`) and none in base 3 (every base-3 digit pair sums to at most 2).
- The walk builds the six allowed base-3 digit pairs from the top down with numpy arrays and prunes with the same proved bit-length cut as `census.py`, then keeps the pairs with `x AND y = 0`.
- The row count `row(n) = card{k : gcd(binomial(n, k), 6) = 1}` walks the submasks of `n` in binary and keeps those with every base-3 digit at most that of `n`.
- The box count `B(N) = card{(i, j) in [0, N)^2 : gcd(binomial(i + j, i), 6) = 1}` adds one index at a time.
- `kummer.gp` is the PARI route: `gp -q kummer.gp` from this folder prints `K(3^m)` at `m = 0..7` and `row(n)` at `n = 0..299` from `gcd(binomial(n, k), 6) = 1` alone, about 30 s, equal to `kummer.py` terms `0 7` and `rowlist 299`.
- Direct controls use `math.comb` and `gcd` only; a PARI `gp -q` count of `gcd(binomial(i + j, i), 6) == 1` is the independent second generator.

## KUMMER VERBS

- `terms LO HI` prints `K(3^m)`, the ratio to the level below and `log_3 K(3^m) / m`. Runtime `0.1` s at `m = 10`, `30` s at `m = 13`; `m = 14` was stopped after two minutes and is not run.
- `control M` rebuilds `K(3^m)` for `m = 0..M` by the walk and by testing every pair with `gcd(binomial(i + j, i), 6)`, and prints whether they agree.
- `rows M` prints `row(n)` beside the direct `gcd` count for `n = 0..M`; `rowlist M` and `bfile M` print `n row(n)`.
- `ident M` prints `K(3^m)` beside the sum of `row(n)` over `n < 3^m` for `m = 0..M`.
- `boxcheck N` builds `B(n)` for `n <= N`, compares it with the `gcd` count for `n <= 60`, and prints `K(3^m)` beside `B(3^m)`; `boxlist N` prints `n B(n)`.
- `budget` prints `log_2 3 + log_3 6 - 2`, the exponent of two independent gaskets.

```
uv run python research/cobham/lab/py/two-base-gasket/kummer.py terms 0 13
uv run python research/cobham/lab/py/two-base-gasket/kummer.py control 6
uv run python research/cobham/lab/py/two-base-gasket/kummer.py ident 10
uv run python research/cobham/lab/py/two-base-gasket/kummer.py rows 300
uv run python research/cobham/lab/py/two-base-gasket/kummer.py boxcheck 2200
uv run python research/cobham/lab/py/two-base-gasket/kummer.py bfile 10000
uv run python research/cobham/lab/py/two-base-gasket/kummer.py budget
```

## RUN

- From the repo root.

```
uv run python research/cobham/lab/py/two-base-gasket/census.py budget
uv run python research/cobham/lab/py/two-base-gasket/census.py control 6
uv run python research/cobham/lab/py/two-base-gasket/census.py terms 0 24
uv run python research/cobham/lab/py/two-base-gasket/census.py hankel 24 12
uv run python research/cobham/lab/py/two-base-gasket/census.py selftest
```

## WITNESSES

- `gasket-three-base-and-two-frequencies.md:12` the budget `0.584963`
- `gasket-three-base-and-two-frequencies.md:13` the twenty-five terms `C(3^m)` at `m = 0..24` and the control that rebuilds them to `m = 6`
- `gasket-three-base-and-two-frequencies.md:14` the axis lower bound `C(3^m) >= 2^(m+1) - 1`
- `gasket-three-base-and-two-frequencies.md:15` the excess `log_3 2 >= 0.630929` over `0.584963`
- `gasket-three-base-and-two-frequencies.md:17` the falling readings `log_3 C(3^m) / m` at `m = 14..24`
