# sides-holding-a-number

- The numbers of Object E on [Two bases](../../../notes/object-e-and-the-schanuel-wall.md): the even-digit design `E_N` at odd side `N`, its integers `Z_N = 2 K_N`, and which sides hold one number.
- Exact arithmetic throughout: integers, `Fraction`, and numpy integer arrays; the only floats are the printed ratios, rounded outward where they are bounds.
- Membership of `p/q` is the orbit rule `N^j p mod 2q in {0, ..., q}`; its control is an independent walk on the digits of `p/q` by `Fraction`.

## VERBS

- `period Q`: for every `p/q` in `(0, 1)` with `q <= Q`, the residues `R(p/q)` of the odd sides holding it mod `2q`; checks the orbit rule against the digit walk at every `p/q` in `[0, 1]` with `q <= Q`, over odd sides `3..6q+1` at `q <= 40` and one period `3..2q+1` above, then the symmetry `p -> q - p`, the least period `2q`, the bounds `2/q` and `(q+1)/(2q)` or `1/2`, the rule mod `q` at odd `q`, the cyclic-subgroup count of the unit residues, the prime formula, and the invariance `r -> q - r` at even `q`; prints the largest shares and the odd primes with `share(1/q) > 2/q`. Runtime `19.6` s at `Q = 200`.
- `eisenstein Q`: at every odd prime `q <= Q` and every `p`, the parity of the number of odd sides `3 <= N < q` whose first digit of `p/q` is odd, against the Legendre symbol of the odd one of `p, q - p` by Euler's criterion; also the exact sum of `floor(pN/q)` over the odd `N < 2q` against `(2p - 1)(q - 1)/2 + p`, and how often its parity matches the symbol. Runtime `0.03` s at `Q = 200`.
- `count K`: `c(k) = card{odd N : 3 <= N <= 2k, 2k in Z_N}` at every `k <= K`, the small sides by enumerating `K_N` and the sides above `sqrt(2K)` by the parity of `floor(2k/N)` on a difference array; checks every `k <= 3000` by the direct digit test; prints the extremes of `(c(k) - (1 - log 2) k)/sqrt(k)` and the largest error against the bound `sqrt(2k) + 1`. Runtime `2.1` s at `K = 10^6`.
- `second E`: `c(k)` by an `O(sqrt k)` block count, checked against the direct test at every `k <= 3000`, at `60` random `k` in each `[10^e, 2 10^e)`, `e = 6..E`, seed `376`; prints the mean and spread of the error over `sqrt(k)` and of the small-side part against `kappa` and `beta`. Runtime `7.6` s at `E = 11`.
- `family`: the integers below `3000` held by every odd side; the rejection of every `p/q`, `q <= 60`, at side `2q - 1`; `E_N` inside `E_(N^e)` on integers and on fractions, and its strictness at `e >= 2` by `N + 1` and `(N + 1)/N^e`; Kummer's reading of `K_p` at odd primes `p <= 23` by `math.comb`; the carry reading of `K_(p^a)`; and the least witnesses at sides `9` and `15`. Runtime `1.2` s.
- `inter X [FILE]`: the integers below `X` held by each of nine side sets, by a walk over `K_3` from the top digit that cuts a branch when the least member of another `K_N` at or above the branch's least value lies past its largest; checks the walk against the direct digit test below `k = 10^6` at four side sets, and, when `FILE` is a copy of the [A030979](https://oeis.org/A030979) record from `https://raw.githubusercontent.com/oeis/oeisdata/main/seq/A030/A030979.seq`, the sides `{3, 5, 7}` against its terms; without `FILE` that comparison is skipped and the URL printed. Runtime `1.2` s at `X = 10^12`.
- `deep D`: the same walk for the sides `{3, 5, 7, 11}` below `10^D`. Runtime `14.6` s at `D = 1000`.
- `weyl M L`: for five irrationals from `90`-digit truncations, the share of the odd sides `3 <= N <= 2M + 1` holding each to level `L`, times `2^L`. Runtime `0.2` s at `M = 10^5`, `L = 6`.

## RUN

```
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py period 200
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py eisenstein 200
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py count 1000000
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py second 11
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py family
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py inter 1e12 A030979.seq
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py deep 1000
uv run python research/cobham/lab/py/sides-holding-a-number/sides.py weyl 100000 6
```

## WITNESSES

- `object-e-and-the-schanuel-wall.md`, Object E: the orbit rule against the digit walk, `1661996` checks, the least period `2q`, the bounds and the maximum `2/3`, the rule mod `q`, the coset count and the prime formula at all `12231` fractions with `q <= 200`, and the shares printed there (`period 200`).
- The half-period parity law and the full-period identity at all `4180` pairs with `q` an odd prime `<= 200` (`eisenstein 200`).
- The readings `0.9821` to `1.0496` of the share to level `L` times `2^L` (`weyl 100000 6`).
- The error extremes `0.547191`, `0.357533`, `[-0.348228, -0.082525]` and `c(10^6) = 306665` (`count 1000000`).
- The means `-0.220703` to `-0.215000` against `kappa = 0.213864` and the small-side readings against `beta = 0.139689` (`second 11`).
- The family checks, the strictness at `e >= 2`, and the witnesses `20, 6, 10, 2, 3, 15` (`family`).
- The counts `10072` and `50`, the `17` integers equal to twice the A030979 terms, and the lists at sides adding `9, 11, 13, 15` (`inter 1e12`); `0, 2, 6320` alone below `10^1000` (`deep 1000`).
