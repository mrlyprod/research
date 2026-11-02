# Chaos Game

- The chaos game of a design read as a shift register on digits, in [beneath](../../../notes/chaos-game.md).
- A step `x -> (x + a)/b` pushes the draw `a` on as the new leading base-`b` digit, so the level-`k` cell of `x_n` is the last `k` draws, newest first.
- Every law of the chapter is checked here by exact enumeration in `Fraction` arithmetic, never by sampling, except the two sampled sanity columns named below.

## VERBS

- `register`: the closed form `x_n = sum_(j=1..n) a_(n+1-j) b^(-j) + b^(-n) x_0` and the cell of `x_n` at every level `k <= n` against the reversed last `k` draws, every draw word of length `0..4` at `(base, dim) = (2,1), (2,2), (3,2), (2,3)`, four starts in `[0,1)^dim`, two of them off every attractor.
- `forget`: law 1; `27` designs, every non-empty digit set at `(2,1)`, `(2,2)`, `(3,1)`, the base-3 carpet and the octahedron flake `bang dim 3, code 126`, uniform and weighted draws `p(a) proportional to 1, 2, 3, ...`; the law of the level-`k` cell at every `n = 0..k+1` by summing over every draw word, against the product law `p^k` on cells.
- `paint`: law 2; a de Bruijn schedule from the Lyndon-word construction, played on every design of `forget` with `card F >= 2`; the cells of steps `k..m^k+k-1` and of a second lap; the cyclic perfect paintings counted by brute-force Euler circuits against `(m!)^(m^(k-1))/m^k`.
- `cover`: the expected random painting time by two exact routes that share no code: the alternating sum over restricted games of their counting series at `1/m`, and a recursion over the set of cells already painted; a seeded sample of `20000` runs per case is the sanity column.
- `census`: law 3; every one-step rule on the four corners of the base-2 square, `65536` codes in `8548` classes under the square's `8` symmetries, with `rho` from the squarefree part of the exact characteristic polynomial, irreducibility, products, the rank-one index of the move matrix, the closed cells of each game's limit set against all `15` design sets, and the named games.
- `families`: the magic, mosaic and slice games of `### Across the families`.

## CONVENTIONS

- A corner is `c = x + 2y`, digit vector `(c & 1, c >> 1)`, the convention of the memory dial.
- A game rule is a set of time-ordered pairs `(previous, next)`, code bit `4 previous + next`; the dial rule of its limit set is the reversed rule, since the newest draw is the leading digit.
- The cyclic order of the square is `0, 1, 3, 2`, anticlockwise with `y` up.
- A rule is irreducible when its move graph is strongly connected on the corners it uses.
- The move matrix `P` draws the next corner uniformly among the allowed successors; its rank-one index is the least `m` with every row of `P^m` equal, `None` when there is none.
- The closed cells of a limit set at level `L` are computed by the extremal-tail rule: a word `v` reaches the neighbour cell `int(v) + e`, `e` in `{-1, 0, 1}^2`, iff an infinite continuation has every coordinate `i` digit `0` where `e_i = -1` and `1` where `e_i = +1`; the rule is checked against the cells of exact eventually periodic points.

## RUN

- `uv run python -I research/beneath/lab/py/chaos-game/game.py` from the repository root, about `75` seconds, `cover` and `census` taking most of it; each verb runs alone.
- Prints only, writes nothing; the samples are seeded.

## WITNESSES

- `register`: `241768` checks, `0` mismatches.
- `forget`: `n >= k`: `1328` of `1328` cases exactly the product law, largest distance `0`; `n < k`: `872` of `872` at positive distance, and `436` of `436` uniform cases at distance at least `1 - card F^(n-k)`.
- `paint`: `54` of `54` (design, `k`) schedules paint every level-`k` cell exactly once; perfect cyclic paintings `1, 1, 2, 16` at `m = 2`, `k = 1..4`, `2, 24` at `m = 3`, `6, 20736` at `m = 4`, each equal to `(m!)^(m^(k-1))/m^k`.
- `cover`: `E[T] = 3, 11/2, 25/3` at `k = 1`, `m = 2, 3, 4`, the coupon collector; `19/2`, `18850259/711620 = 26.489220`, `89259/3640 = 24.521703`, `55.068831`, `58.632878` at `(m, k) = (2,2), (3,2), (2,3), (4,2), (2,4)`, both routes equal as rationals, against the coupon value `(k-1) + m^k H(m^k)` `9.333333, 26.460714, 23.742857, 55.091664, 57.091664`.
- `census`: `8548` classes, Burnside `8548`; `3421` irreducible, `5` products, codes `1, 51, 1632, 1911, 65535`; distinct `rho` `185` over all classes and `184` over the irreducible ones; characteristic polynomials `333`; irreducible classes by corners used `1, 6, 76, 3338`; `376` irreducible non-products with integer `rho`.
- `census`: rank-one index `1` on `5` classes, `2` on `16`, `3` on `18`, none on `3382`; the cell law exactly stationary from `n = k + m` on and not before in `68` of `68` cases, never stationary in `80` of `80` sampled cases.
- `census`: `3416` of `3416` irreducible non-products separated from all `15` design sets by their closed cells, at least level `1, 2, 3, 4` in `33, 2059, 1289, 35` classes; `510` (set, level) pairs of the extremal-tail rule against exact points, `0` mismatches.
- `census`, named: no repeat `31710`, not one place anticlockwise `48765`, not the opposite corner `60855`, each `card W = 12`, `rho = 3`, `kappa = 0.207519`, never exact; not a neighbour `38505`, reducible; only a neighbour `27030`, `rho = 2`; only anticlockwise `16770`, `rho = 1`; the shift rule `50115`, `rho = 2`, rank-one index `2`; the gasket product `1911`; no neighbour after a repeat, width `3`, `card W = 56`, `rho = 2 + sqrt 2 = 3.414213562`, dimension `1.771553`, `kappa = 0.164232`.
- `families`, magic: `27` of `27` laws exact on the reversed window, `18` of `18` early laws not; perfect laps `12` and `72` at sizes `3, 2` and `k = 2, 4`, `54` at sizes `3, 3, 2` and `k = 3`, none at the unbalanced `k`.
- `families`, mosaic: substitution matrix `[[2, 1], [1, 1]]`, `rho = 2.618033989`, dimension `1.388484`, filled cells `3, 8, 21, 55, 144, 377, 987, 2584`; the graph game visits `377` of `377` level-6 cells, `0` outside.
- `families`, slices at `x + y = 2/3`: code `7`, states `1/3, 2/3`, `rho = sqrt 2`, dimension `1/2`; code `11`, states `1/3, 2/3, 4/3, 5/3`, the path matrix, `rho = 1.618033989`, dimension `0.694242`, counts `2, 3, 5, 8, ..., 377`; the carry-state game visits `8` of `8` and `21` of `21` level-6 cylinders, `0` outside.
