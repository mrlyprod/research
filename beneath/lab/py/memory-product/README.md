# Memory Product

- The coupling rigidity of the memory dial: a live width-`k` rule has `kappa = 0` exactly when it is a product `F^k`, over any finite alphabet.
- The proof is three steps, in [beneath](../../../notes/memory-dial.md), The coupling: Fekete gives `N_W(n) >= rho^n` for every `n`, `kappa = 0` then pins `N_W(2k) = #W^2`, and that splice closure forces `W = F^k`.
- It also re-reads the five preset rules of the dial from a second generator, verb `named`.
- This study is the falsification: the equivalence `N_W(2k) = #W^2` iff `W` is a product is checked on every non-empty rule of eight small cases, on every single-window perturbation of every product at five larger ones, and on seeded random rules.

## THE TEST

- A rule is a set `W` of windows over `a` letters; `N_W(2k)` is counted as the number of pairs `(u, v)` in `W x W` whose concatenation `uv` has every window in `W`, which is every accepted word of length `2k` because the first and second windows determine it.
- A rule is closed when `N_W(2k) = #W^2` and a product when `#W = (number of letters used)^k`.
- `mismatch` counts non-empty rules where closed and product disagree; the theorem says `0`, dead rules included.
- `rho` is the largest modulus of `numpy.linalg.eigvals` of the transfer matrix, a rule being live when it exceeds `1/2`, since a nonnegative integer matrix has `rho = 0` or `rho >= 1`.
- `kappa = log_2(#W)/k - log_2 rho`, compared against `log_2(#W^2 / N_W(2k)) / (2k)` and against the uniform gap `log_2(#W^2 / (#W^2 - 1)) / (2k)` with tolerance `1e-6`; both must sit below `kappa` on every live non-product.
- The generator is cross-checked: `N_W(n)` by brute-force words, by `1^T A^(n-k+1) 1` and by the pair count agree on every non-empty rule at `(a, k) = (2,2), (2,3), (3,2)`.

## RUN

- `uv run python -I research/beneath/lab/py/memory-product/product.py` from the repository root, `12` seconds; the verbs `golden`, `named`, `generator`, `fekete` and `splice` run alone.
- Prints only, writes nothing; the random sweep is seeded.

## WITNESSES

- [beneath](../../../notes/memory-dial.md), The coupling - the golden rule, code `7` at `(1,2)`: `N(4) = 8 < 9`, the first broken pair `01`, `10` with the window `11` outside `W`, `kappa = 0.098239` against the pair bound `0.042481`: verb `golden`.
- The rule allowing `00` and `11` only: `N(4) = 2 < 4`, the broken pair `00`, `11` with the window `01` outside `W`, `kappa = 0.5`: verb `golden`.
- `N(n) >= rho^n` for `n = k-1..14` on every live rule at `(2,2)`, `(2,3)`, `(3,2)`, `0` below one: verb `fekete`.
- The equivalence closed iff product, `0` mismatches, on every non-empty rule at `(a, k) = (2,1), (2,2), (2,3), (2,4), (3,1), (3,2), (4,1), (4,2)`, `65535` rules each at `(2,4)` and `(4,2)`: verb `splice`, lines `exhaustive`.
- The same on all single-window perturbations of every product at `(2,5), (2,6), (3,3), (3,4), (4,3)`, `93, 189, 183, 561, 950` distinct rules, none a product, none closed: verb `splice`, lines `near-product`.
- The least `kappa` over live non-products falls with `k`: `0.098239, 0.056639, 0.029945` exhaustively at `(2,2), (2,3), (2,4)` and `0.015614, 0.008105` on the perturbations at `(2,5), (2,6)`, so the gap of the theorem is a gap at each `(#W, k)` and not a constant: verb `splice`.
- `named`: the golden, supergolden, plastic and tribonacci rules at codes `7`, `23`, `54`, `127` and the no-repeat rule at `(alphabet, k) = (4,2)` code `31710`, each with its counts at levels `1` to `8` (`1` to `4` for the no-repeat rule) by brute-force words and by the transfer matrix, `rho` to nine digits and `kappa` to six: Fibonacci `2, 3, 5, 8, 13, 21, 34, 55` at `1.618033989` and `0.098239`; `2, 4, 4, 6, 9, 13, 19, 28` at `1.465571232` and `0.115204`; `2, 4, 4, 5, 7, 9, 12, 16` at `1.324717957` and `0.260981`; `2, 4, 7, 13, 24, 44, 81, 149` at `1.839286755` and `0.056639`; `4, 12, 36, 108` at `3.000000000` and `0.207519`.
