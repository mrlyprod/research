# Nyman-Beurling

- The Nyman-Beurling distance restricted to a digit design: `d_N^2(S) = inf over real c of norm(chi - sum_(n in S, n <= N) c_n rho_n)^2` in `L^2(0, infinity)`, `rho_n(t) = {1/(nt)}`, `chi = 1_(0,1]`, for the full set `S = N` and for `S = S_F` at base 3 `{0,1}`, `N = 3^level`.
- The Gram entries `<rho_h, rho_k>` for coprime `h, k` are Vasyunin's formula as printed by Bettin and Conrey, `(log 2 pi - gamma)/2 (1/h + 1/k) + ((k - h)/(2hk)) log(h/k) - (pi/(2hk)) (V(h/k) + V(k/h))` with `V(h/k) = sum_(m=1)^(k-1) {mh/k} cot(pi m/k)`, and `<rho_(gh), rho_(gk)> = <rho_h, rho_k>/g`; the target is `<chi, rho_n> = (log n + 1 - gamma)/n` and `d_N^2 = 1 - b^T G^(-1) b`.
- The `(0,1)` variant drops `int_1^infinity`, so its Gram matrix is `G - v v^T` with `v_n = 1/n`.
- Every Gram system is solved in PARI at 60 digits and again at 90, and the 1-norm condition number `norm(G)_1 norm(G^(-1))_1` is printed beside each distance.
- The local floor `delta_K(S)`: the least value of `sigma^2 + sum_(k=1)^(K-1) int_k^(k+1) (sigma u - 1 - sum_(n in S, n <= k) c_n floor(k/n))^2 du/u^2` over real `sigma` and `c`, a finite least squares; `d_N^2(S) >= delta_K(S)` for every `N` and `K`, and `delta_K(S) > 0` exactly when a squarefree integer below `K` is missing from `S`. `q_0(S)` is the least squarefree integer outside `S`.

## RUN

- `uv run python research/zeta/lab/py/nyman-beurling/nyman_beurling.py gram`: Vasyunin's formula against exact quadrature, `(1/(hk)) int_0^1 {hw}{kw} zetahurwitz(2, w) dw` on the pieces between the breakpoints, at eleven pairs, and `d_N^2` at four small sets against the quadrature of the residual itself, with a Hurwitz tail. About 2 seconds.
- `uv run python research/zeta/lab/py/nyman-beurling/nyman_beurling.py table`: the table, full set at `N = 3^2..3^6` and the design at `N = 3^1..3^7`, both precisions. About 165 seconds.
- `uv run python research/zeta/lab/py/nyman-beurling/nyman_beurling.py floor`: `q_0` for eight designs, `delta_K` for base 3 `{0,1}` at `K = 2..40`, the full set at `K = 2..12` and seven other designs at `K = 40`. About 6 seconds; `floor 100` runs the same to `K = 100` in about 125 seconds and reads `0.100877` at `K = 60` and `0.101926` at `K = 100`.
- `uv run python research/zeta/lab/py/nyman-beurling/nyman_beurling.py zeros`: `2 sum_(0 < gamma < T) 1/(1/4 + gamma^2)` plus the smooth tail `2 (log(T/2 pi) + 1)/(2 pi T)` against `2 + gamma - log(4 pi)`, `T = 1000`. About 5 seconds; `zeros 3000` takes about 215 seconds and reads `0.04619150` against `0.04619142`.
- Prints only; the PARI script it runs is written to the study's own scratch folder under a gitignored `data/` at the repo root.

## WITNESSES

- meter-and-closure.md, THE CLOSURE PROBLEM, the formula row: the eleven pairs agree to `1e-57` at 40 digits, `norm(rho_1)^2 = log 2 pi - gamma`, and the four residual quadratures agree with `1 - b^T G^(-1) b` to every printed digit.
- meter-and-closure.md, THE CLOSURE PROBLEM, the table row: `0.024015, 0.015314, 0.010949, 0.008254, 0.006738` on the full set, `0.116043, 0.106008, 0.104951, 0.104689, 0.104610` on the design, rounded up, the condition numbers, `d_N^2 log N` falling through `0.046191`.
- meter-and-closure.md, THE CLOSURE PROBLEM, the floor rows: `delta_3 = 0.070873`, `delta_40 = 0.098689`, `delta_100 = 0.101926` against `d_2187^2 = 0.104579`, the full set at `1e-75`, and the seven designs at `K = 40`.
- meter-and-closure.md, THE CLOSURE PROBLEM, the constant: `0.0461922` from 649 zeros against `0.0461914`.
