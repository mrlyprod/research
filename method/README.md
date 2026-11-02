---
title: Method
lead: How the results here are produced and checked, worked through on the odd-side fill polynomial.
figure: research-method
slug: method
date: 2026-08-28
---

# Method

How the results here are produced and checked, worked through on the odd-side fill polynomial.

- [notes/method.md](notes/method.md): Method
- [notes/baselines.md](notes/baselines.md): Against the baselines
- [claims/coin.md](claims/coin.md): The coin
- [claims/odd-side-fills.md](claims/odd-side-fills.md): Odd-side fills
- [refs.md](refs.md): the references this topic cites
- [lab/py/walsh-fill/](lab/py/walsh-fill/): Computes the Walsh coefficients `hat f(S) = 2^-dim sum_(c in F) (-1)^(|S cap c|)` of every base-2 design at `dim 1..4` by one Hadamard product, and their level sums `W_j = sum_(|S| = j) hat f(S)`, kept as the integers `2^dim W_j`.
- [lean/Method.lean](lean/Method.lean): Proves by kernel `decide` that the cell-by-cell render `solidCells` of `bang dim 3, code 23` fills `4*k^3 - 3*k^2` cells at every odd side `2k - 1` of `1..12`.
- [figures/research-baselines.ts](figures/research-baselines.ts): research-baselines
- [figures/research-baselines.json](figures/research-baselines.json): research-baselines census
- [figures/research-method.ts](figures/research-method.ts): research-method
