# eisenstein-farey

- `new(n)`: the distinct positive distances from the origin to points of `(1/n) Z[w]` in the closed disc of radius `sqrt 3` that no scale `m < n` holds.
- `new.gp`: the Mobius sum `sum_(d | n) mu(d) L(3 n^2 / d^2)` over a Loeschian sieve (`qfrep`); prints `new(1..n)` for the `n` set on its fourth line (200); `gp -q new.gp`, about 1 s.
- `union.py`: stdlib only; norms `a^2 - a b + b^2` of every lattice point, each reduced to lowest terms as `k / n^2` per scale, with the new ones counted against the union of the scales below; neither the Loeschian criterion nor the Mobius sum.
- Run: `uv run --project <repo> python -I union.py 200`, about 1 s; the argument is the top `n`, default 200.
- The two outputs agree at `n = 1..200` and equal the first 30 terms of `spin-census hex`.
- Raising that `n` and the `union.py` argument to `400` keeps the two equal, about 9 s and 8 s.
