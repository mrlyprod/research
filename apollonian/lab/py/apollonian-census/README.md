# apollonian-census

- Counts the circles of curvature at most `10^n` in the Apollonian packing with root curvatures `-1, 2, 2, 3`, a second route beside the crate verb `census` of `../../rs/apollonian`.
- `census.py` is a chain-collapse Descartes recursion in the standard library: a chain of circles tangent to two fixed circles has quadratic curvatures, and its childless tail is counted in closed form.
- `census.gp` is the plain Descartes recursion in PARI, one call per circle.
- Run: `uv run --project ../../../../.. python -I census.py 6` takes 7 s; the argument is the top `n`, and 7 takes about 2.5 minutes.
- Run: `gp -q < census.gp` takes under a second.
- `census.py` prints the counts with the root quadruple excluded (`5, 165, 3325, 67163, 1359167, 27463391` at `n = 1..6`) and with every circle counted (`1, 9, 169, 3329, 67167, 1359171, 27463395` at `n = 0..6`, the `n = 0` entry a constant), then the local exponents, the ratios to `10^(n delta)` and the cut residues mod 24.
- `census.gp` prints the same two lines for `n = 1..4` and `n = 0..4`.
- Both agree with the crate verb on every count.
