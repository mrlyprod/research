# flake-band-gap

- Builds the base-2 code-23 flake at any `level`: `4^level` cells, face adjacency by plus-one lookup on each axis, rooted by breadth-first layers from the origin cell.
- Asserts the face graph is a tree: `4^level` nodes, `4^level - 1` edges, every node reached.
- Counts eigenvalues below a shift by Sylvester's law on the leaf-to-root elimination of `Lap - shift*I`, never forming a matrix, and bisects both band edges: `lo(level)` in `(1, 3)` against the target `3*4^(level-1)`, and the top of the spectrum `hi(level)` in `(4, 8)` against `4^level`.
- Repeats the elimination over `Fraction` at shifts `2`, `4`, `4 - 1e-9`, `4 + 1e-9`, then back-substitutes the null vector at `4`, clears it to integers, and takes `(Lap - 4I)v` in exact integer arithmetic.
- One generator. Float bisection reaches `level 11`; the exact counts reach `level 6` and the integer null vector `level 4`.

## RUN

```
uv run python research/complexity/lab/py/flake-band-gap/flake.py
uv run python research/complexity/lab/py/flake-band-gap/flake.py 11 6 4
```

- Both commands are the same run; the arguments are `top exact_top vector_top`. Under ten seconds on one core.
- The source ran the float sweep to `level 11` and the page table stops at `level 10`; this study prints both. At `level 12` the bisection has spent its double precision and the ratio column stops rising, so `11` is the top it certifies.

## WITNESSES

- `README.md:51` lower edge `1.000000, 1.827520, 1.975680, 1.996862, 1.999605, 1.999950`, upper edge exactly 4 and simple, `3*4^(level-1)` below 2, `c` near `12.9868`.
- the claims line the same edge list, and the top of the spectrum `3 + sqrt(5) = 5.2360679775` at `level 2` climbing to `5.7090316570` at `level 6`.
- `complexity.md:451-456` the edge table `1.000000, 1.827520, 1.975680, 1.996862, 1.999605, 1.999950`, upper edge 4, zero strictly inside.
- `complexity.md:465-468` 4 is not the top: the spectrum climbs to about `5.7090`, split `3*4^(level-1)` below 2 and `4^(level-1)` at or above 4.
- `complexity.md:480-489` the defect table to `level 10`, ratios `5.7978, 7.0920, 7.7494, 7.9384, 7.9850, 7.9963, 7.9991, 7.9998, 7.9999`, scaled defect `12.986750`.
- `complexity.md:494` `c` near `12.9868`: the run prints `12.986773` at `level 11`.
- `complexity.md:621` upper edge exactly 4: root pivot `0`, no earlier zero pivot, residual `0`, jump of one at 4.
- The source's stored float values `4.000000000000001` and `3.9999999999999893` are a dense solver's noise; this study prints the exact upper edge instead and carries no such constant.
- `complexity.md:639` names the source generators and carries no number.
