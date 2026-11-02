# tree-matching

- Rebuilds the level-`level` Sierpinski triangle of design code 7 at base 2 straight from `i AND j = 0`, no Kronecker step, with the axis-nearest-neighbour graph, `level = 1..10`.
- Asserts the graph is a tree: `3^level - 1` edges, one component; and the recursion the proof uses: each of the three blocks equals the level below, shifted, and exactly two link edges cross, `(0, h-1) - (0, h)` and `(h-1, 0) - (h, 0)`.
- Computes the matching number `nu` two ways, leaf stripping and a Hopcroft-Karp bipartite matching on the parity classes of `i + j`, and again with the corner `(0, 0)` removed, which must lower it by exactly one.
- Checks `n - 2 nu = 3^(level-1)`; at `level = 1..8` also the nullity of `A` and the multiplicity of eigenvalue 1 of `I - D^-1/2 A D^-1/2` by dense diagonalisation, tolerance `1e-9`.
- At `level = 1..4` checks `det(xI - A) = sum (-1)^k m_k x^(n-2k)` in exact integers, the characteristic polynomial by Faddeev-LeVerrier against the matching counts `m_k` by tree dynamic programming.
- Exits 1 and names the level if any check fails.

## RUN

`uv run python -I research/complexity/lab/py/tree-matching/matching.py`

Needs numpy and scipy. One pass, about 75 seconds, exits 0.

## WITNESSES

- complexity.md `1, 3, 9, 27, 81, 243, 729, 2187`, the multiplicity of eigenvalue 1, here the `n - 2 nu`, `null_A` and `mult_1` columns
- complexity.md `6561` and `19683`, the multiplicity at levels 9 and 10, the `n - 2 nu` column
- complexity.md `59049 nodes`, the top level of the tree and matching check
- complexity.md `at level 1..4 the exact integer characteristic polynomial equals the matching expansion`
