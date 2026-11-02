# tree-walks

- The trees of the walk census: `bang dim 2, base 3, code 79` and `bang dim 2, code 7` (base 2, 3 of 4), the L designs with a full row and a full column, and every design of bases 2, 3, 4 that is a tree at all levels. Cell `(a, c)` of a tile is bit `base*a + c`; the render is the Kronecker power built cell by cell, with no code shared with the crate.
- `render` builds both designs at levels 1..6 (base 3) and 1..9 (base 2), prints nodes, edges, and two FNV-1a digests (the sorted cell list, the sorted edge list over flat cell indices), and asserts every level is a tree (edges = nodes - 1, connected). `../../rs/tree-walks` prints the same lines from the crate's own render and graph; the two outputs diff to nothing.
- `resist` computes the corner resistances `R(A, B)`, `R(O, A)`, `R(O, B)` by exact star-mesh (Schur) elimination over `Fraction`, which never assumes a tree, and checks them against the BFS distance, the closed forms `2 b^n - 2` and `b^n - 1`, and a float sparse solve. From level 2 it also measures the graph distance between the two facing diagonal cells of the copies at tile cells `(0, 1)` and `(1, 0)` and asserts `4 b^(n-1) - 2`. At 130 nodes or fewer it solves the hitting times by exact Gaussian elimination and asserts the commute time `E_A T_B + E_B T_A = 2 (nodes - 1) R(A, B)`.
- `volume` builds the level `L + 1` graph and takes the centres in its corner copy of level `L` (`L = 5` for code 79, `L = 8` for code 7), so a ball of radius up to `b^L - 2` is exact in the infinite union. All-pairs distances by unweighted Dijkstra. It asserts `V(x, 2 (b^k - 1)) >= fill^k`, `V(x, b^k - 2) <= 7 fill^k` and that a ball of radius `b^k - 2` meets at most 7 blocks of side `b^k`, prints the maximum observed (3), the range of `V(x, r)/r^alpha` per window `[b^k, b^(k+1))`, and the largest `d(O, y)/max(y0, y1)`.
- `scale` reproduces the spectral reading of the walk census: the four lowest nonzero eigenvalues of the combinatorial Laplacian (dense `eigh`) at two level pairs per design, `d_w = log(lambda_k(n)/lambda_k(n+1))/log b`, against the exact `log(fill b)/log b`; and the mean exit time from the corner at radius `b^k` by a sparse solve.
- `exact` computes that exit time `E_O tau_(O, b^k)` in exact rationals by leaf-to-root elimination on the tree (`k <= 5` for code 79, `k <= 8` for code 7) and asserts it equals `b^k (fill^k + 1)/2`.
- `siblings` takes every code of dimension 2 and base 2, 3, 4 (15, 511, 65535), applies the tile test (the tile graph is a tree and the two opposite faces of each present direction share exactly one cell), asserts it agrees with "the render is a tree" at levels 1..4, 1..3, 1..2, then for each class under the 8 symmetries of the square builds the port-distance recursion `d(n+1) = M d(n) + c`, asserts it equals the BFS port distances on the rendered graph at levels 1..3, and sorts the class as uniform (every active pair grows at the one Perron rate `lambda > 1`, with no chain of two classes at that rate) or not. Prints `fill`, `lambda`, `d_w = log(fill lambda)/log lambda` in the graph metric and `log(fill lambda)/log b`.

## RUN

- `uv run --project ../mrlyprod python -I research/walks/lab/py/tree-walks/trees.py` runs every verb in about 22 s: `render` 0.2 s, `resist` 2 s, `volume` 13 s, `scale` 3 s, `exact` 0.2 s, `siblings` 3 s.
- One verb by name: `... trees.py volume`.
- Pin to the crate: `cargo run --release --manifest-path research/walks/lab/rs/tree-walks/Cargo.toml` prints the same 15 lines as `render`.
- Needs numpy and scipy; prints only, writes nothing.

## WITNESSES

- trees.md, section "Two designs are trees": the 15 digest lines of `render` equal to the crate's; nodes `3^n` and `5^n`, edges one fewer.
- trees.md, section "Resistance and distance": `R(A, B) = 4, 16, 52, 160, 484, 1456` at levels 1..6 for code 79 and `2, 6, 14, 30, 62, 126, 254, 510, 1022` at levels 1..9 for code 7; `R(O, A)` half of each; the commute times `32, 768, 12896` and `8, 96, 728, 4800`.
- trees.md, section "Volume growth in the graph metric": the per-window ranges of `V/r^alpha`, the block counts, the ratio `d(O, y)/max(y0, y1)` at most 1.9961.
- trees.md, section "What the readings say": spectral `d_w` 2.4662 (code 79, levels 4 to 5) and 2.5865 (code 7, levels 6 to 7) against 2.464974 and 2.584963; the exit times `9, 117, 1701, 25353, 379809` and `4, 20, 112, 656, 3904, 23360, 140032, 839936`.
- trees.md, section "The siblings": 8, 15 and 508 tree codes in bases 2, 3, 4 in 2, 5 and 69 classes, uniform in 2, 5 and 33; the table by `fill` and `lambda`.
