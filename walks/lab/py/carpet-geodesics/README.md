# carpet-geodesics

- Shortest paths in `bang dim 2, code 7` at odd side `N` and level `L`: the side-`N` tile voids a cell iff both digits are odd, the render is its `L`-th Kronecker power in the unit square, and a path runs in the union of the closed filled cells.
- `code` prints the side-`N` tile as a base-`N` code, bit `N*i + j`, for `N = 3..11`, and asserts `495` at `N = 3`.
- `corner` prints `D`, `N (D - sqrt(2))` and `2 - D`. It computes the corner distance `D(N, L) = d((0,0), (1,1))` exactly: Dijkstra with the straight-line heuristic on the visibility graph of void corners, a segment blocked iff it meets the interior of a void (Liang-Barsky on integer coordinates), voids and corners pruned by the ellipse `|p| + |p - (1,1)| <= bound`, the bound widened until a path is found under it, which makes the pruning exact. It asserts `sqrt(2) <= D <= 2 - (2 - sqrt(2)) (1 - 1/N)^L`, `D` nondecreasing in `L`, and `D(3, L) = 2 sqrt(5)/3` to `1e-12`.
- `level1` runs `corner` at `L = 1` for every odd `N` from 3 to 41 against `sqrt(2) + (2 sqrt(5) - 3 sqrt(2))/N`.
- `map` iterates the homogenisation map one-sidedly. The street grid has period 2 and void `(1, 2)^2`; the edges are the visible segments between translates of the four void corners within a window of 3 periods (244 edges) and of 6 periods (484 edges). For each of 720 directions `p`, bisection on the 24 simple node cycles of the four corners finds the cycle of largest ratio `p.v` over cost, and its displacement over its cost is a point of the true unit ball. The convex hull of those points is inside the true ball, so its gauge is an upper bound on `nu_L`; that gauge prices the edges of the next level. It prints the upper bounds at four angles and the lower bound they give on the gap to the regular octagon, at `L = 1..4, 5, 10, 15, 20` for both windows, and the explicit path of displacement `(4, 2)`. A wider window only lowers the bounds, and both windows read the same at `L = 1`.
- `bridge` compares the exact level-1 distance from `(0,0)` to `(1, (N-1)/(2N))` with `nu_1` of the same vector at `N = 11, 21, 31, 41`. Here `nu_1` comes from the support function on 1440 directions, by max-plus Floyd-Warshall on the 3-period window.
- `hull` counts the vertices, edges and faces of the convex hull of the 18 unit vectors along the axes and the face diagonals of the cube, and asserts `18, 48, 32`.

## RUN

- `uv run python research/walks/lab/py/carpet-geodesics/geodesics.py` runs every verb in about 7 s: `code` 0 s, `corner` 5 s, `level1` 0.3 s, `map` 1.7 s, `bridge` 0.2 s, `hull` 0.1 s.
- One verb by name: `uv run python research/walks/lab/py/carpet-geodesics/geodesics.py corner 5,3 11,2` runs the two long cases, 52 s and 31 s.
- Needs numpy and scipy; prints only, writes nothing.

## WITNESSES

- shortest-paths.md, section "Shortest paths: two limits that do not commute": the codes `495` and `33226431`; every `D(N, L)` and `N (D(N, L) - sqrt(2))` printed there; the level-1 formula at `N = 3..41` to `4.4e-16`; the upper bounds on `nu_L` at `22.5` degrees, the lower bounds on the gap to the octagon and the explicit path; the bridge readings; the counts `18, 48, 32`, a check of the proof on the page.
