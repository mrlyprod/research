# carpet-nodal

- Builds the Sierpinski carpet as the Kronecker power at `level` of the 3 x 3 tile with the centre void, `8^level` cells, asserted at 512 for `level 3` and 4096 for `level 4`.
- Takes the cell graph on 4-neighbour adjacency between filled cells, checks it is one component, and diagonalises the combinatorial Laplacian `D - A`, degree minus adjacency, densely with `numpy.linalg.eigh`.
- For the `k`-th returned eigenvector counts strong nodal domains `nu_k`: connected components of the positive cells plus connected components of the negative cells, a cell with `|v| < 1e-9` belonging to neither.
- Clusters the spectrum at gap tolerance `1e-8` into eigenvalue classes, reads the multiplicity `r_k` and the last index of each class, and checks every `k` against Courant `nu_k <= k` and against the discrete nodal domain theorem `nu_k <= k + r_k - 1`.
- Prints per graph: `nu_k`, `r_k` and the zero-cell count for `k = 1..12`, the maxima, the Courant and discrete-theorem fractions, the class and degeneracy counts, the five-bin histogram of `nu_k / k`, the top eigenvector's sign pattern, and the first double eigenvalue where the two returned vectors, their sum and their difference disagree in `nu`.
- Controls: the `22 x 22` and `64 x 64` grid graphs, once in numpy's basis and once in the separable cosine basis ordered by eigenvalue then by the mode pair `(p, q)`, whose count is checked against `(p + 1)(q + 1)` derived from the one-dimensional sign changes.
- On a degenerate eigenvalue the nodal count depends on the basis, so every `nu_k` is a fact about the returned basis only; the theorem bound `k + r_k - 1` is the only basis-free statement, and the study claims nothing about the eigenspace.
- Domain: `level 3, 4`, zero tolerance `1e-9`, multiplicity tolerance `1e-8`; `level 5` has 32768 cells and is not attempted.

## RUN

`uv run python research/automata/lab/py/carpet-nodal/nodal.py`

Needs numpy and scipy. One pass, about 30 seconds, exits 0.

## WITNESSES

- automata.md section 13, the nodal domain lines: `nu_k` for `k = 1..12` at `level 3` and `level 4`, the Courant fraction 1 on all 512 and all 4096 returned eigenvectors, the degenerate index counts 262 of 512 and 2062 of 4096, the maximal multiplicities 4 and 20, and the double eigenvalue at `k = 6, 7` whose returned vectors have 4 domains while their sum and difference have 6.
- The control lines: the separable count `(p + 1)(q + 1)` on all 484 and 4096 grid eigenvectors, and the mean of `nu_k / k` over `k >= 2`, `0.605` and `0.536` on the carpet against `0.508` and `0.465` on the separable grid.
- The top eigenvector line: no edge joins two nonzero cells of one sign on any of the six graphs, and at `level 4` twelve cells of that vector fall below the zero tolerance, so its printed count is 4084.
