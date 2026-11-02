# JSR Schedules

- Regenerates the joint spectral radius section of [word order, closed and bounded](../../../notes/word-order-closure-and-spectral-radius.md): the observable frame, the six frame images, the triangularity of the component cocycle, the cross-polytope certificate, the joint and lower spectral radius of every subfamily, and the bracket on two named pairs.
- The question is the standard one for a noncommuting matrix family: `JSR(F) = lim_L max_(|w| = L) ||M_w||^(1/L)`, and whether the maximum is attained by a periodic word.
- A word is an ordered list of base-2, dim 2 codes folded by the Kronecker product, first letter outermost; bit `i` of a code is residue corner `i` in row-major order, so code 3 is the top row and code 5 the left column. The empty word is the one filled cell.
- The component representation is rebuilt from scratch rather than imported: Hankel-basis elimination in exact rational arithmetic over the 241 suffixes of length at most 2, prefixes in breadth-first order, with the basis words, the matrices, `lambda` and `gamma` all outputs of the elimination.
- The four observables are drawn from cells: `comp`, 4-connected components; `H` and `V`, the maximal horizontal and vertical runs; `fill`. Each is read on the basis words to give a frame vector, and the frame matrix's determinant is checked to be a unit.
- The frame images `M_c gamma`, `M_c h`, `M_c v`, `M_c phi` are printed in frame coordinates for all 15 codes, asserted nonnegative, and their column sums asserted equal to `(comp, nonempty rows, nonempty columns, fill)` of the letter's own tile.
- Triangularity is asserted, not assumed: every frame matrix is checked lower triangular with nonnegative entries, and its characteristic polynomial is computed as an integer determinant at five points and matched against the product over its diagonal, so the spectral radius is the fill by proof rather than by a solver.
- The transfer laws are checked against drawn cells on 3856 words, printed as the 3616 of length at most 3 and the 240 seeded of length 4 to 7, in both directions: the observable tuple times the frame matrix against the next tuple, and `lambda M_w` against all four observables.
- The certificate is the cross-polytope `P = conv{+/- gamma, +/- h, +/- v, +/- phi}`. The study prints the exact integer residual `k_c` minus the frame `l1` norm of each image and asserts it is nonnegative, which is the inclusion `M_c P subset k_c P` in exact arithmetic.
- The bracket is printed for `{3, 6}` and `{3, 7}`: an exhaustive scan over all `2^L` words to `L = 16`, the largest spectral radius and its word, the largest frame 1-norm, both taken as exact integers and printed as roots truncated down for the lower bound and rounded up for the upper. Both are asserted equal to `(max fill)^L` at every one of the 16 lengths, and the smallest spectral radius to `(min fill)^L`.
- The lifting bound is exact rather than iterated. Kronecker powers of triangular matrices are triangular, so `rho(M_a^(x)k + M_b^(x)k)` for the pair's two letters `a` and `b` is the largest sum of two matched diagonal products, which is `k_a^k + k_b^k`. Both ends are printed at `k = 2, 4, 6`, the upper `(k_a^k + k_b^k)^(1/k)` rounded up and the lower `2^(-1/k)` times it truncated down, with the lower asserted at most the joint spectral radius and the joint spectral radius at most the upper.
- The telescope pair is read separately: the leading `2 x 2` blocks, the common invariant line, and an exhaustive check that every word over `{3, 6}` of length 1, 2, 3, 6 and 10 has spectral radius exactly `2^L`.
- The dead routes are printed with their witnesses: the common right eigenvector that makes every subfamily reducible, and the polytope algorithm's orbit, normalised by the letter fill `k_c`, closing on one vertex and spanning a line.
- The study exits nonzero if the representation, a transfer law, a triangularity, a characteristic polynomial, a certificate residual or a bracket order fails.

## RUN

- `uv run python research/connectivity/lab/py/jsr-schedules/jsr.py`
- About nine seconds; prints only, writes nothing, and holds one `128 x 128` raster at a time.

## WITNESSES

- word-order-closure-and-spectral-radius.md the frame: `gamma = (1,1,1,1)`, `h = (1,1,2,2)`, `v = (1,2,1,2)`, `phi = (1,2,2,4)`, determinant `-1`.
- word-order-closure-and-spectral-radius.md the frame images: the six rows of the transfer table and their column sums `(1,1,1,1)`, `(1,1,2,2)`, `(1,2,1,2)`, `(2,2,2,2)`, `(1,2,2,3)`, `(1,2,2,4)`.
- word-order-closure-and-spectral-radius.md the diagonals: `(0,0,0,1)`, `(0,1,0,2)`, `(0,0,1,2)`, `(0,0,0,2)`, `(1,1,1,3)`, `(1,2,2,4)`, and `rho(M_c) = k_c` at all six classes.
- word-order-closure-and-spectral-radius.md the certificate: residuals `[0,0,0,0]`, `[1,1,0,0]`, `[1,0,1,0]`, `[0,0,0,0]`, `[2,1,1,0]`, `[3,2,2,0]` at codes 1, 3, 5, 6, 7, 15.
- word-order-closure-and-spectral-radius.md the spectral table: the 15 class pairs with their joint and lower spectral radius, and 78 of the 105 letter pairs carrying two different fills.
- word-order-closure-and-spectral-radius.md the bracket: the scan alone gives `[2.000000000, 2.000000000]` on `{3, 6}` and `[3.000000000, 3.000000000]` on `{3, 7}` at `L = 16`; the lifting alone gives `[2.000000000, 2.244924097]` and `[2.710444581, 3.042371177]` at `k = 6`.
- word-order-closure-and-spectral-radius.md the lifting: upper `2.828427125, 2.378414231, 2.244924097` and lower `2.000000000, 2.000000000, 2.000000000` on `{3, 6}`; upper `3.605551276, 3.138288993, 3.042371177` and lower `2.549509756, 2.638975964, 2.710444581` on `{3, 7}`, at `k = 2, 4, 6`.
- word-order-closure-and-spectral-radius.md the telescope: `A = [[0,1],[-2,3]]`, `B = [[2,0],[4,0]]`, `p = (1,2)`, and spectral radius exactly `2^L` on every word of length 1, 2, 3, 6, 10.
- word-order-closure-and-spectral-radius.md the blindness witness: `max |entry(M_3^L)| = 2^(L+2) - 2` reading 6, 14, 62, 1022, 262142, 17179869182 against `comp(A_(3^L)) = 1`.
