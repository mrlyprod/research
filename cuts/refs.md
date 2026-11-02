# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A001316 | Gould's sequence: number of odd entries in row n of Pascal's triangle; a(n) = 2^A000120(n) | https://oeis.org/A001316 |
| A047999 | Sierpinski's triangle (or gasket): Pascal's triangle mod 2 | https://oeis.org/A047999 |
| A048883 | a(n) = 3^wt(n); number of odd values in the n-th layer of Pascal's tetrahedron | https://oeis.org/A048883 |
| A268240 | Pascal's tetrahedron of trinomial coefficients (A046816) read mod 2 | https://oeis.org/A268240 |
| A001045 | Jacobsthal sequence `a(n) = (2^n - (-1)^n)/3`, `0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341, 683, 1365, ...`, read from the OEIS dump; the corners of the tent rank law | https://oeis.org/A001045 |
| A299916 | Name is `a(n) = A299914(2n+1)` only, offset 0, from Sahin and Tan's divisibility paper; the six-pointed-star hole count is a comment, not the definition. Holes of the `n`-th size, not tiles: this tree's slice census is `A299916(n+1)` | https://oeis.org/A299916 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Nakajima and Watanabe 2026 | Topology of slices through the Sierpinski tetrahedron, arXiv:2603.06004 | https://arxiv.org/abs/2603.06004 |
| Nakajima and Watanabe 2026 journal | same paper, Chaos, Solitons & Fractals 209, 118353 | https://doi.org/10.1016/j.chaos.2026.118353 |
| Diaconis and Fulman 2009 | Carries, shuffling, and an amazing matrix, Amer. Math. Monthly 116, 788-803; read in the arXiv version: the introduction recalls Holte's carries chain on `{0, ..., n-1}` for `n` summands with uniform digits in base `b`, its transition matrix (H1) `P(i, j) = b^(-n) sum_r (-1)^r C(n+1, r) C(n - 1 - i + (j + 1 - r) b, n)` and the base-2 form (H2); the eigenvalue and eigenvector properties were not read here | https://arxiv.org/abs/0806.3583 |
| Holte 1997 | Carries, combinatorics, and an amazing matrix, Amer. Math. Monthly 104, 138-149; carried from Diaconis and Fulman 2009, not read at its source | https://doi.org/10.2307/2974981 |
| Glaisher 1899 | On the residue of a binomial-theorem coefficient with respect to a prime modulus, Quart. J. Pure Appl. Math. 30, 150-156 | https://babel.hathitrust.org/cgi/pt?id=hvd.32044102924578 |
| slice-recurrence-order | Menger Diagonal Slices: A Recurrence of Order ceil(D/2); the shelf paper of the order law, the product formula over 3-adic angle towers and the pinning of rho_dim | papers/slice-sign-law.md |
| slice-sign-even-half | The Even Half of the Slice Sign Law; the shelf paper of the even half at bases 3 and 5, the certificate machines, the transient constant ln(R)/4, the tent rank law and the layer-2 window law | papers/slice-sign-law.md |

## REFERENCE PAGES

| ref | title | url |
|---|---|---|
| n-flake article | Wikipedia: octahedron flake of dimension log(6)/log(2), and the Cantor cube projecting to a hexaflake | https://en.wikipedia.org/wiki/N-flake |
| Kummer's theorem | Wikipedia: the 2-adic valuation of C(i+j,i) counts base-2 carries, so C(i+j,i) is odd iff i AND j = 0 | https://en.wikipedia.org/wiki/Kummer%27s_theorem |

## PRIOR ART ON THE BASE-3 SLICE

| ref | title | url |
|---|---|---|
| Hocking, Bridges 2023 | Three-Dimensional Diagonal Cross-Sections of Four-Dimensional Menger Sponges - generalises the cut along dimension `n` and along a hole-iness parameter `k`, with "the base three expansion" fixed throughout | https://archive.bridgesmathart.org/2023/bridges2023-291.pdf |

## DATASETS

| ref | title | url |
|---|---|---|
| b001316 | OEIS b-file for A001316, 50001 terms | https://oeis.org/A001316/b001316.txt |
| b047999 | OEIS b-file for A047999, 10585 terms (rows 0..144) | https://oeis.org/A047999/b047999.txt |
