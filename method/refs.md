# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A000290 | The squares, n^2; the odd-side fills of the low-corner design in the plane | https://oeis.org/A000290 |
| A000384 | Hexagonal numbers, n(2n-1); the odd-side fills of the flat tree | https://oeis.org/A000384 |
| A000567 | Octagonal numbers, n(3n-2); the odd-side fills of the flat carpet | https://oeis.org/A000567 |
| A001844 | Centered square numbers, 2n(n+1)+1; the odd-side fills of the flat void | https://oeis.org/A001844 |
| A016754 | Odd squares, (2n+1)^2; the odd-side fills of the solid square, also the centered octagonal numbers | https://oeis.org/A016754 |
| A000578 | The cubes: a(n) = n^3 | https://oeis.org/A000578 |
| A000616 | a(-1)=1 by convention; for n >= 0, a(n) = number of irreducible Boolean functions of n variables | https://oeis.org/A000616 |
| A001024 | Powers of 15: a(n) = 15^n | https://oeis.org/A001024 |
| A003215 | Hex (or centered hexagonal) numbers: 3*n*(n+1)+1; also the lattice lines down the space diagonal of a cube of side n+1, hence the sponge's diagonal shadow at n = 3^level - 1 | https://oeis.org/A003215 |
| A009971 | Powers of 27 | https://oeis.org/A009971 |
| A016755 | Odd cubes: a(n) = (2*n + 1)^3 | https://oeis.org/A016755 |
| A129824 | a(n) = Product_{k=0..n} (1 + binomial(n,k)) | https://oeis.org/A129824 |
| A395241 | a(n) = n^2*(4*n + 3) - this tree's own submission, void subcubes of the odd sponge tile | https://oeis.org/A395241 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Odlyzko and te Riele 1985 | Disproof of the Mertens conjecture, J. reine angew. Math. 357, 138-160 - the exhaustive-verification lesson in [method](notes/method.md) | https://doi.org/10.1515/crll.1985.357.138 |
| Griffin, Ono, Rolen and Zagier 2019 | Jensen polynomials for the Riemann zeta function and other sequences, PNAS 116(23), 11103-11110 - cited by [method](notes/method.md) as a false friend, never as support | https://doi.org/10.1073/pnas.1902572116 |
| Collatz 1942 | Einschliessungssatz fuer die charakteristischen Zahlen von Matrizen, Mathematische Zeitschrift 48(1), 221-226; the first half of the Collatz-Wielandt name, the min-max characterization of the Perron root and not the 3n + 1 map. The original is paywalled at Springer; the statement is read through the Wikipedia Perron-Frobenius article and the metadata verified in the Crossref record of the DOI, author Collatz, volume 48, pages 221-226, year 1942 | https://doi.org/10.1007/BF01180013 |
| Wielandt 1950 | Unzerlegbare, nicht negative Matrizen, Mathematische Zeitschrift 52(1), 642-648; the second half of the Collatz-Wielandt name, the irreducible case behind the brackets this tree certifies with. The original is paywalled at Springer; the statement is read through the Wikipedia Perron-Frobenius article and the metadata verified in the Crossref record of the DOI, author Wielandt, volume 52, pages 642-648, year 1950 | https://doi.org/10.1007/BF02230720 |
| O'Donnell 2014 | Analysis of Boolean Functions, Cambridge University Press, read at source in arXiv:2105.10386v1: Definition 1.2 and Proposition 1.8 the Fourier coefficient, Definition 1.19 the weight at degree k, Definitions 2.42 and 2.46, Proposition 2.47, Fact 2.48 and Theorem 2.49 noise stability and the noise operator, Exercise 1.8 odd functions, Exercise 5.28 the Kravchuk polynomials, section 8.4 the p-biased mean E[f^(p)] = T_mu f(1, ..., 1) and the Margulis-Russo formula (8.9) | https://arxiv.org/abs/2105.10386 |

## REFERENCE PAGES

| ref | title | url |
|---|---|---|
| Perron-Frobenius theorem | Wikipedia: the Collatz-Wielandt formula, the min-max over positive vectors behind the even-half certificates of the `slice-sign-even-half` shelf paper | https://en.wikipedia.org/wiki/Perron%E2%80%93Frobenius_theorem |

## DATASETS

| ref | title | url |
|---|---|---|
| OEIS dump | `stripped.gz`, the full-sequence dump the novelty searches on [sequences](../shared/sequences.md) run against; 398556 lines | https://oeis.org/stripped.gz |
