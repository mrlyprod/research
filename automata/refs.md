# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A000370 | Number of NPN-equivalence classes of Boolean functions of n or fewer variables | https://oeis.org/A000370 |
| A000616 | a(-1)=1 by convention; for n >= 0, a(n) = number of irreducible Boolean functions of n variables | https://oeis.org/A000616 |
| A001316 | Gould's sequence: number of odd entries in row n of Pascal's triangle; a(n) = 2^A000120(n) | https://oeis.org/A001316 |
| A047999 | Sierpinski's triangle (or gasket): Pascal's triangle mod 2 | https://oeis.org/A047999 |
| A071053 | Number of ON cells at generation n of elementary rule 150 started from a single ON cell | https://oeis.org/A071053 |
| A087206 | a(n) = 2*a(n-1) + 4*a(n-2); with a(0)=1, a(1)=4; equals 2^n Fibonacci(n+2), the ON cells of rule 150 over its first 2^n rows | https://oeis.org/A087206 |
| A160239 | Number of ON cells at generation n of Fredkin's Replicator, the odd-rule automaton on the eight-cell Moore neighbourhood | https://oeis.org/A160239 |
| A246035 | Number of odd terms in f^n where f = (1/x+1+x)*(1/y+1+y), the nine-cell Moore odd-rule automaton; equals A071053(n)^2, so a separable planar kernel is the tensor square of rule 150 | https://oeis.org/A246035 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Wolfram 1983 | Statistical mechanics of cellular automata, Rev. Mod. Phys. 55(3), 601-644; the table partitioning the 256 elementary rules into 88 classes under reflection and state complementation | https://doi.org/10.1103/RevModPhys.55.601 |
| Li and Packard 1990 | The structure of the elementary cellular automata rule space, Complex Systems 4(3), 281-297; the derivation of the 88 count. The only full text reached is a scan, so its statements here are carried from Schaller and Svozil's restatement | https://www.complex-systems.com/abstracts/v04_i03_a03/ |
| Schaller and Svozil 2025 | Irreducible rules and equivalence classes of one-dimensional cellular automata, arXiv:2512.08117; state permutations and lattice isometries generate the symmetry group, extended by neighbourhood scaling, and free permutation of the three neighbourhood coordinates is not among them | https://arxiv.org/abs/2512.08117 |
| Harrison 1963 | The number of transitivity sets of Boolean functions, J. Soc. Indust. Appl. Math. 11(3), 806-828, doi:10.1137/0111059, the JSTOR page beside answering automated fetches with a bot wall; the cycle-index formula OEIS credits for A000616, the 22 classes of three-variable Boolean functions under the order-48 hyperoctahedral group | https://www.jstor.org/stable/2946322 |
| Harrison 1963 IEEE | The number of equivalence classes of Boolean functions under groups containing negation, IEEE Trans. Electron. Comput. 12, 559-561; the source OEIS credits for A000370, the 14 NPN classes | https://doi.org/10.1109/PGEC.1963.263656 |
| Langton 1990 | Computation at the edge of chaos: phase transitions and emergent computation, Physica D 42(1-3), 12-37; lambda is the fraction of rule-table entries mapping to a non-quiescent state. Every copy reached is a scan, so the definition is carried from restatements | https://doi.org/10.1016/0167-2789(90)90064-V |
| Mitchell, Hraber and Crutchfield 1993 | Revisiting the edge of chaos: evolving cellular automata to perform computations, Complex Systems 7(2), 89-130; the abstract states the earlier interpretation of lambda's role is not correct | https://www.complex-systems.com/abstracts/v07_i02_a01/ |
| Martin, Odlyzko and Wolfram 1984 | Algebraic properties of cellular automata, Comm. Math. Phys. 93(2), 219-258; rule 90 is the Laurent polynomial x + x^-1 over GF(2) and rule 150 is x + 1 + x^-1. It observes rule 90's single-cell pattern has fractal dimension log2 3; it does not prove it and never says Pascal | https://doi.org/10.1007/BF01223745 |
| Willson 1987 | Computing fractal dimensions for additive cellular automata, Physica D 24(1-3), 190-206; the abstract states the largest eigenvalue of an integer matrix built from the automaton gives the fractal dimension. Paywalled; the body was not read, so no rule-150 number is drawn from it | https://doi.org/10.1016/0167-2789(87)90074-1 |
| Cook 2004 | Universality in elementary cellular automata, Complex Systems 15(1), 1-40; rule 110 is Turing universal | https://www.complex-systems.com/abstracts/v15_i01_a01/ |
| Culik and Yu 1988 | Undecidability of CA classification schemes, Complex Systems 2(2), 177-190; the Wolfram classes are made formal and membership is undecidable | https://www.complex-systems.com/abstracts/v02_i02_a02/ |
| Ollinger 2009 | Intrinsically universal cellular automata, EPTCS 1, 199-204; intrinsic universality is undecidable though recursively enumerable, rule 110's Turing universality does not give intrinsic universality, and whether rule 110 is intrinsically universal is open | https://doi.org/10.4204/EPTCS.1.19 |
| Delorme, Mazoyer, Ollinger and Theyssier 2011a | Bulking I: an abstract theory of bulking, Theoret. Comput. Sci. 412(30), 3866-3880 | https://doi.org/10.1016/j.tcs.2011.02.023 |
| Delorme, Mazoyer, Ollinger and Theyssier 2011b | Bulking II: classifications of cellular automata, Theoret. Comput. Sci. 412(30), 3881-3905; bulking is a quasi-order comparing space-time diagrams up to rescaling, not a tensor operation on neighbourhood masks | https://doi.org/10.1016/j.tcs.2011.02.024 |
| Amoroso and Patt 1972 | Decision procedures for surjectivity and injectivity of parallel maps for tessellation structures, J. Comput. System Sci. 6(5), 448-464; both properties are decidable in one dimension | https://doi.org/10.1016/S0022-0000(72)80013-8 |
| Sutner 1991 | De Bruijn graphs and linear cellular automata, Complex Systems 5(1), 19-30; quadratic-time tests for reversibility and surjectivity | https://www.complex-systems.com/abstracts/v05_i01_a03/ |
| Sloane 2015 | On the number of ON cells in cellular automata, arXiv:1503.01168; for an odd-rule automaton with neighbourhood F the state at generation n is F^n and the ON count is the number of its nonzero terms. Fredkin's Replicator is the eight-cell Moore neighbourhood, (1/x+1+x)(1/y+1+y) - 1, which is not separable | https://arxiv.org/abs/1503.01168 |
| Evans 2001 | Larger than Life: Digital Creatures in a Family of Two-Dimensional Cellular Automata, DMTCS Proceedings AA; Life rescaled to box neighbourhoods of range r with rescaled gliders and spaceships | https://doi.org/10.46298/dmtcs.2288 |
| Pivato 2007 | RealLife: the continuum limit of Larger Than Life cellular automata, Theoretical Computer Science 372; Larger than Life at large radius has life forms virtually identical after rescaling | https://arxiv.org/abs/math/0503504 |
| Rafler 2011 | Generalization of Conway's Game of Life to a continuous domain - SmoothLife; ring neighbourhood and disk cell, typically outer radius three times the inner | https://arxiv.org/abs/1111.1567 |
| Chan 2019 | Lenia - Biology of Artificial Life, Complex Systems 28(3); scale invariance through large neighbourhoods with exponential, polynomial and rectangular kernel cores | https://arxiv.org/abs/1812.05433 |
| Chladni 1787 | Entdeckungen ueber die Theorie des Klanges; the sand figures of vibrating plates, digitised by the Deutsches Textarchiv | https://www.deutschestextarchiv.de/book/show/chladni_klang_1787 |
| Davies, Gladwell, Leydold and Stadler 2001 | Discrete nodal domain theorems, Linear Algebra and its Applications 336, 51-60; an eigenvector of index k and multiplicity r on a connected graph has at most k + r - 1 strong and k weak nodal domains | https://arxiv.org/abs/math/0009120 |
| Morrison 2004 | Cosine Products, Fourier Transforms, and Random Sums; the infinite cosine product whose Fourier transform is the Cantor measure | https://arxiv.org/abs/math/0411380 |
| Tsujii 2012 | On the Fourier transforms of self-similar measures; the transform of a self-similar measure as an infinite product | https://arxiv.org/abs/1212.1553 |
| Liu, Wang and Zheng 2026 | Spectrality of product-form self-similar measures and tiles; the digit-set transform as the product of the mask polynomial at scaled frequencies | https://arxiv.org/abs/2603.21031 |
| Liu-Wang-Zheng 2026 | Spectrality of product-form self-similar measures and tiles, unrefereed arXiv preprint, read the abstract, Section 1.1 with Theorem 1.2 quoting An-Wang Theorem 1.7, Theorems 1.3 and 2.3, the proof of Theorem 3.5 and the references | https://arxiv.org/abs/2603.21031 |

## REFERENCE PAGES

| ref | title | url |
|---|---|---|
| NKS note, rule 150 | A New Kind of Science note for page 885: there are 2^m Fibonacci[m+2] black cells up to step 2^m, so the fractal dimension is Log[2, 1 + Sqrt[5]] | https://www.wolframscience.com/nks/notes-3-2--rule-150/ |
| NKS note, surjectivity | A New Kind of Science note for page 959: 30 surjective elementary rules, and in two dimensions such properties are in general undecidable | https://www.wolframscience.com/nks/notes-6-7--surjectivity-and-injectivity-of-cellular-automaton-maps/ |
| NKS page 436 | A New Kind of Science: of the 256 elementary rules only six are reversible; the page does not name them | https://www.wolframscience.com/nks/p436--the-notion-of-reversibility/ |
