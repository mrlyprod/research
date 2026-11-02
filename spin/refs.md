# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A001481 | Numbers that are the sum of 2 squares; the ring radii squared of a spun square-lattice picture | https://oeis.org/A001481 |
| A003136 | Loeschian numbers, the norms x^2 + xy + y^2; the ring radii squared of a spun hexagonal picture | https://oeis.org/A003136 |
| A004016 | Theta series of the planar hexagonal lattice; the weight of each hexagonal ring | https://oeis.org/A004016 |
| A064533 | Decimal expansion of the Landau-Ramanujan constant K = 0.7642236..., the constant in `#{k <= X : k = a^2 + b^2} ~ K X / sqrt(ln X)` | https://oeis.org/A064533 |
| A004018 | Theta series of the square lattice, r2(n); the weight of each square ring | https://oeis.org/A004018 |
| A003215 | Hex (or centered hexagonal) numbers: 3*n*(n+1)+1; also the lattice lines down the space diagonal of a cube of side n+1, hence the sponge's diagonal shadow at n = 3^level - 1 | https://oeis.org/A003215 |
| A220978 | a(n) = 3^(2n+1) - 3^(n+1) + 1, the left Aurifeuillian factor of 3^(6n+3) + 1; the same diagonal shadow indexed by the level | https://oeis.org/A220978 |
| A301429 | Decimal expansion of the Loeschian analog of the Landau-Ramanujan constant, K_3 = 0.6389094054..., the constant in `L(x) ~ K_3 x / sqrt(ln x)` | https://oeis.org/A301429 |
| A400227 | The Gaussian Farey `new(n)`: distances of `(1/n) Z^2` in the disc of radius `sqrt 2` new at scale `n`, the square-lattice twin of the hexagonal count | https://oeis.org/A400227 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Landau 1908 | Über die Einteilung der positiven ganzen Zahlen in vier Klassen nach der Mindestzahl der zu ihrer additiven Zusammensetzung erforderlichen Quadrate, Archiv der Mathematik und Physik (3) 13, 305-312; read at source in the scan: equation (2) proves that `lim B(x) sqrt(ln x) / x = K` exists and is positive, by a contour around the cut of the branch point `s = 1` of the Dirichlet series of the indicator of A001481, and the introduction says it asks for no estimate sharper than the leading term | https://archive.org/details/archivdermathem37unkngoog |
| Landau 1909 | Handbuch der Lehre von der Verteilung der Primzahlen, vol. 2, section 176, pp. 643-645; read at source in the scan: restates `B(x) ~ b x / sqrt(ln x)` for sums of two squares and sends the proof and the value of `b` to sections 177-183 (headings read, proof not); no error term is stated in section 176 | https://archive.org/details/handbuchderlehre02landuoft |
| Moree 2011 | Counting numbers in multiplicative sets: Landau versus Ramanujan, arXiv:1110.0708; read at the LaTeX source: Theorem 1 (Wirsing 1961, Math. Ann. 143, 75-102, unread, in the form of Finch, Martin and Sebah) gives `S(x) ~ C_0(S) x (ln x)^(delta-1)` for a multiplicative set whose primes satisfy `pi_S(x) ~ delta x / ln x`, with `C_0(S) = Gamma(delta)^-1 lim prod_(p < P) (sum_j chi_S(p^j)/p^j) (1 - 1/p)^delta`; the section on races states `S_H(x) ~ C_0(S_H) x / sqrt(ln x)` for the Loeschian numbers `S_H`, `C_0(S_H) = 0.639...` | https://arxiv.org/abs/1110.0708 |
| Moree and Osburn 2006 | Two-dimensional lattices with few distances, Enseign. Math. 52, 361-380, arXiv:math/0604163; read at the LaTeX source: Bernays' theorem `B_f(x) = C(f) x / sqrt(ln x) + O(x / (ln x)^(1/2 + delta))`, `delta < min(1/h, 1/4)`, for a primitive positive definite binary form `f` of negative discriminant, `B_f` counting the integers `f` represents; Table 1 gives `C(X^2 + 3 Y^2) = 0.6389094054453438822549426747`, the form `X^2 + 3 Y^2` representing exactly the Loeschian numbers | https://arxiv.org/abs/math/0604163 |
| Bernays 1912 | Über die Darstellung von positiven, ganzen Zahlen durch die primitiven, binären quadratischen Formen einer nicht-quadratischen Diskriminante, Dissertation, Göttingen; unread, cited through Moree and Osburn 2006 | http://www.math.uni-bielefeld.de/~rehmann/DML/ |
| Glass 1969 | Moire effect from random dots, Nature 223, 578-580 | https://doi.org/10.1038/223578a0 |
| Hardy 1915 | On the expression of a number as the sum of two squares, Quarterly Journal of Mathematics 46, 263-283; the Bessel series for the circle-problem error | https://zbmath.org/?q=an:45.1253.01 |
| Cherny, Anitas, Osipov and Kuklin 2011 | Deterministic fractals: extracting additional information from small-angle scattering data, Phys. Rev. E 84, 036203 | https://doi.org/10.1103/PhysRevE.84.036203 |
| Mattila 1987 | Spherical averages of Fourier transforms of measures with finite energy; dimension of intersections and distance sets, Mathematika 34, 207-228 | https://doi.org/10.1112/S0025579300013462 |
| Falconer, Fraser and Jin 2015 | Projections of self-similar and related fractals: a survey of recent developments, Fractal Geometry and Stochastics V | https://doi.org/10.1007/978-3-319-18660-3_4 |
| Patterson 1944 | Ambiguities in the X-ray analysis of crystal structures, Phys. Rev. 65, 195-201 | https://doi.org/10.1103/PhysRev.65.195 |
| Rosenblatt and Seymour 1982 | The structure of homometric sets, SIAM J. Algebraic Discrete Methods 3, 343-350 | https://doi.org/10.1137/0603035 |
| Kac 1966 | Can one hear the shape of a drum?, Amer. Math. Monthly 73, no. 4, 1-23 | https://doi.org/10.2307/2313748 |
| Gordon, Webb and Wolpert 1992 | One cannot hear the shape of a drum, Bull. Amer. Math. Soc. 27, 134-138 | https://doi.org/10.1090/S0273-0979-1992-00289-6 |
| Hu 1962 | Visual pattern recognition by moment invariants, IRE Trans. Inform. Theory 8, 179-187 | https://doi.org/10.1109/TIT.1962.1057692 |
| Teague 1980 | Image analysis via the general theory of moments, J. Opt. Soc. Amer. 70, 920-930 | https://doi.org/10.1364/JOSA.70.000920 |

## REFERENCE PAGES

| ref | title | url |
|---|---|---|
| Landau-Ramanujan constant | Wikipedia: the count of sums of two squares below X is asymptotic to K X / sqrt(ln X), Landau 1908 | https://en.wikipedia.org/wiki/Landau%E2%80%93Ramanujan_constant |
