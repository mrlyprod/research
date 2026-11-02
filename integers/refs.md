# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A000070 | Partition sums: a(n) = Sum_{k=0..n} p(k); the qualifying-signature count of the divisor-avatar material, now the `divisor-avatars` lane | https://oeis.org/A000070 |
| A000567 | Octagonal numbers, n(3n-2); the odd-side fills of the flat carpet | https://oeis.org/A000567 |
| A043635 | Numbers whose base-9 representation has exactly 6 runs; its 30 listed terms lie wholly inside the integer census's miss set | https://oeis.org/A043635 |
| A049537 | Values of k for which A075059(k) = A003418(k) + 1 is prime; the first record to carry a 4-term window of the integer census's miss set, at offset 417 | https://oeis.org/A049537 |
| A100290 | Numbers divisible by the smallest number with the same binary weight, so A038573(a(n)) divides a(n); shares a 12-term window with the integer census's champion set, then gives 21 where the census gives 20 | https://oeis.org/A100290 |
| A103532 | Number of divisors of 240^n; the odd bisection of A011934 | https://oeis.org/A103532 |
| A112820 | Numbers k with lcm(1,2,...,k)/17 equal to the denominator of the k-th harmonic number; carries 20 consecutive integers of the census's miss set | https://oeis.org/A112820 |
| A118471 | a(0)=1, a(n) = a(n-1)*(n+1) if n is in the sequence and a(n-1)+1 otherwise; carries 20 consecutive integers of the census's miss set | https://oeis.org/A118471 |
| A129824 | a(n) = Product_{k=0..n} (1 + binomial(n,k)) | https://oeis.org/A129824 |
| A209631 | Square array of the exponential transform applied n times to the identity function; the only record carrying the census's written-per-decade run `9, 90, 859`, and it continues 6689 where the census gives 5452 | https://oeis.org/A209631 |
| A336231 | Integers whose binary expansion has an even number of 0's between any two consecutive 1's; shares a 12-term window with the census's champion set, then gives 19 where the census gives 20 | https://oeis.org/A336231 |
| A361796 | Primes preceded by two consecutive products of four distinct primes; the longest record lying wholly inside the census's miss set, 41 terms | https://oeis.org/A361796 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Hunter 1957 | The minimum discriminants of quintic fields | https://doi.org/10.1017/S2040618500033463 |
| Pohst 1982 | On the computation of number fields of small discriminants including the minimum discriminants of sixth degree fields | https://doi.org/10.1016/0022-314X(82)90061-0 |
| Jones and Roberts 2014 | A database of number fields | https://arxiv.org/abs/1404.0266 |
| divisor avatars | Divisor Avatars: Which Parity Designs Count the Divisors of a Power | claims/divisor-avatars.md |
| Tan 2018 | Degree bound of Polya Positivstellensatz, arXiv 1802.02752; the univariate statement, P(x) > 0 for x >= 0 gives (1 + x)^m P(x) with positive coefficients for m large, and the Powers-Reznick exponent m > (d^2 - d) L(P) / (2 lambda(P)) - d with L(P) = max_j abs(a_j) / C(d, j) and lambda(P) = inf over x >= 0 of P(x) / (1 + x)^d, read at source | https://arxiv.org/abs/1802.02752 |
| Powers and Reznick 2001 | A new bound for Polya's Theorem with applications to polynomials positive on polyhedra, Journal of Pure and Applied Algebra 164, 221-229; the explicit exponent in Polya's theorem. The original is paywalled at Elsevier; the bound is read through Tan 2018 and the metadata verified in the Crossref record of the DOI, authors Powers and Reznick, volume 164, pages 221-229, year 2001 | https://doi.org/10.1016/S0022-4049(00)00155-9 |

## REFERENCE PAGES

| ref | title | url |
|---|---|---|
| A001792 | a(n) = (n+2)*2^(n-1); the level-1 cell count of the Menger sponge in every dimension | https://oeis.org/A001792 |
| A005408 | The odd numbers; the odd-side fills of the dimension-one solid | https://oeis.org/A005408 |
| Polya 1928 | Ueber positive Darstellung von Polynomen, Vierteljahrsschrift der Naturforschenden Gesellschaft in Zuerich 73, 141-145; a form positive on the simplex has a multiple (x_1 + ... + x_n)^m with nonnegative coefficients. The original is unread; the statement is read in the Wikipedia article Positive polynomial, which cites it, and in the univariate restatement of Tan 2018 | https://en.wikipedia.org/wiki/Positive_polynomial |
| positive polynomial article | Wikipedia: Polya's theorem, a homogeneous polynomial positive on the nonnegative orthant away from the origin has (x_1 + ... + x_n)^m p with nonnegative coefficients for some m, cited to Polya 1928 | https://en.wikipedia.org/wiki/Positive_polynomial |

## DATASETS

| ref | title | url |
|---|---|---|
| OEIS dump | `stripped.gz`, the full-sequence dump the novelty searches on [sequences](../shared/sequences.md) run against; 398556 lines | https://oeis.org/stripped.gz |
| LMFDB 2.0.3.1 | LMFDB number field 2.0.3.1 | https://www.lmfdb.org/NumberField/2.0.3.1 |
| LMFDB 2.2.5.1 | LMFDB number field 2.2.5.1 | https://www.lmfdb.org/NumberField/2.2.5.1 |
| LMFDB 3.1.23.1 | LMFDB number field 3.1.23.1 | https://www.lmfdb.org/NumberField/3.1.23.1 |
| LMFDB 3.3.49.1 | LMFDB number field 3.3.49.1 | https://www.lmfdb.org/NumberField/3.3.49.1 |
| LMFDB 4.0.117.1 | LMFDB number field 4.0.117.1 | https://www.lmfdb.org/NumberField/4.0.117.1 |
| LMFDB 4.0.1225.1 | LMFDB number field 4.0.1225.1 | https://www.lmfdb.org/NumberField/4.0.1225.1 |
| LMFDB 4.4.725.1 | LMFDB number field 4.4.725.1 | https://www.lmfdb.org/NumberField/4.4.725.1 |
| LMFDB 5.1.4429.1 | LMFDB number field 5.1.4429.1 | https://www.lmfdb.org/NumberField/5.1.4429.1 |
| LMFDB 5.3.4511.1 | LMFDB number field 5.3.4511.1 | https://www.lmfdb.org/NumberField/5.3.4511.1 |
