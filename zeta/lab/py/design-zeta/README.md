# design-zeta

- The Dirichlet series of a digit design, `zeta_F(s) = sum_(n in S_F) n^(-s)` over the base-`base` integers whose digits all lie in a digit set `F`, continued to the whole plane and read for its zeros.
- Abscissa `alpha = log_base fill`, `fill = card F`; the object is Burnol 2026's `K(s)`, and the continuation is his Proposition 4.1 in a peeled form.
- Peeling: with `E_j(w)` the sum over the exactly-`j`-digit elements and `G_P(w) = sum_(j >= P) E_j(w)`, the digit map gives `(1 - fill base^(-w)) G_P(w) = E_P(w) + sum_(l >= 1) binom(-w, l) base^(-w-l) gamma_l G_P(w+l)` with `gamma_l = sum_(a in F) a^l`, and `zeta_F(w) = D_(P-1)(w) + G_P(w)` for the finite Dirichlet polynomial `D_(P-1)` over `S_F` below `base^(P-1)`.
- The `l`-series has ratio `max F/base^P`, so the peel depth buys the convergence; the small quantity `G_P` is carried directly and never as a difference of two large ones.
- Truncation is proved and propagated: `abs(G_P(w)) <= k_1 (fill base^(-Re w))^(P-1)/(1 - fill base^(-Re w))` bounds the base of the ladder, a binomial majorant with a geometric remainder bounds the `l`-cut, and the same recursion that carries the value carries the bound.
- Every printed value carries its bound, the census included: the ladder raises rather than returns when the requested tolerance is not reached, and the census prints the largest bound it met anywhere on its contour beside the largest phase step.
- Residues come off the same ladder: at `s_(m,j) = alpha - m + 2 pi i j/log base` the factor `1 - fill base^(-s)` has derivative `log base`, so the residue is the ladder numerator over `log base`.

## THE CENSUS

- The census reads the Lyndon cofactor `Z(s) = zeta_F(s)(1 - fill base^(-s))` and not `zeta_F`. The level-`m` denominator of the ladder is `1 - fill base^(-(s+m))`, so `Z` is analytic on `Re s > alpha - 1`, every zero of `zeta_F` right of that line is a zero of one analytic function, and no pole-free strip is needed.
- One strip, `alpha - 0.92 < Re s < alpha + 3.02`, `0.02 < Im s < 60`, split at `Re s = alpha` exactly. The two halves are therefore right of the abscissa and left of it with nothing unscanned between them; the stated strip is the computed strip and its edges are printed with the count.
- Below `Re s = alpha - 0.92` the ladder cannot reach tolerance against the cofactor's own pole line `Re s = alpha - 1`, so that is the left edge and the left-hand count is a count on the strip actually scanned.
- The phase is accumulated along an adaptively refined contour that halves any segment whose principal-branch step exceeds one radian; a box of winding one is bisected by quadrant and polished by Muller.
- The count is Verified and not Proved. The largest surviving phase step is printed, but nothing here bounds `zeta_F'/zeta_F` on the contour, so the argument principle is resolved and not certified and a zero pair closer than the surviving spacing would stay invisible.

## RUN

- `uv run python research/zeta/lab/py/design-zeta/design_zeta.py`
- `uv run python research/zeta/lab/py/design-zeta/design_zeta.py --full` runs the census to `Im s = 60` and adds the base-10 columns; `--y10=40` caps the base-10 height alone.
- Prints only, writes nothing.

## CONTROLS

- Base 2 on the full digit set is the whole integer line, so `zeta_F` is `zeta` there. The engine matches mpmath's `zeta` to the printed bound at `s = 2`, at the first zero, at `0.3 + 40i`, at `-1 + 2i` and at `0.5 + 100i`, gives residue exactly `1` at `alpha = 1`, and gives residue zero at `alpha + 2 pi i/log 2`, where the design family has a pole and the integers do not.
- That column is censused on BOTH sides of its abscissa `alpha = 1`, and each side names the object it counts. Zeros of `zeta` right of the abscissa, `alpha + 0.02 < Re s < alpha + 3.02`: winding zero, which is what the Euler product forbids, now computed and not quoted. Zeros of `zeta` left of it, `alpha - 0.98 < Re s < alpha - 0.02`: exactly the first thirteen below `Im s = 60`, six of them below `Im s = 40`.
- The cofactor is a different object on that column and is counted separately. `1 - 2 base^(-s)` has teeth exactly ON `Re s = alpha`, `floor(T log base/2 pi)` of them below height `T`, so the one-strip cofactor count reads `19 = 13 + 6` at `Im s < 60` and never puts a zero of `zeta` right of `Re s = 1`.
- The residues at `base = 3` reproduce the certified interval enclosures of `dimensions/lab/py/burnol-residue` digit for digit on both `F = {0,1}` and `F = {0,2}`, by a route that shares no code with them.
