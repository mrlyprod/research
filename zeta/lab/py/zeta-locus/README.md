# zeta-locus

- The locus of the zeros of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)`, with `F` the digits of the code and `S_F` the integers they write: where the zeros sit across bases and digit sets, and which law the real part obeys.
- The ladder, the contour engine and the truncation bounds are imported from `../design-zeta` and never copied; this study adds the cofactor, the Laurent data at the poles, the derived shadow law and the falsification sweep.

## THE COFACTOR

- `Z(s) = zeta_F(s) (1 - fill base^(-s))` is the Lyndon cofactor: `1/(1 - fill base^(-s))` is the zeta of the free monoid on `F` with norm `base^(length)`, is zero free, and carries the whole pole lattice.
- `Z` is analytic on `Re s > alpha - 1`, because `1 - fill base^(-s)` cancels exactly the `m = 0` line of the ladder's poles and no other; its own poles are the `s_(m,j)` with `m >= 1` at which `zeta_F` has a nonvanishing residue, the nearest to the strip being `Re s = alpha - 1` with residue `-s_(1,j) gamma_1 r_j/fill`, and on a full digit set there are none at all, `Z` being entire.
- `Z(s) -> a_min^(-s)` as `Re s -> +infinity` with `a_min` the least nonzero digit, so every zero of `zeta_F` right of `alpha - 1` is a zero of one analytic function.
- The transfer runs one way without exception and the other way with one: a zero of `zeta_F` is always a zero of `Z`; a zero of `Z` is a zero of `zeta_F` except at a pole `s_(0,j)` where the residue vanishes, and there `Z` vanishes while `zeta_F` is regular.
- The one-level digit recursion gives it in closed form: `Z(s) = E_1(s) + sum_(l >= 1) binom(-s, l) base^(-s-l) gamma_l zeta_F(s+l)` with `E_1(s) = sum_(a in F, a != 0) a^(-s)` and `gamma_l = sum_(a in F) a^l`; the peeled ladder of `../design-zeta` evaluates it as `(1 - fill base^(-s)) D_(P-1)(s)` plus the ladder numerator, so the printed VALUE of `Z` carries the same propagated bound. The Laurent coefficients do not: `Z_1` and `Z_2` are central differences at step `1e-5` carrying about `1e-10` of truncation, far above the ladder's own bound, enough for a prediction compared at `1e-3` and not a certificate.

## THE DERIVED SHADOW

- At a pole `s_(0,j) = alpha + 2 pi i j/log base` one has `fill base^(-s_(0,j)) = 1` exactly for every `j`, so `1 - fill base^(-s) = 1 - base^(-u)` in `u = s - s_(0,j)` with no `j` dependence at all.
- With `Z(s) = Z_0 + Z_1 u + Z_2 u^2 + ...`, the identity `1/(1 - e^(-x)) = 1/x + 1/2 + x/12 - x^3/720 + ...` at `x = u log base` gives the Laurent expansion of `zeta_F` at the pole term by term.
- Residue `r_j = Z_0/(log base)`, regular part `R_j = Z_1/(log base) + Z_0/2`, its derivative `R'_j = Z_2/(log base) + Z_1/2 + Z_0 (log base)/12`. The residues are Burnol's `lambda_(0,j)`, certified by `../burnol-residue`; on the full digit set `r_0 = 1` and `R_0` is Euler-Mascheroni, which is what `zeta` has at `s = 1`.
- A zero near the pole solves `u (R_j + R'_j u + ...) = -r_j`. First order `u_1 = -r_j/R_j`; second order the root of `R'_j u^2 + R_j u + r_j = 0` nearest `u_1`. Nothing is fitted: the prediction is built from the residue and the regular part alone and is then compared with the polished zero.
- The pole index `j` is the only label; the disc count, not the size of the prediction, decides whether a pole carries a zero at all.

## THE SWEEP

- Every printed zero carries `Re s`, `Im s`, `Im s log base/(2 pi)`, its fractional part, `alpha`, `fill/base`, the pole index and the first-order prediction, so a curve law, a comb law and a family law each have a column to die in.
- Equal-`alpha` pairs at different bases are the falsifier: base 4 `{0,1}`, base 9 `{0,1,2}` and base 16 `{0,1,2,3}` all have `alpha = 1/2`, and base 4 `{0,1}` and base 16 `{0,1,2,3}` also share `fill/base = 1/2`.
- Equal-`alpha` pairs at one base are the cheaper falsifier: all fill 2 sets at base 4 share `alpha` and `fill/base` and differ only in `F`.
- Scaling classes are quotiented first: `zeta_(aF)(s) = a^(-s) zeta_F(s)` has the same zeros, so one representative per class is censused.
- The base 2 full digit set is the line control: there `zeta_F = zeta` and the locus is the critical line `Re s = 1/2 = alpha/2`.

## THE CENSUS

- The strip is `alpha - 0.92 < Re s < alpha + 3.02`, one box, no pole inside and no blind sliver, the winding of `Z` on it counted by the argument principle with the largest surviving phase step printed as the certificate.
- Each pole gets its own contour, a circle of radius `0.45` about `s_(0,j)`, an assignment radius fixed by the discs not overlapping and not by the tooth law, so every count is conditional on it; that ONE radius serves both the count and the tooth: the winding gives the number of zeros of `Z` in the disc, the residue-null centre is subtracted to give the number of zeros of `zeta_F`, and exactly that many are then located by a polar grid inside the disc.
- The teeth are located without the prediction: the grid seeds the polish and the prediction is compared afterwards, so no tooth is selected by the law it tests, and a pole whose zero the search fails to pin is reported as such.
- `N_F(T)` is the running winding at the box boundaries, printed per design, and the zeros per period is `N_F(T) 2 pi/(T log base)`.
- A zero is either matched to a tooth or tagged second family; the second family is what remains when the pole lattice is stripped, and on the base 2 full digit set it is the critical line.

## THE ROUCHE CERTIFICATE

- The disc count is an argument principle on a resolved contour, which is a measurement. The certificate replaces it by an inequality, and the inequality is what a proof needs.
- Split `Z(s_0 + u) = P(u) + T(u)` with `P(u) = (1 - base^(-u)) D_(P-1)(s_0+u) + E_P(s_0+u)` and `T` the `l >= 1` part of the ladder numerator. The split matters: a modulus bound on the Dirichlet polynomial `D_(P-1)` throws away its cancellation and runs about fifty times its true size, so `P` is never bounded, only expanded.
- `P` is entire with Taylor coefficients in closed form, `d_m = sum_(n) n^(-s_0)(-log n)^m/m!` over the two string pools convolved with the coefficients of `1 - e^(-u log base)`, computed exactly to `m = 60` with a remainder bounded by `sum_n n^(rho - Re s_0)(rho log n)^(M/2+1)/(M/2+1)!` and its mirror on the other factor.
- `T` is bounded on `abs(u) <= R_2` by the ladder's own majorant, `sum_(l >= 1) binom(abs(s_0)+R_2+l-1, l) base^(-sigma-l) gamma_l G(sigma+l)` with `sigma = Re s_0 - R_2`, `G` bounded by summing two peel levels exactly before the geometric tail, worth about a factor three over the raw tail. The `l` sum is closed by a majorant ratio and not by an observed one: `gamma_(l+1)/gamma_l <= a_max` and `G(sigma+l+1)/G(sigma+l) <= base^(-(P-1))` because every string in the pools is at least `base^(P-1)`, so the term ratio is at most `R_l = ((abs(s_0)+R_2+l)/(l+1)) a_max base^(-P)`, which decreases in `l` once `abs(s_0)+R_2 >= 1` and is below `a_max base^(-P)` otherwise; stopping at the first `l` with `R_l < 1` and adding `term_l R_l/(1-R_l)` is a proof, where the term ratio itself is not monotone.
- `Z_0` comes from the ladder at `s_0` with its propagated bound. `Z_1` is NOT the central difference: it is the first Fourier mode of `T` on a circle of radius `R`, whose aliasing is `(B_T/R_2)(R/R_2)^N/(1 - (R/R_2)^N)` by Cauchy, plus the exact `p_1`, so no step of the certificate uses a differenced quantity.
- Rouche against the model `Z_0 + Z_1 u`, whose zero count in the disc is exactly one when `abs(Z_0/Z_1) < rho`: certified when `margin = (abs(Z_1)_low rho - abs(Z_0)_up) - sum_(m >= 2) abs(P_m) rho^m - B_T tau^2/(1 - tau) > 0`, `tau = rho/R_2`.
- A positive margin proves exactly one zero of `Z`, hence of `zeta_F` when the residue does not vanish, in `abs(u) < rho`. `rho` and `R_2` are free parameters of the proof and the verb reports the pair it used.
- The peel depth is a free parameter of the proof, not a constant: the verb raises `P` at a pole until the certificate fires or the string pool passes its cap, and prints the depth each certified pole used. A pole that fails at the depth `build` picks automatically may certify two levels deeper.

## RUN

- `uv run python research/zeta/lab/py/zeta-locus/zeta_locus.py shadow 40 all` - the derived law at every design: residue, regular part, prediction, disc counts of `Z` and of `zeta_F`, the located tooth and both misses.
- `uv run python research/zeta/lab/py/zeta-locus/zeta_locus.py rouche 40 rou` - the certificate at every pole of the nine designs the landed census counts: margin, the `rho` and `R_2` achieving it, `B_T`, the bounds on `abs(Z_0)` and `abs(Z_1)`, the sample count, and the argument principle count in the same disc for comparison.
- `uv run python research/zeta/lab/py/zeta-locus/zeta_locus.py census 40 all` - the same plus the strip census, `N_F(T)`, one row per zero and the falsification lines.
- The third argument selects a family: `q3`, `q4`, `q5`, `half`, `wide`, `ctl`, `rou` the nine designs of the Rouche census, `all`. The second is the height.
- Prints only, writes nothing. The full census to height 40 is about fifteen minutes; base 10 costs about nine times base 3 per evaluation.
