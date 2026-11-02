# zeta-family

- The SECOND FAMILY of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)`: the zeros that survive when every comb of the pole lattice is stripped at a stated radius, and the law they obey.
- On the full digit set the second family is exactly the critical-line zeros of `zeta`, so it is the half of the zero set that carries the Riemann hypothesis; on a design it has no law.
- The ladder, the contour engine and the truncation bounds come from `../design-zeta`, the cofactor, the Laurent data and the residue comb from `../zeta-locus`; neither is copied. This study adds the split into second family and comb, the level-one comb, and three derived tests.

## THE SPLIT

- The pole lattice of `zeta_F` is `s_(i,j) = alpha - i + 2 pi i j/log q`, one vertical line per level `i = 0, 1, 2, ...`, and `Z(s) = zeta_F(s)(1 - k q^(-s))` kills the level-zero line and is analytic on `Re s > alpha - 1`.
- Every zero of `Z` in the census strip is a zero of `zeta_F`, except at a lattice point where the residue vanishes: there `Z` vanishes and `zeta_F` does not.
- THE ASSIGNMENT. A zero within `abs(u) < 0.45` of a pole with nonvanishing residue is a tooth of that pole's comb; a zero sitting at a pole with vanishing residue is a zero of the cofactor alone. Everything else is SECOND FAMILY.
- The next cofactor `Z_m(s) = zeta_F(s) prod_(i <= m) (1 - k q^(-(s+i)))` has exactly the zeros of `Z` inside the strip, since the extra factors vanish only on `Re s = alpha - 1, ..., alpha - m`: the survivors are the same set under every comb. What `Z_m` adds is the level-one comb, whose teeth reach into the strip.
- THE LEVEL-ONE RESIDUE, derived and not measured. `Z(s) = E_1(s) + sum_(l >= 1) binom(-s,l) q^(-s-l) gamma_l zeta_F(s+l)` is singular at `s_(1,j)` through its `l = 1` term alone; `q^(-s_(1,j)-1) = q^(-s_(0,j)) = 1/k` and `1 - k q^(-s_(1,j)) = 1 - q`, so `zeta_F` has residue `r_(1,j) = s_(1,j) gamma_1 r_(0,j)/(k(q - 1))` at the level-one pole.
- Two readings for free: the level-one comb is empty wherever the level-zero comb is, and at the full digit set `s_(1,0) = alpha - 1 = 0` makes `r_(1,0)` vanish, which is `zeta` having no pole at `s = 0`.
- A level-one tooth lies within `abs(u) < 0.45` of `Re s = alpha - 1`, so it is visible in the strip only in the band `alpha - 0.92 < Re s < alpha - 0.55`, and there the assignment is applied against both lattices.
- THE TRANSPORT READING. A zero of `zeta_F` at `Re rho > alpha` forces the design's own Mobius `nu_F`, the Dirichlet inverse of `1_(S_F)`, to have Mertens exponent at least `Re rho`. So the rightmost second-family real part is a proved lower bound on that exponent, printed per design.

## THE THREE TESTS

- The expectations are derived first, from the full digit set and from the classical formulas, and never from the design sweep.
- SYMMETRY. On the full set the functional equation puts a zero at `(1 - sigma) + i t` for every zero at `sigma + i t`, so the real parts are symmetric about `Re s = 1/2 = alpha/2` and each zero is its own partner. The test reads `c_F` as the midpoint of the real parts of the two second-family zeros of least `Im s`, then asks of every remaining zero whether `(2 c_F - sigma) + i t` is again a second-family zero, and compares `c_F` against `alpha/2`, `alpha` and `alpha - 1/2`.
- The design has no functional equation: the reflection moves the kernel and not the design, so the expectation is failure at the third zero, and the falsifier is a design where the pairing holds.
- COUNT. On the full set `N(T) = (T/2 pi) log(T/2 pi e) + O(log T)`, hence `N(T) = c T log T + d T` with `c = 1/(2 pi) = 0.15915494` and `d = -(1 + log 2 pi)/(2 pi) = -0.451662163`. Reading `c` and `d` from two heights rather than the closed form costs the control about two percent, which is the accuracy of the test; the question is whether `c_F` is a function of `alpha`, of the fill, or of the design alone.
- The count needs no zero located: `N_2(T)` is the winding of `Z` on the box less the zeros of `Z` inside every pole disc, both read by the argument principle.
- THE FULL-SET LIMIT. As `alpha -> 1` through designs missing one digit and then two, the second family's real parts contract to `alpha/2` or they do not. A contraction is the rule that makes RH a special case of MrlyMath; a non-contraction says the critical line is a property of the Euler product alone. The statistic is the mean and the spread of `Re s - alpha/2` per design, ordered by `k/q`.
- The controls are the base 2, 3 and 4 full digit sets, where the answers are known: `c_F = 1/2`, `N(T)` classical, spread zero.

## THE BLIND BAND

- The ladder does not reach tolerance in `alpha - 1 < Re s < alpha - 0.92`, a band against the cofactor's own pole line, so every census here runs on `alpha - 0.92 < Re s < alpha + 3.02` and every count is a lower bound for the strip to `alpha - 1`.
- The box height is nudged up to the midpoint between two poles whenever a pole sits within `0.5` of it, so no pole disc straddles the top edge; the height actually used is printed with every count.

## WHAT THE SWEEP FOUND

- Twenty-two designs: the twenty of `../zeta-locus`, plus base 5 `{0,1,2,3}` and base 10 missing two digits as the rungs that carry `alpha` toward 1. Each is censused to its own printed height, `40` except the four base 3 designs at `42.894`, base 9 `{0,1,2}` at `41.464` and base 10 missing two at `25.923`, every height nudged off the pole lattice.
- 377 zeros wound on the strips, 351 located, 171 teeth of which 9 are level-one teeth, 19 cofactor-only zeros at null-residue poles, and 161 second family, all at `rho = 0.45`. Where the located count falls short of the winding the printed `N_2` is a lower bound, base 4 `{2,3}` at 12 of 18 being the worst.
- THE RADIUS IS AN UNDEFENDED CONSTANT AND THERE IS NO GAP AT IT. `rho = 0.45` is fixed only by the pole discs not overlapping, while the tooth law is accurate only inside `abs(u) < 0.2` and says nothing past `0.3`. The nearest live pole to a second-family zero is `0.45510938` away at base 4 `{2,3}`, `0.45909168` at base 3 `{0,1}`, `0.48696667` at base 4 `{0,1,2}` and `0.50481072` at base 4 `{1,3}`, five of base 4 `{2,3}`'s eight sitting inside `0.45 < abs(u) < 0.6`, so the cut lands in a pile-up and every count is conditional on `rho`.
- Both radii are printed. `N_2` reads `8, 7, 13, 9, 14` at `rho = 0.45` and `7, 6, 12, 8, 7` at `rho = 0.6` on base 3 `{0,1}` and the four base 4 two-digit designs. The full digit set is where the gap is real: at base 2 the nearest live pole is `14.143566` away and no radius below `0.9` moves a count.
- Base 3 `{0,1}` reads `N_2 = 7` below `Im 40` and `8` to its nudged height `42.894`; the seventh zero `-0.273079611 + 39.262315320 i` sits `0.0160` right of the strip edge, inside the censused box and outside the blind band.
- THE LEVEL-ONE COMB BITES. Nine of the twenty-two designs carry a zero the level-zero comb calls second family and the level-one comb claims, all in the band against the cofactor's pole line; without it the second family would be overcounted by nine.
- THE CONTROLS. Base 2, 3 and 4 full digit sets give `N_2 = 6, 7, 6` with every real part `0.5` to `1.4e-22`, `1.1e-21` and `1.8e-22`, the comb empty and every tooth a cofactor-only zero.
- THE TRANSPORT READING. Nineteen of the twenty-two designs carry a zero right of `alpha`, twelve of them from the second family alone, and the seventeen of those whose digit set contains `1` get a proved lower bound on the Mertens exponent of their own Mobius `nu_F`; base 4 `{2,3}` and `{0,2,3}` omit the digit `1` and have no `nu_F`.
- The strongest is base 10 missing two digits, `alpha = 0.903089987`, with a zero at `1.00151438765 + 2.77402670058 i`: the Mertens exponent of `nu_F` exceeds 1 at a second base-10 column. Base 4 `{1,2}` at `alpha = 1/2` has a SECOND-FAMILY zero at `0.940012431696 + 13.0678968771 i`, an exponent of `0.94` against a design mass exponent of `0.5`.

## WHAT DIED

- NO SYMMETRY, AND NOTHING PAIRS AT ALL. No second-family zero in any design has a reflection partner: the branch needs two second-family zeros within `0.05` in `Im s`, and the smallest ordinate gap inside a design is far above that on every design tested, so it cannot fire. Every recorded pair is a self-pair, a real part landing within `0.05` of `c_F`, and self-pairs run BELOW chance: 8 of 47 with 0 partners over the ten designs recensused, a rate of `0.170213` against the `0.229904` that a uniform draw from each design's own band predicts, and 22 of 117 over the full sweep. The three full-set controls pair 13 of 13 at `c_F = 1/2` to `1e-22`, where the functional equation makes every zero its own partner.
- `c_F` is not a quantity: `c_F - alpha/2` runs `-0.28413232` to `+0.47788515` and `c_F - 1/2` runs `-0.78413232` to `+0.28664994`. The reflection moves the kernel and not the design, and the census says so zero by zero.
- NO COUNTING LAW IN `alpha` OR FILL, AT EITHER RADIUS. The four base 4 two-digit designs share `alpha = 1/2` and `fill/base = 1/2` exactly and read `N_2(40) = 7, 13, 9, 14` at `rho = 0.45` and `6, 12, 8, 7` at `rho = 0.6`, with `N_2(80) = 20, 30, 22, 29` and `17, 26, 21, 20`: a factor of two at both radii, so the refutation is radius-robust even though the integers are not.
- The spread is carried by COMB OCCUPANCY and not by the second family: base 4 `{0,1}` and `{2,3}` differ 29 percent in total winding, 14 against 18, and a factor of two in `N_2` because 7 of 8 poles are occupied against 4 of 8, and by the pile-up above an unoccupied pole is a pole whose zero drifted past `rho`. Every winding is the nearest integer to an integrated phase whose largest surviving step runs `0.9205` to `0.9998` against a cap of `1`, so these counts are Verified and not Proved.
- Read as `N_2(T) = c_F T log T + d_F T` from the two heights, `c_F` is `0.10820213, 0.072134752, 0.072134752, 0.018033688` at `alpha = 1/2`, a spread of `0.09016844`, while the two full-set controls read `0.15486803` and `0.16230319` against the classical `1/(2 pi) = 0.15915494` and spread `0.0074351582`. The equal-`alpha` spreads are `0.071807313, 0.09016844, 0.013413767, 0.09016844` against that control spread.
- NO CONTRACTION OF THE REAL PART. The refuted law is `max abs(Re s - alpha/2) -> 0` as `alpha -> 1`. Undivided the statistic stays FLAT along the ladder, reading `0.2275679549, 0.5549589411, 0.5885444877, 0.4233198337, 0.5365616661, 0.3151426744` at `alpha = 0.430676558, 0.5, 0.630929754, 0.792481250, 0.861353116, 0.903089987`, and then collapses to `1.43e-22`, `1.10e-21` and `1.76e-22` at the base 2, 3 and 4 full sets. At base 10 missing 9, `alpha = 0.954242509`, the real parts read `0.216084781875` to `0.70401657869` about `alpha/2 = 0.477121255`, a band of `0.488` against `1 - alpha = 0.0458`.
- Divided by `1 - alpha` the same statistic runs `0.39971647` to `5.7047812` with NO monotone in `alpha`, falling from `3.8699872` to `3.2519104` on the last two rungs that carry `alpha` toward 1, so the refutation rests on the undivided spread and the ratio is not quoted for it.
- So the critical line is not reached by contraction. It appears at `alpha = 1` as a jump from a band of width `0.49` to `1e-22`, and that is the Euler product speaking, not the digit set.

## WHAT THE ORDINATES DO INSTEAD

- The heights converge even though the real parts do not. Against the derived null - a quarter of the mean gap between consecutive `zeta` ordinates in the range, exact for an equally spaced set of the same density and about 18 percent conservative for one with gap variance - the mean distance from a second-family ordinate to the nearest `zeta` ordinate divided by that null reads `2.0495374` at `alpha = 0.430676558`, `1.8953371` at `0.5`, `0.75419266` at `0.630929754`, `0.51648744` at `0.792481250`, `0.32356636` at `0.861353116`, `0.090501352` at `0.903089987` and `1.0429899e-23` at `alpha = 1`.
- Monotone in `alpha` over seven rungs with a control that reads zero exactly, and the base and fill confounds are dead: the fall is monotone at FIXED BASE, `2.0495374` to `0.32356636` inside base 5 and `1.8953371` to `0.51648744` inside base 4, and at FIXED fill 2 across bases, `2.0495374, 1.8953371, 0.75419266, 1.0429899e-23`. The nulls move only `1.0425839` to `1.3595166` across the ladder while the raw mean distance falls `2.4032315` to `0.12303809`, so the denominator is not driving it.
- `alpha` is a trend and not a function: the four base 4 two-digit designs at one `alpha = 1/2` spread `0.79050661` to `2.8404536`. Over the ladder `mean abs(Re s - 1/2)` reads `0.36482392, 0.39426128, 0.3901396, 0.25540269, 0.31452367, 0.20473972` and `2.4065966e-23` and does not fall monotonically at all.
- Reading: A FILLING DESIGN LEARNS WHERE `zeta`'S ZEROS ARE BEFORE IT LEARNS THEY LIE ON A LINE. At `alpha = 0.903` the heights are pinned to `2.3%` of the mean gap while the real parts are still off `1/2` by `0.20`. The Riemann hypothesis is the last thing to appear and it appears only at the full digit set.
- One caveat printed with the rows: the matching is nearest-ordinate and not injective, `3` distinct `zeta` ordinates for `4` design zeros at base 10 missing two and `5` for `6` at base 5 missing one, so the design carries a few zeros with no `zeta` partner.

## WHAT STANDS AND WHAT IS OWED

- The second family is comb-independent and censused at a stated radius; it is not symmetric, it has no counting law in `alpha` or fill at either radius, and its real parts do not contract to `alpha/2`.
- What does converge is the ordinate set, and that is the open object: the rate of the ordinate convergence in `1 - alpha`, and whether the extra zeros with no `zeta` partner are the design's own.
- Owed: a defended assignment radius, which the tooth law does not supply past `abs(u) = 0.3`; the shadow row at base 10 missing 9, one command; a Rouche certificate for the disc counts; and the blind band `alpha - 1 < Re s < alpha - 0.92`, which `Z_1` removes as soon as the ladder can stop one level higher.

## RUN

- `uv run python research/zeta/lab/py/zeta-family/zeta_family.py family 40 all` - the second family split from the comb: per design the counts `N_2` at both radii, the real parts, the rightmost real part with its Mertens reading, and per zero the residue of `Im s` against `2 pi/log q`.
- `uv run python research/zeta/lab/py/zeta-family/zeta_family.py symmetry 40 all` - the family sweep plus the reflection test about `c_F`.
- `uv run python research/zeta/lab/py/zeta-family/zeta_family.py count 40 all` - the winding count at two heights and the `c_F T log T + d_F T` read.
- `uv run python research/zeta/lab/py/zeta-family/zeta_family.py limit 40 ladder` - the same on the rungs that carry `alpha` toward 1, followed by the ordinate shadow against the zeros of `zeta` and its derived null.
- The second argument is the height, the third a family: `q3`, `q4`, `q5`, `half`, `wide`, `ctl`, `near`, `ladder`, `rungs`, `q34`, `all`; `tests` runs symmetry, limit, shadow and transport together. Prints only, writes nothing.
- The full sweep to height 40 is about thirty-five minutes and the ladder about seven; base 10 costs about nine times base 3 per evaluation.
