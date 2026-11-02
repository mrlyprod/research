# transport-census

- The transport census: for every design, the zeros of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)` that lie right of the abscissa `alpha = log k/log q`, and the lower bound on the Mertens exponent of the design's own Mobius `nu_F` each one forces.
- The transport theorem of [mrly-pairing](../../../../mobius/lab/py/mrly-pairing/README.md) is the engine: a zero `rho` of `zeta_F` with `Re rho > alpha` gives `sigma_c(N_F) >= Re rho`, so `sum_(n <= x) nu_F(n)` is not `O(x^(Re rho - eps))` for any `eps > 0`. This study turns that one-zero statement into a table over every design the locus sweep censuses.
- The ladder, the contour engine, the truncation bounds and the box machinery are imported from [design-zeta](../design-zeta/README.md), the cofactor and the polisher from [zeta-locus](../zeta-locus/README.md); nothing is copied and no sweep is repeated.

## THE STRIP IS CLOSED ON THE RIGHT

- A census of the zeros right of `alpha` needs no arbitrary right edge. `zeta_F(s) a_min^s -> 1` as `Re s -> +infinity` with `a_min` the least nonzero digit, and the coefficients are nonnegative, so `abs(zeta_F(s) a_min^s - 1) <= a_min^sigma zeta_F(sigma) - 1` at `sigma = Re s`.
- Hence `sigma_1`, the least point of the grid `alpha + 0.05 n` at which `a_min^sigma (zeta_F(sigma) + bound) < 2`, is a PROVED zero free edge: for `sigma_1 > alpha`, which the grid forces, `zeta_F` has no zero with `Re s >= sigma_1`, the ladder's own error bound carried through.
- So the census box needs no hand-chosen right edge, and its count is EXACT ON THE BOX: uncounted are the sliver `alpha < Re s <= alpha + 1e-6`, the band `0 < Im s < 0.02`, everything above the census height and the conjugate half plane, so the same number is a LOWER BOUND for the half plane `Re s > alpha`.

## THE COUNT

- The winding is taken on the Lyndon cofactor `Z(s) = zeta_F(s)(1 - k q^(-s))`, which is analytic on `Re s > alpha - 1` and has no pole in the box, so the argument principle counts zeros with no pole correction.
- Right of `alpha` a zero of `Z` is a zero of `zeta_F` without exception: the one place the transfer fails is a pole `s_(0,j)` with vanishing residue, and every such point sits ON the line `Re s = alpha`, outside the box.
- The left edge is `alpha + 1e-6`, which is what keeps the full digit sets' residue null teeth outside the contour.
- The box is cut into unit sub-boxes in `Im s`; each sub-box of nonzero winding has its zeros located by a coarse grid and polished by Muller on the cofactor, and the count located against the count wound is printed as the completeness flag.

## THE BOX

- The rightmost zero of each design is boxed, resolved and not certified, by a winding box on `zeta_F` itself, half-width `5e-5` in `Re s` and in `Im s`, the argument principle on the same engine as [mrly-pairing](../../../../mobius/lab/py/mrly-pairing/README.md) verb `box`.
- A winding of `1` resolves exactly one zero inside the rectangle, sampled and not certified, since nothing bounds the derivative between samples, so `Re rho` is pinned to the box edges and the transport bound is the LEFT edge, truncated down.
- A zero right of `alpha` may be a tooth of the level-zero comb and so sit near a pole of `zeta_F`, so every box prints its distance to the nearest pole of the lattice `s_(i,j) = alpha - i + 2 pi i j/log q`. That distance is positive by construction, since the lattice lies on `Re s <= alpha` and the box lies right of `alpha`; what the column carries is how far, against the box's own half-width.
- The printed contour minimum is a minimum over the SAMPLED points of the contour and not over the contour, so it reads as evidence beside the winding and not as a certificate of its own.

## THE HEIGHT

- A rightmost real part is a statement below a cut and nothing more. The teeth of the level-zero comb drift right with the pole index, so the rightmost zero of a design moves right as the census height rises, and every rightmost printed here carries the height it was read below.
- Twenty-three designs are censused to `Im s = 40` and base 50 missing one digit to `Im s = 4` alone, so its column is a height-`4` maximum of a drifting shape and is not on the same clock as the others.

## THE FALSIFICATION

- The base 2, 3 and 4 full digit sets are the controls: there `zeta_F = zeta`, no zero lies right of `alpha = 1`, and the census must return zero.
- At a full set the boxes around the cofactor only teeth `s_(0,j)`, `j >= 1`, separate the two functions: `Z` winds `1` there and `zeta_F` winds `0`, which is the transfer exception seen directly.
- The two boxes, resolved not certified, of [mrly-pairing](../../../../mobius/lab/py/mrly-pairing/README.md) are re-run at their own edges, with the control rectangle.

## THE HYPOTHESES

- `nu_F` is the Dirichlet inverse of `1_(S_F)` and exists only when `1 in S_F`, that is when `1 in F`. A design missing the digit `1` carries its zeros and carries no transport bound; the table says so rather than quoting the exponent.
- Every bound is `theta(nu_F) >= x_0` with `x_0` the left edge of a winding-`1` box, resolved not certified, truncated down; nothing is read off a polished root or a small residual.

## RUN

- `uv run python research/zeta/lab/py/transport-census/transport_census.py census 40 locus` - the per design table: `sigma_1`, the strip, the winding, the zeros right of `alpha`, the box on the rightmost, resolved not certified, and the bound it gives.
- `uv run python research/zeta/lab/py/transport-census/transport_census.py root 40 locus` - the crossing: the root of `a_min^sigma zeta_F(sigma) = 2` bisected against the ladder bound and rounded up, beside the grid edge `sigma_1` and `a_min^sigma zeta_F` at it; about four seconds.
- `uv run python research/zeta/lab/py/transport-census/transport_census.py law 40 locus` - the same, followed by the columns in `k/q` order and in `alpha` order, the equal-key ties, the designs with `Re rho > 1`, the corollary line per design, the full set null boxes and the reproduction of the landed boxes.
- The second argument is the height, the third a family: `locus` the twenty of the locus sweep, `rungs` the two of the family sweep, `b20` and `b50` the two high rungs, `ctl` the three full digit sets, `all` the twenty-two. Prints only, writes nothing.
- `census 40 locus` runs in about eleven minutes, `census 40 rungs` in three and `census 40 b20` in six; base 50 is censused to height `4` in two and a half. Peak resident memory is under `0.2` GB throughout.

## WITNESSES

- the zero free edge is one real evaluation: `sigma_1` reads `0.5` to `1.75` over the twenty-four designs, with `a_min^sigma zeta_F(sigma)` between `1.86` and `1.999` against the threshold `2` and ladder bounds `1e-11` to `1e-47`, so `zeta_F` has no zero at all in `Re s >= 1.75` on any design censused and the `alpha + 3.02` right edge of [zeta-locus](../zeta-locus/README.md) is more than twice as wide as the zeros need
- the winding counts `157` zeros of `zeta_F` right of `alpha` over the twenty-three designs censused to `Im s = 40`, all `157` located at the seek grid `15 x 11` with `2 want + 8` seeds, and `2` more at base 50 missing one digit censused to `Im s = 4`; the count is exact on the box and a lower bound for the half plane `Re s > alpha`
- twenty-one of the twenty-four carry a zero right of `alpha`, nineteen of the twenty-two the locus and family sweeps censused; the three that do not are the base 2, 3 and 4 full digit sets, whose windings on the strip read `-1.97e-33`, `1.73e-33` and `1.53e-33`, and that reproduces the nineteen of [zeta-family](../zeta-family/README.md) from a construction with no assignment radius in it
- the rightmost real parts below the census height run `0.441505537191` at base 5 `{0,1}` to `1.002685494780` at base 20 missing one digit, both below `Im s = 40`, and the bounds they prove run `theta(nu_F) >= 0.4414555` to `theta(nu_F) >= 1.0026354`, each the left edge of a winding-`1` box truncated down; the size of a rightmost is a statement below its own cut, base 20 missing one digit reading `1.000285484146`, `1.000549674321` and `1.002685494779` at `Im s = 2.0988`, `4.1971` and `14.6920`
- every box returns winding `1` with a sampled contour minimum of `1.2e-4` to `6.1e-3` against engine bounds of `1e-13` to `1e-33`, at least eight orders of magnitude at every box, and the tightest pole distances are `0.00517845` at base 50 missing one digit and `0.0223021` at base 20 missing one digit, four hundred times the box's own half-width
- nineteen of the twenty-four designs meet every hypothesis of the transport bound, and seventeen of the twenty-two the locus and family sweeps censused: base 4 `{2,3}` and base 4 `{0,2,3}` have zeros right of `alpha` and do not contain the digit `1`, so `1` is outside `S_F`, the indicator vanishes at `1` and the Dirichlet inverse `nu_F` does not exist
- the gain `Re rho - alpha` is not a function of `alpha` and `k/q`: at `alpha = 1/2`, `k/q = 1/2` the four base 4 two-digit designs read `0.0853043873`, `0.4400124317`, `0.2706238545`, `0.3439264581`, a spread of `0.35470804`, and the other three equal-key families spread `0.37474232`, `0.17605693` and `0.060972003`
- four designs have `Re rho > 1`, so their `nu_F` outruns the count of integers below `x`: base 10 missing two, base 10 missing `9`, base 20 missing one and base 50 missing one, rightmost `1.001514387650`, `1.001589275290`, `1.002685494780` and `1.000061474950` at `alpha = 0.9030900, 0.9542425, 0.9828779, 0.9948357`, the first three below `Im s = 40` and the last below `Im s = 4`, so only the SIGN of `Re rho - 1` is read and never its size; the `k/q` reading dies on base 5 `{0,1,2,3}`, the same `k/q = 0.8` as base 10 missing two, whose rightmost is `0.989748105861`, below `1`; base 20 and base 50 missing one digit are two further rungs the base 10 pair does not fix, and both land above `1`
- the least rightmost real part at each `alpha` reads `0.4485242462`, `0.4415055372`, `0.5853043873`, `0.7207876015`, `0.9126562295`, `0.9897481059`, `1.0015143877`, `1.0015892753`, `1.0026854948`, `1.0000614750` up the ladder `alpha = 0, 0.4307, 0.5, 0.6309, 0.7925, 0.8614, 0.9031, 0.9542, 0.9829, 0.9948`, rising at every step but the first and the last, and the last is where the census height drops from `40` to `4`; ten rungs, one design each above `alpha = 0.86` against six at `alpha = 0.5`, and no fit taken
- the cofactor only teeth separate the two functions directly: around `s_(0,1)` at `1 + 9.06472028365 i`, `1 + 5.71920173476 i` and `1 + 4.53236014183 i` the winding of `Z` is `1` and the winding of `zeta_F` is `-9.42e-32`, `-6.28e-32` and `0`, with `min abs(Z)` about `5e-5` and `min abs(zeta_F)` `1.351`, `0.8747` and `0.7295`
- the root of `zeta_F(sigma) = 2` reads `1.1790285003` at base 3 `{0,1}` and `1.69933660058` at base 10 missing `9`, below their grid edges `alpha + 0.55` and `alpha + 0.75`, where `zeta_F` reads `1.995208619` and `1.991335591`, and the three full sets read `1.728647239`, the root of `zeta(s) = 2`
- the two boxes, resolved not certified, of [mrly-pairing](../../../../mobius/lab/py/mrly-pairing/README.md) reproduce at their own edges, winding `1` and `1` with contour minima `8.298e-4` and `6.865e-4`, and the control rectangle returns winding `0` with contour minimum `1.543e-2`

## SOURCES

- [DLMF 25.2](https://dlmf.nist.gov/25.2) - the Dirichlet series and the abscissa of convergence, the classical `A(x) = O(x^theta)` implies `sigma_c <= theta` the transport theorem closes on.
- [Baker and Harman 1991](https://doi.org/10.1112/jlms/s2-43.2.193) - the exponential sum bound the study [mrly-pairing](../../../../mobius/lab/py/mrly-pairing/README.md) costs, cited there and not used here.
