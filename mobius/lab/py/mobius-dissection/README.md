# mobius-dissection

- The numbers of the unconditional dissection on [mobius](../../../notes/unconditional-dissection.md): the wall its region A sets, the block split that carries it off the powers of the base, and a falsification of every region's bookkeeping on the whole grid of three small designs.
- `wall`: the least base with `PB < (base - 1) base^(-4/5)` at one excluded digit, for the step 3 constant `PB_base(1)` and for the chord form `PB'_base(1, e_0)`, the latter at the worst digit and at `e_0 in {0, base - 1}`; the up-set flag and the `held` count over the scan, the least float gap above each wall, the gap at the wall and one below at 40 digits, and `1/5 - alpha_1` at each wall; then readings of the unshifted mass ratio `c_k/c_(k-1)` at base 33 and base 17 against `2(base - 1)`.
- `blocks`: the split of `S_F` below `x` into blocks `P base^k + D_k`, checked against the exact `M_F(x)` at 400 random `x` below `2 * 10^6` and at every `base^e - 1`, in five one-missing-digit families, with the count per scale against `fill + 1` and the mass floor `A_F(x) >= fill^(L-1) - 1`.
- `regions`: every grid point `a mod y` of base 10 missing 5 and missing 0 at level 6 and base 5 missing 2 at level 9, assigned its Dirichlet fraction at `Q = y^(3/5)` by the last convergent, and cut into regions A, B, C1, C2 at the printed `Z`; then the region sums against the exact sum of `mu` over the strings, the two counts, the C2 algebra, the second approximation at every point of B with `h >= 1`, and the hybrid `l^1` bound at every class meeting the two conditions its proof uses, `V_1 V_2 <= y` and `16 D H <= y`, a wider set than its hypothesis `16 base^2 D (D + H) <= y`, with `alpha_1` from the step 3 constant, printed past the class of `a = 0`, whose ratio is a closed form; and the Gallagher step of the proof at every such class with `V_1 <= 10^5`, the fractions' largest `|hat F_(i_1)|` at their residues summed against `4 D^2 ||f||_1 + ||f'||_1`, both norms read on `4 V_1` points, with every fraction's residues checked distinct mod `V_2`. The perturbed Lemma A' is then checked at seeded grid points of level 30 in base 10 and level 45 in base 5, where its perturbation hypothesis holds, with the phases `base^i a mod y` in exact integers.

## THE ROUNDING

- The walls are float scans with the least gap above each wall printed; that gap sits at the wall, above `10^-6`, and the wall and the base below it are re-read in `mpmath` at 40 digits, so the sign at each is certified far past float error.
- `1/5 - alpha_1` is printed from the 40-digit gap through `log(1 + gap/PB)/log base`, never by differencing two numbers of size `1`, its mantissa floored at four decimals in `mpmath`.
- The hybrid ratio and the `alpha_1` of the small designs add `5 * 10^-7` before a six-digit print, so each printed value is an upper bound.
- The minor-arc bound carries an unstated constant, so `regions` prints its readings on the grid as readings; they bound nothing.

## RUN

- `uv run python research/mobius/lab/py/mobius-dissection/dissection.py wall` in under a second.
- `uv run python research/mobius/lab/py/mobius-dissection/dissection.py blocks` in 18 seconds.
- `uv run python research/mobius/lab/py/mobius-dissection/dissection.py regions` in five seconds, peak resident memory 0.59 GB.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The wall: `39363` over every excluded digit and `28352` at `e_0 in {0, base - 1}` with the chord form, `92317` with the step 3 constant, each an up-set of its scan to `2 * 10^5`; gaps `<= -2.3195 * 10^-6` at `39362` and `>= 2.3677 * 10^-5` at `39363`, `<= -1.3539 * 10^-5` at `28351` and `>= 1.8837 * 10^-5` at `28352`; `1/5 - alpha_1 >= 2.6965 * 10^-7` at `39363`, `2.3643 * 10^-7` at `28352` and `1.8712 * 10^-8` at `92317`.
- The unshifted mass: `c_4/c_3 = 74.1654` at base 33 missing 16 and `70.5661` missing 0, `c_5/c_4 = 35.5234` at base 17 missing 8, readings above `2(base - 1)`.
- The block split: exact at every tested `x` in all five families, at most `fill` blocks at one scale.
- The regions: the four sums meet the exact total at all three designs; `74778` and `128942` second approximations in range; `73` hybrid classes at each base-10 design and `86` at base 5 under the proof's two conditions, largest ratio past the class of `a = 0` at most `0.001356`, the Gallagher step at most `0.110814`; the perturbed Lemma A' bound exceeds `|hat F_k(a/y)|` by at least `27.9` in the logarithm at all `6000` draws.
