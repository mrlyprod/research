# design-dimers

- Domino tilings of a plane design: a code at base `base` is the mask of cells `(i, j)` with bit `base*i + j`, level `n` its `n`-th Kronecker power, and `T(n)` the number of perfect matchings of the cell graph.
- Counts are `abs(det K)` with `K` the black-by-white Kasteleyn matrix: sign `(-1)^c` on a vertical edge in column `c`, sign `(-1)^h` on a horizontal edge `(r, c)(r, c+1)` with `h` the empty cells below it in column `c`. Determinants are exact integers from PARI (`gp -q`, `matdet`).
- The control is a brute-force transfer count over a broken row profile, exact integers in Python, and an augmenting-path matching test.
- `imbalance` checks the black-minus-white count of every code at base 2 through level 6 and base 3 through level 4 against `s^n` at odd base and `fill^(n-1) s` at even base, and counts the codes with `s = 0`.
- `census` finds the first tileable level of every code with `s = 0` at base 2 and 3 through level 4 and at base 4 through level 3, counts the untileable ones by orbit and by a sealed unbalanced component of `D_m` for `m = 1, 2, 3`, tries the base-4 codes left at level 4 and by the run certificate for codes whose blocks meet in one direction only, and checks nine late codes at bases 5 and 7: untileable at level 1, unsealed, `T(2)` by determinant and, at base 5, by brute force.
- `search` draws dense codes from a fixed seed at bases 5, 6 and 7, each cell filled with probability `0.7`, keeps the distinct codes, and counts those that tile first at level 2, and at level 3 among those with `fill^3 <= 16000` whose level-2 deficiency falls below `fill` times the level-1 deficiency; then it walks every base-5 code within Hamming distance 4 of the late orbits, by orbit, and tries each code untileable at levels 1 and 2 at level 3.
- `count` runs the control on 300 random cell sets, then prints `T(n)` with its 2-adic split for code 15 at base 2 through level 6 (against the rectangle product formula), code 63 at base 3 through level 4 (against `F(3^n+1)^(2^(n-1))`), the carpet 495 and its sibling 255 through level 4, brute force at side 9 and below. It caches the counts in its data folder for `cross` and `growth`.
- `cross` weights the edges between the big blocks by `x`, evaluates the determinant at `x = 0..E` and interpolates in PARI: the number of tilings by crossing count for the carpet and its sibling at levels 2 and 3 and for code 15 at base 2 at levels 2 to 4, the parity, the mean, the share that respects the blocks, and the square test.
- `growth` prints `log T(n) / fill^n`, its increments and their ratios, the Hadamard upper bound from the 16-state exposure recursion at levels up to 400, with the idempotence check, the bracket, the geometric-tail closures and the exact constants.

## RUN

- `uv run python research/dimers/lab/py/design-dimers/dimers.py` runs every verb; one verb by name, e.g. `... dimers.py count`.
- Runtimes: `imbalance` 0.1 s, `census` about 5 s, `search` 134 to 160 s, `count` 49 s (the three 4096-cell determinants), `cross` 3 s, `growth` under 1 s with the cache; about 3.5 minutes in all.
- Needs `gp` on the path; standard library only.

## WITNESSES

- dimers.md, "The colour imbalance": the law at every code, the counts 5 and 125 of codes with `s = 0` (`imbalance`).
- dimers.md, "The census of tileable codes": the base-2 and base-3 split 97 and 28, the 28 codes and their 5 orbits, the base-4 counts 5699, 7170, 945, 2253, 4709, 68, the 72 codes in 10 orbits of the run certificate, the 2803 one-direction codes that tile level 1, and the 68 open codes in 9 orbits, the nine late codes with their `T(2)` (`census`); the sample counts at bases 5, 6 and 7 (`search`).
- dimers.md, "Kasteleyn on a design": the 300 agreements, the 291 of the column signs alone, the 64 tileable sets, the zero determinant on the carpet without the hole signs (`count`).
- dimers.md, "The counts": the table, the control against A004003 and the product formula, the code-63 formula, the 2-adic splits (`count`).
- dimers.md, "What crosses": the parity at odd and even base, `P_2(x)`, the means, ratios and shares (`cross`).
- dimers.md, "The growth constant": the lower bounds, the Hadamard limits, the brackets, the increment ratios and the tail closures (`growth`).
