# Boolean Measures

- Computes seven Boolean complexity measures, `s`, `bs`, `C`, `dt`, `deg`, `dnf`, `cnf`, on every hyperoctahedral class of designs at `dim 3` (22 classes) and `dim 4` (402 classes), from the truth-table definitions.
- Diffs the result cell by cell against `measures_d3.csv` and `measures_d4.csv`, which travel as data; fits the fill polynomial from rendered odd-side grids and compares it with the closed form.
- Groups the catalog by genus, `GF(2)` degree, popcount and fill polynomial, and reports which keys determine which measures, with the smallest witness pair at `dim 4`.
- Checks the pin family, the inequality web, total influence in exact rationals at `dim 1..4`, and the `(s, bs, C)` profile of every design at `dim 4` plus seeded samples at `dim 5`.

## RUN

- `uv run python research/complexity/lab/py/boolean-measures/measures.py`
- Domain: every design at `dim 1..4`; at `dim 5` a seeded sample of 50000 uniform designs, plus 20000 uniform and 20000 thinned designs; about three seconds.

## WITNESSES

- complexity.md:61 286 cells at `dim 3` and 5226 at `dim 4`, zero mismatches.
- complexity.md:66-87 the `dim 3` table of 22 classes, printed row by row.
- complexity.md:90-91 genus `iso 10, axis 2, comp 10`; `GF(2)` degree histogram `-1: 1, 0: 1, 1: 3, 2: 9, 3: 8`.
- complexity.md:94-96 the three classes `bang dim 3, code 7`, `bang dim 3, code 31`, `bang dim 3, code 23` with `bs = 2`, `dt = 3`.
- complexity.md:104-118 the pin family, six measures equal `r`, `dnf = 1`, `cnf = r`, at `dim 3` and `dim 4`.
- complexity.md:141-150 fill polynomial 21 of 22 at `dim 3`, one collision `4k^3 - 4k^2 + k` at odd side `2k - 1`; at `dim 4` 183 of 402, 94 shared, 92 groups, 81 split, 279 measure-splits.
- complexity.md:154-164 the witness rows `bang dim 4, code 27` and `bang dim 4, code 281`, six of seven measures split, their ANFs, one of 14 size-two groups splitting `bs`.
- complexity.md:194-209 `C = bs` on all 424 classes and on 50000 samples at `dim 5`; `bang dim 4, code 7128` the sole `s != bs` class, orbit 24; 216 samples at `s = 3, bs = 4`; `deg = s^2` at `bang dim 4, code 855` and `bang dim 4, code 1911`.
- complexity.md:248 `P[s(f) = dim]` is `1/2, 5/8, 69/128, 18253/32768`.
- complexity.md:269-278 no design with one bichromatic edge at `dim 3`, realised counts `0, 3, 4, 5, 6, 7, 8, 9, 12`; influence mean `1/2, 1, 3/2, 2`, variance `1/4, 1/4, 3/16, 1/8`.
- README.md:50 the `dim 4` witness bullet.
