# pascal-shear

- The level-`level` cell set of the 2D design with code `7`, tile `[[1, 1], [1, 0]]`, is `{(i, j) in [0, 2^level)^2 : i AND j = 0}`.
- Builds that set by literal Kronecker powers for `level = 1..9` and counts `3^level` cells, then recounts `3^level` by the digit sum `Sum_i 2^(level - popcount(i))` for `level = 1..14`.
- Checks Kummer directly: `v2(C(i+j, i))` against the base-2 carry count, on exact binomials for `0 <= i, j < 128`.
- Rebuilds Pascal mod 2 from the additive recurrence alone, rows `0..1023`, and compares entry `k` of row `n` against `k AND (n - k) = 0`.
- Checks the shear `(i, j) -> (i, i + j)`, determinant `1`, as a set bijection onto the odd Pascal entries at `level = 6`.
- Confirms row sums equal `2^popcount(n)` to `n = 1023`, and that inside a level-8 square the antidiagonal count is `2^popcount(n)` below `2^level` and strictly smaller above.
- Compares both OEIS b-files against the recomputed triangle, index column included, `A047999` from the copy kept here and both from the live files when the network answers.

## RUN

- `uv run python research/cuts/lab/py/pascal-shear/pascal_shear.py`
- Full domain, no arguments; about two seconds plus the two b-file reads.
- `a047999.txt` is the whole b-file, `10585` terms. `b001316` has `50001` terms and is too large to keep here, so that count needs the network; without one the run still checks `A047999` and exits `0`.

## WITNESSES

- cuts.md:181-182: level sets to `level = 9` with `3^level` cells each, `3, 9, 27, ..., 19683`, exact binomials for `i, j < 128`, `16384` of them with `0` faults.
- cuts.md:182-183: Pascal mod 2 from the additive recurrence to row `1023`, `524800` entries, `0` mismatched cells.
- cuts.md:183-184: the shear checked as a bijection, `729` cells onto `729` odd entries at `level = 6`.
- cuts.md The same page, one dimension down: `50001` terms of `A001316` against `2^popcount(n)`, `0` differences, only when the live b-file read succeeds.
- cuts.md:185-186: `A047999` is checked against the triangle with `0` differences, but at `10585` terms, rows `0..144`, not the `8256` the line says; `8256` is a `128`-row prefix, and REFS.md:176 carries the same `10585`.
- REFS.md:175: `b001316`, `50001` terms. REFS.md:176: `b047999`, `10585` terms, rows `0..144`.
- README.md:49: `A047999` with antidiagonal populations `A001316`, both b-files checked term for term.
