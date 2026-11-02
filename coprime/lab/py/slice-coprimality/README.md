# slice-coprimality

- Coprimality of the points of a 3D parity design on a diagonal slice `x + y + z = s`.
- Checks `A_s = Sum_{d | s, d squarefree} mu(d) N_s^(d)` at every height, the base-prime peel `N_s^(base)(level) = [000 in P][base | s] N_{s/base}(level-1)`, and the prime-slice visibility bound.
- Measures the local price of a prime aggregated over the heights it divides, against `1/p^2`, and at `p = 2` against the parity-walk formula.
- Counts the central slice `s* = 3(base^level - 1)/2`: its visible density, its per-prime locals, the independence product, its peel, the `(9, -12)` recurrence and the peel ratio.
- Designs are the parity codes carpet `{000,100,010,001}` (`fill 20` at `base 3`), net `{110,101,011,111}` (`fill 7`), tree `{000,001}` (`fill 12`).

## RUN

- `uv run python research/coprime/lab/py/slice-coprimality/slices.py`
- One core, about three seconds.

## DOMAIN

- Enumerated: carpet `base 3` to `level 4`, `base 4` to `level 3`, `base 5` to `level 3`; net `base 3` to `level 7`; tree `base 3` to `level 7`.
- The `level 6` aggregated locals come from a residue transfer on `(Z/p)^3` and the central slice from a meet-in-the-middle over half levels; the transfer is checked against enumeration at `level 3`, the meet-in-the-middle against the height recursion at `level 2..7`.
- The central count and its peel run to `level 14` on the one-dimensional height recursion.

## WITNESSES

- `visibility-parity-and-slices.md:99,101` - zero Mobius mismatches at every height over the four design-base pairs, worst hidden count on a prime slice 3, peel exact at `base 3, 5`, no net point with `3 | gcd` on all `7^7` points.
- `visibility-parity-and-slices.md:103` - `0.040902` against `1/25`, `0.020446` against `1/49`, `p = 11, 13` still converging.
- `visibility-parity-and-slices.md:104` - `0.2850378 = 9121792/32002048`, walk formula and count equal on the integer.
- `visibility-parity-and-slices.md:105` - tree dichotomy, zero visible on even heights and zero even gcds on odd, 1.49 million points each side.
- `visibility-parity-and-slices.md:107,108` - `0.89216, 0.89776` against `0.57143, 0.61067, 0.65218`; `R_7 = 1093`; `base 5` reading `0.345, 0.492, 0.560`; independence product `0.64780, 0.89764, 0.55741` against measured `0.65218, 0.89776, 0.56006`.
- `visibility-parity-and-slices.md:109,110` - `6, 42, 306, 2250, 16578, 122202, 900882`, the peel `3, 27, 207, 1539, 11367, 83835, 618111`, both on the `(9, -12)` recurrence to `level 14`, ratio `0.093070331` against `(sqrt(33) - 5)/8 = 0.0930703308`.
- the claims line - the same numbers, plus `2^1092 = 1 mod 1093^2` and the tree dichotomy carried to `level 7`.
- `slices.md:415` - the pointer: `1/p^2` on the slice against `1/p^3` in the solid.
