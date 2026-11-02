# odd-base-slice-grammar

- Computes the centroid diagonal slice of the odd-base parity solid, the solid whose cells have at most one odd coordinate, as a tile census with no raster.
- The plane `x + y + z = 3*base^level/2` meets exactly three integer layers, so a cell is cut as a hexagon at layer offset `1` and as a triangle at offsets `0` and `2`.
- Peeling the top digit sends a window to one of three windows only, and complementing digits pairs the outer two, so the census closes on two tile symbols.
- Prints the `2x2` substitution matrix on `(hexagons, triangles)`, its two-term recurrence, the dominant root, `log_base` of that root, and `log_base(fill) - 1`.
- Cross-checks the tile recursion against a direct layer census, the full digit-sum distribution built by convolution with no tile reduction.
- Also runs the middle-digit solid, the cells with at most one coordinate equal to `(base-1)/2`, which is the other reading of the base-3 sponge.
- Domain: matrices and dimensions at base 3, 5, 7, 9; closed forms at odd base 3..21; the `mod 4` side test at odd base 3..401 by exact rational comparison of the dominant root against `fill/base`; direct layer cross-check to level 6 at base 3, 5 at base 5, 4 at base 7, 9.

## RUN

- `uv run python research/spectra/lab/py/odd-base-slice-grammar/slice_grammar.py`

## WITNESSES

- `spectra.md:26` the four rules `x9 -12`, `x11 +62`, `x42 -288`, `x28 +693`
- `spectra.md:26` dimensions `1.8184 / 1.6869 / 1.8026 / 1.7204`
- `spectra.md:26` `log_base(fill) - 1 = 1.7268 / 1.7304 / 1.7430 / 1.7544`
- `spectra.md:28` the base-3 target `[[6,1],[6,3]]` and hexagons `1, 6, 42, 306, 2250, 16578, 122202`
- `spectra.md:29` `[[7,3],[30,4]]`, `[[30,3],[24,12]]`, `[[19,9],[96,9]]`
- `spectra.md:30` both closed forms reproduce every census matrix at base 3..21
- `spectra.md:37` `112` of `125`, `dim_slice = 1.960651` against `1.931768`, excess `+2.888e-02`
- `spectra.md:37` `81` of `125`, `1.6869` against `1.7304`
- `spectra.md:38` the `mod 4` side holds at every odd base up to `401`
