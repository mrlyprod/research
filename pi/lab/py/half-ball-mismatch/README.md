# Half Ball Mismatch

- The half-disk chord constant split into an integer and an area, and the two half-ball families whose parity blocks any repeat of that coincidence.
- `zerr.py`: the inner integral `Integral_(-1)^(1) (-a u + sqrt(1 - a^2 + a^2 u^2))^3 du = 2` for every `a`, symbolically and at 50 digits on 50 values of `a`, then `I_diam = 4` and `P = I_diam/(3 Area(H)^2) = 16/(3 Pi^2)`.
- `version_l.py`: the two-point family from its Wallis recurrences in exact rationals to ambient dimension `d = 11`, with an independent `10^7`-sample random-point check at `d = 2..7`.
- `version_h.py`: the `d`-point hyperplane family, derived here by an unoriented-normal Blaschke-Petkantschin reduction and integrated in closed form at `d = 2..7`, cross-checked by 60- and 80-digit quadrature and by a `10^8`-sample random-point estimate at `d = 3, 4, 5`.
- `mismatch.py`: every base-2 and base-3 design at dim 2 against every Version L value to `d = 24`.
- Beyond the page: Version H is exact here at `d = 6` and `d = 7`, and the dim 2 sweep is a record to `d = 24` rather than to `d = 11`.

## RUN

- `uv run python research/pi/lab/py/half-ball-mismatch/zerr.py`
- `uv run python research/pi/lab/py/half-ball-mismatch/version_l.py`
- `uv run python research/pi/lab/py/half-ball-mismatch/version_h.py`
- `uv run python research/pi/lab/py/half-ball-mismatch/mismatch.py`
- Under four minutes together, almost all of it the `10^8`-sample check in `version_h.py`. Prints only, writes nothing.

## WITNESSES

- pi.md:165 and pi.md:167 `16/(3*Pi^2) = 0.5403796460924681` and `1 - 16/(3*Pi^2) = 0.4596203539075319`.
- pi.md:169 and research/claims/:185 the inner integral is `2` for every `a`, `I_diam = 4`, `Area(H)^2 = Pi^2/4`, `P = 16/(3*Pi^2)`, three symbolic residuals exactly zero, 50 digits at 50 values of `a`.
- pi.md:171 and research/claims/:350 the parity law and the uniqueness at dim 2 and `d = 2`: 11 base-2 and 502 base-3 designs, base-2 numerators `4, 16/3, 6, 8`, exactly one match, `16/3` at `d = 2`, carried by 3 designs.
- pi.md:173-180 Version L `16/(3*Pi^2) = 0.5403796`, `3/8`, `128/(45*Pi^2) = 0.2882025`, `15/64`, `1024/(525*Pi^2) = 0.1976246`, `175/1024`.
- pi.md:182-187 Version H `1 - 16/(3*Pi^2)`, `4 - 19845*Pi/16384 = 0.1947689081`, `4 - 549978112/(14189175*Pi^2) = 0.0727502984`, `16 - 178919214166875*Pi/35184372088832 = 0.0244047160`.
- pi.md:189 and research/claims/:351 exact rationals to `d = 11`, and the 60- against 80-digit agreements `2.3e-62`, `7.2e-64`, `1.5e-63` at `d = 3, 4, 5`.
- Not regenerated: the three `10^8`-sample deviations `2.1e-05`, `9.7e-06`, `4.7e-06` at pi.md:189 and research/claims/:351 belong to one sample. This study draws its own and reports the deviation against its own one sigma.
