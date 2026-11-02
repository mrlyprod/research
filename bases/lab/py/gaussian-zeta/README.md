# gaussian-zeta

- Counts the fill of the diagonal design, corners `(0,0)` and `(1,1)`, cell by cell on the `side x side` grid and checks the fluctuation `rho(side) - 1/2` in exact rationals: `0` at even side, `1/(2 side^2)` at odd side, `-1/(2 side^2)` for the orbit mate.
- Reduces `Z(s) = (1/2) lambda(s+2)` exactly to `1/192` and `1/1920` at `s = 2, 4`, then matches `Z(2) = pi^4/192` and `Z(4) = pi^6/1920` at 90 digits with pi by Machin against the library pi and `lambda` by Euler-Maclaurin over the odd integers against `(1 - 2^-s) zeta(s)`.
- Sums the series straight off the counted fills, no closed form, and prints the gap at `side <= 49, 199, 999`.
- Checks the three parity shell identities on `Z^2` from shell counts, evaluates `Z_c(2) / (pi^2 G)` for all fifteen nonempty base-2 designs, and checks each against a truncated lattice sum with the `pi/(4R^2)` tail.

## RUN

`uv run python research/bases/lab/py/gaussian-zeta/gaussian_zeta.py`

About 3 s. Domain is the source domain: fluctuation exact to side 80, series to side 999, shells to norm `1000000`, lattice radius `6000`.

## WITNESSES

- pi.md:127-128 `Z(2) = pi^4/192 = 0.50733901580...` and `Z(4) = pi^6/1920 = 0.50072353832...`, both matched at 90 digits
- pi.md:132-136 exact to side 80, reductions `1/192` and `1/1920`, two pi routines agreeing to the last digit, two lambda routes at `1.3e-82` and `5.2e-87`, gap `8.3e-11` at `side <= 999`
- pi.md:139-142 orbit mate code 6, fluctuation `-1/(2 side^2)` on odd `side <= 40`
- bases.md:108-114 the three shell identities and the truncated lattice sums for all fifteen nonempty designs, a second pass over a wider domain than shell-energies, which is where the page number on those lines comes from
- bases.md:120 `Z_c(2) / (pi^2 G) = (a_ee + 3 a_oo + 6 a_eo + 6 a_oe) / 24`
- bases.md:124 fifteen rationals `1/24` to `2/3`
- bases.md:126 `13 pi^2 G / 24 = 4.8967847822`

## NOTE

- The `9.4e-10` at bases.md:114 is the shell-energies truncation at norm `200000`; at radius `6000` with the tail added the worst gap here is `1.2e-13`, so the two agree on the identity and differ only in domain.
- The 4D grading at bases.md:137-144 is not computed here.
