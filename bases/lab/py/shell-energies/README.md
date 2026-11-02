# Shell Energies

- The three exact shell relations behind the base-2 design energy `Z_c(s)`, as integer identities on the Dirichlet coefficients of `r2`.
- Every parity class of `Z^2` counted by direct lattice enumeration to norm `n = 200000`, `r2` checked against `4 (d1 - d3)`.
- Doubling gives `r2_ee(n) = r2(n/4)`; the 45 degree rotation gives `r2_oo(n) = r2_mix(n/2)`; the coordinate swap gives `r2_eo = r2_oe`.
- The rotation also run as a map: every odd-odd point of every even norm to `n = 6000` sent to `((i+j)/2, (i-j)/2)`, the image required to be mixed parity of half the norm and the map required bijective onto that class.
- The four linear relations solved symbolically to recover `Q_c(t) = a_ee t^2 + a_oo t (1-t) + (a_eo + a_oe)(1-t)/2`.
- The closed form `Z_c(s) = 4 zeta(s) beta(s) Q_c(2^-s)` against truncated lattice sums at `s = 2, 3, 4` for all fifteen nonempty designs, the `s = 2` comparison adding the leading tail `k pi / (4N)`.
- `beta(s)` taken from the Hurwitz zeta and self-tested against Catalan and `pi^3/32`.

## RUN

- `uv run python research/bases/lab/py/shell-energies/shell.py`
- Under two seconds on one core; prints only, writes nothing.

## WITNESSES

- bases.md:109 `S_ee = t^2 * S` as `r2_ee(n) = r2(n/4)`, 0 mismatches to `n = 200000`.
- bases.md:110-111 `S_oo = t * S_mix` as `r2_oo(n) = r2_mix(n/2)`, 0 mismatches to `n = 200000`, with the rotation bijective on 4716 odd-odd points to `n = 6000` and 0 faults.
- bases.md:112 `S_eo = S_oe` as `r2_eo(n) = r2_oe(n)`, 0 mismatches to `n = 200000`.
- bases.md:113-115 the same `n = 200000` domain, and the worst gap `9.4e-10` at `s = 2`; beyond the page, `3.9e-11` at `s = 3` and `1.3e-16` at `s = 4`.
