# shifted-sums

- Walks the shifted sums `T_s = sum_{w in B_s} mu(6w + 1)` of [the Franel note](../../../notes/threshold-surrogate-repunits.md)'s repunit section, `B_s` the `2^s` integers below `3^s` whose base 3 digits lie in `{0,1}`, to `s = 26`.
- The walk takes `w` in ascending order, so `T_s` is the prefix sum `P_n` at `n = 2^s`; it records every sign change of the last nonzero sign of `P_n` for `n <= 2^26`.
- The control is the whole class: `M(x; 6, 1) = sum_{n <= x, n = 1 mod 6} mu(n)` at `x = 3^(s+1) - 2`, the top of the image at level `s`, to `3^16 - 2`, beside the drift `3/4 - Psi_(2,3)(x)`, `Psi_(2,3)` the count of 3-smooth numbers up to `x`: the term at `s = 0` of its Dirichlet series `(1/2)(1/(zeta(s)(1 - 2^-s)(1 - 3^-s)) + 1/(L(s, chi_-3)(1 + 2^-s)))` when each `M(x/(2^a 3^b))` is read as its constant `-2` and `1/(L(0, chi_-3)(1 + 2^0)) = 3/2`.

## METHOD

- PARI `moebius`, batched from Python into 8 `gp -q` processes; the walk is cut into chunks of `2^16` consecutive indices, chunk `j` holding `w = w(j) 3^16 + w(r)` with `w(i)` the base 3 reading of the binary digits of `i`.
- Pass one returns each chunk's sum and the minimum and maximum of its local prefix; a chunk can hold a sign change only if its prefix range reaches the side opposite the last nonzero sign, and pass two rewalks exactly those chunks from their true offsets.
- Self-checks, each from an independent source: `T_1` to `T_16` against the list printed on the note; `T_(t-1) - T_(t-2) + mu(3^t + 1)` against the note's endpoint readings `M_F(x_t; R_t)` at `t = 2` to `24`, which [restricted-franel](../restricted-franel/) walks over a different set; a plain sequential walk to `2^18` repeats `T_18`, the flip count and the last flip.

## RUN

- `uv run python research/franel/lab/py/shifted-sums/shifted_sums.py` from the repository root; append `walk` or `control` for one verb.
- `walk` 83 s on 8 cores (449 s of CPU), `control` 7.7 s on one; `gp` on the path.

## WITNESSES

- [The Franel note](../../../notes/threshold-surrogate-repunits.md), the repunits: `T_s` at `s = 17` to `26` reading `-88, 121, 665, 1523, 1857, 1450, 4479, 4939, 1417, 12641`, negative at `s = 2` to `4` and `6` to `17`, zero at `s = 1, 5`, positive at `s = 18` to `26`, `abs T_s / 2^(s/2)` largest at `1.8562` at `s = 13`.
- The prefix walk: 132 sign changes to `n = 2^26`, the last at `n = 238418`, element `6w + 1 = 1128943015`, never negative from there to the cut and `0` at `n = 238419`.
- The control: `M(x; 6, 1)` negative at `s = 2` to `13`, `8` at `s = 14`, `265` at `s = 15`, 3996 sign changes below `3^16`; the drift reads `-173.25` at `s = 13` against `-169` and `-223.25` at `s = 15` against `265`.
