# design-averages

- Averages `abs(M_S(x))^2`, `M_S(x) = sum_(n <= x, n in S) mu(n)`, over the three flat families of [the Franel note](../../../notes/strict-set-and-average.md)'s section The average over designs: the translates of the level set mod `base^level`, the digit shifts, and the digit subsets.
- Computes each average three ways: member by member, through the pair form `sum mu(n) mu(m) P(n, m)`, and through the Fourier form, exactly in `Z[e(1/3)]` for the digit families at base 3 and by FFT for the translates.

## METHOD

- `check`: at base 3 with `fill = 2`, every `level <= 8` at `x = 3^level` and at `x = 100, 1000, 2000, 5000` (levels `5, 7, 8, 8`, the cell `2000` one above `level(2000) = 7`), the direct average equals the pair form and, for the digit families, the exact Fourier form as exact rationals, the translates' FFT form to relative `1e-9`; the mod `9` lower bound with weight `1/16` is tested exactly on the digit families and with `kappa_2` in floating point on the translates; the flat mean is counted at every `n <= 3^level`; the two lower bounds of the note (frequency `0`, and the `base` frequencies with weight `q = 1/4`) are tested as exact rationals; base 5 `{0,1}` to `level = 4` and base 4 `{0,2}` at `level = 3` repeat the pair form in all three families. The expectations are independent of each other: the direct average enumerates members, the pair form never forms a member, the Fourier form never forms a pair.
- `read`: the half plane edge `1 - log 2/(2 log 3) = 0.684535123`; then at base 3 with `fill = 2` and `x = 3^level`, `level = 4` to `13`, the ratio of the mean square to the mean mass `(2/3)^level x`, the share of the frequency `0`, and the share of the frequencies of level at most `floor((1 - alpha) level)`, for the digit subsets (exact) and the translates of `{0,1}` (FFT).

## RUN

- `uv run python research/franel/lab/py/design-averages/design_averages.py` from the repository root; append `check` or `read` for one verb.
- `check` 2.4 s, `read` 0.6 s, one core.

## WITNESSES

- The expansion, verified: `930830/6561` for the translates of `{0,1}` and `337520/2187` for the digit shifts at `x = 3^8`, `M(3^8) = -13`; 0 failures over every cell.
- The flat mean: 0 counts off `2^level`; the unwrapped translates hold `n = 1` and `n = 8` in `2` and `4` of `9` at `level = 2`; one digit set at every position holds `n = 1` in `2` of `3` designs and `n = 15` in none.
- The converse readings: low-level share `0.270590` at `level = 4` and `0.013826` at `level = 13` (digit subsets), `0.264550` and `0.015085` (translates), least `0.010210` and `0.013292` at `level = 8`; mean square over mean mass between `0.465892` and `0.617104` over `level = 6` to `13`.
