# staircase-constants

- Computes the fill ratio of every base-2 design along the `row` word of odd sides `3, 5, ..., 2L+1`, against the closed form `R_L = (w/2^dim)^L prod_i Gamma(L+2-r_i)/Gamma(2-r_i) / (Gamma(L+3/2)/Gamma(3/2))^dim` over the roots `r_i` of the fill polynomial `P_F(n) = sum_j a_j n^(dim-j) (n-1)^j`, `a_j` the corners with `j` odd coordinates.
- Corner `i` is the binary digits of `i`, most significant first, as `mrly_math::bang::factory::code_to_corners` reads a code; the fill depends on the profile `a_j` alone, so 77 profiles cover all 273 nonempty codes at `dim 1, 2, 3`, `3 + 15 + 255`.
- `section_render` counts filled cells of every code at `dim 1..3` on odd sides `3..9` and even sides `2..8`, and the row word by `np.kron` at `L = 2..5`, `2..3`, `2` in `dim 1, 2, 3`, against `P_F` and `w/2^dim`.
- `section_exact` checks the Gamma identity against the exact integer product at `L = 1, 2, 3, 10, 100, 1000` at 50 digits.
- `section_asymptotic` sums `log1p` of the exact per-letter excess in float to `L = 10^6` and compares with `drift log L + log C + c1/L`, `drift = dim/2 - mean`, `c1 = dim/8 + drift - var/2`, mean and variance of the odd count over the corners.
- `constant` asserts `C` real and positive before taking its real part.
- `section_named` compares `C` with its closed form on named codes; `section_parity` compares the Gamma form at roots of unity for the two parity designs at `dim 2..6` with the independent zeta log series `exp(-+ sum_m (+-1)^m (lambda(dim m) - 1)/m)` at 50 digits, with the closed form where one is printed, and with the direct product; `section_family` takes the roots of the at-most-one-odd design from its fill polynomial at `dim 1..8`, so the 50-digit match tests only the Gamma recurrence and the direct sum at ten digits is the independent check.
- Values print at 12 significant digits and rounded to nine decimals, the form the notes copy.
- `section_pairing` checks `P(1-n) = (-1)^dim P(n)` iff the profile is a palindrome on 776 profiles at `dim 1..4`, the mirror identity `C_F C_F' = prod sin(pi r)/(4 r (1-r))` on all 77, and counts the profiles whose roots outside `0, 1/2, 1` pair under `r -> 1-r`.
- `section_fence` multiplies the renormalised factors along sides `3^k` and along the odd primes; `section_staircase` evaluates the stacked-prefix staircase dimension by Barnes G at the roots, against the direct fill sum, and the gap to `dim` against `log(2^dim/w)/(log 2n - 3/2)` to `n = 2 x 10^5`.

## RUN

```
uv run python research/magic/lab/py/staircase-constants/staircase.py
```

- From the repo root; one core, about 7 seconds: render 0.1s, roots 0.1s, exact 0.3s, asymptotic 1.0s, pairing 0.3s, staircase 2.8s, the rest under 0.2s each.

## WITNESSES

- `magic.md` The staircase: the Gamma identity on 77 profiles, worst log gap `4.4e-47`; the direct sum to `L = 10^6` against `drift log L + log C + c1/L`, worst residual `1.13e-12`, residual times `L^2` flat at `1.12` from `L = 10^4`.
- `magic.md` The staircase: renormalised limits along `3^k`, `1.56493401857`, `1.31484053105`, `1.87429848245` at `dim 1 code 1`, `dim 2 code 7`, `dim 3 code 23`; the five staircase dimensions `1.892789261, 1.892315261, 1.893034267, 1.894190425, 1.895495742` from Barnes G; Barnes G against the direct sum on 77 profiles, gap `2.1e-50`; `n` times the rate error `-0.75, -0.71, -0.70, -0.69, -0.69` at `n = 10^2 .. 2 x 10^5`, `dim 2 code 7`, and `-0.41 .. -0.36` on the row word.
- `pi.md` Pi on the staircase: the parity constants `0.785398163`, `1.254589239`, `0.948815486`, `1.052420668`, `0.985352084`, `1.014708181`, `0.995477949`, `1.004525474`, `0.998553027`, `1.001447181` at `dim 2..6`, each against the zeta log series at 50 digits, gap below `1e-50`, and its closed form where one is printed, gap below `1e-50`.
- `pi.md` Pi on the staircase: the named constants and the mirror products `9 sqrt3 pi/64`, `sqrt2 pi^2/24`, `pi/4`; 18 flip-closed codes against 46 palindromic codes at `dim 2, 3`; 37 of 77 profiles reduce by reflection, 13 of them palindromic.
