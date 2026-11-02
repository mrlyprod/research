# slice-at-large-side

- Computes the central hexagon `x + y + z = 3(N-1)/2` of the dim 3 parity designs at large odd side and on words of growing odd sides, from the carry block `S_b` of `notes/spectra.md`.
- A word `(b_1, .., b_L)` lists its letters coarsest first; its count is `(S_(b_L) .. S_(b_1))[0, 0]` and its ink is the count over the `(3N^2 + 1)/4` lattice points of the hexagon.
- `words` counts the slice by brute force over every lattice point of the plane, by the `3x3` carry automaton and by the `2x2` block, on twelve words for code 23 and four words for each of seven other codes, and checks that no carry leaves `abs(c) <= 1`.
- `letter` derives the block as polynomials in `b` per class of `b mod 4` by inclusion-exclusion on boxes, checks them against the digit polynomial at odd `b = 3..101` and against the closed forms of spectra, splits `S_b/b^2 = A + B/b + B_2/b^2`, checks the three-carry limit `f(c') q(c + c' + g)`, matches `A` to `[[3u/4, v/8], [3v/2, u/4]]`, and prints the Perron roots, the exact `rho_b`, the first-order coefficient `mu` and the side expansion `log(2 lambda) + (mu - 3/2)/b`.
- `drift` derives the drift exponents `1/4`, `1/4 + 1/(2 sqrt(7))` and `(2 + sqrt(2))/4` by first-order perturbation of the limit letter and of the limit pair in both orders, then runs four words to length `32000` in floats and prints the local exponent at `L -> 2L` for `L = 4000, 8000, 16000` and one Richardson step.
- `constants` runs the same four words at 40 digits to `L = 16384` and `16385`, divides the ink by `lambda^L L^gamma`, extrapolates by Neville in `1/L` at two depths and prints their agreement, then the blink ratio against its derived value `sqrt((5 sqrt(2) - 1)/3)`.
- `closed` checks `S_(4k-1) = (k/2)(k K_1 + K_0)`, the polynomial identity behind `W = (1 - 24z + 96z^2)^(-3/4) (1 - 6z, 12z)`, the series against the integer product at `L = 1..60`, and `Gamma(3/4)(1 + sqrt(3)) / (3 (sqrt(3) - 1)^(3/4))` against the extrapolated constant, printed to 25 digits.
- `designs` runs all 255 nonempty codes through their 63 weight classes: the symbolic block per class against the digit polynomial at odd `b` up to `61`, the limit letter against `[[3u/4, v/8], [3v/2, u/4]]`, the real roots of the two sign polynomials and of the case guard `s_01 s_10` (`s_00` where `s_01` vanishes), the codes whose block has a zero diagonal on a class, an exact sign check at every odd base up to `15`, the law `sgn(o - e)` at `3 mod 4` and `sgn(e - o)` at `1 mod 4`, and the automaton row sums of the `e = o` designs at odd `b = 3..41`.
- Exact arithmetic everywhere except `drift`, in floats, and `constants`, in mpmath at 40 digits. Runtimes: `words` 0.1s, `letter` 0.2s, `drift` 0.2s, `constants` 1.9s, `closed` 0.3s, `designs` 0.7s; all six 3.7s.

## RUN

```
uv run python research/spectra/lab/py/slice-at-large-side/large.py
uv run python research/spectra/lab/py/slice-at-large-side/large.py words
uv run python research/spectra/lab/py/slice-at-large-side/large.py letter
uv run python research/spectra/lab/py/slice-at-large-side/large.py drift
uv run python research/spectra/lab/py/slice-at-large-side/large.py constants
uv run python research/spectra/lab/py/slice-at-large-side/large.py closed
uv run python research/spectra/lab/py/slice-at-large-side/large.py designs
```

## WITNESSES

- `spectra.md` THE SLICE AT LARGE SIDE, The word product, the counts `60, 72, 2412, 2688, 300`: verb `words`, no mismatch.
- `spectra.md` THE SLICE AT LARGE SIDE, The letter at infinite side, the limits, the Perron roots, the exact roots and the side expansion: verb `letter`.
- `spectra.md` THE SLICE AT LARGE SIDE, The drift, the exponents and the local fits: verb `drift`.
- `spectra.md` THE SLICE AT LARGE SIDE, The constants, the closed form at sides `3 mod 4`: verb `closed`; the blink and the extrapolated constants: verb `constants`.
- `spectra.md` THE SLICE AT LARGE SIDE, Every parity design, the letter of every design, the ties and the sign law: verb `designs`, no mismatch and no violation.
