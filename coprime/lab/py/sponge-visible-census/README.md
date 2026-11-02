# sponge-visible-census

- Exact count `A(level)` of visible (gcd 1) nonzero points of the Menger sponge at that level, without enumerating the `20^level` points.
- Split at a cutoff `G`: points with `1 < gcd <= G` cancel exactly inside `Sum_{d <= G} mu(d) (T_d(n) - 1)`, where `T_d(level)` counts points divisible by `d` from one transfer pass on `(Z/d)^3`.
- Points with `gcd > G` are `g*y` with `y` primitive in a box of side `(3^level - 1)/g`, tested for membership digit by digit and weighted by `Sum_{d | g, d <= G} mu(d)`.
- Cross-checks: direct enumeration at level 1..6 matching `12, 270, 5916, 123504, 2538447, 51497040`, and `A(9)` recomputed at `G = 100` and `G = 150`, which split head and tail differently.
- Also prints the second-order gap `(delta*20^level - A(level)) / 12^level` with `12 = fill*lambda`.
- The bracket is derived, not quoted: `19` of the `20` digits are nonzero, times `27/26` for the prime `3`, gives `513/520`, so `delta = (513/520)/zeta(3) = 0.8207086195`.

## RUN

```
uv run python research/coprime/lab/py/sponge-visible-census/census.py 9
```

- From the repo root. Argument is the top level; default 9.
- Whole run takes about 1 min 24 s on one core, `A(9)` alone 23 s.
- Domain run is the full source domain: level 1..9, coordinates `0..19682`, cutoff ladder `G = round(3^(level/2))` capped at 120.

## WITNESSES

- the claims line - `A(7) = 1038074187`, `A(8) = 20860210527`, `A(9) = 418429711224`, and `G = 100` and `G = 150` agreeing on `A(9)`.
- the claims line - the gaps `0.347, 0.349, 0.344` in units of `12^level`.
- `visibility-parity-and-slices.md:84` - the base 3 carpet bracket `513/520`.
- `visibility-parity-and-slices.md:115` - the same three counts under "Counting without enumerating".
- `visibility-parity-and-slices.md:118` - the same three gaps.
