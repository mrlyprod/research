# slice-ladder-controls

- Checks the anti-diagonal profile identity `P_{A (x) B}(t) = P_A(t^side_B) P_B(t)` on every base-2 word of length 2 and 3, counting mismatches.
- Checks that the level-1 central diagonal slice of the base-3 Menger analog equals the vertex count of the cube's central cross-section, `dim 2..14`.
- The slice side is enumerated over all `3^dim` digit vectors and the vertex side is binomial, so the two sides share no formula.
- Prints the `dim 4` central diagonal census at levels 1..6, the order-2 recurrence fitted on its first four terms and checked on all six, its dominant root and its slice dimension.
- Prints the `3, 5, 7, 9, 11` staircase dimensions at rung `n = 1..5`, assuming fill `base^2 - ((base-1)/2)^2` for each base.
- Domains run are the page's own: words of length 2 and 3 over the 15 nonempty 2x2 tiles, `dim 2..14`, levels 1..6, rungs `n = 1..5`. Whole run is about three seconds.

## RUN

```
uv run python research/cuts/lab/py/slice-ladder-controls/controls.py
```

## WITNESSES

- `slices.md:118` - the profile identity, zero mismatches over all words of length 2 and 3.
- the claims line - the same identity, zero mismatches.
- the claims line - level-1 slice counts `2, 6, 6, 30, 20, 140, 70` at `dim 2..8`.
- `cuts.md:324` - the same counts, and the level-1 slice identified with the vertex set.
- `cuts.md:330` - the identity at `dim 2..8` and the `dim 2..14` ladder.
- `shared/sequences.md:96` - `6, 132, 1848, 29040, 441408, 6772128`, `a(n) = 11a(n-1) + 66a(n-2)`, root `(11 + sqrt(385))/2`, dimension `2.483635500`.
- the claims line - the `dim 4` characteristic polynomial `x^2 - 11x - 66` and `rho_4 = 15.310708`.
- `shared/sequences.md:97` - `2, 6, 6, 30, 20, 140, 70, 630, 252, 2772, 924, 12012, 3432`.
- `dimensions.md:296-297` - `1.892789261`, `1.892315261`, `1.893034267`, `1.894190425`, `1.895495742`.
- The run states its carpet fill assumption before any staircase number.
