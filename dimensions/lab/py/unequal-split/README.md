# Unequal Split

- The patch: one child of side `1/2` in the corner of the unit square and five of side `1/3`, the letter drawn on the 6-grid; its open set condition as exact rationals, its `dimension` solving `2^(-s) + 5*3^(-s) = 1` by PARI, the level render's fill `(9 + 4k)^level` on side `6^level` to level 4, and the census of cells by size `2^-a 3^-b` against `C(a+b, a) k^b` to word length 7, and the level-2 render against the Kronecker square of the tile.
- Its complex dimensions: the zeros of `1 - 2^(-s) - k 3^(-s)` for `k = 1..27` in the box `Re [D_l - 1/4, D + 1/4]`, `Im [-60, 60]`, seeded by PARI `polroots` of the lattice approximant `1 - x^41 - k x^65` and polished by Newton, the count certified by the winding number of the boundary with a Lipschitz margin; the first ones at `k = 5`; the printed root of the 2-3 nonlattice equation against its lattice approximant; the roots at `D + 2 pi i q/ln 2` for the convergents `p/q` of `log2 3` at 60 digits against the law `P Q (ln 2)^2 theta^2/(2 f'(D)^3)`.
- The count of cells `N(r) = #{cells of side >= r}` as exact integers to `r = e^(-300)`, the sizes sorted by exact integer comparison, the renewal identity at every breakpoint, `N(r) r^D` against its limit `1/(D f'(D))` over windows with its swing, the carpet as the lattice control, the folded variance at `ln 2`, `ln 3`, `ln 6`, and the periodogram of `ln N(e^(-u)) - D u` against the complex dimensions and their residues.
- The spin mass `mu(B(c, r))` of the natural measure about three corner fixed points, enclosed by cells: the identity `M(r) = ratio^D M(r/ratio)`, the ripple's swing against its enclosure, and the folded variance at `ln 2`, `ln 3`, `ln 6`.

## RUN

- `uv run python research/dimensions/lab/py/unequal-split/unequal_split.py patch`, about 1 second.
- `uv run python research/dimensions/lab/py/unequal-split/unequal_split.py poles`, about 5 seconds.
- `uv run python research/dimensions/lab/py/unequal-split/unequal_split.py count`, about 1 second.
- `uv run python research/dimensions/lab/py/unequal-split/unequal_split.py ripple`, about 1 second.
- Every verb prints only and writes nothing; PARI runs through `gp -q`.

## WITNESSES

- staircase-schedules.md The unequal split, The patch: the open set condition, `dimension = 1.778602507`, the render's fill `29^level` and `log(29)/log(6) = 1.879323585`, the census of cells (`patch`).
- staircase-schedules.md The unequal split, Its complex dimensions: 21 zeros and winding number 21 at every `k = 1..27`, the margins of at least 161, the gap `0.0116` from `D_l`, the floored offsets `0.470, 0.058, 0.461`, the first eleven at `k = 5`, the root of the 2-3 nonlattice equation, the convergent law to 7 digits at `q >= 665` (`poles`).
- staircase-schedules.md The unequal split, The count of cells: 59448 sizes, the renewal identity, `1/(D f'(D)) = 0.573459971`, swing times `sqrt(U)` between `1.486` and `1.546` on the six printed windows, folds `0.003, 0.005, 0.004` against the carpet's `0.999` at `ln 3`, the ten peaks within `0.001` of a complex dimension and within `1.3%` of its residue amplitude (`count`).
- staircase-schedules.md The unequal split, What the detectors see: the identity at 96 shifted radii per centre, swings `0.23906`, `0.15388`, `0.17843` against enclosures at most `8.2e-5`, at least 2661-fold, folds `1.000` at the centre's own period (`ripple`).
