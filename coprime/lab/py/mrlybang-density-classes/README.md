# mrlybang-density-classes

- A mrlybang code is a set `P` of parity corners in `{0,1}^3`; at base `base` it fills the design `F_base(P) = { v in {0..base-1}^3 : v mod 2 in P }`.
- Computes the exact rational part `delta * zeta(3)` of the coprime density of every one of the 255 nonempty codes at every base, as a Mobius bracket times Euler factors.
- Checks the even-base band `(8/7)(1 - t/|P|)` at even `base <= 40` and the same band in `dim 2`, and the odd-base self-similarity `fill_e(base) = fill_P(base/e)` at odd `base <= 75`, both domains the widest the pages claim.
- Sorts every code by its mod-2 difference span `H`, groups the codes by weight enumerator, takes the design lattice rank and index at `base 2, 3, 4, 6, 8, 10`, and finds the codes whose even value equals their odd limit.
- Measures the visible fraction by exact enumeration of all `fill^level` points at that level, against the limit and against the exact finite-level walk factor `1 - (1/8) Sum_t lambda_t^level`; every enumerated case is the level the page measured.

## RUN

```
uv run python research/coprime/lab/py/mrlybang-density-classes/classes.py
```

Whole run is under half a minute.

## WITNESSES

- `visibility-parity-and-slices.md:60` `visibility-parity-and-slices.md:63` - lattice index in `{1, 2, 4, 8}` at even `base 4, 6, 8, 10`; `151` codes full rank at `base 2`, index `1` or `2`, the smallest of size `4`.
- `visibility-parity-and-slices.md:61` `visibility-parity-and-slices.md:62` - nine band values `0, 4/7, 16/21, 6/7, 32/35, 20/21, 48/49, 1, 8/7`, zero exceptions to `(8/7)(1 - t/|P|)`; carpet `0.712853` and `0.709137`, net `0.951771`, against `0.713063` and `0.950751`.
- `visibility-parity-and-slices.md:64` - the `dim 2` parity carpet at every even base, `delta * zeta(2) = 8/9`, `delta = 0.5403796460924681 = 16/(3 pi^2)`.
- `visibility-parity-and-slices.md:66` `visibility-parity-and-slices.md:67` - zero exceptions to `fill_e(base) = fill_P(base/e)`; `63` weight classes, `6` of them mixing regimes, `4` spanning against not and `2` corrected against no-limit.
- `visibility-parity-and-slices.md:68` - twins `0.707501` and `0.711167` at `base 4`, splitting at `base 3` into `0.698204` against `153/182` and `0.783436` against `51/52`.
- `visibility-parity-and-slices.md:70` - `7` codes of lattice rank `<= 1` at `base 3`, `37` codes with a coordinate pinned odd, `{110,100}` measuring `0.972222 0.994084 0.999979` against the trichotomy's `0.987319`.
- `visibility-parity-and-slices.md:71` `visibility-parity-and-slices.md:72` `visibility-parity-and-slices.md:74` `visibility-parity-and-slices.md:78` - class counts `149`, `43`, `63`; the subgroup codes `1 3 5 9 15 17 33 51 65 85 105 129 153 165 195 255`, exactly the parity-stable ones; `239` codes with two values, where the page says `240` and counts the empty code.
- `visibility-parity-and-slices.md:73` - void `base 3, level 8` measures `0.381949` against limit `0.438808` and finite level `0.380043`; tree `0.451821` against `0.452521`, tree `base 5` `0.468340` against `0.468560`.
- `visibility-parity-and-slices.md:75` - `{111}` at `base 5` is exactly `0` at every even level and `0.925781 0.952881 0.957167` against `0.958419`; axes `0.987338 / 0.739563` against `0.987319 / 0.740489`.
- `visibility-parity-and-slices.md:83` to `visibility-parity-and-slices.md:92` - the family table `513/520 27/26 99/182 48/91` at `base 3` up to `200981/201096 1331/1330 9559/16758 16456/28861` at `base 11`; carpet `0.98654 0.99562 0.99810 0.98959 0.99943` with `fill_P(base) = 20, 81, 208, 425, 756`; net `B(9) = 297/304 = 1 - 7/304` and no odd `base <= 81` below `1`, the page's `20/425` and `1/208` being the carpet's fill, not the net's.
- the claims lines - the band, the lattice indices, the `dim 2` constant, the enumerations `0.712853, 0.951771, 0.709137`; the trichotomy counts at odd `base <= 75`, the axes pair, the `{111}` zeros.
