# Design Census

- The three design censuses in this tree, rerun in full.
- Coprime: every design with `fill >= 2` at base 2 in dimensions 2 and 3 and at base 3 in dimension 2, the Menger sponge and two base-6 samples; the measured coprime density at the deepest level inside a budget of 200000 points against the predicted `B(F) / zeta(dim) * prod (1 - p^-dim)^-1` for `F` the filled digits, the exact base-local identity at `level 4`, the spanning index, the split of the spanning lines by `log_base(fill)` against one, and how many lines widen their gap to the prediction at some level.
- Fill: one row per cube-group orbit at bases 2 and 3 in dimensions 2 and 3, base 3 dimension 3 by the four named designs only; level-1 fill and void at sides 1..12, the fill quasi-polynomial, orbit size, popcount and three algebraic degrees, every row rendered cell by cell and checked against the closed form and the level law.
- Fill classes: the parity-rule generator against the popcount-profile closed form over every base-2 design at `dim 2, 3`, the distinct fill sequences at `dim 1..4` against `Prod_k (1 + C(dim, k))`, the 16 published terms of A129824, and the per-dimension counts table with A000616 by Burnside to `dim 6`.
- The base-2 cube group counted two ways at `dim 1..4`: orbit walk, Burnside average.
- Boundary: the base-7 coset design `v_1 + v_2 = 1 mod 3`, index 3, to `level 7` against its level-corrected constant; the base-2 `dim 4` coprime census at `level 4` over the least code of each of the 402 orbits; the first terms of one design by code.

## RUN

- `cargo run --release --manifest-path research/coprime/lab/rs/design-census/Cargo.toml`
- Half a minute on 8 cores; writes `sequences.csv` and `counts.csv` beside this file.
- `... -- limits` runs the base-7 coset design in three seconds, `... -- dim4` the `dim 4` census in a second, `... -- terms 2 2 14` the terms of one design to twenty million points.

## WITNESSES

- coprime.md:23 the base-7 coset design: limits `0.698166, 0.698166, 0.465444` along `level mod 3 = 1, 2, 0`, measured `0.698172, 0.465147, 0.698172` at `level 5, 6, 7`, gaps to the corrected constant `5.3e-06, -3.0e-04, 5.4e-06`.
- coprime.md:48 the or-triangle, code 14, terms `3, 6, 22, 58, 200, 576` and `0.811236` at `level 15` against `8/Pi^2`.
- research/claims/ the base-2 `dim 4` row: 65536 designs, 402 orbits, 384 group elements, 400 at `fill >= 2`, 336 spanning with `fill > 2`, 189 distinct `(A(1)..A(4))`, 87 collisions covering 298 classes, 10 classes on `(4, 16, 88, 436)` and on `(5, 25, 165, 985)`.
- the-master-density.md:14 the exact identity at `level 4` on all 763 lines, zero failures.
- the-master-density.md:20 code 34376528265 at base 6 with `fill = 8`, `fill_2 = 3`, `fill_3 = 2`, `fill_6 = 1`, bracket `1/2`.
- coprime.md:35-44 the family table 11/5/3, 247/149/27, 502/365/175, 1/1/1, 2/2/2, total 763 designs, 522 spanning, zero flagged; the sponge at `level 4` measuring 0.7719 against 0.8207; `fill 3` reaching `level 11`.
- coprime.md:45-46 219 distinct term-vectors among the 502 base-3 designs and 175 spanning; 490 spanning lines above dimension one, 32 at exactly one, four more at index 3, none below; 426 lines whose gap to the prediction widens at some level.
- core.md:47-51 fill at `side 3` of 8 of 9, 20 of 27 and 12 of 27, with `20^2 = 400` at level 2.
- core.md:81-86 and core.md:156-163 carpet and net share one orbit at base 2, rep code 7 in dimension 2 and 23 in dimension 3 so only the net label survives, and the six 2D classes with orbit sizes 1, 4, 4, 2, 4, 1 and degrees -1, 2, 1, 1, 2, 0.
- core.md:96-106 base-2 orbit counts 3, 6, 22, 402 agreeing three ways at `dim 1..4`.
- core.md:243-248 the census this study writes: representatives, orbit sizes, fill polynomials, three degrees, per-dimension counts.
- method.md:109-112 58 designs across bases 2 and 3 and dimensions 2 and 3, rendered cell by cell, orbit counts against Burnside, a 59-line csv.
- method.md:113-117 763 designs measured against the predicted density, 522 spanning, 241 degenerate, none flagged.
- method.md:322-330 `bang dim 3, code 23` filling `4*k^3 - 3*k^2` at odd side `2k - 1` and reading 20 at `side 3`, `bang dim 3, code 15` and `bang dim 3, code 27` sharing one polynomial, 21 distinct polynomials over 22 classes and 64 over 256 designs.
- bijection.md:88-99 group orders 2, 8, 48, 384, orbit counts 22, 402, 1228158, 400507806843728, and 6 canonical designs in dimension 2, 22 in dimension 3; bijection.md:184 the base-3 `dim 3` Burnside 111618.
- shared/sequences.md:11 and shared/sequences.md:105-108 fill classes counted by A129824 with zero profile collisions, A000616 `3, 6, 22, 402, 1228158, 400507806843728`, all 16 published A129824 terms, 144 and 2304 generator checks with zero mismatches, and the A129824 column to `dim 8`.
