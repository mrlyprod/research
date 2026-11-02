# Life Census

- Reads the `2^18` outer-totalistic life-like rules on the Moore neighborhood as designs on the parity cube `{0,1}^9`, and computes fill, Langton `lambda`, `GF(2)` degree and genus for all of them; the Walsh level sums are computed for `B3/S23` alone.
- A rule is a birth set `B` and a survive set `S` inside `{0..8}`, and as a design it is the dim 9 Boolean function `f(c, n) = [c = 0][|n| in B] + [c = 1][|n| in S]` on the centre `c` and the eight outer bits `n`.
- The code is the 512-bit corner bitmask in the crate's row-major corner order, corner `i` carrying coordinate `j` equal to `(i / 2^(8-j)) mod 2`, with the centre as the first coordinate; the code under the 3x3 block order, centre in the middle, is printed beside it, the two being a permutation of variables apart and so equal in fill, degree and genus.
- Dimension is the core page's `log(fill) / log(side)` at level 1, so at base 2 it is `log2(fill) = dim + log2(lambda)` with `fill` the number of filled corners of `{0,1}^dim`.
- The Moore mask is rebuilt from the crate's definition, ones on a `3^dim` grid with the centre index cleared, and compared with the level-1 side-3 tile of the design "void iff every coordinate is odd" at dim 1, 2, 3.
- That comparison is a theorem in every dimension: at side 3 the coordinates are `0, 1, 2`, the only odd value is `1`, so the only cell with every coordinate odd is the centre `(1, ..., 1)`; hence the level-1 side-3 tile of "void iff every coordinate is odd" is the full `3^dim` block minus its centre, which is `moore(dim)`.
- Genus follows the core page: a design is isotropic when some member of its `B_dim` orbit is a level set of the popcount, which by the orbit law means `f(x)` depends only on `|x xor t|` for one of the 512 flips `t`, and axial when its support is a subcube. Both tests are orbit-wide, so a per-rule verdict is well defined.
- The fill histogram is computed twice: as the self-convolution of the subset-sum distribution of the binomials `{1, 8, 28, 56, 70, 56, 28, 8, 1}`, and by brute force over all `2^18` pairs.
- The axial test is computed twice as well: a closed form in `B` and `S` from the fill and the edge counts, and a brute subcube test on the support of every one of the `2^18` truth tables.
- The degree sweep runs a Mobius transform over every one of the `2^18` truth tables, and is cross-checked against the composition law `deg(B, S) = deg(B)` when `B = S` and `max(deg(B), 1 + deg(B xor S))` otherwise, which follows from `f = g_B + c (g_B xor g_S)` and `g_B xor g_S = [|n| in B xor S]`; the exception is real, since the empty rule `B/S` has degree `-1` where the unconditional form gives `0`.
- The Walsh level sums `Sigma_k`, the sums of `W(S) = sum_x (-1)^(S.x) f(x)` over the subsets `S` of size `k`, are checked against the generating identity `sum_k Sigma_k t^k = sum_x f(x) (1+t)^(9-|x|) (1-t)^|x|`, which is derived from the support alone and not from the transform it checks, and against Parseval, `sum_S W(S)^2 = 512 fill` for a `0/1` valued `f`.
- The parity rule `B1357/S02468` is the sum mod 2 of all nine Moore cells, generating polynomial `(1 + x + x^2)(1 + y + y^2) = k(x) k(y)` with `k` the rule 150 kernel; the study evolves it from one seed to `t = 64` on a side `2t+3` grid with constant-0 boundary and compares every slice cell for cell with the outer product of two rule 150 rows.
- The named replicator `B1357/S1357` is the sum mod 2 of the eight outer cells only, kernel `k(x) k(y) - xy`, so at `t = 2^j` it is the outer product with the centre copy removed; that is checked too, and the two population sequences are checked against 20 published OEIS terms each.
- Structural laws are asserted and the study exits nonzero if one fails; headline counts are printed.

## RUN

- `uv run python research/automata/lab/py/life-census/life.py`
- Domain: every one of the `2^18` life-like rules, all `2^18` axial checks against the closed form, evolution to `t = 64` on a `131 x 131` grid; about three seconds, prints only, writes nothing.

## WITNESSES

- In every dimension the level-1 side-3 tile of the design "void iff every coordinate is odd" is the `3^dim` block minus its centre, which is the Moore mask, because `1` is the only odd coordinate value at side 3. (Proved.)
- The tile equals the crate's `moore(dim)` array cell for cell at dim 1, 2, 3, fills 2 of 3, 8 of 9 and 26 of 27; in the plane that design is `bang dim 2, code 7`, dimension `log(8)/log(3) = 1.892789`. (Verified.)
- `B3/S23` has fill `C(8,3) + C(8,2) + C(8,3) = 140` of 512, so `lambda = 140/512`. (Proved; Verified.)
- At base 2 a design's dimension is `log2(fill) = dim + log2(lambda)`, so at fixed dim it is a strictly increasing function of `lambda`, and it is not a function of `lambda` across dimensions; for `B3/S23`, dim 9 and `log2(140) = 7.129283`. (Proved; Verified.)
- `B3/S23` has `GF(2)` degree 8 with 184 monomials, and Walsh level sums `140, 308, -224, -896, -168, 840, 448, -224, -196, -28`, whose squares over the 512 subsets sum to `71680 = 512 * 140`. (Verified.)
- `B3/S23` is not a level set of the popcount under any of the 512 flips and its fill 140 is not a power of two, so it is neither isotropic nor axial, hence compound. (Verified.)
- Outer-totalistic means a level set on the eight outer axes for each value of the centre, so the family has `2^18` members and the totalistic rules, level sets of all nine, number `2^10`; all 1024 of them are isotropic. (Proved; Verified.)
- The fill histogram over the `2^18` rules is mirror symmetric about 256, takes 479 of the 513 values, peaks at fill 256 with 3270 rules, and is 1 at fill 0 and at fill 512. (Verified.)
- The 34 unreachable fills are `5, 6, 7, 13, 14, 15, 21, 22, 23, 41, 42, 43, 49, 50, 51, 69, 77` and their mirrors, so `lambda` has gaps near both ends. (Verified.)
- 165 rules share the fill 140 of `B3/S23`. (Verified.)
- Genus over the `2^18`: isotropic 2044, axial only 4, compound 260096. (Verified.)
- The degree histogram over the `2^18` is `-1:1, 0:1, 1:6, 2:24, 3:96, 4:384, 5:1536, 6:6144, 7:24576, 8:98304, 9:131072`. (Verified.)
- The composition law is `deg(B, S) = deg(B)` when `B = S`, and `max(deg(B), 1 + deg(B xor S))` otherwise; it holds on all `2^18`, and the empty rule `B/S`, of degree `-1`, is the case the exception is there for. (Proved; Verified.)
- A life-like rule is affine over `GF(2)` iff `f = alpha c + beta (n1 + ... + n8) + gamma`, so `B` and `S` are each empty, full, the odds or the evens and share `beta`; exactly 8 rules qualify, the four degenerate ones and `B1357/S1357`, `B1357/S02468`, `B02468/S1357`, `B02468/S02468`. (Proved; Verified by the degree of all `2^18`.)
- Count sets of degree at most `d` number `2^(d+1)`, and with the composition law the rules of degree at most `d` number `2 * 4^d` for `d = 1..8`. (Verified; the product form Proved.)
- `B1357/S02468` is `k(x) k(y)` with `k` the rule 150 kernel, so every slice is the outer product of two rule 150 rows and the population is the rule 150 population squared, and at `t = 2^j` the pattern is nine copies of the seed at spacing `2^j`, checked at `t = 1, 2, 4, 8, 16, 32, 64` and cell for cell to `t = 64`. (Proved; Verified.)
- `B1357/S1357` is the same product minus the centre, so at `t = 2^j` it is the outer product with the centre copy removed, since the rule 150 row has a live centre at `t = 2^j`, eight copies, checked at `t = 1, 2, 4, 8, 16, 32, 64`. (Proved; Verified.)
- The three population sequences are OEIS `A071053` for rule 150, `A246035` for `B1357/S02468`, and `A160239` for `B1357/S1357`, 20 terms each. (Verified.)
