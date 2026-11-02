# evasive-saturation

- The line-evasive set of [Fuhrer and Taranchuk](https://arxiv.org/abs/2403.18611) at `t = 3`, `q = 7`, built explicitly in `F_7^6`, and the test of the question in [planar-block](../../../notes/planar-block.md), section Four-term progressions, bullet The saturation of the evasive set: can a point outside it join without closing a four-term progression? None can.
- No dependency: one file, `std` only, `8` threads.

## THE METHOD

- `F_49 = F_7(sqrt 3)`, `F_343 = F_7(cbrt 2)`: `3` is a non-residue and `2` a non-cube mod `7`, so both polynomials are irreducible; `N_2(a + b sqrt 3) = a^2 - 3b^2`, `N_3(a + b r + c r^2) = a^3 + 2b^3 + 4c^3 - 6abc` with `r^3 = 2`, checked against the determinant of the multiplication matrix, multiplicative, zero only at `0`, onto `F_7`.
- A point `v = (v_1, v_2, v_3)` of `F_7 x F_49 x F_343 = F_7^6`, coordinates `c_0 = v_1`, `(c_1, c_2) = v_2`, `(c_3, c_4, c_5) = v_3`; `P(v) = v_1 + N_2(v_2) + N_3(v_3)`; the set is the level set `S_u = {P = u}`, `u = 0` by default (the level sets are translates of each other along `v_1`).
- Points are packed four bits a coordinate in a `u32`; addition is nibble-wise with one mask fix, negation likewise; membership is a bit set over `2^24`.
- `lines`: one direction per class (first nonzero coordinate `1`, `19608` classes), one base point per line (that coordinate `0`), `329551656` lines, the histogram of `card(line n S)`.
- `saturate`: for every `b` outside `S`, the number of `d != 0` with `b + d, b + 2d, b + 3d` in `S` (`b` an end, each progression through `b` as an end counted once) and with `b - d, b + d, b + 2d` in `S` (`b` an inner term, likewise once); the point is addable exactly when both are `0`.
- `check`: an independent test that `S` holds no four-term progression with distinct terms, by digit arrays and base-`7` indices instead of packed nibbles and the bit set, over every `x` in `S` and every `d != 0`; then the negative control, `S` plus one outside point, which must hold some.

## THE CONTROLS

- `cargo test`: packed addition and negation against digit arithmetic on `20000` random pairs; the two norms are multiplicative, anisotropic and onto; the largest subset of `Z_7` free of four-term progressions has `4` points; `200000` random lines meet the level set in at most `3` points.
- The run asserts `card S = 16807`, the line count, no progression in `S`, and at least one in the spoiled set.

## RUN

- `cargo run --release --manifest-path research/progressions/lab/rs/evasive-saturation/Cargo.toml -- all` from the repo root, `9` seconds on `8` threads; `-- lines`, `-- saturate`, `-- check` run one part; a second argument picks the level `u`.
- Prints only, writes nothing.

## WITNESSES

- `card S_0 = 16807 = 7^5`; the `329551656` lines meet it in `0, 1, 2, 3` points `107699256, 147683109, 40639326, 33529965` times and never in `4`.
- Every one of the `100842` outside points is an end of at least `372` and an inner term of at least `372` four-term progressions whose other three terms lie in `S_0`; the counts take the values `372, 393, 407, 414, 456` (`16758, 33516, 16758, 33516, 294` points each), the same histogram for ends and for inner terms; addable points: `0`.
- The independent check finds `0` progressions in `S_0` and `1614` ordered ones after adding the outside point `(3, 2, 2, 2, 2, 2)`.
