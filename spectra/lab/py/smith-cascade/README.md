# smith-cascade

- Builds the even-convention carry matrix `M_even` at base 3, the odd carry parameter `D = 2R + 1`, size `n = R + 1`, from the digit polynomial `P = (sum_k C(D-1,k) x^(2k))(1 + D x + x^2)`, entry `P[c + D - 3c'] + P[-c + D - 3c']`.
- Computes its 2-adic Smith divisor valuations `a_i` by min-valuation pivoting modulo `2^256`, at every odd `D = 3..511`.
- Reads off `v_2 = sum a_i`, the layers `L_j = #{a_i >= j}`, the nullity `L_1` against the Jacobsthal tent `min_t |D - t|/2 + 1`, `t in {2J(k)+1, 2J(k)+3}`, the excess `X = v_2 - L_1`, and `a_max`.
- Tabulates per octave `[2^k, 2^(k+1))` the maxima of `L_2, L_3, L_4, L_5, a_max, X` and the slack `v_2 - ceil(n/3)`, and the min-of-cones shape of layers 1, 2, 3.
- Cross-checks `v_2` against an exact Bareiss determinant at `D <= 61`, the three largest rows at precision 1024, and the direct pencil `det(fill I - 3 M_even)` at `D = 7`.

## RUN

`uv run python research/spectra/lab/py/smith-cascade/smith_cascade.py`

About three minutes. Domain is the page domain, odd `D = 5..511`, 254 rows.

## WITNESSES

- research/claims/:225 254 rows; octave maxima `L_2 = J(k-2)`, `L_3 = J(k-4)`, `a_max = floor(log_2 D) + 4`; `L_4 <= 1`; cones at layers 1, 2, 3; `X` octave maxima `6, 6, 7, 8, 10, 17, 28`; `v_2 <= ceil(n/3) + 9` with slack `5, 5, 5, 7, 7, 9, 9`, extremal `255, 257`; `v_2 <= n` at odd `D >= 9`, equality `9, 15`, violations `5, 7`; class `D = 1 mod 6` 84 rows `13..511`, `v_2 <= n - 3`, ratio `7/18` at `D = 19`; `95 < 99` at `D = 511`
- research/claims/:281 tent law 255/255 at odd `D = 3..511`; `max a_i <= 9` first fails at `D = 127`, 12 by 511; `#{a_i >= 3} = 1` first fails at `D = 175`, reaches 5; `#{a_i >= 2} <= 5` first fails at `D = 183`, reaches 21; `D = 7` profile `{0,0,0,7}`, `D = 5` profile `{0,0,4}`
- research/claims/:210 `D = 7` has `v_2(det) = 7 >= 6` and `det(fill I - 3M) != 0`

## NOTE

- The octave-2 rows `D = 5, 7` break both octave laws: `max L_2 = 1` against `J(0) = 0` and `a_max = 7` against `6`; the laws hold on octaves 3..8.
- Every non-spike divisor has `a_i <= 3`, so `L_4 <= 1` is the spike alone and `L_5 = #{a_i >= 5}` is 1, not 0, at 233 of 254 rows; `L_5 = 0` only if the spike is excluded.
- `v_2/(D-1)` at `D = 13` is `1/3`, not `7/18`; the slack 9 is also attained at `D = 511`.
