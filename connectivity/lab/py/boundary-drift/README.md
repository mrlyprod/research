# boundary-drift

- Regenerates the boundary law of [structure against noise](../../../notes/connectivity.md): the closed form of the boundary per cell `(2*dim*N - 2*E)/N` of a self-similar design at every level, from the tile alone.
- The tile gives `f` (cells), `A` (face-adjacent tile pairs) and `c` (positions on a face where the high face and the low face of the tile are both occupied, asserted equal on every axis); the closed form is `E_k = A (f^k - c^k)/(f - c)`, computed in exact rationals.
- Two builds share nothing: a Kronecker power with shifted-AND edge counts, and the digit rule (a point is in level `k` iff every base-`b` digit tuple is a tile cell) with set-lookup edge counts.
- The digit build also counts the pairs that cross copies at the top digit and asserts the count equals `A c^(k-1)`, the recursion step of the lemma, at every level from 2.
- Seven designs: gasket, diagonal, seven-of-eight, carpet, sponge, the base-3 `L` of code 79 and the base-4 `L`; each boundary per cell is asserted equal as an exact fraction to the closed form, and to the value the note prints (2.0082, 2.0027, 2.0009, 2.0003, 4.0000, 1.8542, 0.8633, 2.2560).
- The gasket is rebuilt with scipy component labelling: one piece, `E = N - 1`, perimeter `2N + 2` at levels 1 to 8; the carpet at levels 1 to 4 is one piece with perimeter strictly below `2N + 2`.
- The words section draws every word of length 1 to 4 over the 15 codes of base 2 (54240) and asserts `E(w) = sum_i f_1 ... f_(i-1) A_i c_(i+1) ... c_n` per axis, then counts the multisets whose orderings differ in `E`.
- The census renders every code of dimension 2 at base 2 (levels 1 to 4) and base 3 (levels 1 to 3), keeps those that are one piece with `E = N - 1` at every level, asserts there are 8 and 15, and asserts `2 + 2/f^k` on each.

## DOMAIN

- Kronecker levels: gasket, diagonal, seven-of-eight and carpet 1 to 8, sponge 1 to 5, code 79 1 to 6, base-4 `L` 1 to 5.
- Digit-rule levels: gasket and diagonal 1 to 8, seven-of-eight and carpet 1 to 5, sponge 1 to 4, code 79 1 to 5, base-4 `L` 1 to 4.
- Words: every word of length 1 to 4 over the 15 codes; multiset counts at length 2 over the 15 codes and the ten-code library, length 3 over the library.
- Nothing is sampled; the program prints this domain first and cuts nothing.

## RUN

- `uv run python -I research/connectivity/lab/py/boundary-drift/drift.py`
- About four seconds; prints only, writes nothing, exits nonzero on any failed check. Needs numpy and scipy.

## WITNESSES

- connectivity.md, section "The boundary law": gasket `E_k = 3^k - 1` and boundary `2 + 2/3^k`, reading 2.666667, 2.222222, 2.074074, 2.024691, 2.008230, 2.002743, 2.000914, 2.000305 at levels 1 to 8; perimeters 8, 20, 56, 164, 488, 1460, 4376, 13124.
- connectivity.md, section "The boundary law": carpet `E_k = (8/5)(8^k - 3^k)`, edges 8, 88, 776, 6424, 52040, 418264, 3351944, 26833048 and boundary 2.000000, 1.250000, 0.968750, 0.863281, 0.823730, 0.808899, 0.803337, 0.801251 at levels 1 to 8.
- connectivity.md, section "The boundary law": seven-of-eight edges 9, 90, 711, 5220, 37269, 263070, 1848051, 12956040 and sponge edges 24, 672, 14976, 311808, 6334464, read against `(9/4)(7^k - 3^k)` and `2 (20^k - 8^k)`.
- connectivity.md, section "The boundary law": code 79 `2 + 2/5^k` at levels 1 to 6 and the base-4 `L` `2 + 2/7^k` at levels 1 to 5; 8 tree codes at base 2 with fills 2 and 3, 15 at base 3 with fills 3 and 5.
- connectivity.md, section "The boundary law": 54240 words with zero mismatches; orderings differ in `E` for 78 of 105 multisets (length 2, 15 codes), 32 of 45 (length 2, library) and 188 of 210 (length 3, library).
