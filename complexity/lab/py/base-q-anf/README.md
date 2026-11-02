# base-q-anf

- Builds the base normal form of a design: the inverse Vandermonde on `0..base-1` applied axis by axis over `GF(base)`, and the integer Mobius transform along the product of chains.
- Proves both round-trip on every value table by the matrix identities `T E = E T = I` over `GF(base)` and `B M = I` over `Z`, which settles base 3 at `dim 3` without touching its `2^27` designs.
- Sweeps every base-2 design at `dim 2, 3, 4` against the classical XOR ANF and the signed subset real ANF, recording the variable reversal between the two labelings.
- Sweeps all 512 base-3 designs at `dim 2` for both round-trips and both degree histograms.
- Evaluates the four historical rules at base 3, `dim 2, 3`, and the Sierpinski carpet and Menger sponge as actual shapes.
- Shows the Vandermonde is never invertible for a composite `base`.
- Domain run: exhaustive at every `(base, dim)` named above; base 3 at `dim 3` rests on the identity alone. Whole run is 10 seconds on one core.

## RUN

```
uv run python research/complexity/lab/py/base-q-anf/anf.py
```

## WITNESSES

- `complexity.md:296` - at `base 2` the inverse Vandermonde is `[[1,0],[1,1]]`.
- `complexity.md:297-299` - the Vandermonde is not invertible mod `base` at `base 4, 6, 8, 9, 10, 12`, and is at every prime through 13.
- `complexity.md:305-306` - `T E = E T = I` at `base 3, dim 3` covers all `3^27` value tables, so all `2^27` designs.
- `complexity.md:308-311` - all `16, 256, 65536` base-2 designs at `dim 2, 3, 4` round-trip with `0` failures and match the XOR ANF with `0` differences after variable reversal.
- `complexity.md:312-315` - `512` base-3 designs at `dim 2`, both round-trips `0` failures, `GF(3)` degree histogram `-1: 1, 0: 1, 2: 24, 3: 144, 4: 342`, no degree 1.
- `complexity.md:326-329` - `GF(3)` degrees void `2, 4`, tree `2, 4`, carpet `3, 6`, net `4, 6`.
- `complexity.md:332-336` - Sierpinski carpet `8/9` cells, `GF(3)` degree `4`; the carpet row keeps `3/9`; the `dim 3` row keeps `4/27`; Menger sponge `20/27`, `GF(3)` degree `6`, the ceiling `dim (base-1)`.
- `complexity.md:341-342` - at base 3, `dim 2`, the tree has `GF(3)` degree `2` and integer degree `1`; the void `GF(3)` degree `2` and integer degree `4`.
- `complexity.md:646-647` - the provenance lines this study replaces.
