# Code Factorisation

- Regenerates the composite and monoid material of [magic](../../../notes/magic.md): which tiles are Kronecker products, how a tile can factor in two shape-distinct ways, what the irreducibles are, and when two letters render one tile at one side.
- A tile is a 0/1 square array. A plane code at base `base` is the tile at side `base` whose bit `i` is the cell `(i / base, i mod base)`, so base 2 gives the 15 non-empty codes of [core](../../../../core/notes/core.md) and base 3 the 511 non-empty codes.
- The monoid is all non-empty tiles under the Kronecker product, graded by side, with the one filled cell as its unit. A word is a factorisation; a composite is a product.
- The whole study is exact integer comparison. Nothing here is a rearrangement singular value, a fit or a sample: the nearest-Kronecker-product literature solves an approximation problem, which is the wrong problem and the wrong arithmetic for a decision.
- The block test is the engine: `C` of side `N` cuts at `d | N` when every non-zero `d`-block of `C` is one tile `B`, and then `C = A (x) B` with `A` the 0/1 indicator of the non-zero blocks. It decides factorability at a named shape in `O(N^2)` and names both factors.
- Census universes are stated before they are counted and de-duplicated once: side 6 is the two shape images `15 x 511` and `511 x 15`, side 8 the two images `15 x 65535` and `65535 x 15`, side 9 the `511 x 511` base-3 pairs, side 12 the three ordered plane-code shapes.
- Counting at prime-power side runs the irreducible series `I = T/(1+T)` over the grading in `num-bigint`, and the same machinery is cross-validated in one dimension against exhaustive brute force at `N = 4, 8, 16, 9`.
- One dimension is not decoration: the diagonal embedding `D -> {(x,x)}` carries cut sets and irreducibility, so a line sweep is a tile sweep, and the exhaustive line sweep runs every non-empty subset of `{0..N-1}` at `N = 1..20`.
- Both readings of the one-cell statistic are printed, because a one-cell letter anywhere and a one-cell outer factor are different counts that happen to give 48 on the same 50 tiles; the study checks the three 48-sets are equal rather than assuming it.
- Every geometric figure is cross-checked against `mrly_math::bang::factory::create` and `mrly_core::Tensor::kron` cell for cell, on all 526 letters and all 7665 shape-(2,3) products, so two independently written renderers stand behind every count.
- Structural laws are asserted, not the headline counts: the study exits non-zero if the block test and the product census disagree, if `121` is not the rectangle set, if `3375` is not the triple-product set, if the criterion mismatches, if gcd closure fails, or if a witness stops being a witness.

## RUN

- `cargo run --release --manifest-path research/magic/lab/rs/code-factorisation/Cargo.toml`
- About four seconds; prints only, writes nothing, and holds two 983025-key sets at a time.

## WITNESSES

- magic.md the block test: 225 ordered base-2 pairs give 225 distinct side-4 tiles with largest preimage 1, and the block test over all 65535 non-empty masks returns the same set.
- magic.md the side-4 line: reducible 225, irreducible 65310, reducible share `0.343328%`.
- magic.md the side-6 census: images 7665 and 7665 with zero internal collisions, cross-shape 171, reducible `7665 + 7665 - 171 = 15159` of 68719476735, irreducible 68719461576, share `0.0000221%`.
- magic.md the anatomy of the 171: axis-separable 121 against 50, commutations 11 against 160 rewritings, not separable and not commuting 48, not separable and commuting 2.
- magic.md the fill tables: `1:36 2:64 3:32 4:16 6:14 12:8 36:1` over the 171 and `2:16 3:32 6:2` over the 50, outer-fill signature `(1,1):24 (1,2):8 (1,3):16 (2,3):2`.
- magic.md the one-cell readings: `0:23 1:8 2:140` and `0:2 2:48` for a letter anywhere, `0:23 1:60 2:88` and `0:2 1:24 2:24` for an outer factor, with the three 48-sets equal.
- magic.md the two-radix lines: the 11 subsets `{0} {1} {0,1} {2} {0,2} {3} {4} {5} {3,5} {4,5} {0..5}` and `121 = 11 x 11` by set equality.
- magic.md the commuting pairs: `(1,1) (2,4) (3,7) (4,64) (5,73) (6,84) (8,256) (9,273) (10,292) (12,448) (15,511)`, with base-2 codes 7, 11, 13, 14 unpartnered.
- magic.md the commutation criterion: `gcd(m-1,n-1)+1` singleton pairs and `gcd(m-1,n-1)+2` commuting pairs in one dimension at nine side pairs, exceeded at `(3,9)` at 7 against 4, and the side-15 witness `[3]{(1,1)} (x) [5]{(2,2)} = [15]{(7,7)} = [5]{(2,2)} (x) [3]{(1,1)}`.
- magic.md the side-8 line: 983025 and 983025 with intersection 3375, equal as a set to the triple products, reducible 1962675.
- magic.md the side-9 line: 261121 ordered base-3 pairs, 261121 distinct, zero collisions.
- magic.md the prime-power counts: 225, 1962675, 261121, 553402322215537199175 and `(2^25 - 1)^2 = 1125899839733761`, with the one-dimensional cross-check 9, 63, 1431, 49.
- magic.md the two witnesses: `[6]{(0,0),(2,2)}` factoring as `c1 (x) c257.q3` and `c17.q3 (x) c1`, and `[12]{(0,0),(3,3)}` factoring at profiles `(2 x 2 x 3)` and `(3 x 4)` with the side-4 letter irreducible and cut set `{1,2,3,4,12}`.
- magic.md the word census at side 12: 114975 words and 114975 composites per shape, pairwise 2565, 2565, 483, triple 483, union 339795, with `2565 = 15 x 171` and the side-12 witness in one image only.
- magic.md the cut-set sweeps: zero gcd-closure failures over every line to `N = 20` and over the 339795 side-12 composites, 132 and 2376 lcm-closure failures, first at `N = 12` with cut set `{1,2,3,4,12}`, and zero mismatches of the incomparable-divisor criterion.
- magic.md the render collisions: `2:480 3:15 4:1 5:1 6:1 7:1 8:1 9:1 12:1 18:1`, the carpet's unique side-3 partner `c495` of fill 8, and the side-9 fills 65, 72, 64 on three pairwise distinct tiles.
- `mrly_math::bang::factory::create` and `mrly_core::Tensor::kron`, the crate paths both renderers are checked against.
