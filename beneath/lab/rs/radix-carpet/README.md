# radix-carpet

- The Sierpinski carpet as a radix design at base `3` on `Z[i]`, digits the box `{a + c i : 0 <= a, c < 3}` minus its centre `1 + i`, built by `mrly_math::radix::tile(3, 495)` of [beneath](../../../notes/radix-dial.md).
- Prints the digit count, the fill and the distinct-point count at level `2`, the similarity dimension, whether the digits are the canonical representatives, and the code of the same digit set read two ways.
- The box code `495` reads bit `r m + c` as the cell at row `r` and column `c`; the class code from `Radix::code` reads bit `i` as the class of the `i`-th canonical residue, and the canonical system at base `3` holds `-1` and `-i` where the box holds `2` and `2i`.
- Asserted inside the run: size `8`, fill `64`, distinct `64`, dimension `1.892789`, digits not canonical, class code `479`, and `479 != 495`.

## RUN

```
cargo run --release --manifest-path research/beneath/lab/rs/radix-carpet/Cargo.toml
```

- Prints only, `0.00` seconds.

## WITNESSES

- The carpet rows of [the radix dial claims](../../../claims/radix-dial.md): the box digits fill `64` words at level `2` with similarity dimension `1.892789`, and the same carpet codes `479` over the canonical classes where it codes `495` in box row-major order.
