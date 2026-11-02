# Crop Counts

- Regenerates the level tables and the radius sweep on [crop](../../../notes/crop.md): exact shape crops over the carpet and the sponge.
- The remaining crop.md numbers are pinned by `mrly_math::shape` tests or derived from the printed lines.
- The designs are code 7 at dim 2 (the carpet, `mrly_math::bang::factory::create(7, 3, 2, 2, level)`) and code 23 at dim 3 (the sponge, `create(23, 3, 3, 2, level)`), with level 0 the single filled cell.
- The shapes are `mrly_math::shape::named` ball and diamond, centered at one half on every axis, classified by `mrly_math::shape::census` in exact integer arithmetic, no floats anywhere.
- A line prints `filled_in` and `filled_cut`, the filled cells fully inside and crossing the boundary, and `exposed_after`, the exposed unit faces of `crop(types, shape, true)` - a face counts as exposed when its neighbour is empty or off the grid.
- The level series runs the inscribed radius `1/2` at levels `0..5` in 2D and `0..4` in 3D; the sweep runs radii `1/24 .. 24/24` at carpet level 4 (side 81) and sponge level 3 (side 27).
- Every line asserts the partition and complement identities before printing: the three filled tallies sum to the design's fill, `crop` at each `keep_cut` matches its census columns, and the crop and the `Shape::Anti` crop with the complementary cut rule partition the filled set exactly, both ways.

## RUN

- `cargo run --release --manifest-path research/crop/lab/rs/crop-counts/Cargo.toml`
- Under five seconds; prints only, writes nothing.

## WITNESSES

- crop.md the carpet level series: cut cells `1, 8, 28, 76, 204, 580` and the full ball and diamond tables at level 0..5.
- crop.md the sponge level series: cut cells `1, 20, 216, 1224, 7968` and the full ball and diamond tables at level 0..4.
- crop.md the dead zones: carpet crops empty below `r = 1/6` for ball and diamond alike, sponge diamond below `r = 1/3`; the sponge ball sweep is empty through `r = 5/24` and first cuts at `r = 6/24`, the first grid radius past the exact contact `sqrt(2)/6`.
- crop.md the saturation radii: the carpet ball crop holds all 4096 cells from `r = 17/24`, the sponge ball crop all 8000 from `r = 7/8`, and the sponge diamond never saturates in the sweep.
- crop.md the partition and complement identities, asserted on all 118 printed configurations.
