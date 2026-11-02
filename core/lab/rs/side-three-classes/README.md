# side-three-classes

- The shape classes of the base-2 designs drawn at side 3, and the fills that vary inside one cube class at every side from 2 to 7.
- Shapes: all `2^(2^dim)` codes at `dim 1..4` drawn cell by cell from the parity rule on the `3^dim` grid, counted distinct, then reduced to the least pattern under the `2^dim * dim!` signed axis permutations of the grid.
- The reflection `x -> 2 - x` fixes every parity, so only axis permutations move a pattern: the shape classes are the designs up to axis permutation, `4, 12, 80, 3984`, against `3, 6, 22, 402` cube classes of codes.
- Burnside: `(1/dim!) * sum over axis permutations of 2^c`, for `c` the cycles on the `2^dim` corners, to `dim 6`, the second builder.
- Fills: the drawn cell count against the closed form `fill_from_corners` for every code, and per side the cube classes whose members fill differently.
- Class of `bang dim 2, code 7`: the codes `7, 11, 13, 14`, their cell counts at side 3 and `log(cells)/log(3)`.

## RUN

- `cargo run --release --manifest-path research/core/lab/rs/side-three-classes/Cargo.toml`
- About a second on 8 cores; prints every table.

## WITNESSES

- core.md:102 the four terms `4, 12, 80, 3984` with `37333248, 25626412338274304` from Burnside, the class `7, 11, 13, 14` with `8, 7, 7, 5` cells, mixed-fill classes `4, 20, 400` of `6, 22, 402` at side 3 and none at sides 2, 4, 6.
- method.md:103 the same class of `bang dim 2, code 7`.
