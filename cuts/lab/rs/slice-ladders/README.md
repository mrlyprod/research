# Slice Ladders

- Prints the central diagonal slice census of the base-3 middle-digit solid, the `dim`-dimensional Menger sponge, at `dim = 2..6`: the cells of level `L` whose digit vectors never put the middle digit in two coordinates at one position, with coordinate sum `dim(3^L - 1)/2`.
- Three routes: `mrly_math::counts::ladder::ladder`, the carry automaton's even block to the `i128` limit; the product of the digit polynomial at the scales `3^j` to level 9; a direct enumeration of the kept cells while their number stays under `3 * 10^7`.
- The recurrence is the characteristic polynomial of the even block, checked on every automaton term, with the Hankel determinants of order `ceil(dim/2)` from level 0 and from level 1.

## RUN

- `cargo run --release --manifest-path research/cuts/lab/rs/slice-ladders/Cargo.toml -- <verb>`, verb one of `check`, `ladders`, or `all` for both; each prints only and writes nothing.
- `check`: about one second. The three routes against each other and the recurrence on every term, at `dim = 2..6`.
- `ladders`: under one second. The even block, its determinant, the signature, the two Hankel determinants and the terms from level 0, at `dim = 2..6`.

## PARI

- `check.gp`: the second route for `dim = 4, 5, 6`, `gp -q check.gp` from this folder, under one second: builds the digit polynomial, the even block and its characteristic polynomial, prints the signature and the first 8 terms, and checks the digit-polynomial product against the block iteration.
