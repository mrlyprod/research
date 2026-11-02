# Encyclopedia

- Runs every basic recipe of mrly life on the infinite plane and reads it: the recipe grammar, the two-mode engine, the readings, the cache and the timings behind the encyclopedia census.
- A recipe is one name that regenerates its run bit for bit: `carpet.3.2 x3 island mask:net.5.1 birth:carpet_fills_zeros survive:carpet_voids x:1024 cap:301`; a wallpaper drops the cap and tiles once.
- The grammar (`Grammar::standard`, every list a parameter): seeds the plane classes `code1 code3 code6 code7 code15` and the ten catalog designs at sides 3, 5, 7 and levels 1, 2, equal designs run once; tilings 1, 3, 5 as islands and one wallpaper; masks the same designs popped at sides 3, 5, 7, 9 and levels 1, 2, equal masks run once, plus `copy` and `copyinv` of the tiled seed when they differ from every design mask; rules the 23 sources with the four rolls on both sides, `B3/S23`, and the eight family pairs with the four rolls, fills for birth and voids for survive.
- Island: the seed alone on the plane, the canvas grows with the pattern, the run stops when a live cell comes within the mask's reach of the cap wall or at X; every frame before the stop is the infinite plane.
- Wallpaper: the seed is the torus, which is the seed-tiled plane, for X generations.
- A zero in birth runs as the rule pair on a flipping background, the frames and readings being the pattern against its background; the `flips` and `background` fields say so.
- A repeated frame stops the run early as still, loop or mover with its period and displacement; a mover reports the generation by which it reaches the wall.
- Counts come from the mask's Kronecker factors, one dilated pass per factor, equal to `next_grid` on both boundaries.
- Readings: the stop and its generation, the settled generation, the population and box curves, the box growth slope, entropy, churn, symmetry, the last frame as hex rows with its origin, the Kronecker cuts of the frame, the heatmap's visited count, peak and median level set with its cuts; readings are never filters.
- The cache is one JSON line per recipe under the repo's `data/research/automata/lab/rs/encyclopedia/chunks/`, keyed by the name, so a rerun of a chunk skips every cached recipe; the census JSON of the note is written from these lines by a later bin.
- The note of the catalogue is [the encyclopedia](../../../notes/encyclopedia.md): the object, the two modes, the census as run, the laws, the sequences, the library and what the catalogue does not contain.
- The test-only torus configuration (seed centred in a wrapped canvas) replays the retired mini preset: 1272 runs on canvas 27 for 64 generations agree rule for rule with the old atlas run today; the old number-rule line of the claims no longer replays on its own witness either, since the `random_0` source has changed.

## RUN

- `cargo run --release --manifest-path research/automata/lab/rs/encyclopedia/Cargo.toml -- count --x 1024 --cap 301` prints the grammar's exact counts.
- `... -- list FROM TO` prints names; `... -- run NAME` prints one cache line.
- `... -- chunk FROM TO --x X --cap CAP --threads 4` runs a slice of the recipe list into the cache, resumable.
- `... -- time --threads 4 --per 8 --cap 301` prints the cost per recipe of every class at X in 64, 256 and 1024.
- `... -- big --x 256 --cap 151 --threads 4 --minutes 20 --chunk 50000 --floor 15` runs the catalogue in chunks of 50000 in the grammar's order on a wall clock, launches no chunk past the minutes or under the disk floor, lets the running chunk finish, and prints the cut: the last index and recipe, the rows, the bytes a recipe, the thread-seconds.
- `... -- laws --oeis PATH` reads every chunk line and prints the laws of the catalogue: the cut and the classes complete or cut, the stop tally by class, the wall law by mask kind, the budget law, island against wallpaper, the odd-rule periods against the Frobenius order, the ten population sequences with their windows looked up in the OEIS dump at PATH; it writes `tally.json` to the data dir for the census bin `research-automata-encyclopedia`.
- `... -- library` writes `library.json` to the data dir: one shelf cell per stop, class and seed family, the first 8 recipes of each in grammar order, readings only.
- `... -- mini` prints the mini control table.
- The cache line of a recipe is sparse: the stop, generation, period, displacement, wall_by, settled, growth slope, entropy, churn, symmetry, the Kronecker cuts and tilings when any, the heat visited count and peak, flips and background when set, and the population and box curves sampled as the first 16 values, every 16th after, the last, the max and its generation; `run NAME --full` prints the full line with the last frame and the heat level set.
- `cargo test --manifest-path research/automata/lab/rs/encyclopedia/Cargo.toml` runs the hand checks: blinker, glider, R-pentomino, island against a large canvas, wallpaper against the torus, the sublattice lemma, the flipping background, the mini control.

## LAWS

- The odd-rule law, proved: a rule whose birth and survive sides are both the odd counts (`odds_ones`, and `odds_zeros_ones`, which draws the same counts) sets the next cell to the neighbour count mod 2, so the next frame is `x_{t+1} = A x_t` over `GF(2)` with `A = sum over the mask offsets o of the shift sigma^o`, the centre being popped. On a torus of odd side `n` every shift has order `n` and `A^{2^k} = sum of sigma^{2^k o}` by the Frobenius map of characteristic 2; with `k = ord_n(2)` this is `A` again, so `x_{t + 2^k - 1} = A^{2^k - 1} x_t` and, applied to `x_1`, `x_{2^k} = A^{2^k} x_0 = A x_0 = x_1`: from generation 1 the run is periodic with a period dividing `2^k - 1`. Tori of side 3, 5, 7, 9, 25 and 49 have `k = 2, 4, 3, 6, 20, 21`, so the periods divide 3, 15, 7, 63, 1048575 and 2097151; the last two lie past `X = 256`, so the catalogue reads those runs as limits.
- The affine rule `evens_zeros` on both sides (even counts, zero included) sets the next cell to `1 + count mod 2`, which is `x_{t+1} = A x_t + 1` with the same `A`; the catalogue prints its periods, and the Frobenius bound is not claimed for it.

## WITNESSES

- `count --x 256 --cap 151`: 2744500 recipes, 66 seeds, 82 design masks, 125 rules, 38500 copy-mask recipes (Verified).
- `big --x 256 --cap 151 --threads 4 --minutes 20 --chunk 50000 --floor 15`: the whole grammar in 55 chunks, 919 s wall, 3663 thread-seconds, 944266606 bytes, 344 bytes a recipe, no cut (Verified).
- `laws --x 256 --cap 151 --oeis PATH`: the stop tally, the wall law, the budget law, island against wallpaper, the odd-rule periods, the ten population sequences; every number of `claims/encyclopedia.md` (Verified, the Frobenius bound Proved above).
- `library`: 207 shelf cells, 1648 rows, 721011 bytes (Verified).
- The census bin `research-automata-encyclopedia` of `research/shared/census` reads `tally.json` from the data dir and writes the figure's JSON: the stop shares dead, periodic and open by mask kind against rule class.
