# Arc Loops

- Counts the closed loops and the strands of a dimension-two design drawn in quarter-circle arcs: every cell of level `n` carries two arcs, a filled cell around its lower-left and upper-right corners, a deleted cell around the other two, and arcs meeting at an edge midpoint join into curves.
- The cells come from `mrly_math::bang::factory::create`, at `bang dim 2, base b, code c` with bit `b y + x` for the mask cell in column `x` and row `y`, row 0 at the bottom, or at a parity code drawn at a side number; the study prints both names of every parity design, so the carpet reads `bang dim 2, base 3, code 495 = bang dim 2, code 7 at side number 3`.
- Three counts: union-find over the edge midpoints, a component being a loop when every midpoint in it has degree 2; `c(G) - 2 side - 1` with `c(G)` the components of the mirror graph, isolated lattice points included, one diagonal per cell; and the block recursion, which glues the strand matchings of `base^2` blocks of the level below and never builds the cells. The first two share the cells built by `mrly-math` and one union-find structure; the block recursion shares neither.
- The block recursion holds the matching of a level as one array of `4 side` ports and takes the void block by its formula `B_t - R_(N-1-t)`, `L_t - T_(N-1-t)`, so memory is linear in `side`.
- Fits are Berlekamp-Massey over the prime `2^61 - 1`, lifted to the integers and checked exactly on every term, tried from each start level 0 to 7, and printed only with at least two terms of margin; the integer roots are split off by synthetic division and any factor left over is printed as a factor with no integer root.
- With `OEIS_STRIPPED` set to a local copy of the [OEIS stripped dump](https://oeis.org/stripped.gz), every nonzero sequence is looked up by its first seven terms from the first nonzero one.

## RUN

- `cargo run --release --manifest-path research/arcs/lab/rs/arc-loops/Cargo.toml -- <verb>`, verb one of `check`, `two`, `carpet`, `gains`, `lengths`, `bases`, `census`, or `all` for every one; each prints only and writes nothing.
- `check`: under one second. Every code at base 2 to level 7 and at base 3 to level 4, the three counts and the strand count against each other, the three images of each code under the half turn and the two diagonal reflections, the glued void block against its formula, and the parity route against the full code.
- `two`: under one second. The base 2 lemmas of codes 7, 11 and 9 against the glued matching at levels 1 to 14, and their gains.
- `carpet`: about one second. Union-find and cycle rank to level 7, the block recursion to level 15, the law, the matching lemma at levels 0 to 10 and the gain split `10 |P_n| + 3n`.
- `gains`: about two seconds. The loop sequence and the gains `J(n) = L(n + 1) - k L(n)` of the named designs, each tested against its stated gain law at every gluing and its stated loop law at every level.
- `lengths`: under one second. The new loops of each gluing by their number of block strands, void-block strands counted, and the longest, for the carpet and codes 43, 171, 175 and 181 at base 3 to level 12.
- `bases`: about fifty seconds. The parity codes and the one-cell deletions at base 4 to level 12 and at base 5 to level 10, with the distinct nonzero sequences they give and their OEIS hits.
- `census`: about two minutes and 1.1 GB. Every code at base 2 to level 24 and every class at base 3 to level 14, the classes without a fit rerun to level 17.

- `cargo test --manifest-path research/arcs/lab/rs/arc-loops/Cargo.toml`: three tests, the three counts at base 2 to level 5, the carpet law to level 8, and the glued void block.

## WITNESSES

- The study is the generator behind every section of [arcs.md](../../../notes/arcs.md).
- `check`: 2688 levels, 1584 images, 27 void sides and 32 parity levels, zero mismatches; the strand law, the cycle-rank theorem, the symmetry, the void block and the block recursion.
- `two`: the base 2 table, every lemma at 14 of 14 levels and every gain at 13 of 13 gluings.
- `carpet`: the carpet law at 16 of 16 levels, the lemma at 11 of 11 levels, the gain split at 15 of 15 gluings.
- `gains`: the base 2 and carpet laws at every gluing printed, the parity table at side number 3 at gluings 0 to 13, and the Fibonacci gains `F_(2n-3)` of code 13 at gluings 1 to 13.
- `lengths`: the loop of 6032 strands that code 43 closes at the gluing from level 11 to 12, and the carpet's longest new loop, 8 strands, over gluings 0 to 11.
- `bases` and `census`: every count, fit and OEIS line of the census sections.
