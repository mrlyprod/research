# Schedules

- The schedule study: which schedules of design rounds, with every stamp, generate a primitive group, read off local conditions and checked against the exact block census.
- The machine is the one of the round census, rebuilt here: a grid of `N = n^2` bits, bit `n i + j` the cell `(i, j)`, a split `A` written from its complement `B` by a symmetric mess `[count in S]` over a base-3 level-1 neighbourhood (`carpet 495`, `net 186`, `void 341`), a spread `id`, `shift`, `drop`, `slide`, `turn` or `flip`, and the stampless round `h = L . F`.
- A schedule is a list of rounds; the group is `G = <T, h_1, ..., h_k>` with `T` all `2^N` stamps. A pair is the split and its complement written in turn with the identity spread; a single is one split with one spread.
- The block census: for each cell the closure under derivatives `h(x + u) + h(x)` of every unit vector and of the all-ones vector, then of every nonzero vector when those are all the whole space; a proper closure is re-verified as a block. A cell with no block is primitive.
- The long prime: a primitive cell is searched for a word of stamped rounds with a cycle of prime length `p`, `M/2 < p <= M - 3`, which alone proves the group primitive and alternating.
- The conditions, each computed exactly: `R` the reads digraph is strongly connected; `S` the second differences along each read cell span its written neighbourhood; `W` every live window rule is neither affine nor self-complementary; `I` every live read cell is isolated; `G` the second-difference graph is strongly connected. The theorem's hypotheses are `G` and `I`.
- The laws, for imprimitive cells: `closed` (no `R`), `mess` (a proper spread-invariant module of the mess classes), `complement`, `parity` (a kept parity of a torus partition), else `other`.
- `trap` counts the splits that contain a spread orbit and checks the formula, then checks at `3 x 3` that an orbit inside the written side is exactly its own block.
- `lemma` lists the linear structures of every symmetric function of `m = 2..10` variables.

## RUN

- `cargo run --release --manifest-path research/random/lab/rs/schedules/Cargo.toml -- <verb>`, single-threaded, prints only, writes nothing.
- The verbs: `lemma`, `trap` (2 s), `pairs` (10 s, with the cells it covers per grid and in total), `verbose` (the pairs with every nonlinear cell printed), `singles` (about 130 s), `onto` (6 s, the singles whose spread carries the split onto its complement), `family` (the checkerboard and stripes pairs at `n = 4, 6, 8, 10`), `five <seconds>` (the `5 x 5` line pairs; about 36 s to its long prime), `all` (all but `five`).

## WITNESSES

- `lemma`: at every `m = 2..10` no non-affine symmetric function has a linear structure other than the all-ones vector, and the all-ones vector is one exactly when `g(t) + g(m - t)` is constant; at `m = 10`, 2044 non-affine functions, 60 with it.
- `trap`: trapped splits at `3 x 3` id 510, shift, drop and slide 168, turn 286, flip 484 of 510; at `4 x 4` id 65534, shift, drop, slide and turn 14910, flip 64806 of 65534; at `5 x 5` id 33554430, shift, drop and slide 4925280, turn 22163806, flip 33495382 of 33554430; each equal to `2^N - 1 - prod (2^|o| - 1)`. Trapped fraction of all masks at `n = 16`: shift 0.000244, turn 0.983925, flip 1.000000. 123930 orbit closures at `3 x 3`, 0 faults.
- `pairs`: `3 x 3`, 12 splits, 144 pair cells (12 splits x 3 neighbourhoods x 4 Sources), 108 nonlinear, 12 alternating, each by a long prime (257 to 499); imprimitive 96 by first law closed 62, mess 32, other 2. `4 x 4`, 23 splits, 276 pair cells, 207 nonlinear, 10 alternating by a long prime; imprimitive 197: closed 179, mess 13, parity of the torus 3, other 2. Both catalogs: 420 pair cells, 315 nonlinear. Theorem hypotheses met in no pair cell; `R S W` met in 6 cells at `4 x 4`, all alternating.
- `pairs` other: `3 x 3 0x004e net binary` and `0x001b void binary` keep a one-dimensional block although every unit and the all-ones closure is the whole space; `4 x 4 0x134c` and `0x125a carpet primes` keep the single parities `0xcc33` and `0xa5a5`.
- `singles`: 2856 cells, 2142 nonlinear: imprimitive by first law closed 1802, mess 113, complement 8, parity 18; primitive 201, 195 of them alternating by a long prime. Theorem hypotheses met in 59 cells, all primitive, 58 alternating by a long prime. A strongly connected second-difference graph with a proper unit closure: 0 cells. Every unit closure is the whole space in 238 nonlinear cells, and every imprimitive one of those has a proper all-ones closure.
- `singles` falsifiers of `R S W`: `0x0055 net primes slide` keeps the parities of the cosets of a subgroup of order 2, `0x0505 net primes turn` those of order 4.
- `onto`: 84 cells, 63 nonlinear; `R S W` met in 14, all alternating; 49 imprimitive, closed 23, mess 19, parity 7.
- `family`: the checkerboard pair with `carpet` or `net` and `binary` or `fibonacci`, and the stripes pair with `carpet` and `primes` or `binary`, meet the hypotheses at `n = 6, 8, 10`, and not at `n = 4`; the second-difference graph of the checkerboard `carpet binary` pair at `n = 10` has depths `1, 4, 8, 12, 16, 18, 16, 12, 8, 4, 1` from cell 0, the histogram the figure redraws.
- `five`: 81 nonlinear `5 x 5` line pair cells, 7 meet the hypotheses; for the first, columns `{0, 2}` with `carpet primes`, word 1 (two stamped pair rounds) has a lone prime cycle 1229311 and word 13 (two) a cycle of prime length 31896037 in `(2^24, 2^25 - 3]`.
