# Walk Dimension

- The walk dimension `d_w` of the nine census designs by two readings: Laplacian eigenvalue level ratios `lambda_k(level)/lambda_k(level+1) = base^d_w` on the giant component, and the mean squared displacement of 20000 blind ants, `MSD(t) = t^(2/d_w)`, fitted from `t = 32` below `(side/6)^2`.
- Three anchors first: the solid grid and the path drawn by code 7 against the closed form `2 - 2 cos(pi/n)` at `n` nodes and `d_w = 2`, and the corner-glued Sierpinski gasket against `tau = 5`, `d_w = log 5 / log 2`.
- The base-3 class census: 26 classes of the 72-element wreath group with orbits summing to 512, the giant share of every representative at levels 4 and 5, and over all 511 codes which giant components touch all four walls at level 5 and how many classes they fall in.
- Code 127 read on four modes at two level pairs, and the weight of its slowest mode inside the block hanging from tile cell (2,0) at level 4.
- Level pairs 4 to 5 in base 3 with 3 to 4 as the drift bar, 6 to 7 for the base-2 design, 3 to 4 for the sponge (160000 nodes); walkers at level 5, 8 and 4, gasket at level 8.
- The two-panel census figure: `d_w` by both readings per design, and `d_s = 2 dimension / d_w` against `dimension` with the gasket starred.
- Dense eigenvalues through faer up to 2000 nodes, above that a block Krylov Rayleigh-Ritz projection through a projected conjugate gradient; the walker stream is seeded but its digits are not the recorded ones, so the walker column agrees to about 0.01 and the gasket gate to 1.6%.

## RUN

- `cargo run --release --manifest-path research/walks/lab/rs/walk-dimension/Cargo.toml`
- About six minutes; writes `walks-fig.png` to `data/research/walks/lab/rs/walk-dimension/`.

## WITNESSES

- walks.md:32-35 level pairs 4 to 5 and 3 to 4 with 160000 nodes, modes `k = 1..4`, 20000 walkers.
- walks.md:43-45 solid `d_w = 1.99990`, path `d_w = 1.99999`, machine-exact `lambda_2`; gasket ratios `4.9973, 4.9973, 4.9945, 4.9972` at levels 7 to 8; the walker anchors 2.0020, 2.0002 and 2.3486 (1.15%) come back as 2.0035, 2.0094 and 2.3590 (1.60%) on this stream.
- walks.md:57-58 rep 15 giant share 0.42 then 0.32, rep 31 0.22 then 0.13; walks.md:70 rep 238 with 1556 components and share 0.0312.
- walks.md:67-71 95 of the 511 codes have a giant component touching all four walls at level 5, in seven classes with 238 among them; only 83 codes in six classes are a single component, and 245, 350, 371, 413 carry 1556 components with giant share 0.0624.
- walks.md:94-102 the census table: spectral 2.466, 2.543, 2.640, 2.190, 2.097, 2.000, 2.586, 2.164 and the walker column 2.487, 2.583, 2.257, 2.643, 2.170, 2.123, 2.002, 2.743, 2.182 to about 0.01.
- walks.md:159-162 code 127 modes 1 and 2 near 2.53, modes 3 and 4 at 2.20 to 2.29 at both level pairs, walkers 2.245 with the spectral drift 0.004.
- walks.md:199 the `lambda_2` readings 2.4662 and 2.5433 for reps 79 and 95.
- README.md:57 `dimension = 1.7712` shared by 127 and 239, walk dimensions 2.25 and 2.64, the carpet 2.097 to 2.124.
- the figure `walks-fig`, from a rerun.
