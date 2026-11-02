# spectral-spacings

- Computes the unfolded nearest-neighbour level spacings of graph Laplacians, both combinatorial `D - A` and normalised `I - D^-1/2 A D^-1/2` for the degree matrix `D` and the adjacency matrix `A`, for the carpet cell graph at level 3, 4 (512 and 4096 nodes), the Menger sponge cell graph at level 2 (400), the sponge's diagonal slice on the plane `x + y + z = const` at level 2, 3 (306 and 2250) and the Sierpinski gasket at level 5, 6 (243 and 729), after two controls, a seeded random graph on 400 nodes at `p = 0.1` and the `20 x 20` square lattice.
- Two unfolders on every spectrum: a degree-12 Chebyshev least-squares fit of the counting staircase, and a local window of 21 levels; negative steps of the polynomial are clamped to zero and counted.
- Per row: `P(s < 0.5)` against GOE's `1 - exp(-pi/16)` and Poisson's `1 - exp(-1/2)`, the Kolmogorov-Smirnov distances to both laws with their asymptotic p-values, and at tolerance `1e-9` the distinct eigenvalue count, the zero spacings `n - distinct`, the eigenvalues lying in repeated classes and the largest multiplicity.
- Builds the level-4 slice as a graph and prints its node and edge counts without diagonalising it.
- Dense eigenvalues through faer; every design comes from `mrly_math` and the cell graphs from `mrly_num`.

## RUN

- `cargo run --release --manifest-path research/complexity/lab/rs/spectral-spacings/Cargo.toml`
- About ten seconds after the build, exits 0.

## WITNESSES

- complexity.md:510-512 the objects and their sizes 512, 4096, 400, 306, 2250; complexity.md:514-518 the random control reads GOE with 100% distinct eigenvalues, KS not rejecting GOE (p = 0.24 and 0.99) while rejecting Poisson (p = 0.000), and the square lattice reads clustered (`P(s < 0.5)` 0.54 to 0.59).
- complexity.md:523 the GOE prediction 0.1783; complexity.md:522-524 the band 0.44 to 0.57 and the interpolating band 0.32 to 0.36 are not reproduced here: the polynomial unfolder reads 0.50 to 0.69 and the window unfolder 0.49 to 0.67 across the five fractal rows, on both Laplacians, with GOE excluded everywhere.
- complexity.md:528-530 the KS distance is nearer Poisson than GOE in every row; the square lattice ties to four decimals under the window unfolder, 0.4987 against 0.4987, both rejected.
- complexity.md:534-537 the random control has no repeated eigenvalue; the slice on the combinatorial Laplacian is 62.09% and 63.16% distinct at `1e-9` with largest multiplicity 12 at level 2 and 48 at level 3; complexity.md:538-539 the level-4 slice has 16578 nodes and 21546 edges.
- complexity.md:623 clustering established up to 4096 nodes.
- research/claims/:294 the sponge at 400 nodes, normalised Laplacian, degree-12 unfolder: `P(s < 0.5) = 0.6892`, Sierpinski level 6 0.8297, the square lattice 0.5940, the sponge's zero spacings 245 of 400 (61.25%) forcing `P(s < 0.5) >= 0.61`; GOE's constant prints as 0.17828, the page's 0.17826 is a rounding slip.
- README.md:59 the negative result, GOE and GUE excluded up to 4096 nodes.
