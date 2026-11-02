# percolation-race

- Races five self-similar designs against random cell sets of exactly the same size on the same grid, on three metrics: face-connected components, largest-component fraction, and boundary per cell `(2*dim*N - 2*E)/N`.
- Pass A builds each design as a Kronecker power of its level-1 tile, finds components by union-find, and draws the control with `numpy.random.default_rng` on PCG64.
- Pass B builds each design by substitution over coordinate sets, finds components by breadth-first search on a padded grid, and draws the control with `random.Random` on the Mersenne Twister.
- The dispersing tile is the antidiagonal in pass A and the diagonal in pass B, so the two passes race reflected copies of it.
- Pass B also checks the substitution route against the digit rule for level 1, 2, 3 on all five designs, and counts the per-draw ties the means hide.
- The passes share no build, no search and no generator; the last block prints every overlapping comparison as a multiple of the sample standard deviation and exits nonzero if any check fails.

## DOMAIN

- Pass A: 2D grids to `81x81`, 3D to `81^3`, 400 seeds thinning to 25 at the largest.
- Pass B: 2D grids to `256x256`, 3D to `81^3`, 400 seeds thinning to 20 at the largest.
- The program prints this domain first and runs it whole in about 9 s, cutting nothing.

## RUN

```
uv run python research/connectivity/lab/py/percolation-race/race.py
```

## WITNESSES

- connectivity.md, section "The race" - the two passes agree on every comparison to within a standard deviation; the widest gap printed is 0.2061 sd.
- connectivity.md, section "Connectivity" - the components table, 1 against `134.39 +/- 7.34` at the `32 x 32` gasket up to 1 against `31576.40 +/- 152.71` at the `81^3` sponge, with every density, and the largest fraction 1.0000 against `0.0411 +/- 0.0098`, 0.0011 and 0.0169.
- connectivity.md, section "Boundary" - the boundary table, including the diagonal's exact 4.0000 against `3.8808 +/- 0.0897`, and the gasket's measured drift 2.0082, 2.0027, 2.0009, 2.0003, which `../boundary-drift` proves exact.
- connectivity.md, section "How decisive, honestly" - the per-draw ties 0 of 400 and 0 of 100 at the gasket, 73 of 400 then 0 of 200 at seven-of-eight, 61 of 400 then 27 of 200 at the diagonal.
- connectivity.md, section "How decisive, honestly" - the dispersing extreme, random 30.10 against 32 at `32 x 32` and 125.93 against 128 at `128 x 128`.
- README.md:56 - the matched-noise bullet: `5257.23 +/- 32.94`, `31576.40 +/- 152.71`, 0 of 400 at the `32 x 32` gasket, boundary 2.0003 against `3.6006 +/- 0.0103`.
