---
title: Complex dimensions
lead: Complex dimensions: every design in the lattice class, the carpet not Minkowski measurable by earlier theorems (Deniz et al., Lapidus, Pearse and Winter) with an explicit profile added here, the sponge's non-measurability found in none of the sources read, and the arithmetic pole at `s_0 + 2 pi i/log 3` certified genuine.
figure: research-dimensions
slug: dimensions
---

[The core page](../../core/notes/core.md) gives every design one real number, [the dimension](https://en.wikipedia.org/wiki/Fractal_dimension) `log(fill)/log(base)`. Fractal-string theory promotes that number to the real part of an infinite family - the *complex dimensions*, the poles of a zeta function attached to the set's gaps. This page computes them for the 1D designs, watches the imaginary parts surface as an oscillation in the box count, and follows the theory to its structural consequence, which needs careful qualifying.

`lab/py/complex-dimensions` regenerates every number below except those of the arithmetic-pole section, which `lab/py/burnol-residue` prints; both print only and keep no log.

## Chapters

- [The string of a design and its arithmetic pole](string-to-composition.md): Computes the complex dimensions of the 1D designs, watches them surface as an oscillation in the box count, certifies the arithmetic pole, and shows composition multiplies the base.
- [Measurability and the sponge](measurability-and-sponge.md): States the measurability criterion with its hypotheses, proves the carpet case, and follows the sponge through the criterion without pluriphase, including the tube past 1/6 and what would open the door it shuts.
- [Staircase schedules and the unequal split](staircase-schedules.md): Builds the cheapest non-stationary object, a staircase schedule, and works the unequal split: its complex dimensions, cell count, detectors, and what survives.

## Where the numbers live

`lab/py/complex-dimensions` is the one pass behind every number on this page - the poles, the box-count periodogram and folding, the composition, the two-ratio control and the tube - and it prints only. `cuts/lab/py/slice-ladder-controls` prints the staircase dimensions, and `lab/py/burnol-residue` the arithmetic-pole section: the certified residues, the band, the controls and the column. `lab/py/unequal-split` prints the unequal split: the patch, its complex dimensions and their approach to the line, the count of cells with its periodogram, and the ripples. The dimension formula this page extends, and the designs it names, are [the core page](../../core/notes/core.md); the spectral side of the same fractals is [the complexity page](../../complexity/notes/complexity.md).
