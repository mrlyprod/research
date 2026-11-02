---
title: Complexity
lead: Boolean complexity of the catalog, and the Laplacian spectra of the fractals it builds.
figure: research-complexity
slug: complexity
date: 2026-08-28
---

# Complexity

Boolean complexity of the catalog, and the Laplacian spectra of the fractals it builds.

- [notes/complexity.md](notes/complexity.md): Complexity
- [claims/boolean-complexity.md](claims/boolean-complexity.md): Boolean complexity
- [claims/flake-band-gap.md](claims/flake-band-gap.md): Flake band gap
- [claims/laplacian-degeneracy.md](claims/laplacian-degeneracy.md): Laplacian degeneracy
- [claims/spectral-spacings.md](claims/spectral-spacings.md): Spectral spacings
- [refs.md](refs.md): the references this topic cites
- [lab/py/base-q-anf/](lab/py/base-q-anf/): Builds the base normal form of a design
- [lab/py/boolean-measures/](lab/py/boolean-measures/): Computes seven Boolean complexity measures, `s`, `bs`, `C`, `dt`, `deg`, `dnf`, `cnf`, on every hyperoctahedral class of designs at `dim 3` (22 classes) and `dim 4` (402 classes), from the truth-table definitions.
- [lab/py/flake-band-gap/](lab/py/flake-band-gap/): Builds the base-2 code-23 flake at any `level`
- [lab/py/laplacian-degeneracy/](lab/py/laplacian-degeneracy/): Builds the level-`level` Sierpinski triangle as the design-code-7 base-2 tile Kronecker-substituted into itself, `level = 1..8`, and checks the array against the `i AND j = 0` characterisation instead of assuming it.
- [lab/py/tree-matching/](lab/py/tree-matching/): Rebuilds the level-`level` Sierpinski triangle of design code 7 at base 2 straight from `i AND j = 0`, no Kronecker step, with the axis-nearest-neighbour graph, `level = 1..10`.
- [lab/rs/spectral-spacings/](lab/rs/spectral-spacings/): Computes the unfolded nearest-neighbour level spacings of graph Laplacians, both combinatorial `D - A` and normalised `I - D^-1/2 A D^-1/2` for the degree matrix `D` and the adjacency matrix `A`, for the carpet cell graph at level 3, 4 (512 and 4096 nodes), the Menger sponge cell graph at level 2 (400), the sponge's diagonal slice on the plane `x + y + z = const` at level 2, 3 (306 and 2250) and the Sierpinski gasket at level 5, 6 (243 and 729), after two controls, a seeded random graph on 400 nodes at `p = 0.1` and the `20 x 20` square lattice.
- [figures/research-complexity.ts](figures/research-complexity.ts): research-complexity

## Lost

- Graph-Laplacian spacings as a Hilbert-Polya test: every fractal clusters, and nothing in it concerns the Riemann hypothesis (Spacings: the spectra cluster).
