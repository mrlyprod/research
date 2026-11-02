---
title: Complex dimensions
lead: Complex dimensions: every design in the lattice class, the carpet not Minkowski measurable by earlier theorems (Deniz et al., Lapidus, Pearse and Winter) with an explicit profile added here, the sponge's non-measurability found in none of the sources read, and the arithmetic pole at `s_0 + 2 pi i/log 3` certified genuine.
figure: research-dimensions
slug: dimensions
date: 2026-08-28
---

# Complex dimensions

Complex dimensions: every design in the lattice class, the carpet not Minkowski measurable by earlier theorems (Deniz et al., Lapidus, Pearse and Winter) with an explicit profile added here, the sponge's non-measurability found in none of the sources read, and the arithmetic pole at `s_0 + 2 pi i/log 3` certified genuine.

- [notes/dimensions.md](notes/dimensions.md): Complex dimensions
- [notes/string-to-composition.md](notes/string-to-composition.md): The string of a design and its arithmetic pole
- [notes/measurability-and-sponge.md](notes/measurability-and-sponge.md): Measurability and the sponge
- [notes/staircase-schedules.md](notes/staircase-schedules.md): Staircase schedules and the unequal split
- [claims/arithmetic-pole.md](claims/arithmetic-pole.md): The arithmetic pole
- [claims/complex-dimensions.md](claims/complex-dimensions.md): Complex dimensions
- [claims/sponge-measurability.md](claims/sponge-measurability.md): Sponge measurability
- [claims/unequal-split.md](claims/unequal-split.md): The unequal split
- [papers/sponge-measurability.md](papers/sponge-measurability.md): The Menger Sponge Is Not Minkowski Measurable
- [refs.md](refs.md): the references this topic cites
- [lab/py/burnol-residue/](lab/py/burnol-residue/): The residue of `K(s) = sum n^(-s)` over the positive integers whose base-3 digits lie in `{0,1}` at the off-real poles `s_(0,k) = log_3 2 + 2 pi i k/log 3`, `k = 0..10`, enclosed in a complex rectangle by interval arithmetic (`mpmath.iv`, 128 bits) with every truncation bounded by a proved tail
- [lab/py/complex-dimensions/](lab/py/complex-dimensions/): The complex dimensions of the four one-base designs, base 3 `{0,2}`, base 5 `{0,2,4}`, base 15 `{0,4,10,14}` and base 15 `{0,2,4,10,12,14}`
- [lab/py/sponge-tube/](lab/py/sponge-tube/): The Menger sponge through Kombrink, Pearse and Winter 2016, Corollary 3.2, with the open unit cube as the feasible open set
- [lab/py/unequal-split/](lab/py/unequal-split/): The patch
- [lab/rs/sponge-window/](lab/rs/sponge-window/): Falsifies the tube of the Menger sponge's plus past `1/6` and the corner volume `Deep` at every radius, as written in measurability and the sponge, subsection "The tube past `1/6`, and `Deep` in closed form", against the crate's exact distance `mrly_math::three::sponge::distance`.
- [figures/paper-sponge-measurability.ts](figures/paper-sponge-measurability.ts): paper-sponge-measurability
- [figures/research-dimensions.ts](figures/research-dimensions.ts): research-dimensions
- [figures/research-dimensions.json](figures/research-dimensions.json): research-dimensions census

## Lost

- The Lapidus-Maier inverse spectral problem on a design: every one-base design is lattice, so the machinery is empty here ([notes/measurability-and-sponge.md](notes/measurability-and-sponge.md)).
- A nonlattice design or periodic composition: none exists, by block reduction (Composition multiplies the base).
- A single ripple about a fixed point as a nonlattice detector: it keeps that map's period, so one point at a time never sees the unequal split; two ripples about different fixed points do (What the detectors see).
