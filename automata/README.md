---
title: The automata
lead: An elementary cellular automaton is a three-dimensional design read as a rule and Life a nine-dimensional one: what the identity buys, and what it cannot.
figure: research-automata
slug: automata
date: 2026-09-02
---

# The automata

An elementary cellular automaton is a three-dimensional design read as a rule and Life a nine-dimensional one: what the identity buys, and what it cannot.

- [notes/automata.md](notes/automata.md): The automata
- [notes/encyclopedia.md](notes/encyclopedia.md): The encyclopedia
- [notes/life-and-the-atlas.md](notes/life-and-the-atlas.md): Life, the siblings and the atlas
- [claims/atlas.md](claims/atlas.md): Atlas
- [claims/automata.md](claims/automata.md): Automata
- [claims/encyclopedia.md](claims/encyclopedia.md): Encyclopedia
- [claims/waves.md](claims/waves.md): Waves
- [refs.md](refs.md): the references this topic cites
- [lab/py/carpet-nodal/](lab/py/carpet-nodal/): Builds the Sierpinski carpet as the Kronecker power at `level` of the 3 x 3 tile with the centre void, `8^level` cells, asserted at 512 for `level 3` and 4096 for `level 4`.
- [lab/py/kernel-waves/](lab/py/kernel-waves/): Asks whether the wavelength of a Larger-than-Life still is set by the mask or by the rule, and whether it is readable off the mask's Fourier transform alone.
- [lab/py/life-census/](lab/py/life-census/): Reads the `2^18` outer-totalistic life-like rules on the Moore neighborhood as designs on the parity cube `{0,1}^9`, and computes fill, Langton `lambda`, `GF(2)` degree and genus for all of them
- [lab/py/sibling-census/](lab/py/sibling-census/): Reads a MrlyMath automaton as a triple `(dim, mask, kind)`
- [lab/rs/automata-census/](lab/rs/automata-census/): Reads the 256 elementary cellular automata as the 256 base-2 designs of `bang dim 3`, then measures how much of the design tree survives as dynamics.
- [lab/rs/encyclopedia/](lab/rs/encyclopedia/): Runs every basic recipe of mrly life on the infinite plane and reads it
- [figures/research-automata.ts](figures/research-automata.ts): research-automata

## Lost

- Chladni figures in Larger-than-Life: the search found none and its code left the library; the Chladni figures of a design are the nodal domains of section 13, The waves.
- The waves as a route to the Riemann hypothesis: they are not one; what stays open is the lobe conjecture.
- The atlas and its mini, medi and maxi presets: a padded finite torus is neither the infinite plane nor the seed-tiled plane; the encyclopedia replaced it and replays mini as a control.
- The mini number line 106, 24, 56, 132 (dead, still, loop, timeout): the `random_0` source changed and the six number rules now end 99, 18, 65, 136 on the same preset, on the atlas as on the encyclopedia.
- That a tessellated seed's still peaks on the seed comb: the comb hits read on the retired atlas were lattice-mask tilings and four-cell frames, and the Moore stills read off the comb; no live generator replays it.
- The lobe law on the sequence rules: the peak ring of a drawn-mask frame sat in the mask's first negative lobe no more often than the cut itself allowed, on a still-hunt generator that retired with mrlygame.
- The seed-comb product law: peak rings on comb multiples came at the rate of chance on the same still-hunt frames.
- The atlas closure census of index-1 frames cutting properly in place: with the changed `random_0`, the count was read again on the retired atlas and gave index-1 `random_0` frames that cut properly in place; no live generator replays it.
- Three atlas rows the retirement dropped: the quest replay outside the engine, the ring cut over the default quests and the capped-canvas copy-path seed; each was read on the retired atlas, the ring cut and the seed on a still-hunt generator that had already retired, and no live generator replays them.
