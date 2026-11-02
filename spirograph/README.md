---
title: The spirograph
lead: The roulette a design draws: the coincidence law that counts its distinct curves, the node count `2ab C(k,2) + k a(b-1)` inside the window `min(1, A)` with the plane cut into `nodes + 2` regions, the exact crossing staircase above the loop threshold whose thresholds are algebraic numbers obeying a reciprocity, the alignment reaches where four of the carpet's curves run through one point, and the walls, whose cover no closed form tried survives.
figure: research-spirograph
slug: spirograph
date: 2026-09-10
---

# The spirograph

The roulette a design draws: the coincidence law that counts its distinct curves, the node count `2ab C(k,2) + k a(b-1)` inside the window `min(1, A)` with the plane cut into `nodes + 2` regions, the exact crossing staircase above the loop threshold whose thresholds are algebraic numbers obeying a reciprocity, the alignment reaches where four of the carpet's curves run through one point, and the walls, whose cover no closed form tried survives.

- [notes/spirograph.md](notes/spirograph.md): The spirograph
- [claims/spirograph-loops.md](claims/spirograph-loops.md): The spirograph loops
- [claims/spirograph-nodes.md](claims/spirograph-nodes.md): The spirograph nodes
- [claims/spirograph-reaches.md](claims/spirograph-reaches.md): The spirograph reaches
- [claims/spirograph-walls.md](claims/spirograph-walls.md): The spirograph walls
- [refs.md](refs.md): the references this topic cites
- [lab/rs/roulette-cover/](lab/rs/roulette-cover/): Measures the shape a roulette walls off inside its own disc, the question the spirograph page asks, and prints only.
- [lab/rs/roulette-loops/](lab/rs/roulette-loops/): Reads the self crossing count of one trochoid as an exact step function of the pencil's reach `|p|`, past the loop threshold `|p| = 1`, and the crossing count of two trochoids drawn from one wheel.
- [lab/rs/roulette-nodes/](lab/rs/roulette-nodes/): Counts the nodes of a design's roulette
- [lab/rs/roulette-reaches/](lab/rs/roulette-reaches/): Finds the reaches at which three or four of a wheel's trochoids run through one point on a circle track, and reads them back off `mrly_math::spirograph::point` and `mrly_math::spirograph::trace`.
- [figures/research-spirograph.ts](figures/research-spirograph.ts): research-spirograph
- [figures/research-spirograph.json](figures/research-spirograph.json): research-spirograph census

## Lost

- The loop threshold as the edge of a design's node laws: inside a track with `a < 2b` the seat threshold `(a - b)/b` bites first; the Refuted row in [claims/spirograph-nodes.md](claims/spirograph-nodes.md).
- Three closed forms for the cover between the walls, the full annulus, the summed areas of the distinct curves and the single simple curve's form: each fails; the Refuted rows in [claims/spirograph-walls.md](claims/spirograph-walls.md).
