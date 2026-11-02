---
title: Design masks
lead: A Toeplitz attention mask whose offsets form a digit design reaches every offset in two hops exactly when its digit sumset covers the digits; the sparsest ties the strided mask in density and loses to it on tiny Shakespeare, which reads the nearest offsets.
figure: research-masks
slug: masks
date: 2026-10-19
---

# Design masks

A Toeplitz attention mask whose offsets form a digit design reaches every offset in two hops exactly when its digit sumset covers the digits; the sparsest ties the strided mask in density and loses to it on tiny Shakespeare, which reads the nearest offsets.

- [notes/masks.md](notes/masks.md): Design masks
- [claims/design-masks.md](claims/design-masks.md): Design masks
- [claims/offset-masks.md](claims/offset-masks.md): Offset masks
- [lab/py/design-masks/](lab/py/design-masks/): The mask census across the five tile groups, the hex projections and two base-3 digit designs, and the paired sweep against dense, window and random masks at equal density.
- [figures/research-masks.ts](figures/research-masks.ts): research-masks
- [figures/research-masks.json](figures/research-masks.json): research-masks census

## Lost

- Two-dimensional designs laid over the query-key square: at side 81 and 300 steps, nineteen masks from every family are no better than their random twins at equal density on copy and a little worse on recall; a lag task reads only the mask's offset profile, so the route moved to offset masks.
