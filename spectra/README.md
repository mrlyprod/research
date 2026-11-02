---
title: Spectra
lead: The tile grammar of the diagonal slice at every odd base: the two-tile claim, the closed forms, and the mod-4 split.
figure: research-spectra
slug: spectra
date: 2026-08-28
---

# Spectra

The tile grammar of the diagonal slice at every odd base: the two-tile claim, the closed forms, and the mod-4 split.

- [notes/spectra.md](notes/spectra.md): Spectra
- [claims/2-adic-smith-cascade.md](claims/2-adic-smith-cascade.md): The 2-adic Smith cascade
- [claims/diagonal-cut.md](claims/diagonal-cut.md): The diagonal cut
- [claims/hexagram-provenance.md](claims/hexagram-provenance.md): Hexagram provenance
- [claims/menger-slice-sources.md](claims/menger-slice-sources.md): Menger slice sources
- [claims/odd-base-slice-grammar.md](claims/odd-base-slice-grammar.md): Odd-base slice grammar
- [claims/slice-at-large-side.md](claims/slice-at-large-side.md): Slice at large side
- [claims/tent-rank-law.md](claims/tent-rank-law.md): The tent rank law
- [refs.md](refs.md): the references this topic cites
- [lab/py/odd-base-slice-grammar/](lab/py/odd-base-slice-grammar/): Computes the centroid diagonal slice of the odd-base parity solid, the solid whose cells have at most one odd coordinate, as a tile census with no raster.
- [lab/py/slice-at-large-side/](lab/py/slice-at-large-side/): Computes the central hexagon `x + y + z = 3(N-1)/2` of the dim 3 parity designs at large odd side and on words of growing odd sides, from the carry block `S_b` of `notes/spectra.md`.
- [lab/py/smith-cascade/](lab/py/smith-cascade/): Builds the even-convention carry matrix `M_even` at base 3, the odd carry parameter `D = 2R + 1`, size `n = R + 1`, from the digit polynomial `P = (sum_k C(D-1,k) x^(2k))(1 + D x + x^2)`, entry `P[c + D - 3c'] + P[-c + D - 3c']`.
- [lab/py/smith-window/](lab/py/smith-window/): Builds the even carry core at base 3, at the odd carry parameter `D = 2R + 1`, through its symbol rather than its entries
- [lab/py/spectra-from-cuts/](lab/py/spectra-from-cuts/): Composes the two-tile claim of `notes/spectra.md` with the carry automaton of `notes/cuts.md` run at an odd base `b` instead of base 3.
- [figures/paper-slice-recurrence-order.ts](figures/paper-slice-recurrence-order.ts): paper-slice-recurrence-order
- [figures/paper-slice-recurrence-order.json](figures/paper-slice-recurrence-order.json): paper-slice-recurrence-order census
- [figures/paper-slice-sign-even-half.ts](figures/paper-slice-sign-even-half.ts): paper-slice-sign-even-half
- [figures/paper-slice-sign-even-half.json](figures/paper-slice-sign-even-half.json): paper-slice-sign-even-half census
- [figures/research-spectra.ts](figures/research-spectra.ts): research-spectra
