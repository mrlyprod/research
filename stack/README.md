---
title: The algebra of the stack
lead: The algebra of the stack: stacking is Dirichlet convolution and the group knows nothing new, the exact spun stack is a Gaussian Farey whose Fourier `L^2` discrepancy obeys Franel's identity one field up and is equivalent to the Riemann hypothesis for `zeta(s) L(s, chi_-4)`, its hexagonal twin puts `L(2, chi_-3)` in a proved node constant, the layers are a dilation system whose symbol is `zeta(1 + s)` so no reweighting moves the line, and the lean is the operation outside the Dirichlet group that still lands a closed form.
figure: research-stack
slug: stack
date: 2026-09-09
---

# The algebra of the stack

The algebra of the stack: stacking is Dirichlet convolution and the group knows nothing new, the exact spun stack is a Gaussian Farey whose Fourier `L^2` discrepancy obeys Franel's identity one field up and is equivalent to the Riemann hypothesis for `zeta(s) L(s, chi_-4)`, its hexagonal twin puts `L(2, chi_-3)` in a proved node constant, the layers are a dilation system whose symbol is `zeta(1 + s)` so no reweighting moves the line, and the lean is the operation outside the Dirichlet group that still lands a closed form.

- [notes/stack.md](notes/stack.md): The algebra of the stack
- [notes/stacking-and-spinning.md](notes/stacking-and-spinning.md): Stacking, selecting and spinning
- [notes/franel-and-eisenstein.md](notes/franel-and-eisenstein.md): Franel and the Eisenstein twin
- [notes/levels-dilations-nodes.md](notes/levels-dilations-nodes.md): Levels, dilations and nodes
- [notes/lean-and-meter.md](notes/lean-and-meter.md): The lean and the novelty meter
- [notes/horocycle-and-moire.md](notes/horocycle-and-moire.md): The design horocycle and the moire limit
- [claims/design-horocycle.md](claims/design-horocycle.md): The design horocycle
- [claims/eisenstein-stack.md](claims/eisenstein-stack.md): The Eisenstein stack
- [claims/franel-one-field-up.md](claims/franel-one-field-up.md): Franel one field up
- [claims/layers-as-a-dilation-system.md](claims/layers-as-a-dilation-system.md): The layers as a dilation system
- [claims/leaning-stack.md](claims/leaning-stack.md): The leaning stack
- [claims/levels-and-designs.md](claims/levels-and-designs.md): Levels and designs
- [claims/moire-local-limit.md](claims/moire-local-limit.md): Moire local limit
- [claims/novelty-meter.md](claims/novelty-meter.md): The novelty meter
- [claims/parity-carpet-stack-spectrum.md](claims/parity-carpet-stack-spectrum.md): The parity-carpet stack spectrum
- [claims/parity-fill.md](claims/parity-fill.md): The parity fill
- [claims/spun-stack.md](claims/spun-stack.md): The spun stack
- [claims/stack-algebra.md](claims/stack-algebra.md): The stack algebra
- [papers/novelty-meter.md](papers/novelty-meter.md): The Stack Hears the Zeros
- [refs.md](refs.md): the references this topic cites
- [lab/py/design-horocycle/](lab/py/design-horocycle/): The design horocycle of stack
- [lab/py/eisenstein-stack/](lab/py/eisenstein-stack/): The hexagonal twin of `lab/py/spun-stack`
- [lab/py/gaussian-franel/](lab/py/gaussian-franel/): Franel's identity one field up
- [lab/py/leaning-stack/](lab/py/leaning-stack/): Leans the line stack
- [lab/py/moire-local-limit/](lab/py/moire-local-limit/): The local limit of the overlay of two `dim 2` designs at sides `N` and `N + h`, both on the unit square, under the difference filter (a point is on when exactly one design inks it), computed exactly from window integrals.
- [lab/py/parity-fill/](lab/py/parity-fill/): Computes the inked fraction of the unit square under the parity blend of the odd carpet stack
- [lab/py/smoothed-novelty/](lab/py/smoothed-novelty/): The novelty meter
- [lab/py/spin-render/](lab/py/spin-render/): Renders the odd parity-carpet stack with each layer spun about the centre, and measures what survives a change of raster.
- [lab/py/spun-stack/](lab/py/spun-stack/): What the line stack does when every layer is rotated
- [lab/py/stack-algebra/](lab/py/stack-algebra/): Stacking a stack is Dirichlet convolution on the layer weights, and this study checks that by literal stacking, then reads off the closed form of every selected stack.
- [lab/py/stack-dilations/](lab/py/stack-dilations/): The parity stack's layers read as a system of dilated functions in `L^2(0,1)`, and the Dirichlet series that is its symbol.
- [lab/py/stack-levels/](lab/py/stack-levels/): The correlation law of two stacked layers when the layers are different designs, or the same design at different levels, computed exactly.
- [lab/py/totient-constant/](lab/py/totient-constant/): The node count of the spun stacks is a totient sum over ideals, and this study proves its constant instead of fitting it.
- [lab/rs/node-error/](lab/rs/node-error/): The error of the spun stacks' node count, recounted to norm bound 8388607 in Rust and read over dyadic windows, Gaussian and Eisenstein.
- [figures/moire-local-limit.ts](figures/moire-local-limit.ts): moire-local-limit
- [figures/moire-local-limit.json](figures/moire-local-limit.json): moire-local-limit census
- [figures/paper-moire-correlation-laws.ts](figures/paper-moire-correlation-laws.ts): paper-moire-correlation-laws
- [figures/paper-novelty-meter.ts](figures/paper-novelty-meter.ts): paper-novelty-meter
- [figures/paper-novelty-meter.json](figures/paper-novelty-meter.json): paper-novelty-meter census
- [figures/research-stack.ts](figures/research-stack.ts): research-stack
- [figures/research-stack.json](figures/research-stack.json): research-stack census

## Lost

- The smoothed novelty meter as a route to the Riemann hypothesis: on a finite range it renders the first zeros, so it can neither certify nor refute the hypothesis; the Refuted row in [claims/novelty-meter.md](claims/novelty-meter.md).
- The design horocycle as a route: no semigroup orbit carries a density-zero design, the horocycle sampled on the design's Cantor set puts the digits on the position and not on the cusp index, and the meter gives back only holomorphy of `Z_F` (The design horocycle).
- Reweighting the layers: they form a dilation system with symbol `zeta(1 + s)`, so no reweighting moves the line (The layers as a dilation system).
- The `N log N` size of the spun stacks' node-count error, the integers' totient error carried over: it climbs where `N^(5/4)` is flat, the circle problem's exponent on the `b = 1` term; the error is `O(N^(4/3))` Proved, `O(N^(5/4))` in mean square Proved, and reads `N^(5/4)` (the final rows of [claims/spun-stack.md](claims/spun-stack.md)).
