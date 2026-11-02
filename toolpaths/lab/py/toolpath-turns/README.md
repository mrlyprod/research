# Toolpath Turns

- Replacement curves on the triangular lattice: a generator of `n` unit steps from `0` to a chord `c` with `N(c) = n`, and a flag `F`, `R`, `M` or `MR` per step; `M` and `MR` only when `c/conj(c)` is a unit.
- `census [orders] [stop]`: every gentle generator walk (turns at most 60 degrees, point-avoiding), every flag assignment gentle by the junction closure and point-avoiding at level 2, expanded to the first crossing, up to level 4; mirrors where the chord allows them; prints the (generator, flags) pairs and, when the order is exhausted, their classes up to reversal and mirror; `stop` cuts after that many walks.
- `plain [orders] [stop]`: the same with `F, R` only.
- `gosper`: the Gosper curve's turn counts, axis counts and `abs(h2)^2` at levels 1 to 6, asserted against `7^(k-1)`, `4*7^(k-1) - 1`, `2*7^(k-1)` and `7^k`; its mirror image at level 2 asserted against the first 49 terms of A229214.
- `bead`: the second harmonic `h2` of Hilbert's curve at levels 1 to 6 and Peano's at levels 1 to 4.
- `rate [orders]`: every curve avoiding points through level 4 at an order, one per class up to reversal and mirror, with its 120-degree turn counts at levels 1 to 5.
- Asserts exit nonzero; nothing is written.

## RUN

- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py census 3,4,7,9,13`: 4 s.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py plain 12`: 1 s.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py census 12 379`: 84 s, order 12 with mirrors cut after walk 379 of 2646.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py plain 16 97426`: 60 s, order 16 with `F, R` cut after walk 97426 of 122964.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py gosper`: under 1 s.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py bead`: under 1 s.
- `uv run python research/toolpaths/lab/py/toolpath-turns/toolpath.py rate 7`: under 1 s; order 9 with mirrors passes 300 s.

## WITNESSES

- `toolpaths.md`, The census: orders 3, 4, 7, 9 (mirrors), 12 (`F, R`) and 13 exhausted, pairs and classes, no gentle curve avoids points at level 3 (`census 3,4,7,9,13`, `plain 12`, Verified); the cuts at order 12 with mirrors and order 16 (`census 12 379`, `plain 16 97426`, Conjecture).
- `toolpaths.md`, The Gosper count: the turn triple at levels 1 to 6 and the mirror image against A229214 (`gosper`, Verified).
- `toolpaths.md`, The bead axis: `abs(h2)^2 = 7^k` for Gosper, `h2 = -1` for Hilbert, `h2 = (9^k - 1)/2` for Peano (`gosper`, `bead`, Verified).
- `toolpaths.md`, The printer's line: the two order-7 classes and their sharp counts `2*7^(k-1)` and `3*7^(k-1)` (`rate 7`, Verified).
