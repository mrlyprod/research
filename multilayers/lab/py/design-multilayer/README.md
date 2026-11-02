# design-multilayer

- Light at normal incidence through a digit-design multilayer: layer A (`nA`) on the level-`k` design cells (base `b`, digits `D`), host B (`nB = 1.45`) elsewhere and on both sides, every cell of equal optical thickness `delta`.
- Builds the level-`k` characteristic matrix by the level recursion `M_(k+1) = W_D(M_k, M_B(b^k delta))`, `b - 1` multiplications a level; a stack is dark once an entry passes `1e100`, and a dark stack's matrix is rescaled by a positive number each level so the phase of `t_k` stays exact.
- Mahler measures: `m(P_D)` by Jensen on the roots, `m(Q_D)` by Jensen in `x` on `2^15` midpoints in `y`; the Boyd class of `Q_D` (monomial times cyclotomic polynomials in monomials) by exact factoring in sympy.
- One generator for every number in `notes/multilayers.md`.

## RUN

```
uv run python research/multilayers/lab/py/design-multilayer/multilayer.py
uv run python research/multilayers/lab/py/design-multilayer/multilayer.py check born
```

- No argument runs every verb; about 46 seconds on one core. Needs numpy, sympy, mpmath and PARI (`gp`).
- `check` (2.4 s): the recursion against a 40-digit layer-by-layer product on five designs at `nA = 1.6, 2.3, 3.5`, asserting a relative gap under `1e-12`; the zero-contrast derivative against `prod P_D(e^(2i b^j delta))` at `nA = nB(1 + 1e-8)` in `abs(gap)/max(1, abs(prod))`, asserting `1e-5`.
- `born` (2.1 s): per-level mean drift of `ln abs(r/t)` at `nA = 1.455` against `m(P_D)` on six designs, with its standard error and `abs(mean - m(P_D))`, and the largest gap over the six.
- `smyth` (0.3 s): `m(Q_(0,1,3))` against `3 sqrt(3) L(chi_(-3), 2)/(4 pi)` from PARI, asserting `1e-8`.
- `deep` (27 s): on nine designs at `nA = 1.6, 2.3, 3.5`, the deep drift (mean step where `abs(r_k/t_k) < 0.02`) with per-level standard errors, the tail mean of the last four levels with `se` the mean of their standard errors (an upper bound for the tail mean's standard error) and its distance from `m(Q_D)` in `se`, the deep lift and its largest tail gap from the drift, the step-law median gap where `abs(r_k/t_k) < 1e-3` on 50 or more frequencies (`n/a` otherwise), the coupling `abs(mean e^(2i(alpha - beta)))` of block and spacer phases, the torus mean of `ln abs(Q_D(u^2, v^2))` over every frequency, dark or open, and the dark fraction; then the step-law cell count and the zero-set band over the seven table designs.
- `refute` (3.8 s): the three designs with cyclotomic `P_D` and non-Boyd `Q_D`, asserting both, with `f_k sqrt(k)` to level 12 at `nA = 2.3`, beside two Boyd-cyclotomic controls.
- `census` (10.9 s): every design at bases 3 to 6 with `0` in `D` and `2 <= abs(D) < b`, up to shift and mirror (39), its Boyd class and the flatness `f_12 sqrt(12)/(f_6 sqrt(6))` at `nA = 2.3`, verdict critical above `0.9`, and the flatness band of each class.

## WITNESSES

- `notes/multilayers.md` "The level recursion": the `1e-12` agreement (`check`).
- "The Born drift": the derivative identity (`check`) and the weak-contrast drifts `0.3832`, `0.4425`, `0.2818`, within `0.001` of `m(P_D)` (`born`).
- "The block-spacer lift": the step law below `3e-6` on 25 of 27 cells (`deep`).
- "Cyclotomic is not enough": the three counterexamples and the controls `1.097`, `0.901` (`refute`).
- "The deep drift": the table, the refutation at `nA = 1.6`, the zero-set band and the `{0,1,3}` rows (`deep`, `smyth`).
- "Two regimes": the 39-design census, 26 critical and 13 decaying, matching the Boyd class on all (`census`).
