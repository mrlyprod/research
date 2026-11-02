---
title: Dimers on a design
lead: Domino tilings of a plane design: the colour imbalance read off the code, tilings that never switch off, codes that tile late, and the carpet's growth constant bracketed.
figure: research-dimers
slug: dimers
date: 2026-10-03
---

# Dimers on a design

Domino tilings of a plane design: the colour imbalance read off the code, tilings that never switch off, codes that tile late, and the carpet's growth constant bracketed.

- [notes/dimers.md](notes/dimers.md): Dimers on a design
- [claims/dimers-on-a-design.md](claims/dimers-on-a-design.md): Dimers on a design
- [refs.md](refs.md): the references this topic cites
- [lab/py/design-dimers/](lab/py/design-dimers/): Domino tilings of a plane design
- [lean/Dimers.lean](lean/Dimers.lean): Proves by `native_decide` that the colour-imbalance law matches the coordinate count of `level` at every code of base 2 through level 6 and base 3 through level 4, and that `bruteCount` is 0 at every code with `s != 0` at base 2 levels 1..4 and base 3 levels 1..2; base 2 through level 3 also by kernel `decide`.
- [figures/research-dimers.ts](figures/research-dimers.ts): research-dimers
- [figures/research-dimers.json](figures/research-dimers.json): research-dimers census

## Lost

- Level 1 deciding whether a plane code tiles some level: base 5 code `19920882` has no tiling at level 1 and `252236412223488` at level 2; the Refuted row in [claims/dimers-on-a-design.md](claims/dimers-on-a-design.md).
