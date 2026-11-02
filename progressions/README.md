---
title: Progression-free digit designs
lead: A digit design with `0` among its digits is free of three-term progressions exactly when its digit set is and the signed set `F + F - 2F` closes no carry walk, carries living in `{-1, 0, 1}`, and three levels decide; the no-carry case is the digit lift of a progression-free subset of `Z_b`; the census to base `40` by the criterion against brute force finds the half-base construction beaten at bases `19`, `35`, `37` and `39`, every winner a cyclic lift and a carry that leaves `0` never paying, and the dimension `log 2/log 3` of the Stanley sequence first passed at base `37` by ten digits; the planar block of Elsholtz, Hunter, Proske and Sauermann is run as a discrete block problem whose potential is a peeling order, every exact grid optimum bounding every planar block from above, so no planar block has area past `19/49` and the route cannot bring the constant below `2.3381`; the grid search tracks `7/24 + 1/m` to `m = 48` and every grid set found above `7/24` dies under dilation, while peelable blocks lifted by a slice give `r_3(Z_m^n)` the constants `4.362, 5.196, 5.477, 6.633, 7.483, 7.745` at `m = 7, 8, 9, 11, 12, 13`, above every published constant read, and a `36`-point block in `Z_5^3` gives `3.301`, above the `35^(1/3)` read for `F_5^n`.
figure: research-progressions
slug: progressions
date: 2026-10-17
---

# Progression-free digit designs

A digit design with `0` among its digits is free of three-term progressions exactly when its digit set is and the signed set `F + F - 2F` closes no carry walk, carries living in `{-1, 0, 1}`, and three levels decide; the no-carry case is the digit lift of a progression-free subset of `Z_b`; the census to base `40` by the criterion against brute force finds the half-base construction beaten at bases `19`, `35`, `37` and `39`, every winner a cyclic lift and a carry that leaves `0` never paying, and the dimension `log 2/log 3` of the Stanley sequence first passed at base `37` by ten digits; the planar block of Elsholtz, Hunter, Proske and Sauermann is run as a discrete block problem whose potential is a peeling order, every exact grid optimum bounding every planar block from above, so no planar block has area past `19/49` and the route cannot bring the constant below `2.3381`; the grid search tracks `7/24 + 1/m` to `m = 48` and every grid set found above `7/24` dies under dilation, while peelable blocks lifted by a slice give `r_3(Z_m^n)` the constants `4.362, 5.196, 5.477, 6.633, 7.483, 7.745` at `m = 7, 8, 9, 11, 12, 13`, above every published constant read, and a `36`-point block in `Z_5^3` gives `3.301`, above the `35^(1/3)` read for `F_5^n`.

- [notes/progressions.md](notes/progressions.md): Progression-free digit designs
- [notes/carry-criterion-census.md](notes/carry-criterion-census.md): The carry criterion and census
- [notes/planar-block.md](notes/planar-block.md): The door: the planar block
- [claims/carry-criterion.md](claims/carry-criterion.md): The carry criterion
- [claims/block-area.md](claims/block-area.md): The block area
- [claims/four-baselines.md](claims/four-baselines.md): Against four baselines
- [claims/four-term.md](claims/four-term.md): Four-term progressions
- [refs.md](refs.md): the references this topic cites
- [lab/rs/carry-census/](lab/rs/carry-census/): The numbers of progressions
- [lab/py/block-search/](lab/py/block-search/): The discrete block problem
- [lab/rs/block-exact/](lab/rs/block-exact/): The exact grid optimum
- [lab/py/four-baselines/](lab/py/four-baselines/): Against four baselines
- [lab/rs/evasive-saturation/](lab/rs/evasive-saturation/): The saturation of the evasive set
- [lab/rs/layer-spectra/](lab/rs/layer-spectra/): Graphs and spectra
- [figures/research-progressions.ts](figures/research-progressions.ts): research-progressions

## Lost

- A digit design as the largest known progression-free subset of `[0, N)`: the base-37 design passes greedy, Behrend and the Salem-Spencer family at `N = 37^7` and loses to the sets of Wroblewski cited at A065825 by `1.21` to `1.31`; the rows are in [claims/four-baselines.md](claims/four-baselines.md).
- The planar climb at `m = 21, 23, 25` as a lower-bound constant for `r_3(Z_m^n)`: the published block at its best shift holds `136, 164, 191` points against the climb's `132, 156, 182`; the rows are in [claims/block-area.md](claims/block-area.md).
- A planar 4-peel block as a record constant for `r_4(Z_7^n)`: the exact planar optimum is `21` points, constant `sqrt 21 = 4.5826`, below Fuhrer and Taranchuk's `7^(5/6) = 5.0611`; and augmenting their evasive set fails, since every outside point ends a 4-AP with three terms in it; the rows are in [claims/four-term.md](claims/four-term.md).
