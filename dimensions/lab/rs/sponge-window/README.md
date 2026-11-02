# Sponge Window

- Falsifies the tube of the Menger sponge's plus past `1/6` and the corner volume `Deep` at every radius, as written in [measurability and the sponge](../../../notes/measurability-and-sponge.md), subsection "The tube past `1/6`, and `Deep` in closed form", against the crate's exact distance `mrly_math::three::sponge::distance`.
- Thresholds: `delta_m = sqrt(1/36 + s_m^2/4)`, `s_m = 3^(-m-1)`, for `m = 1..5`, with `sqrt(10)/18` and `sqrt(2)/6`.
- Levels: at 15 radii, eleven spread over the sub-intervals of `(1/6, sqrt(2)/6]` down to `(delta_5, delta_4)`, then `1/6`, `0.148`, `1/8` and `1/12`, it enumerates the crossing holes of the line `across = min(radius, 1/6)` to level 12, integrates each level's pair integral by its own adaptive Simpson rule in `a`, prints the sum against the crate's `deep`, which sums the same integrals by Gauss-Legendre in the angle, and the largest gap over all radii, and rasters the crate's distance on every crossing row to level 4 inside the box the lemma confines it to, `40 x 40` across and `800` along, printing the midpoint volume and the bracket of cells beyond the radius plus or minus the half diagonal; it asserts every bracket contains its level integral.
- Outside the boxes: on a `72^3` raster of one arm it prints the largest distance at a cell centre outside every confinement box, asserted at most the radius, for every radius from `1/6` on.
- Centre: on `[1/6, sqrt(2)/6]` it rasters the centre cube near its middle, `120^3` cells on `[1/3 + c, 1/2]^3` with `c = sqrt(radius^2 - 1/36)`, times 8, against `7/27 - T - 24 Deep`, asserting the bracket.
- Period: `p(eps) = eps^(D-3) (20/27 + sum_(l = 0..40) (27/20)^l T(eps/3^l))` with `T` the crate's `exact`, on 3000 and 6000 steps of `u = ln(1/eps)` over one period from `sqrt(2)/6`; golden-section search for the maximum and the minimum, the swing, the logarithmic mean, the range on `(1/6, sqrt(2)/6]` and whether `p` falls there at every step, and `p` at five phases beside the crate's upper-edge `profile`.

## RUN

- `cargo run --release --manifest-path research/dimensions/lab/rs/sponge-window/Cargo.toml`
- About 8 seconds after the build; prints only, writes nothing.

## WITNESSES

- measurability-and-sponge.md, subsection "The tube past `1/6`, and `Deep` in closed form": the thresholds `delta_2 = 0.167692`, `delta_3 = 0.166781`, `delta_4 = 0.166679`; the level integrals against the crate's `deep`, gap at most `3.2e-21` at all fifteen radii; the row rasters within `6.12e-3` relative of their level integrals on `[1/6, sqrt(2)/6]`, down to `3.988024e-13` at level 4, every bracket containing them; the centre rasters within `2.31e-3` of the closed uncovered part; the largest distance outside the boxes, below the radius; `Deep(1/6) = 1.894504e-6`, `Deep(1/8) = 4.827978e-8`, `Deep(1/12) = 0`, `T(1/6) = 0.2570649366`.
- The same subsection: the maximum `p = 2.139869` at `eps = 0.148577`, the minimum `2.122663` at `0.081380`, the swing `0.8106 %`, the logarithmic mean `2.130510`, `p(1/6) = 2.136706`, `p(1/8) = 2.135739`, and `p` falling across `(1/6, sqrt(2)/6]` to `p(sqrt(2)/6) = 2.122798`.
