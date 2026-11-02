---
title: The hexagon moire
lead: The moire of the stacked diagonal slices: the exact cut-ink laws, the doubling sign law, the quarter-line bands, and the ghost star, whose cell-frame decay is a closed form at every band width, `-1/4` at the arm and `-1/8` in the limit, each layer-pair limit along an affine phase map an exact rational, and no frame-free coefficient.
figure: research-hexagon
slug: hexagon
---

Cut the solid cube of odd side `n` through its centre, perpendicular to the main diagonal, and the section is a regular hexagon; [the slices page](../../slices/notes/slices.md) counts that mesh and the fills that live on it. Do it at every odd side at once, hold the layers up to the light, and the hexagons interfere. This page is the arithmetic of that interference: which lines the stack lights, which pairs of layers agree, how the agreement decays, and which constants survive to the limit. The layers are the four historical families of [slices](../../slices/notes/slices.md) - carpet, net, tree and void - and everything below is either an exact statement about one layer or an exact statement about a pair.

One generator prints every number on this page, `lab/rs/hexagon-moire`, which runs in under a minute and whose sections are the sections below; a number it does not print is not on this page, and the one place a shelf paper is the source instead says so. [The three app](/three/) stacks a cube design into a solid field and cuts it on any plane, the central diagonal cut being this hexagon, and [the slices demo](/three/) draws one layer triangle by triangle.

Three frames run through the page and they are not the same frame. The **ideal frame** is the cut plane's own coordinate square: a point is `(X, Z)` in `[0,1]^2`, the layer's cell is `x = floor(4nX)`, `z = 2*floor(2nZ)`, `y = 6n - 2 - x - z`, the point is on the layer when `0 <= y < 4n`, and every layer is resampled onto one common grid. The **lattice frame** is the rendered raster of the cube's own cells. The **cell frame** is the cube's own cell lattice read one layer at a time, nothing resampled and no raster, every quantity a ratio of exact cell counts. The plane `x + y + z = 6n - 2` sits half a cell off the true centre, so each layer's hexagon is displaced by `1/(2n)` from the last, and that displacement is real geometry rather than an artifact: a scan raster that ignores it misreads the crosshair strengths, and one quantity below - the star's decay coefficient - is a different number in each of the three.

## The cut ink of one layer

**The cut ink laws of all four families are exact closed forms. Proved.** With `chi = (-1)^((3n-1)/2)`, the inked share of the hexagon at odd side `n` is

```
carpet  1/2 + chi/8 + 1/(2n) - chi/(8n^2)
net     1 - carpet
tree    1/4 + (1/3 - chi/12)/n + (1 - chi)/(6n^2)
void    1/4 - chi/(4n) + 1/(2n^2)
```

and the generator checks each against the counted hexagon in exact rational arithmetic: 28 of 28 layers match for all four families at every odd `n <= 55`. The character `chi` is the same `chi_4` the [Walsh spectrometer](../claims/diagonal-slice-stack.md) proves the quasipolynomial by, and the `1/n` and `1/n^2` orders of that law are proved on [slices](../../slices/notes/slices.md) under the slice ink. The wider range `n <= 101` has no generator, so nothing is claimed there.

The centre cell is the first place the families separate. Over the same 28 layers the carpet's centre cell is ink in 14, the net's in 14, the tree's in 14, and the void's in 28 of 28 - the void's centre is ink at every odd `n`, by a two-line [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) argument on the three coordinates. That single cell is the seed of the two star sections below: the carpet's centre agrees with the background and fades, the void's never does. **Verified.**

## Layer pairs and the doubling law

Correlate two layers on the full hexagon, exactly, cell by cell. The seven pairs the generator reads blind are

| pair | `r` | `cov` |
|---|---:|---:|
| `(5,9)` | -0.14179450 | -0.03500978 |
| `(5,7)` | -0.08542646 | -0.01976104 |
| `(9,17)` | -0.12773715 | -0.03098711 |
| `(13,25)` | -0.12529193 | -0.03011364 |
| `(17,51)` | +0.21659159 | +0.05119124 |
| `(5,15)` | +0.21710273 | +0.05156877 |
| `(3,9)` | +0.20367849 | +0.04224750 |

(Verified.) The coprime pairs are not zero, and that is the whole point: on the square stack they are exactly zero, and the reason the hexagon loses it is a theorem of its own, stated under WHERE THE EXACTNESS STOPS.

**The corrected law on the hexagon is dyadic. Proved.** What breaks coprime independence is a hidden half-cell-shifted overtone at doubled frequency, `chi_4(n) s(2nX + 1/2)/8`, that the plane constraint forces into every carpet slice, plus the hexagon's non-product tent marginal; together they couple layer `m` to layers `2m +- 1` and `m +- 2` regardless of [gcd](https://en.wikipedia.org/wiki/Greatest_common_divisor). The sign is forced: `sign r(m, 2m +- 1) = -chi_4(m) chi_4(2m +- 1)`, and the generator holds it on 18 of 18 pairs from `(3, 5)` to `(601, 1201)` across all four residue branches. Two of those rows, the largest computed, read `(501, 1001)` at `r = -0.11732362` and `(601, 1201)` at `r = -0.11729091`.

**The doubling magnitude reads between `0.11711630` and `0.11715991`. Verified.** Richardson extrapolation `r = r_inf + a/m + b/m^2` on sliding triples of `m = [157, 201, 301, 401, 501, 601]` returns `r_inf = -0.11715991, -0.11711630, -0.11714296, -0.11711630`, and the two branch extrapolations in `1/m` land on `0.1171270` at `(501, 601)` and `0.1171274` at `(103, 403)`. That kills the rival `19/162 = 0.11728395`, which sits above every reading in that spread and above both branch extrapolations; the closed form below kills it outright.

**The exact doubling constant is `253/2160`, and the smear integral behind it is closed. Proved.** Write layer `N`'s cut cell at hexagon point `(X, Z)` as `x = floor(4NX)`, `z = 2 floor(2NZ)`, `y = 6N - 2 - x - z`. Its block parities are `s_x = floor(NX) mod 2`, `s_z = floor(NZ) mod 2` and, exactly, `s_y = s_x + s_z + w mod 2` with `w = floor((6N - 2 - p - 2q)/4) mod 2` for `p = floor(4 frac(NX))` and `q = floor(2 frac(NZ))`. That `w` is the half-cell overtone in closed form: it is `1` at exactly the two phase cells `(p, q) = (0, 0)` and `(3, 1)` when `N = 1 mod 4` and its complement when `N = 3 mod 4`, so the plane constraint is the entire coupling and the third coordinate is never free. For `n = 2m + 1` the identity `nX = 2 mX + X` puts layer `n`'s phase at `(2 frac(mX) + X) mod 2`: the two layers' phases lie on one closed geodesic of the torus instead of filling it, which is why coprimality buys nothing on this cut, and as `m` grows `(mX mod 2, mZ mod 2)` equidistributes over the hexagon at rate `O(1/m)`, since the region is a fixed polygon and every nonzero frequency integrates to `O(1/m)`. What is left is a piecewise-constant integral over `(X, Z)` and the phase, with all breakpoints rational, and it evaluates in exact rational arithmetic to covariance `253/9216`. Every layer's variance is `15/64`, so `r = 253/2160` at all four branches, carrying the sign law's sign. The generator reads `-0.11745304` at `(301, 601)` and `-0.11729091` at `(601, 1201)` against `-253/2160 = -0.11712963`, the gap times `m` holding at `-0.0974` and `-0.0969`; those are both `m = 1 mod 4`, and the other class converges to the same magnitude about twice as slowly, `(103, 205)` reading `+0.11914004` and `(203, 405)` reading `+0.11814528` against `+253/2160` with the gap times `m` at `+0.207` and `+0.206`. The rival `19/162` is dead.

## The other pair families

**The breakage is rule-specific in the limit. Proved.** Run the same exact correlation across the four families:

| family | doubling `(201,401)` | adjacent `(199,201)` | adjacent `(249,251)` | echo `(67,201)` | echo `(99,297)` |
|---|---:|---:|---:|---:|---:|
| carpet | -0.11761160 | -0.08115576 | -0.08150224 | +0.21473346 | +0.21476417 |
| tree | +0.00049938 | -0.07016649 | -0.07025454 | +0.14984772 | +0.14929980 |
| void | +0.00107260 | +0.00786627 | +0.00824952 | +0.07721358 | +0.07619141 |

So the persistent doubling coupling is carpet and net only - their doubling covariance is `+-253/9216` and the tree's and the void's are exactly `0` - the tree keeps a neighbour coupling at `(m, m+2)` and no doubling, and the void is nearly independent at both. The gcd echo survives in all four. The controls agree that this is signal and not raster noise: at the coprime pairs `(101, 173)` the three families read `-0.00248023`, `+0.00033407` and `+0.00018183`, and at `(97, 251)` they read `+0.00212976`, `+0.00022204` and `+0.00008585`.

**Each pair family on this page is one integral with a different phase map, and each limit is rational. Proved.** The claim is about a fixed affine phase map `n = a m + c` with `m` growing inside one class mod 4, measured on the mask this page always uses, the full hexagon of the common cut with exact area weighting; a general pair `(m, n)` has no limit theorem here. The map is `n = a m + c`, which sends layer `m`'s phase `alpha = mX mod 2` to layer `n`'s `(a alpha + c X) mod 2`; the doubling pairs are `a = 2`, `c = +-1`, the adjacent pair is `a = 1`, `c = 2`, and the gcd echo is `a = 3`, `c = 0`, where the phase map forgets `(X, Z)` entirely and the hexagon's shape drops out. The same exact integral run on the four families' corner tables gives

| family | doubling `(m, 2m+-1)` | adjacent `(m, m+2)` | echo `(m, 3m)` |
|---|---:|---:|---:|
| carpet | `+-253/2160` | `-11/135` | `+29/135` |
| net | `+-253/2160` | `-11/135` | `+29/135` |
| tree | `0` | `-61/864` | `+4/27` |
| void | `0` | `+7/864` | `+2/27` |

with covariances `+-253/9216`, `-11/576`, `+29/576` for carpet and net, `0`, `-61/4608`, `+1/36` for the tree and `0`, `+7/4608`, `+1/72` for the void. The tree's and the void's doubling zeros are exact and not small: the integral returns their joint ink as `1/16` on the nose against a background of `1/4` each, so the doubling coupling is a carpet and net theorem and nothing wider. The measured rows above sit `O(1/m)` from these: `(249, 251)` reads `-0.08150224` against `-11/135` and `(99, 297)` reads `+0.21476417` against `+29/135`, gaps of `2.1e-5` and `5.1e-5` with `gap * m` at `-0.0052` and `-0.0050`.

## The quarter line and the crosshairs

Stack the carpet's layers in the ideal frame at `N = 55`, 3601 samples per axis, and read one-sided bands of width `0.018` against the hexagon mean `0.543621`:

| line | left | right |
|---|---:|---:|
| `X = 0.2500` | +0.066465 | -0.090881 |
| `X = 0.7500` | -0.090863 | +0.066539 |
| `X = 0.3333` | +0.020969 | -0.058390 |
| `X = 0.2000` | +0.023744 | -0.022475 |
| `X = 0.5000` | -0.001222 | -0.001219 |
| `Z = 0.2500` | -0.010910 | -0.011128 |
| `X-Z = 0.2500` | +0.002907 | -0.067923 |
| `X+Z = 1.2500` | -0.045236 | +0.028988 |

**The quarter-line law. Verified.** The strongest interior lines of the stacked hexagram sit at quarter-cell coordinates `a/4`, generally `a/(4b)` with `b` odd, and the step across such a line is an exact one-sided `+-1/8` that every layer votes for identically, because the overtone's `chi_4` sign meets the layer's own `chi_4` and squares away. Read straight at the line the plateau converges: mean magnitude `0.12210, 0.12340, 0.12412, 0.12446` at `N = 151, 301, 601, 1201` against the limit `1/8`, the one-sided pair at `N = 1201` reading `+0.12320` and `-0.12573`. The odd-fraction crosshairs at `(-1)^a/(4q)`, `q` odd, follow behind. The `X`, `Z` and `X+Z` profiles are numerically identical on the render by the slice's permutation symmetry, so "in five directions but never horizontal" is false, and the missing-`Z`-overtone statement holds only in the rectangle-cell frame.

**The three 60-degree crosshair families obey a limit law. Verified.** The line at coordinate `a/q` carries strength `(-1)^a/(4q)` for odd `q` and nothing for even `q`. In the coarse cut model, 200003 samples along a line, the excesses converge:

| `A` | `N = 55` | `N = 5555` | predicted |
|---|---:|---:|---:|
| `1/3` | -0.093077 | -0.083726 | -0.083333 |
| `2/3` | +0.065218 | +0.082917 | +0.083333 |
| `1/5` | -0.059505 | -0.050272 | -0.050000 |
| `2/5` | +0.032526 | +0.049651 | +0.050000 |
| `1/7` | -0.058065 | -0.036036 | -0.035714 |
| `1/2` | -0.014353 | -0.000217 | +0.000000 |
| `1/4` | -0.023394 | -0.000231 | +0.000000 |

The `N = 55` column is the whole lesson about frames: `1/7` reads `-0.058065` where the `N = 5555` reading is `-0.036036` against a predicted `-0.035714`, and the even denominators read `-0.014353` and `-0.023394` where the prediction is zero. That gap is the per-layer registration drift of `1/(2n)`, and it is what a drifted scan raster misreads. In registration-correct frames the rays are visible in the real render as well, the `X + Z = 1.25` line at `N = 55` reading `-0.045236` and `+0.028988` and the `X = 1/3` one-sided bands reading `+0.020969` and `-0.058390`. A null claiming no rays at all was about the scan geometry and not about the object.

## The void star

**The void slice stack keeps its star forever. Verified.** In the ideal frame at `N = 55`, 3601 samples per axis, centred bands of half-width `0.004`, against a hexagon mean of `0.278670` whose limit is `1/4`:

| line | ink | ratio to mean |
|---|---:|---:|
| `X = 0.50` | 0.513667 | 1.843 |
| `Z = 0.50` | 0.399416 | 1.433 |
| `X+Z = 1.00` | 0.486492 | 1.746 |
| `X-Z = 0.00` | 0.491844 | 1.765 |
| `2X+Z = 1.50` | 0.491525 | 1.764 |

The law behind the table is six central lines at plateau ink `1/2` against a background of `1/4`, ratio 2, every layer voting on all six, with the model-frame `Z` arm weaker at `3/8` - the `Z = 0.50` row - and a centre dot that is ink at every odd `n`, the `28/28` of the first section. Nothing about this decays. The carpet's star, by contrast, is a finite-layer artifact and fades; the two snowflake stacks differ by a theorem, not by a constant.

**Void and carpet have complementary line spectra on the cut. Verified.** The void's lines sit at even-denominator twisted positions `X = a/(2b)` with `b` odd, exactly where the carpet is silent, and the void is silent at the carpet's odd rationals. The tree carries the only untwisted crosshair family plus a permanent ratio-2 line at `2X + Z = 3/2`, and ignores its free axis. The net is the exact pixelwise complement of the carpet, since "at most one odd" and "at least two odd" exhaust the cases - the same partition that makes carpet and net tile one hexagon in [slices](../../slices/notes/slices.md).

## The twist that kills the ray family

**The `chi_4` twist kills the pair-resonance ray family. Verified.** Average `chi_4(n) T(nx)` over odd `n <= N` with `T` the layer's tent, and the average falls like `1/N` at every `x`, rational or not:

| `x` | `N = 55` | `N = 2001` | `N = 40001` | untwisted at `N = 40001` |
|---|---:|---:|---:|---:|
| `0` | +0.00000 | +0.000999 | +0.0000500 | +1.000000 |
| `1/3` | +0.04762 | +0.000333 | +0.0000833 | -0.111111 |
| `2/3` | -0.04762 | -0.000333 | -0.0000833 | +0.111111 |
| `1/5` | +0.02857 | +0.000599 | +0.0000300 | -0.039968 |
| `1/7` | -0.00000 | +0.001284 | +0.0000357 | -0.020363 |
| `1/2` | +0.00000 | +0.000000 | +0.0000000 | +0.000000 |
| `1/4` | +0.00000 | +0.000500 | +0.0000250 | +0.000025 |
| `1/9` | -0.00000 | +0.000333 | +0.0000167 | -0.012294 |
| `0.414214` | +0.09630 | +0.000421 | +0.0000443 | -0.000020 |

The last row is irrational and behaves like the rest. So the snowflake stack has no analogue of the carpet stack's bright main diagonal: in the coarse cut model the `A = C` excess goes `-0.003580` at `N = 55` to `-0.000437` at `N = 555` to `-0.000052` at `N = 5555`, with the background falling `0.550865, 0.507266, 0.500942` to its limit `1/2`. What is left is only the three single-wave crosshair families parallel to the hexagon's edge directions, one per lattice axis, at odd-denominator rational coordinates - the families the previous section measures.

## The ghost star

**The ghost star at the hexagon's centre is a finite-layer artifact. Verified.** Each layer's centre is entirely ink or entirely paper, flipping with `n mod 4`, so 28 layers give exactly `1/2`, which is also the limiting background; the star is the difference between the two, and it goes to zero. On a 1200 by 2399 raster of the carpet's cut layers at odd `n <= 111`, with 2147316 hexagon pixels of which 25440 are star and 1981584 background:

| layers | star | background | star minus background |
|---|---:|---:|---:|
| 5 | 0.57832 | 0.67217 | -0.09385 |
| 28 | 0.51730 | 0.54832 | -0.03102 |
| 56 | 0.50873 | 0.52629 | -0.01756 |

The decay is `(ln L)/L` for `L` the number of layers stacked, and the coefficient is where the frame bites. In the ideal frame, main diagonal `X = Z`, band half-width `0.01`, 2801 samples per axis, `excess * L` runs `-0.7779, -0.9942, -1.1194, -1.2519` at `L = 28, 100, 200, 400`, with slopes against `ln L` of `-0.1807` from 100 to 200 and `-0.1911` from 200 to 400. In the lattice frame, the band `|x - y| <= 0.01` of the cube side read against the exact ink law layer by layer, `excess * L` runs `-0.9337, -1.0212, -1.1387, -1.1923, -1.2791, -1.3645` at `L = 14, 28, 56, 100, 200, 400`, with slopes `-0.1252` from 100 to 200 and `-0.1233` from 200 to 400. The [star app](/star/) stacks the hexagonal cuts one per odd side and draws the arm, whose ink is exactly `1/2 + chi_8(n)/(2n)` and whose decay is exactly `-1/4` on the cube's own cells and a different number in every frame that resamples or widens it.

**The cell frame, defined.** Registration: the cube's own cell lattice, one layer at a time, nothing resampled onto a common grid, so the half-cell displacement `1/(2n)` between consecutive layers is carried rather than averaged away. Band normalisation: the star intensity is the ink share of the exact arm `x = y` of the cut, a diameter of `2n` cells with no width to choose, and the background intensity is the layer's whole-hexagon ink share, the Proved closed form of the first section; the excess is star minus background, and the `L`-layer excess is its mean over the first `L` odd layers.

**The star arm's ink is exactly `1/2 + chi_8(n)/(2n)`. Proved.** The arm `x = y` at odd `n` is the `2n` cells with `z = 6n - 2 - 2x`, and `0 <= z < 4n` pins `n <= x <= 3n - 1`. On the arm the two block parities agree, so "at most one of three odd" collapses to `floor(x/4)` even, and with `f(m) = 4 floor(m/8) + min(m mod 8, 4)` the inked count is `f(3n) - f(n) = n + chi_8(n)`, where `chi_8` is the real character mod 8 of `Q(sqrt 2)`, `+1` at `n = 1, 7` and `-1` at `n = 3, 5`. The generator checks that rational against the counted arm at every odd `n <= 2001`, 1001 of 1001, and the three arms `x = y`, `y = z`, `z = x` return the same rational at every odd `n <= 221`, 111 of 111, by the slice's permutation symmetry.

**In the cell frame at even `L` the decay coefficient is exactly `-1/4`, with an exact constant beside it. Proved.** Writing `excess_L` for the `L`-layer excess, at even `L`

```
L * excess_L = -(ln L)/4 + C + O(1/L^2)
C = ln(1 + sqrt 2)/(2 sqrt 2) - G/8 - gamma/4 - (ln 2)/2 = -0.2937605857
```

| `L` | `excess * L` | `excess * L + (ln L)/4` | residual `* L^2` |
|---|---:|---:|---:|
| 28 | -1.126964 | -0.2939128437 | -0.11937029 |
| 100 | -1.445065 | -0.2937725615 | -0.11975831 |
| 400 | -1.791627 | -0.2937613344 | -0.11978958 |
| 1600 | -2.138200 | -0.2937606325 | -0.11979154 |
| 6400 | -2.484774 | -0.2937605886 | -0.11979188 |

The proof is one subtraction of closed forms. The per-layer excess is `-chi/8 + (chi_8(n) - 1)/(2n) + chi/(8 n^2)`; summed over the first `L` odd `n` the first term is `0` at even `L` and `+1/8` at odd `L`, the character part converges to `L(1, chi_8)/2 = ln(1 + sqrt 2)/(2 sqrt 2)` by the class number formula for `Q(sqrt 2)`, the third converges to `-G/8` with `G` Catalan's constant, and `-sum 1/(2n)` over the odd `n` below `2L` is `-(ln L)/4 - gamma/4 - (ln 2)/2 + O(1/L^2)`, the `1/(4L)` terms of the two harmonic sums cancelling exactly. Every limb of `C` is a closed form that subtraction hands over, so `C` is derived and not a decimal recognised after the fact, and the whole of the ghost star's decay is the `1/(2n)` of the hexagon's ink law read against an arm that has none. The generator prints `L(1, chi_8) = 0.62322524` from its own series against `ln(1 + sqrt 2)/sqrt 2 = 0.62322524`.

**Even `L` is a hypothesis and not decoration. Proved.** The `-chi/8` term is the whole of the difference: at odd `L` it sums to `+1/8`, so the law reads `C + 1/8` in place of `C`, and what is left is `+-1/(4L)` rather than `O(1/L^2)`. The `-1/4` coefficient is untouched by parity. The generator's odd ladder, columns `L * excess_L + (ln L)/4 - C` and that miss less `1/8` times `L`:

| `L` | miss | `(miss - 1/8) * L` |
|---|---:|---:|
| 27 | +0.115715 | -0.250708 |
| 99 | +0.122472 | -0.250244 |
| 401 | +0.125623 | +0.249934 |
| 1601 | +0.125156 | +0.249984 |
| 6399 | +0.124961 | -0.250004 |

**The `1/L^2` term is exact and reads `L` mod 4. Proved.** Three tails carry it: the character tail `-(1/2) sum chi_8(n)/n` over odd `n > 2L` is `-1/(8L^2)` at `L = 0 mod 4` and `+1/(8L^2)` at `L = 2 mod 4`, since the sign pattern `+--+` on the four odd residues starts at `n = 2L + 1` and flips with that residue; the Catalan tail contributes `+1/(64 L^2)` and the harmonic remainder `-1/(96 L^2)`, both blind to the residue. So the coefficient is `-1/8 + 1/64 - 1/96 = -23/192` at `L = 0 mod 4` and `+1/8 + 1/64 - 1/96 = +25/192` at `L = 2 mod 4`. The residual column above is every row at `L = 0 mod 4` and falls by `4.00` at all six doublings from `L = 100` onto `-23/192 = -0.11979167`, the `L = 6400` row overshooting at the floor of double precision where the residual itself is `3e-9`; the generator's `L = 2 mod 4` ladder reads `+0.13020778, +0.13020819, +0.13020825` at `L = 802, 1602, 3202` against `+25/192 = +0.13020833`.

**There is no frame-free decay coefficient. Refuted.** Widen the star from the exact arm to a band of half-width `W` cells about it, change nothing else, and the coefficient walks over a family:

| `W` | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 16 | 20 | 24 | 32 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| coefficient | `-1/4` | `-1/6` | `-1/10` | `-3/28` | `-5/36` | `-3/22` | `-3/26` | `-9/68` | `-5/42` | `-13/100` | `-17/132` | `-33/260` |

Every entry is the exact rational the width family below hands over, and the measurement converges on it: recomputed to `L = 3200`, four times the sweep, `W = 0` through `12` land on `-1/4, -1/6, -1/10, -3/28, -5/36, -3/22, -3/26` within `5.1e-8`, every gap falling by four per doubling of `L` from `L = 800`, where the sweep's own gap is `8.1e-7`. **Proved.** Beside them a band whose half-width is a fixed fraction of the cube side rather than a fixed number of cells reads `-0.1252` and `-0.1233`, and the ideal frame, which resamples every layer onto one 2801-point grid, reads `-0.1807` and `-0.1911` and is still moving at `L = 400`. Resample the layers or do not, and pick the star's width: more than a dozen coefficients come out of one star, so the coefficient is a property of the measurement and not of the stack. What survives the frame change is the `(ln L)/L` order, the sign, and the largest coefficient of them all, the cell frame's `-1/4`, which is the family's first row and the only one a band of zero width can read.

**The width family is exact at every half-width, and the arm is its first row. Proved.** Widen the star from the arm to the band `|x - y| <= W` cells and write `K = floor(W/2)`: `x - y` is even on the cut, so `W` and `W - 1` read one point set at odd `W` and `K` is the only width there is. With `b = 1` when `floor(K/2)` is even and `0` when it is odd, and with `E(K) = #{|j| <= K : j = 3, 4, 5 mod 8} + floor((K + 2)/4) - K`, which is the 8-periodic run `0, -1, -1, 0, 1, 2, 2, 1` at `K = 0..7 mod 8`, the band's excess over the hexagon's ink law at every odd side `n >= K` is exactly

```
excess_W(n) = kappa chi + (m + q chi_8(n))/n + chi/(8 n^2)
kappa = -(-1)^K/(8(2K + 1))
m     = -(K + b)/(2(2K + 1))
q     = (1 - 2 E(K))/(2(2K + 1))
```

with `chi = (-1)^((3n-1)/2)` the ink law's character, `-1` at `n = 1 mod 4`, and `chi_8` the arm's, the real character mod 8 of `Q(sqrt 2)`. Below `n = K` the band is clipped by the hexagon's edge and the identity is false, `W = 6` at `n = 1` missing by `2/7`. At `W = 0` it is `K = 0`, `b = 1`, `kappa = -1/8`, `m = -1/2`, `q = 1/2`, which is the arm's proved excess `-chi/8 + (chi_8(n) - 1)/(2n) + chi/(8 n^2)` term for term, so the arm is not a separate theorem. Summed over the first `L` odd layers the `chi` term vanishes at even `L` and the `1/n` term hands the decay coefficient `m/2 = -(K + b)/(4(2K + 1))`, which at even `W` is the `-(W + 2b)/(8(W + 1))` the sweep predicted, at odd `W` is `-(W - 1 + 2b)/(8W)` - `-1/4` at `W = 1` and `-1/6` at `W = 3`, where the even-`W` formula was false - and tends to `-1/8` as `W` grows, the value the fixed-fraction band of the lattice frame measures.

**The proof is a block-tail sum.** On the cut plane the band's cells are `x = B + j`, `y = B - j`, `z = 6n - 2 - 2B` with `B` running over the `2n` integers of `[n, 3n)` and `|j| <= K`, and nothing is clipped exactly when `K <= n`. The space rule is `s(x) + s(y) + s(z) <= 1` for the block parity `s(c) = floor(c/4) mod 2`, so a row is ink at `2K + 1 - U(B)` cells when `s(z) = 0` and at `P(B)` cells when `s(z) = 1`, where `U(B)` counts the `j` with `s(B + j) = s(B - j) = 1` and `P(B)` the `j` with both `0`; writing `D = 2K + 1 - U - P` for the `j` where the two block parities disagree, the band's ink is `sum_B (2K + 1 - U(B)) - sum_B s(z_B) D(B)`. The second sum collapses to one product. `s(z_B) = 1` exactly on `B = 0, 3 mod 4` when `n = 1 mod 4` and exactly on `B = 1, 2 mod 4` when `n = 3 mod 4`, while `D` takes only two values, `D_a = 2(K - floor(K/4))` on `B = 0, 3 mod 4` and `D_b = 2 floor((K + 2)/4)` on `B = 1, 2 mod 4`, because the disagreement pattern in `j` is the 8-periodic `0, 1, 1, 1, 0, 1, 1, 1` in the first case and `0, 0, 1, 0, 0, 0, 1, 0` in the second. The cut's character and the band's parity are the same partition of `B mod 4`, so `s(z) D` is a constant times an indicator and `sum_B s(z_B) D(B) = (n - 1) D_chi` exactly, with `D_chi = D_a` at `chi = -1` and `D_b` at `chi = +1`. The first sum is a block tail: `U` is 8-periodic in `B`, and for any 8-periodic `h` the tail `sum_{B=n}^{3n-1} h(B) - 2n mean(h)` equals `psi_h(3n mod 8) - psi_h(n mod 8)` for `psi_h` the periodised partial sum, which flips sign between `n = 1` and `n = 3 mod 8` and between `n = 5` and `n = 7`, so it carries no trivial and no `chi_-8` component and is a combination of `chi_4` and `chi_8` alone. Three identities close it. `(D_a + D_b)/2 = A(K)` for `A(K) = sum_{j=1}^{K} w(j mod 4)` with `w = 0, 1, 2, 1`, and `A(K) = K + 1 - b` by the four cases of `K mod 4`, which is `b`'s definition read backwards and is where the `-(K + b)` comes from. `2 mean(U) + A(K) = 2K + 1`, because `w` and the mean pattern `2, 1, 0, 1` sum to `2` at every residue, which is why the constant term is `kappa chi` with no drift. The tail then splits into its two characters by one shift. `s(c + 4) = 1 - s(c)` turns `P` into `U` four steps along, `P(B) = U(B + 4)`, so `U(B) + U(B + 4) = 2K + 1 - D(B)`. Its `chi_4` weight is `(e_1 + e_5)/2` for `e_r = 2(U(r) - mean U)`, while `D_chi` carries `(D_a - D_b)/2`; since `e_1 + e_5 = D_a - D_b = K + (K mod 2)` both weights are `(K + (K mod 2))/2` and cancel identically at every `K`, which is why the `1/n` term has no `chi_4`. Its `chi_8` weight is `(e_1 - e_5)/2 = U(1) - U(5)`, and `U(1)` reads straight off the pattern as the count of `j` in `[-K, K]` with `j = 3, 4, 5 mod 8`, so `U(1) - U(5) = 2 U(1) - (2K + 1) + D_b` and `q = (U(5) - U(1))/(2(2K + 1)) = (1 - 2E(K))/(2(2K + 1))`. `E` is 8-periodic because a block of eight adds `6 + 2 - 8 = 0`, and the generator holds it against the run at `K = 0..200`, 201 of 201. Nothing in the family is fitted.

**The family's constant and both its `1/L^2` branches are closed forms at every width. Proved.** At even `L`

```
L * excess_L = (m/2) ln L + C_W + O(1/L^2)
C_W = m (ln 2 + gamma/2) + q L(1, chi_8) - G/8 + Delta_W
```

with `L(1, chi_8) = ln(1 + sqrt 2)/sqrt 2` by the class number formula for `Q(sqrt 2)`, `G` Catalan's constant, and `Delta_W = sum over odd n < K of (counted excess - the identity)`, the exact rational the finitely many clipped layers contribute, `0` through `W = 5` and `2/7` at `W = 6`. The `chi_4` cancellation is at the `1/n` order only, so `L(1, chi_4) = pi/4` is absent at every width while `L(2, chi_4) = G` sits in every one: `C_W` carries `gamma`, `ln 2`, `L(1, chi_8)` and Catalan's `G`, and widening the band moves only `m` and `q`. Three tails carry the `1/L^2` term. The character tail `-q sum chi_8(n)/n` over odd `n > 2L` is `-q/(4L^2)` at `L = 0 mod 4` and `+q/(4L^2)` at `L = 2 mod 4`, since the sign pattern `+--+` on the four odd residues starts at `n = 2L + 1` and flips with that residue; the Catalan tail of the background's `chi/(8 n^2)` gives `+1/(64 L^2)` at every width, blind to the residue; and the harmonic remainder of `m (H_{2L} - H_L/2)` gives `+m/(48 L^2)`, the `1/(4L)` parts of the two harmonic sums cancelling exactly. So the coefficient is `-q/4 + 1/64 + m/48` at `L = 0 mod 4` and `+q/4 + 1/64 + m/48` at `L = 2 mod 4`, which at `W = 0` is `-23/192` and `+25/192`; at odd `L` the constant is `C_W - kappa` with an `O(1/L)` error, `C + 1/8` at `W = 0`. The sliding-window slope the sweep reads has a residue trap of its own: the window from `L/2` to `L` cancels the `chi` term only at `L = 0 mod 4`, and at `L = 2 mod 4` it converges to `m/2 + kappa/ln 2` instead, `-0.4303` in place of `-1/4` at `W = 0` and `-0.1066` in place of `-1/6` at `W = 2`, so a sweep that never leaves one class of `L` cannot see its own offset. The generator checks the per-layer identity against the counted band in exact rational arithmetic at 14 distinct half-widths from `W = 0` to `W = 64`, every odd `n` from `K` to `201`, and reports the four residue classes `n = 1, 3, 5, 7 mod 8` separately: 1354 of 1354, no class short. The sweep also prints nine odd widths, which repeat their even neighbours cell for cell and are not counted twice. Seven summed ladders confirm the constant to `1e-9` at `L = 1600` and both `1/L^2` branches to `1e-7`; `W = 0` returns `C = -0.2937605857`, `-23/192` and `+25/192` unchanged, and the sliding window reads `-0.24999980` at `L = 1600` against `-0.43078703` at `L = 1602`, whose limit is `-0.43033688`, the `4.5e-4` gap being the odd-`L` error at `L/2 = 801`.

## The constants

**The slice stacks' surviving constants are Leibniz, odd Basel and Catalan. Proved.** Every character series below runs over the exact ink laws of the first section, so the limits are theorems about those closed forms and not fits; the carpet split converges at balanced layer counts, and `N = 55` is 28 layers.

| quantity | closed form | value | measured |
|---|---|---:|---:|
| carpet split `M (I1 - I3 + 1/4)` | `pi/4 + pi^2/32` | 1.0938233009 | 1.0826656664 at `N = 55` |
| void background `M (mean - 1/4)` | `(pi + pi^2)/16` | 0.8131998159 | 0.8065045723 at `N = 55` |
| tree background `M (mean - 1/4 - eps/3)` | `pi/48 + pi^2/48 + G/6` | 0.4237275377 | |
| flat-stack control | `pi^2 ln 2 / (7 zeta(3))` | 0.8130217042 | |

The third and fourth rows are the reason the void's constant has to be stated to eight places: `(pi + pi^2)/16 = 0.8131998159` and `pi^2 ln 2/(7 zeta(3)) = 0.8130217042` collide in the fourth decimal and are different numbers, the second belonging to the flat square stack and not to this one.

**One law prints every row of that table, and the layer count's parity is the only residue it reads. Proved.** Write a family's ink law of the first section as `I(n) = A + B chi + (c + d chi)/n + (e + f chi)/n^2` with `chi = -chi_4(n)`, let `mean` be the average of `I` over the first `M` odd sides `n = 1, 3, ..., 2M - 1`, and let `eps` be the average of `1/n` over the same layers. Subtracting `A` and `c eps` removes exactly the two limbs that do not converge, and what is left is termwise three character sums and one zeta tail, so at every `M`, exactly,

```
M (mean - A - c eps) = -B S - d s_1 + e s_3 - f s_2
S = sum chi_4(n)   s_1 = sum chi_4(n)/n   s_2 = sum chi_4(n)/n^2   s_3 = sum 1/n^2
```

all four sums over those layers, `s_3` the trivial character's tail and the other three `chi_4`'s. `chi_4` alternates on the odd numbers, so `S = 0` at even `M` and `1` at odd `M`, and the other three walk to their L-values, `L(1, chi_4) = pi/4`, `sum 1/n^2 = pi^2/8` over odd `n`, and `L(2, chi_4) = G`:

```
M (mean - A - c eps) -> -B [M odd] - d pi/4 + e pi^2/8 - f G
```

| family | `A` | `B` | `c` | `d` | `e` | `f` | limit at even `M` | at odd `M` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| carpet | `1/2` | `1/8` | `1/2` | `0` | `0` | `-1/8` | `G/8` | `G/8 - 1/8` |
| net | `1/2` | `-1/8` | `-1/2` | `0` | `0` | `1/8` | `-G/8` | `1/8 - G/8` |
| tree | `1/4` | `0` | `1/3` | `-1/12` | `1/6` | `-1/6` | `pi/48 + pi^2/48 + G/6` | the same |
| void | `1/4` | `0` | `0` | `-1/4` | `1/2` | `0` | `(pi + pi^2)/16` | the same |

The leading term is the same phase-map integral the pair limits run, taken at the identity map `a = 1`, `c = 0`: over the hexagon polygon the cut cell law `s_y = s_x + s_z + w` returns `A + B chi` at each residue class of `n mod 4`, the carpet's `3/8` and `5/8` and the tree's and the void's `1/4`, which the generator prints as the `E` fields of every pair limit. So `A` and `B` are the integral's, `B chi` being the half-cell overtone `w` and nothing else, and `c`, `d`, `e` and `f` are the finite-`n` orders no phase integral sees - which is where `pi` and `G` come from and why they cannot reach a pair limit, every one of which is rational.

The hypothesis carries two limbs: the object is the cut ink of a 3D parity design, whose Walsh quasipolynomial carries only `chi_4` and terminates at the `1/n^2` order, proved in [the Walsh spectrometer lane](../claims/diagonal-slice-stack.md), and the quantity is the recentred `M (mean - A - c eps)`. Inside it, `pi` reaches a one-layer constant only through the `1/n` order of the ink law, `G` only through the `1/n^2` order, the residue class only through `B`, and no other character and no other constant can appear - which is why `L(2, chi_-3)` is absent from the table by a theorem and not by a search. Both limbs carry weight: the ghost star's width family above is a one-layer object of this same page whose constant `C_W = m (ln 2 + gamma/2) + q L(1, chi_8) - G/8 + Delta_W` holds `gamma`, `ln 2` and `L(1, chi_8)`, because its character is mod 8 and its `1/n` limb is not subtracted. Since `N = 2M - 1`, even `M` is `N = 3 mod 4` and odd `M` is `N = 1 mod 4`, and the modulus is 4 and no finer: the generator holds the summed identity against the counted hexagons in exact rational arithmetic at every layer count to `N = 55`, all four families, the four classes `n = 1, 3, 5, 7 mod 8` counted apart, 7 of 7 in each. The one row of the table that is not an instance is the carpet split, `M (I1 - I3 + 1/4) = s_1 + s_3/4 -> pi/4 + pi^2/32`, and it needs even `M`: at odd `M` the two classes hold `(M + 1)/2` and `(M - 1)/2` layers and the difference of their means is not a character sum at all.

**The Catalan statement, at both residue classes. Proved.** The carpet is `B = 1/8`, `d = e = 0`, `f = -1/8`, so `M (mean ink - 1/2 - eps/2) -> G/8 = 0.1144956993` along `N = 3 mod 4`, measured `0.1144757884` at `N = 55`, and `-> G/8 - 1/8 = -0.0105043007` along `N = 1 mod 4`, measured `-0.0104828892` at `N = 53`, with `G = 0.9159655942` from its own series. The `1/8` step is the ink law's own `chi` averaged over an odd number of layers, the same parity term the ghost star's even-`L` hypothesis carries, and Catalan enters only as `L(2, chi_4)`, one order below the `pi` the tree and the void collect. Nothing in it is fitted.

**The approach to every one of those limits is a closed form. Proved.** Write `sigma = +1` at even `M` and `-1` at odd `M`, and `a = 2M + 1` for the first layer past the ladder. Each tail solves its own one-line recursion, `T(a) + T(a + 2) = a^-s` twisted and `T(a) - T(a + 2) = a^-s` untwisted, in powers of `1/a`:

```
sum chi_4(n)/n   over odd n >= a = sigma (1/(2a) + 1/(2a^2)) + O(a^-4)    = sigma/(4M) + O(M^-3)
sum chi_4(n)/n^2 over odd n >= a = sigma (1/(2a^2) + 1/a^3) + O(a^-5)     = sigma/(8M^2) + O(M^-4)
sum 1/n^2        over odd n >= a = 1/(2a) + 1/(2a^2) + 1/(3a^3) + O(a^-5) = 1/(4M) + O(M^-3)
```

the `1/M^2` limb of each cancelling on the way from `a` to `M`, so the gap to the limit is `(sigma d - e)/(4M) + sigma f/(8 M^2) + O(1/M^3)`:

| quantity | derived at even `M` | at odd `M` | at `M = 3200` | at `M = 3201` |
|---|---:|---:|---:|---:|
| carpet `gap * M^2` | `-1/64` | `+1/64` | -0.01562500 | +0.01562500 |
| net `gap * M^2` | `+1/64` | `-1/64` | +0.01562500 | -0.01562500 |
| tree `gap * M` | `-1/16 - 1/(48M)` | `-1/48 + 1/(48M)` | -0.06250651 | -0.02082683 |
| void `gap * M` | `-3/16` | `-1/16` | -0.18750000 | -0.06250000 |
| carpet split `gap * M` | `-5/16` | | -0.31249999 | |

The carpet's `1/M` term is absent because `d = e = 0`, which is why its ladder is the sharpest of the five; the tree's second limb is the only place `sigma f/(8 M^2)` is separately visible, and it carries the last three digits of both tree readings. The generator prints all four classes of `M mod 4`, so all four of `N mod 8`, at `M = 400`, `1600` and `3200`, and the two parities are the only split in them.

**Eisenstein is absent from the base-2 slice stack. Verified.** `L(2, chi_-3) = 0.7813024129`, the L-value [the bases page](../../bases/notes/bases.md) finds under base 3, appears nowhere in the list above; the two characters this slice generates are `chi_4`, from the ink law, and the mod-8 character of `Q(sqrt 2)`, from the star arm above, and the hexagonal geometry contributes rational tent integrals rather than an Eisenstein L-function. The hexagon in this page is a shape, not [an Eisenstein lattice](https://en.wikipedia.org/wiki/Eisenstein_integer).

## Where the exactness stops

**The flat analogue is exact and proved on the shelf. Proved.** For the square parity carpets at odd scales, [moire-correlation-laws](../claims/flat-carpet-stack.md) proves that the Pearson correlation of two layers is a closed form in `gcd(m, n)` which vanishes if and only if the scales are coprime, for each of the four fields the plane rule generates. Why none of that reaches this page: the plane rule "ink unless both indices are odd" expands as `3/4 + (sigma_1 + sigma_2)/4 - sigma_1 sigma_2 /4`, and the pair term is the entire mechanism, while the space rule "fill if at most one of three indices is odd" expands as `1/2 + (sigma_1 + sigma_2 + sigma_3)/4 - sigma_1 sigma_2 sigma_3/4`, whose three pairwise coefficients are all exactly zero. A triple product takes the pair term's place, and a planar section of the space rule inherits nothing.

What replaces the exactness is the dyadic overtone of the second section, and the numbers depend on the mask. The paper's own hexagonal remark reports measurements on rendered sections, mean `-0.037` over a range `[-0.205, +0.147]`, and labels them measurements on rendered images rather than theorems, the extreme values mask-dependent. **Conjecture.** The exact full-hexagon computation on this page reads `(5, 9) = -0.14179450` and `(5, 7) = -0.08542646` for the same two pairs, so any number published for a hexagon layer pair must pin the mask convention it was measured under. What is robust across masks is the shape of the failure: coprime pairs correlate, the gcd echo survives, and the sign is the dyadic law's.

## Where the numbers live

`lab/rs/hexagon-moire` prints every number above in one run: the exact cut ink of every odd `n <= 55` against the four closed forms in rational arithmetic, the exact full-hexagon Pearson correlations and the doubling sign law, the pair limit as an exact rational for all four families and all three phase maps with the finite-layer gap beside it, the ideal-frame stack at 3601 samples per axis for the quarter line and the void star, the twisted averages and the coarse crosshair model at 200003 samples, the rendered ghost star on its 1200 by 2399 raster, its decay in the ideal, lattice and cell frames with the arm ink law in exact rationals to `n = 2001`, the cell frame's ladders at `L = 0 mod 4`, `L = 2 mod 4` and odd `L`, the swept band widths with the odd-`W` rows, the recompute to `L = 3200` and the four held-out widths, the width family's per-layer identity in exact rationals at 23 half-widths broken out by `n mod 8`, its seven summed ladders at `L = 0 mod 4`, `L = 2 mod 4` and odd `L` with the sliding window in both even classes, the constants from the `mrly_num` series with the partial character sums at `N = 53` and `N = 55`, and the one-layer law: the six coefficients of each family read off its ink law, the summed identity against the counted hexagons in exact rationals at every layer count to `N = 55` with `n mod 8` broken out, and the approach ladders at `M = 400`, `1600` and `3200` in all four classes of `M mod 4`. The flat stack's ray lines, gcd sums and constants and the xor pair's eternal points print under their own headers. Every domain is the source's own. The mesh those layers are counted on, and the fills that live on it, are [slices](../../slices/notes/slices.md); the ink law's quasipolynomial is proved in [the Walsh spectrometer lane](../claims/diagonal-slice-stack.md); the flat stack this one is measured against is the [moire-correlation-laws](../claims/flat-carpet-stack.md) lane; and the L-value it does not contain is [bases](../../bases/notes/bases.md).
