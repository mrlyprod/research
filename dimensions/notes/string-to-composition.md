---
title: The string of a design and its arithmetic pole
lead: Computes the complex dimensions of the 1D designs, watches them surface as an oscillation in the box count, certifies the arithmetic pole, and shows composition multiplies the base.
figure: research-dimensions
slug: dimensions
---

## The string of a design

A 1D design is a subset `F` of its base's digits, the filled ones; it draws the set of `x` in `[0,1]` whose digits all lie in `F`. That set is the attractor of `fill = |F|` maps `x -> (x + f)/base`, every one with the same contraction ratio `1/base`, and its dimension is `log(fill)/log(base)`. The core page's [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) rule at base 3 fills digits `{0,2}` - the middle-thirds Cantor set, `log(2)/log(3) = 0.630930`.

The complement of the set is a multiset of gaps, and the gap multiset of a one-base design satisfies `G = G_1 + fill` copies of `G/base`. So its geometric zeta function - the sum of `g^s` over all gaps - has the closed form

```
zeta(s) = Z(s) / (1 - fill*base^(-s)), Z(s) = sum of g^s over level-1 gaps,
```

and its poles sit where the complex Moran equation `fill*base^(-s) = 1` holds:

```
s = dimension + 2*pi*i*m/ln(base), m in Z.
```

One vertical line of poles, equally spaced at `omega = 2*pi/ln(base)`. Sets whose complex dimensions line up on such an arithmetic progression are called *lattice*; sets with incommensurable ratios, whose poles spread out, are *nonlattice*. Every one-base design is lattice for the trivial reason that all its ratios are equal.

**Verified** (`lab/py/complex-dimensions`). The progression is derived in closed form above and then evaluated: for the Cantor design and three others (table below), all 81 predicted poles at `m = -40..40` kill the denominator to `5e-14`, and the numerator `Z(s)` stays bounded away from zero at every one of them - minimum `|Z(s)| = 0.500000` for the Cantor design - so no zero is cancelled and every one is a genuine pole. The identity `Z(1) = 1 - fill/base`, the statement that the gaps of a measure-zero set fill the whole interval, holds for all four.

One caveat: the pole set is the zero set of the denominator only where the numerator does not vanish, and a design with no gaps (`fill = base`), with one filled digit (`fill = 1`, a point), or with none has an empty gap multiset, an identically zero zeta function, and no complex dimensions at all.

| object | `dimension` | `omega = 2*pi/ln(base)` |
|---|---|---|
| base 3, digits `{0,2}` | 0.630930 | 5.719202 |
| base 5, digits `{0,2,4}` | 0.682606 | 3.903963 |
| base 15, digits `{0,4,10,14}` | 0.511916 | 2.320188 |
| base 15, digits `{0,2,4,10,12,14}` | 0.661642 | 2.320188 |

## The oscillation in the box count

The imaginary parts are not bookkeeping; they are visible. Let `N(eps)` be the number of `eps`-cells the set meets and detrend it: `g(u) = ln N(exp(-u)) - dimension*u`. A lattice set's `g` oscillates at angular frequency `omega = 2*pi/ln(base)` forever. About a filled corner the oscillation is an identity before it is a measurement: [the corner lemma](../../weights/notes/weights.md) gives `mu(B(r)) = w_0 mu(B(base r))`, so `ln mu(B(r)) - dimension ln r` is exactly `ln base`-periodic at every one-base design.

**Verified** (`lab/py/complex-dimensions`). For all four objects the periodogram peak of `g(u)` lands within 1% of the predicted `2*pi/ln(base)`; for the Cantor design the Blackman periodogram on the direct DFT grid reads 5.7024 against the predicted 5.719202. Four-decimal peak values sit inside one frequency bin of the estimator's resolution and are properties of the estimator, not of the object; what is stable, and what ships, is the sub-1% agreement.

A second reading folds `u` modulo each candidate period and asks how much of the variance of `g` the folded profile explains. At 40 bins on a fixed window:

| object | `ln(3)` | `ln(5)` | `ln(15)` |
|---|---|---|---|
| base 3, `{0,2}` | **0.192** | 0.030 | 0.016 |
| base 5, `{0,2,4}` | 0.018 | **0.509** | 0.034 |
| base 15, `{0,4,10,14}` | 0.007 | 0.022 | **0.667** |
| base 15, `{0,2,4,10,12,14}` | 0.009 | 0.030 | **0.534** |
| two-ratio control | 0.082 | 0.047 | 0.038 |
| aperiodic control | 0.007 | 0.012 | 0.115 |

**Verified** (`lab/py/complex-dimensions`). Each lattice object folds best at its own `ln(base)`, and the two controls fold well at none. But the two lattice signatures are not parallel in strength: the base-15 design's folding explains 66.7% of the variance at its period, the Cantor design's only 19.2% at its own. The Cantor figure is genuinely modest - the box count is a step function and grid alignment injects a large aperiodic component - and it moves between roughly 0.16 and 0.23 as the bin count and window vary (Conjecture; the study prints the 40-bin value only), so it is also partly a property of the estimator. The ordering is stable; the percentages are not constants of the objects.

## The arithmetic pole, certified

The string above lives in `[0,1]`; its arithmetic twin is the set `S` of positive integers whose base-3 digits all lie in `{0,1}` (`1, 3, 4, 9, 10, 12, 13, ...`, the [mobius page](../../mobius/notes/mobius.md)'s `S_F` at base 3), with counting function `A(x)` and Dirichlet series `K(s) = sum n^(-s)`, abscissa `s_0 = log_3 2 = 0.630930`. [Burnol 2026](https://arxiv.org/abs/2602.19727) continues `K` meromorphically to `C` with simple poles among the same lattice `s_(m,k) = s_0 - m + 2 pi i k/log 3` (Proposition 4.1), proves the real pole `s_0` genuine with positive residue, gives every off-real residue as a limit of level sums, `lambda_(0,k) = (log 3)^(-1) lim_j sum_(3^j <= n < 3^(j+1), n in S) n^(-s_(0,k))` (Proposition 5.1), and shows that a vanishing `lambda_(0,k)` empties the whole column `s_(m,k)`, `m >= 1` (Proposition 7.1); the paper states that it does not study the off-real residues further and prints none. Whether a single off-real pole is genuine is therefore open at source, and the [mobius page](../../mobius/notes/mobius.md) accordingly says `among`. This section closes it for `k = 1..10` by a certificate.

**Proved.** Write `G(u) = A(3^u)/2^u`. Splitting an element `n = 3n' + a`, `a in {0,1}`, gives `A(3x) = 2A(x) + 1 - delta(x)` for `x >= 1`, where `delta(x) = A(x) - A(x - 1/3)` is `0` or `1`; so `G(u+1) - G(u) = (1 - delta(3^u))/2^(u+1)` lies in `[0, 2^(-u-1)]`, `G(u+j)` increases in `j`, and `Phi(u) = lim_j G(u+j)` exists, is 1-periodic, satisfies `0 <= Phi(u) - G(u) <= 2^(-u)`, and is continuous: at any point the oscillation of `Phi` is at most the jump of `G(u+j)` there plus twice the uniform distance, `3 * 2^(-u-j)` for every `j`. Hence `A(x) = x^(s_0) Phi(log_3 x) + E(x)` with `|E(x)| <= 1`. The profile is explicit on a third of its period: no element of `S` lies strictly between `11...1 = (3^(level+1) - 1)/2` and `100...0 = 3^(level+1)`, so `A = 2^(level+1) - 1` on `[(3^(level+1) - 1)/2, 3^(level+1))` and `Phi(u) = 2^(1-u)` for `u in [1 - s_0, 1]`, from `Phi(1 - s_0) = 2^(s_0) = 1.548562` down to `Phi(1) = Phi(0) = 1`. The profile is not constant, and that much is elementary.

**Proved.** For `Re s > s_0`, `K(s) = s int_1^inf A(x) x^(-s-1) dx`. Substituting the profile and `x = 3^u`, the `Phi` part is `s log 3 F(s)/(1 - 3^(s_0 - s))` with `F(s) = int_0^1 3^(u(s_0 - s)) Phi(u) du` entire, and the `E` part is holomorphic on `Re s > 0`. So on `Re s > 0` the poles of `K` are among the zeros `s_(0,k)` of `1 - 3^(s_0 - s)`, all simple with derivative `log 3` there, and `Res_(s_(0,k)) K = s_(0,k) c_k` with `c_k = int_0^1 Phi(u) e^(-2 pi i k u) du`. The pole at `s_(0,k)` is genuine exactly when the `k`-th Fourier coefficient of the profile is nonzero; non-constancy says only that some `c_k` is nonzero and nothing about `k = 1`.

**Proved** (computer-assisted, `lab/py/burnol-residue`). Let `Lev_j(w) = sum n^(-w)` over the `2^j` elements of `S` in `[3^j, 3^(j+1))`; Burnol's Proposition 5.1 reads `lambda_(0,k) log 3 = lim_j Lev_j(s_(0,k))`. Three elementary lemmas make the limit computable with a bound. (1) Splitting `n = 3n' + a` and expanding `(1 + a/(3n'))^(-w) = sum_l (-1)^l ((w)_l/l!) (a/(3n'))^l`, absolutely convergent since `a/(3n') <= 1/3`, then summing over the finitely many `n'` of level `j`: `Lev_(j+1)(w) = 3^(-w) sum_(l >= 0) (-1)^l ((w)_l/l!) 3^(-l) gamma_l Lev_j(w + l)` with `gamma_0 = 2` and `gamma_l = 1` for `l >= 1`. (2) `|Lev_j(w)| <= 2^j 3^(-j Re w)` for `Re w >= 0`: `2^j` terms, each of modulus `n^(-Re w) <= 3^(-j Re w)`. (3) At `w = s_(0,k) + i`, cutting the `l`-sum at `T` costs at most `2 * 3^(-ji) sum_(l > T) ((|w|)_l/l!) 3^(-l(j+1))`: termwise `|(w)_l| <= (|w|)_l`, `gamma_l <= 2` (slack by a factor 2, kept for the general engine), and (2) at `w + l` with `2^j 3^(-j s_0) = 1`; the full sum is the binomial series `(1 - 3^(-j-1))^(-|w|)`, increasing in `|w|` so an upper bound of `|w|` may stand in, and the cut is that closed form minus the interval partial sum; the same bound at `i = 0`, cut by `(1 - x)^(-a) - 1 <= a x (1 - x)^(-a-1)` for `0 <= x < 1` and `a >= 0` and summed over `j >= level`, gives the level tail `E_level = (3/2) |s| (1 - 3^(-level-1))^(-|s|-1) 3^(-level-1)`. The generator runs (1) on the vector of shifts `i = 0..I` from `Lev_0 = 1` to level `level = 40` in `mpmath.iv` at 128 bits, `3^(-s_(0,k)) = 1/2` being exact so that `s_(0,k)` and its Pochhammer symbols are the only transcendental inputs, adds the box of (3) at every step and `E_level` at the end, and prints endpoints floored and ceiled from exact fractions. At `k = 1`, `I = 66`: `lambda_(0,1)` lies in `[0.231891517689918, 0.231891517689919] + i [-0.501067414481069, -0.501067414481068]`, width `3.4e-18`, level tail `2.4e-19`, at distance at least `0.552125193` from zero. So the pole of `K` at `s_(0,1) = 0.630930 + 5.719202 i` is genuine, `c_1` is nonzero, and the counting function of `S` carries the frequency `2 pi/log 3` at the fundamental with a certified nonzero amplitude. The claim rests on Burnol's Proposition 5.1, the three lemmas above, and the outward rounding of mpmath's interval arithmetic; every bound used is proved.

**Proved** (computer-assisted, `lab/py/burnol-residue --band`). The same certificate at `k = 0..10`, `level = 40`, `I` from 52 to 164, widths from `2e-18` to `4e-13`, the table printed by the generator:

| `k` | `Re lambda_(0,k)` | `Im lambda_(0,k)` | `abs(lambda_(0,k)) >=` | zero excluded |
|---|---|---|---|---|
| 0 | `[0.799023642655, 0.799023642656]` | `[-0.000000000001, 0.000000000001]` | `0.799023642` | yes |
| 1 | `[0.231891517689, 0.231891517690]` | `[-0.501067414482, -0.501067414481]` | `0.552125193` | yes |
| 2 | `[0.003963750587, 0.003963750588]` | `[-0.033400689772, -0.033400689771]` | `0.033635062` | yes |
| 3 | `[0.329059109066, 0.329059109067]` | `[-0.160535971769, -0.160535971768]` | `0.366130708` | yes |
| 4 | `[-0.039525178566, -0.039525178565]` | `[-0.366998806233, -0.366998806232]` | `0.369121068` | yes |
| 5 | `[-0.149130772496, -0.149130772495]` | `[0.032063501333, 0.032063501334]` | `0.152538701` | yes |
| 6 | `[0.116528526486, 0.116528526487]` | `[0.073732764792, 0.073732764793]` | `0.137896403` | yes |
| 7 | `[0.155123454907, 0.155123454908]` | `[0.001729927879, 0.001729927880]` | `0.155133100` | yes |
| 8 | `[0.252556366907, 0.252556366908]` | `[-0.100290430690, -0.100290430689]` | `0.271740480` | yes |
| 9 | `[0.079103971563, 0.079103971564]` | `[-0.275729544156, -0.275729544155]` | `0.286852261` | yes |
| 10 | `[-0.032664131229, -0.032664131228]` | `[-0.078686575787, -0.078686575785]` | `0.085196964` | yes |

The `k = 0` row is Burnol's positive real residue recovered, its imaginary interval the enclosure's own outward rounding and not a reading, since `s_(0,0)` is real and the residue's imaginary part is exactly `0`. By Proposition 7.3 of the same paper, under its hypothesis `1 < N < base` (here `N = 2` at base 3, so the Pochhammer prefactors relating his `mu_(m,k)` to `lambda_(m,k)` never vanish), a nonzero `lambda_(0,k)` makes `s_(m,k)` a pole exactly when `s_(m,0)` is one, so each certified `k` settles its whole column against the real axis.

**Verified** (`lab/py/burnol-residue`). Two independent floating-point computations land on the certificate. Direct enumeration of `Lev_20(s_(0,1))` over its `2^20` terms gives `0.231891518 - 0.501067414 i` after division by `log 3`, at distance `1.9e-10` from the enclosure against its own bound `7.5e-10`, the successive level differences `1.13e-08, 3.77e-09, 1.26e-09, 4.19e-10` contracting by `1/3`. The counting side, from exact counts `A(x)` on `2^17` midpoints of one period at `J = 30`, gives `c_1 = -0.082138842 - 0.049607499 i` and `s_(0,1) c_1 = 0.231891457 - 0.501067454 i`, at distance `7.2e-08` (`1.0e-06` on `2^15` midpoints: the quadrature error of a rough profile, not a certified bound); the sampled profile runs over `[1.000003, 1.548557]` and matches `2^(1-u)` on `[1 - s_0, 1]` to `7.2e-10`. A second certified enclosure by the functional-equation route, `R_1 = 1 + sum_(m >= 1) (-1)^m ((s)_m/m!) 2^(-1) 3^(-m) K(s + m)`, each `K(s+m)` summed in intervals to level 13 plus its tail `3^(-13m)/(1 - 3^(-m))` and the `m`-series cut at 36 plus its tail, gives `[0.23189098, 0.23189280] + i [-0.50106814, -0.50106632]`, width `1.8e-06`, meeting the first. The engine is tested on two designs with known answers. On `{0,2}`, whose elements are twice those of `{0,1}`, the enclosure `[0.135292875345483, 0.135292875345484] + i [0.329874091244466, 0.329874091244467]` agrees with `2^(-s_(0,1))` times the `{0,1}` enclosure to every printed digit. On the full digit set `{0,1,2}`, where `K` is [the Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function), `k = 0` returns `[0.999999999999999, 1.000000000000001]` around the residue `1` and `k = 1` returns a box of width `1.7e-18` around `0`, the regular point it must be.

**Verified** (`lab/py/burnol-residue`). Burnol's Proposition 7.1 recurrence carried down the `k = 1` column from the certified `lambda_(0,1)`: `lambda_(1,1)` in `[0.6950303416, 0.6950303417] + i [0.3777908610, 0.3777908611]`, `lambda_(2,1)` in `[-0.1945129694, -0.1945129693] + i [0.2161122817, 0.2161122818]`, `lambda_(3,1)` in `[0.0645979040, 0.0645979041] + i [0.1353703406, 0.1353703407]`, all nonzero. The first two meet the closed forms Theorem 7.4 forces, `lambda_(m,k) = (-1)^m (s - 1)...(s - m) lambda_(0,k) [t^m] (1/E(t))` at `s = s_(0,k)`, `E` the moment generating function of the Cantor measure of `{0,1}` with moments `1/4` and `3/32`, so `[t^1] = -1/4` and `[t^2] = 1/64`: `lambda_(1,1) = (s - 1) lambda_(0,1)/4` and `lambda_(2,1) = (s - 1)(s - 2) lambda_(0,1)/64`.

What is Burnol's and what is new. The continuation, the lattice, the positivity at `s_0`, the limit formula and the column proportionality are Burnol 2026, cited; that `A(x) x^(-s_0)` oscillates is elementary (`A(3^level) = 2^level` against `A((3^(level+1) - 1)/2) = 2^(level+1) - 1`, the profile lemma) and nothing here discovers it. New is the certificate at `k = 1..10`, the `k = 0` row recovering Burnol's theorem: the pole at `s_(0,1)`, and at every `s_(0,k)` for `k = 1..10`, is genuine, the residue enclosed to `4e-18` at `k = 1` and `4e-13` at `k = 10`, which is what a complex-dimension claim on the arithmetic side needs. The engine takes any base and digit set `F` containing `0`; only base 3 is certified here, and the zeta check and the `{0,2}` scaling are its tests. No statement is made about `k > 10`, about other bases, or about the real-axis residues `lambda_(m,0)` beyond `m = 0`.

## Composition multiplies the base

Compose two rules by alternating bases across levels: base 3 with digits `{0,2}` at odd levels, base 5 with digits `{0,4}` at even ones.

**Verified** (`lab/py/complex-dimensions`). The alternation produces exactly the one-base design at base 15 with digits `{5*d1 + d2} = {0,4,10,14}` - checked as integer arithmetic and then as geometry, eight alternating levels and four base-15 levels producing the identical 256 intervals as exact fractions. So the composite is a different lattice period, `omega = 2*pi/ln(15) = 2.320188`, not a departure from the lattice class - and the table above shows its signature is the *sharpest* of the family, 66.7% of variance at its own period.

Two qualifications, both load-bearing. First, base 5 with digits `{0,4}` is not a mrly design: the core page's move one fills by parity, and the even digits of base 5 are `{0,2,4}`, three of them - nor is the base-15 composite a parity design. The parity-faithful versions - base 5 `{0,2,4}` and their composite base 15 `{0,2,4,10,12,14}` - run alongside and behave identically, which is the real point: the argument turns on the base, not on which digits the rule picks.

Second, the structural claim that one contraction ratio per level implies lattice is false without a periodicity hypothesis. (Refuted.) An aperiodic control that alternates bases 3 and 5 on a [Thue-Morse](https://en.wikipedia.org/wiki/Thue%E2%80%93Morse_sequence) schedule uses exactly one ratio per level and is not self-similar at all, so it is neither lattice nor nonlattice: its best folding is 0.115 against the composite's 0.667, and its periodogram peak matches no `2*pi/ln(base)`. What survives, and needs no computation, is the statement for mrly designs proper: move two is a [Kronecker](https://en.wikipedia.org/wiki/Kronecker_product) power of one tile, so the schedule is constant, every level subdivides by the same base, and periodic cross-base alternation multiplies into one product base. Neither move can express two ratios inside one level - move one only chooses which cells of a fixed `base^dim` grid survive, and every cell of that grid is the same size. (Verified for the alternation; the one-tile argument is read off the definition.)

"Read off the definition" has a proof behind it, and it covers more than one tile. **Block reduction, Proved:** for any designs `c_1, ..., c_p` and any `level >= 1`, the periodic word `(c_1, ..., c_p)^level` equals `(A_(c_1) (x) ... (x) A_(c_p))^((x) level)`, by associativity of the Kronecker product and nothing else. (The case `p = 2` is the base-15 alternation Verified above; a check on six test cases at periods 2 and 3, lengths to 6, matching the flat word cell for cell against the self-Kronecker power of the composite, is **Verified**, `connectivity/lab/rs/magic-words`.) So every periodic schedule *is* the ordinary self-similar theory of its one-period composite tile, of base `prod_i base_i` and fill `prod_i fill_i`; the base-15 composite above is the case `p = 2`, and its agreement was never in doubt. The corollary is the sharper half. The first genuinely non-stationary behaviour requires an aperiodic word, which is exactly why the Thue-Morse control above is neither lattice nor nonlattice rather than being a third kind of composition. For an aperiodic word at one common base with `fill_i = fill(c_i)`, the scale dimension is `lim_level (sum_(i <= level) log fill_i) / (level log base)` when the limit exists; characterizing which words make it exist, and which leave the dimension fluctuating, is open. (Conjecture, untouched.)

The genuine way out is a two-ratio system: maps of ratio `1/3` and `1/5` mixed *within* one level, outside the mrly family, with dimension `0.518370` solving `3^(-s) + 5^(-s) = 1`. **Verified** (`lab/py/complex-dimensions`). Its 21 complex dimensions in the box `Re` in `[-3,3]`, `Im` in `[-40,40]` - a complete list, by the argument principle: the winding number over the box is 21 and 21 roots are found - have real parts spread from `-0.699926` to `0.518370` and fit no arithmetic progression, the worst offset being 0.43, 0.17 and 0.38 of a step for the three candidate spacings. Nonlattice is a real, different behaviour, and no mrly design or composition exhibits it.
