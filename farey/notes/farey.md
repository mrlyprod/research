---
title: The stack is an RH-observable
lead: The stack's moire is a Farey resonance diagram; each scale `n` adds exactly `phi(n)` bright nodes, the stack is an address rather than a construction (any depth evaluates in closed form, which buys rendering and provably nothing toward RH), and restricting the Farey sequence to a digit design empties a fixed sixth of the line at base 3 `{0,1}` when both coordinates are restricted, while restricting the denominator alone leaves the Franel-Landau shape reading the same.
figure: research-farey
slug: farey
---

Lay the same fractal grid on the unit square at many scales at once - scale `n` puts its cell boundaries at `x = k/n` - drop the opacity and add the layers up. The result is a moire, and a bright point is one that many scales agree on. The question this page answers is what the bright points are, and the answer is not decorative: the lit nodes are the Farey fractions, the amount of new structure each scale contributes is [Euler's totient](https://en.wikipedia.org/wiki/Euler's_totient_function) `phi(n)`, and how evenly those nodes spread is - by a pair of theorems from 1924 - literally equivalent to the Riemann hypothesis.

The [Farey sequence page](https://en.wikipedia.org/wiki/Farey_sequence) builds the stack scale by scale, lights the Farey fractions, and shows `phi(n)` novelty peaking at [the primes](https://en.wikipedia.org/wiki/Prime_number).

## Where the lines land

Stack the scales `n = 1..N`. A point `a/b` in lowest terms receives a grid line from exactly the scales that are multiples of `b`, so over `1..N` its brightness is `floor(N/b)`. **Proved**, and **Verified** by direct simulation at `N = 30`: building the stack node by node and comparing every node's hit count against `floor(30/b)` gives no mismatch anywhere (`lab/rs/farey-discrepancy`).

Brightness therefore falls as one over the denominator, which is the Stern-Brocot ordering of the rationals. The top of the table at `N = 30`:

| node | brightness | `floor(30/b)` |
|---|---|---|
| `0`, `1` | 30 | 30 |
| `1/2` | 15 | 15 |
| `1/3`, `2/3` | 10 | 10 |
| `1/4`, `3/4` | 7 | 7 |
| `1/5` ... | 6 | 6 |

The lit nodes are also exactly the lattice points visible from the origin, since `a/b` is in lowest terms precisely when `gcd(a,b) = 1`. That is the "lighthouse" reading of the picture. **Proved.** The density of visible points is `6/pi^2` - the same constant, and the same base-blindness, discussed in [what base 3 hides](../../bases/notes/bases.md), where it is measured as `0.608042` on a `3000 x 3000` grid. **Verified**, by recounting that grid.

Those same points draw two further pictures of the one node set. The sunburst is the visible lattice points themselves, the `(a, b)` with `gcd(a, b) = 1` and `abs a, abs b <= n`, joined in angular order into a closed polygon, which reads the lit set as directions rather than as points of the line. The resonance diagram is the spike profile, a spike at the node `a/b` of height `1/b`. **The stack is that diagram up to the floor. Proved**, from the brightness law: normalising gives the node `a/b` the height `floor(Q/b)/Q`, which lies in `(1/b - 1/Q, 1/b]` at every depth and equals `1/b` exactly when `b` divides `Q` (`lab/rs/farey-discrepancy`).

## Primes are the maximally novel scales

The nodes scale `n` introduces *for the first time* are the fractions `a/n` with `gcd(a,n) = 1`, since any `a/n` that reduces was already lit by the smaller scale it reduces to. There are exactly `phi(n)` of them. **Proved**, and **Verified** by set difference over the stack for `n = 2..30`:

```
n    2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
new   1 2 2 4 2 6 4 6 4 10 4 12 6 8 8

n   17 18 19 20 21 22 23 24 25 26 27 28 29 30
new  16 6 18 8 12 10 22 8 20 12 18 12 28 8
```

Every count equals `phi(n)`, and the running maxima `1, 2, 4, 6, 10, 12, 16, 18, 22, 28` occur at `n = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29`. The reason is one line: `phi(n) = n - 1` if and only if `n` is prime, because every one of `1..n-1` is coprime to `n` exactly when `n` has no smaller factor. **Proved**.

So primality is readable off the picture. Stack `1..n-1`, then add scale `n`, and count what appeared: `n - 1` new nodes means `n` is prime, fewer means composite. **Proved** (it is the previous claim restated). A composite scale mostly re-lights nodes its own divisors already drew - scale 30 adds only 8 new lines, the rest of its grid falling on lines from 1, 2, 3, 5, 6, 10 and 15.

## Franel and Landau, 1924

Over scales `1..Q` the stack lights exactly the reduced fractions of denominator at most `Q`: the Farey sequence `F_Q`. Its size in `(0,1]` is `m = sum_{k<=Q} phi(k)`. **Proved**, and **Verified** by generating `F_Q` through the next-term recurrence and comparing its length with the totient sum at `Q = 10, 30, 60` (`lab/rs/farey-discrepancy`, which runs the same comparison at `Q = 10, 30, 60, 125`).

Write `rho_1 < ... < rho_m` for those nodes and `delta_j = rho_j - j/m` for how far each one sits from perfect equidistribution. Then:

- Franel (1924) proved that `sum_j delta_j^2 = O(Q^(-1+eps))` for every `eps > 0` is equivalent to the Riemann hypothesis.
- Landau (1924), in a note published immediately after Franel's, proved the same for `sum_j |delta_j| = O(Q^(1/2+eps))`.

**Verified** against the literature: both statements, with the original 1924 citations to the Göttingen Nachrichten, are the standard Franel-Landau formulation, and are reproduced in Edwards, *Riemann's Zeta Function*, chapter 12.

Put the two halves together. The nodes whose discrepancy Franel and Landau are talking about are the nodes the stack draws - not an analogue of them, the same set. So the question "how evenly are the bright points spread?" is not *related to* the Riemann hypothesis; at this level of precision it **is** the Riemann hypothesis. **Proved**, given the identification above, which is what the first two sections establish.

## The meter reads what RH predicts

Both sums are computable. Generating `F_Q` exactly and measuring, with `S2 = sum delta_j^2` and `S1 = sum |delta_j|`:

| `Q` | nodes | `S2*Q` | `S1/sqrt(Q)` | local exponent of `S2` |
|---|---|---|---|---|
| 125 | 4796 | 0.5395 | 0.2040 | - |
| 250 | 19024 | 0.5848 | 0.1942 | -0.884 |
| 500 | 76116 | 0.6241 | 0.1852 | -0.906 |
| 1000 | 304192 | 0.6387 | 0.1634 | -0.967 |
| 2000 | 1216588 | 0.6560 | 0.1512 | -0.961 |
| 4000 | 4863602 | 0.6538 | 0.1314 | -1.005 |
| 8000 | 19455782 | 0.6564 | 0.1123 | -0.994 |

**Verified** by `lab/rs/farey-discrepancy`. `S2*Q` flattens near `0.656` and the local exponent walks to `-1`, which is the Franel condition; `S1` stays under its `Q^(1/2)` envelope and its own local exponent runs between `0.27` and `0.43`, under the Landau threshold of `0.5`. The node count matches `sum phi(k)` exactly at every rung, which is the control that says the object being measured really is the stack's node set.

## The meter on a digit design

A digit design `S_F` is the set of whole numbers whose every digit lies in a digit set `F`, of dimension `alpha = log |F| / log base`. Restricting the meter to one needs a convention, and the strict one is taken here: `F_Q(S_F)` is the set of reduced fractions `a/b` with `0 < a <= b <= Q` and both `a` and `b` in `S_F`. The weaker denominator convention - `b` in `S_F`, `a` free - is measured beside it, and the unrestricted `F_Q` is the control. Neither is a set the stack produces: [Farey order is the stack](#farey-order-is-the-stack-not-the-design) still holds and every design gives the same stack at fixed `Q`, so `F_Q(S_F)` is a filter laid over the stack's nodes by hand, and all that is at stake is what the Franel-Landau functional reads on it. No [Mobius](https://en.wikipedia.org/wiki/M%C3%B6bius_function) sum appears anywhere in this section, so nothing in it bears on the square-root conjecture of [mobius](../../mobius/notes/mobius.md), whose `theta(F)` is a [Mertens](https://en.wikipedia.org/wiki/Mertens_function) exponent and belongs to that page alone; the exponents here are called `e_2` and `e_1` and are this section's own.

Each convention carries a count that never enumerates a fraction. The denominator convention has `card = sum_{b in S_F, b <= Q} phi(b)`. The strict one has `card = sum_{b in S_F, b <= Q} phi_F(b)` with `phi_F(b) = #{a in S_F : a <= b, gcd(a,b) = 1} = sum_{d | b} mu(d) * #{multiples of d in S_F up to b}`, inclusion-exclusion over the divisors of `b`. **Proved.** Both are sieved independently of the enumeration and match it at every rung of every table below, which is the control that says the object measured is the object defined; the largest check is 9538759028 nodes on the base 3 control at `Q = 3^11`. **Verified** by `lab/rs/farey-discrepancy design`.

Rungs are powers of the base, so the design's set is self-similar at every rung. Every exponent below is one ratio between consecutive rungs, `e_2 = ln(S2(Q')/S2(Q)) / ln(Q'/Q)` and `e_1` the same for `S1`, and nothing is fitted.

Write `D_Q = #{b in S_F : b <= Q} ~ Q^alpha` for the denominators the design supplies and `card F_Q(S_F) ~ Q^e` for the nodes. Square-root cancellation in the denominators is a node-count error of order `sqrt(D_Q)`, which puts `e_2` at `alpha - e`, and Cauchy-Schwarz on `S1 <= sqrt(card * S2)` then caps `e_1` at `alpha/2`. Franel and Landau are this pair at `alpha = 1`, `e = 2`: `e_2 <= -1` and `e_1 <= 1/2`, both under RH, both caps on a limsup and never values. The control misses them freely at a single rung - its `e_2` wanders from `-0.463` to `-1.017` and its `e_1` from `+0.236` to `+0.673` across the base 3 ladder - so a lane sitting off a cap at one rung shows nothing by itself. The ladders run to `Q = 3^11 = 177147` and `Q = 10^5`, each with its own control, and the top rung reads:

| set | convention | `alpha` | `e` | `exp card` | `alpha - e` | `e_2` | `alpha/2` | `e_1` |
|---|---|---|---|---|---|---|---|---|
| base 3 `{0,1}` | strict | 0.631 | 1.262 | +1.263 | -0.631 | +1.259 | 0.315 | +1.262 |
| base 3 `{0,1}` | denominator | 0.631 | 1.631 | +1.631 | -1.000 | -0.959 | 0.315 | +0.336 |
| full set to `3^11` | control | 1.000 | 2.000 | +2.000 | -1.000 | -1.017 | 0.500 | +0.236 |
| base 10 without 9 | strict | 0.954 | 1.908 | +1.908 | -0.954 | +1.904 | 0.477 | +1.906 |
| base 10 without 9 | denominator | 0.954 | 1.954 | +1.954 | -1.000 | -0.899 | 0.477 | +0.585 |
| full set to `10^5` | control | 1.000 | 2.000 | +2.000 | -1.000 | -1.000 | 0.500 | +0.347 |

The denominator lanes sit in the control's own band and the strict lanes sit on the mass. **Verified.** Both denominator rows put `e_2` at `-0.959` and `-0.899` against `alpha - e = -1.000` and `e_1` at `+0.336` and `+0.585` against caps of `0.315` and `0.477`, every one of the four inside the range the control walks. Both strict rows instead put `e_2` and `e_1` on `exp card` itself, agreeing with it to two decimals: `+1.259` and `+1.262` against `+1.263`, `+1.904` and `+1.906` against `+1.908`. Sums that ride the mass are sums with no cancellation at all, and the normalised readings settle:

| set | `Q` | `card` | `S1/card` | `S2/card` | widest gap |
|---|---|---|---|---|---|
| base 3 `{0,1}` | 2187 | 4286 | 9.4258e-2 | 1.3307e-2 | 0.16720 |
| base 3 `{0,1}` | 6561 | 17069 | 9.3888e-2 | 1.3100e-2 | 0.16684 |
| base 3 `{0,1}` | 19683 | 67561 | 9.4171e-2 | 1.3170e-2 | 0.16673 |
| base 3 `{0,1}` | 59049 | 269750 | 9.4171e-2 | 1.3147e-2 | 0.16669 |
| base 3 `{0,1}` | 177147 | 1080458 | 9.4057e-2 | 1.3086e-2 | 0.16667 |
| base 10 without 9 | 1000 | 147096 | 5.6111e-3 | 4.1767e-5 | 0.00113 |
| base 10 without 9 | 10000 | 11890654 | 5.2424e-3 | 3.6501e-5 | 0.00011 |
| base 10 without 9 | 100000 | 963170938 | 5.2125e-3 | 3.6159e-5 | 0.00001 |

The strict lane holds those constants from `Q = 2187` at base 3 `{0,1}`, `card` running 4286 to 1080458, and from `Q = 10000` at base 10 without 9, `card` running 11890654 to 963170938: two figures of `S1/card` and of `S2/card` fixed across each span, and no more of either claimed, the base 10 rung below moving the first figure of `S2/card` from `4.1e-5`. **Verified.**

At base 3 `{0,1}` the last column names an exact interval, a sixth of the line. **Proved.** If `b` has its leading base 3 digit at position `level` then `3^level <= b <= (3^(level+1) - 1)/2`. A numerator whose own leading digit sits at `level` gives `a/b >= 2 * 3^level/(3^(level+1) - 1) > 2/3`; one whose leading digit sits at `level - 1` or below gives `a <= (3^level - 1)/2 < b/2`, so `a/b < 1/2`. The closed interval `[1/2, 2/3]` is therefore empty in `F_Q(S_F)` at every `Q`, and a set that misses a fixed interval of positive length does not equidistribute. The widest gap the meter finds contains that interval at every finite `Q` and shrinks onto it from outside, reading `0.16827, 0.16720, 0.16684, 0.16673, 0.16669, 0.16667` at `Q = 3^6 .. 3^11` and starting at `0.49931, 0.49977, 0.49992, 0.49997, 0.49999, 0.50000`.

Base 10 without 9 has no gap to argue from. **Conjecture.** Its widest gap falls like `1/Q`, reading `0.01136, 0.00113, 0.00011, 0.00001` at `Q = 10^2 .. 10^5` against the control's `0.01000, 0.00100, 0.00010, 0.00001` on the same rungs, so the strict set there is as fine-grained as the full Farey sequence and the base 3 argument does not transfer. That it fails to equidistribute rests on the settled constants alone, which is a reading and not a proof.

The denominator convention keeps the transplanted shape. **Conjecture.** Restricting the denominator moves the mass from `Q^2` to `Q^(1+alpha)` and the count error from `sqrt(Q)` to `sqrt(D_Q) = Q^(alpha/2)`, and the two moves cancel in `S2`: `S2*Q` reads `0.8926` at base 3 `{0,1}` and `0.8536` at base 10 without 9 against the control's `0.6782` and `0.6684` at the same top rungs, a constant factor and not a rate. The transplanted Landau reading `S1/Q^(alpha/2)` is flat where the control's `S1/sqrt(Q)` is falling, at `0.243, 0.281, 0.267, 0.268, 0.274` for `Q = 3^7 .. 3^11` and `0.213, 0.222, 0.207, 0.265` for `Q = 10^2 .. 10^5`. So the conjecture is `S2 = O(Q^(-1+eps))` and `S1 = O(Q^(alpha/2+eps))` on the thinner set. Normalising against the node count instead would demand `e_2 = -(1+alpha)/2`, which is the wrong yardstick: it agrees with `alpha - e` on the control only through the accident `sqrt(Q) = card^(1/4)` at `e = 2`, and the denominator lane is not beating a shape by missing it.

The honest cap below covers this section too. Nothing here reaches further than the meter above it. What it adds is negative and clean: digit restriction of the denominator is invisible to the shape, and digit restriction of both coordinates destroys equidistribution outright at base 3 `{0,1}` and appears to at base 10 without 9. The denominator lane's invisibility is measured here and priced in [the restricted Franel identity](../../franel/notes/franel.md), where its exponential sum turns out to have a closed form in Mobius sums.

## The restricted Franel identity

The identity behind the denominator lane has a note of its own, [the restricted Franel identity](../../franel/notes/franel.md): the exponential sums as dilated Mertens sums, Franel's identity in its Fourier and rank forms, the ceiling it gives, the dilate automaton, the converse under (U') and under (SR), the mean-square equivalence, the surrogate `B(Q)`, the strict set's divisor form, and [the repunits](../../franel/notes/threshold-surrogate-repunits.md#the-repunits).

## Weighting the stack by Mobius

Give scale `n` the weight `mu(n)` instead of weight one and the same stack renders a different arithmetic function: the node `a/b` collects `mu` over the scales that are multiples of `b`, so its brightness is `Sum_{k <= N/b} mu(kb) = mu(b) * Sum_{k <= N/b, gcd(k,b) = 1} mu(k)`, a Mertens-type sum over the integers coprime to `b`. **Proved**, by the same divisor count that gives `floor(N/b)` in the unweighted stack. It is not `M(floor(N/b))`: the two agree at only 64 of 200 denominators at `N = 200` (`lab/py/mertens-meter`), and coincide at `b = 1`, where the node reads `M(N)` exactly.

That makes the picture a Mertens meter rather than a Farey one, and the oscillations of `M(x)/sqrt(x)` are where [the nontrivial zeta zeros](https://en.wikipedia.org/wiki/Riemann_zeta_function) live, by the explicit formula. Sampling `M(x)/sqrt(x)` in log-space and taking the power spectrum puts peaks at the first eight zeros:

| known `gamma` | detected | error |
|---|---|---|
| 14.1347 | 13.94 | 0.20 |
| 21.0220 | 20.90 | 0.12 |
| 25.0109 | 24.97 | 0.04 |
| 30.4249 | 30.19 | 0.23 |
| 32.9351 | 32.52 | 0.42 |
| 37.5862 | 37.74 | 0.16 |
| 40.9187 | 40.64 | 0.27 |
| 43.3271 | 42.97 | 0.36 |

**Verified** by `lab/py/mertens-meter`, and weaker than everything above it on this page: `M(x)/sqrt(x)` for `x = 1..50000` from a linear Mobius sieve, resampled uniformly in `log x` on 8192 points, Hann-windowed, the real FFT power spectrum read as `gamma = 2 pi f`, local maxima above three times the band median over `8 < gamma < 55`. The bin width is `0.5806`, so every error in the table sits inside one bin. The honest cap below covers this section too, and covers it harder: the zeros are known to far greater precision than a moire can reach, so what the picture buys is a rendering, not a measurement.

## The honest cap

An observable is not a handle. What [Franel and Landau, 1924](#franel-and-landau-1924) and [the meter](#the-meter-reads-what-rh-predicts) establish, and what [the mean-square equivalence](../../franel/notes/threshold-surrogate-repunits.md#the-threshold-is-a-mean-square) carries one digit design over, is that this picture renders a genuinely RH-equivalent object, which is a real upgrade over the vaguer "fractals and zeta both have self-similar structure" gestures. What it does not do is supply any route to a proof. The Riemann hypothesis is already checked numerically far beyond any range this or any other meter can reach, so the table above can only ever illustrate the expected behaviour - it is consistent with RH, it is not evidence for it, and no amount of extra `Q` changes that. Scored here, the link quality is 6 out of 10 and the meter's tractability 0, and both numbers deserve to be stated together: the connection is exact, and no renderer reaches it - an attack must come through the equidistribution toolkit, never through a picture.

**Two instructions this page hands the rest of the tree.** First: the toolkit flows both ways. The window at dimension one in [coprime](../../coprime/notes/coprime.md) is a discrepancy statement about a discrete arithmetic set and so is Franel-Landau, so the equidistribution methods that attack one are the methods the other needs - that kinship in technique is why an RH equivalence sits on a page of this tree, and it marks the one honest route: theorems, not renders. The verdict above is final for the meter alone; it caps what a picture can claim, never what a proof may attempt. Second: this page renders the Farey set without owning it. The rule whose ORBIT is the Farey set is the mediant, `(a/b, c/d) -> (a+c)/(b+d)`, with the Gauss map `x -> {1/x}` as its continued-fraction twin - simple local rules with emergent complexity, exactly this project's own principle, and they carry the Stern-Brocot and `GL_2(Z)` symmetry that base digit restriction does not. Mayer's theorem lives there: the Selberg zeta function of the modular surface is the Fredholm determinant of the Gauss-Kuzmin-Wirsing transfer operator, a genuine fractal-dynamics-to-zeta bridge. Two cautions travel with it - that is Selberg zeta and not Riemann zeta, and its RH-analogue is known for unrelated reasons; and the alphabet is infinite, so every finite-state tool on this tree needs rebuilding there.

## Farey order is the stack, not the design

- There is no design-specific Farey sequence, and there never was one to find. **Refuted.**
- The stack's lit set at maximum scale `Q` is exactly `{a/b : 1 <= a <= b <= Q, gcd(a,b) = 1}`, because a boundary coordinate `k/n` reduces to `a/b` and reappears at every scale divisible by `b`.
- Farey order is therefore `Q`, the maximum stacked grid scale. Fill count plays no part, and every design gives the same Farey sequence at fixed `Q`. **Proved** from the construction.
- Brightness `hits(a/b) = floor(Q/b)` is checked by literal stacking at `Q = 30` on all 278 lit fractions and up to `Q = 125` (`lab/rs/farey-discrepancy`). **Verified.**
- Under the transparent convention `Q = 3^level`, the geometric side length, the Landau discrepancy `D_Q = sum_i |f_i - i/m|` reads `0.166667, 0.549206, 1.150760, 2.118500, 3.187070` at `Q = 3, 9, 27, 81, 243`, with `m = 4, 28, 230, 2020, 18056`. Both generation routes agree exactly, a Farey next-term recurrence being the independent cross-check; the rows have no generator in `lab/`. **Conjecture.**
- `D_Q/sqrt(Q)` stays inside `[0.0962, 0.2354]` and reads `0.2045` at `Q = 243`. The adjacent log slope falls `1.085, 0.673, 0.556, 0.372`; the all-five log-log fit is `0.660` and the last-three fit `0.464`. Consistent with `O(Q^{1/2+eps})`, discriminating nothing: five nested deterministic points cannot test a statement quantified over every positive epsilon.
- The `Q = 3^level` map is a comparison convention chosen here, not a mapping the tree defines. Mapping `Q` to fill count would be arbitrary and was explicitly rejected.

## The stack is an address, not a construction

Can a stack be created immediately, without stacking? The answer is yes, exactly, and the boundaries of the yes are theorems of their own; every number in this section is regenerated by `lab/py/carpet-stack-address`.

Everything layer `n` does at a rational point `x = (a_1/q, a_2/q)` depends only on `r = n mod 2q`: `n` is odd iff `r` is odd, and `floor(n*a/q)` is odd iff `(r*a) mod 2q >= q`, since `n*a mod 2q = q*(floor(n*a/q) mod 2) + (n*a mod q)`. So the odd-carpet stack's brightness is a residue count with the `N`-dependence in closed form,

```
B_N(x) = ceil(N/2) - Sum_{r in S(x), r <= N} (floor((N-r)/2q) + 1)
```

with `S(x)` the bad residues, and the per-point cost depends on `q` alone, never on `N`. **Proved.** The line-stack's own form is the `floor(N/b)` at the top of this page, `O(1)` per node. **Verified** by two generators sharing no code in `lab/py/carpet-stack-address`, one stacking literally and one forbidden to loop over layers: identical Farey digests at `N = 55` (940 nodes, brightness sum `1540 = N(N+1)/2` landed by count), sha256-identical `512 x 512` renders by three routes, all 48 probes equal at `N = 55` and `5555`, the closed form against literal stacking at `N in {1, 2, 55, 5555, 19945, 19946, 19947, 40001}` with zero mismatches, and a stack of `5*10^17` layers - `N = 10^18` - evaluated in a tenth of a second by both implementations, exactly, values agreeing digit for digit.

The scope is part of the result, each boundary proved. Per-point only: an `R x R` raster costs `R^2` writes no matter what. Exact representations only: on a point supplied as a real oracle the value is undecidable at the discontinuity set `{n*x integer}`, while an irrational with a known continued fraction stays computable by the Ostrowski recursion - the obstruction is representation, not irrationality. Finite `N` only: membership in the infinite-depth limsup set is not decidable. And unweighted only: the Mobius-weighted node of this page carries the Mertens-type sum `Sum_{k <= N/b} mu(kb)`, `M(N)` at `b = 1`, and no polynomial-time algorithm for the Mertens function at binary input is known, the best standing near `x^(2/3)` (Deleglise-Rivat 1996) - the one value on this page without an immediate form, an open computational status and explicitly not a hardness result.

**What immediacy does not buy is the RH question, and the reason is sharp.** The Franel-Landau functional needs each node's rank, and the rank's own closed form is `A(x, Q) = Sum_{d <= Q} mu(d) Sum_{e <= Q/d} floor(x*e)` - classical, **Verified** here at `Q = 12, 25, 40` against brute-force enumeration - so the moment the picture is asked where its nodes sit, Mobius enters the formula. Brightness has a `mu`-free closed form; rank does not; only rank carries the difficulty. Franel's 1924 theorem *is* the symbolic all-`Q` reduction of the discrepancy to Mertens-type sums, so the route "generate every frame at once and read off structure" is not unexplored - it is the proof of the equivalence, and it terminates at Mertens. The Mertens meter's natural global readout collapses outright: `Sum_{n <= N} M(floor(N/n)) = 1` identically (**Proved**, classical Mobius inversion; **Verified** at every `N` through 20000 with zero breaches, `lab/py/mertens-meter`), so the weighted picture aggregates to a constant and informs only where it presupposes `M`. One steelman deserves its named kill so it is closed: the stack's divisibility incidence array is the Redheffer matrix up to its first column, whose entries were always trivial and whose determinant is `M(n)`, RH iff `M(n) = O(n^(1/2+eps))`. Immediate entries, untouchable determinant - the same wall this page's honest cap already describes, stated in the highest shape-adjacency object this tree carries.

**The complexity frontier runs beside this page, not through it.** Deciding a pixel's brightness with every input in binary is in P: the constraint set is a rational polytope in fixed dimension three (the both-even parity branch summed alongside the both-odd), and lattice-point counting in fixed dimension is polynomial (Barvinok 1994). **Verified** against the literature, with the caution that the tree's `O(q)` residue pass is polynomial in `q` and so exponential in bit-length - a unary-input algorithm, the honest name for what runs in `lab/py/carpet-stack-address`. The shared scales are the whole engine: moduli `1..N` give the picture polynomially many faces and closed-form extrema, maximum brightness on the diagonal and `floor(N/b)` at `b = 1`. Destroy the sharing - arbitrary binary moduli, one darkened residue class per layer - and "does any point reach maximum brightness" is Simultaneous Incongruences, NP-complete (Garey and Johnson, SP3); make the ambient dimension part of the input and "is any layer lit at this fixed point" is NP-complete (Lagarias 1985), polynomial at every fixed dimension. **Verified** against the literature, both at source. Evaluation stays easy here exactly because the stack shares its scales; hardness begins where the sharing ends, one structural parameter away. A proved no-shortcut theorem for this stack could therefore never have separated P from NP: the problem it would bound is already in P, and what remains bindable there is fine-grained or expressibility only. Nor does the yes touch RH, for the reasons above - both halves of the question were category errors, and each points at the true theorem beside it.
