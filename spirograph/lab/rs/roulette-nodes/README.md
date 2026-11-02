# Roulette Nodes

- Counts the nodes of a design's roulette: a wheel of radius `r` rolls inside or outside a circle of radius `R` with `R/r = a/b` in lowest terms, a pencil at seat `p` in wheel radii draws a trochoid, and a node is a crossing of one curve with itself or of two curves with each other.
- The counter is `mrly_math::roulette::nodes` over `mrly_math::spirograph::trace`: every pair of polyline segments is tested for a proper crossing by orientation signs on a grid of buckets, each pair judged in the one bucket the two bounding boxes first share, so no crossing is counted twice and no pair is missed.
- The orientation sign is `mrly_math::roulette::side`, exact whenever the two differences it is handed are exact, whatever the size of their products: the determinant is taken by the fused multiply-add identity of Kahan, whose error is at most twice the rounding unit times the determinant itself, so it can neither flip a sign nor invent one. The differences of two `f32` endpoints are exact in `f64` while the picture's coordinates keep their exponents within 29 of one another, as these pictures do.
- Three numbers come back and they are not the same number. `total` adds the crossings, and a node where `n` branches meet contributes `n(n - 1)/2` of them; `points` counts the distinct nodes; `branches` adds `n` over the nodes, which is the edge count of the picture as a plane graph. The three agree as `total = points` and `branches = 2 points` exactly when `crowded` is zero.
- Coincident pencils are dropped before counting by `mrly_math::roulette::spread`: on a circle the seats fall into classes under the rotation group of order `gcd(b, 4)`, and on a line or a polygon every distinct seat draws its own curve, two seats of one radius on a line drawing translates of one shape and never one curve.
- Crossings within `4e-4` of the picture's longer side cluster into one node; every cell with a crowded node is printed, and reprinted at a tenth of that tolerance to separate an alignment from a near miss.
- Every count is three sample counts alike, the last two doubled, to five doublings; a cell that never settles is printed with its count.
- The census sweeps both tracks, every `a/b` in lowest terms with `a <= 20` and `b <= 10`, four designs of one, two, four and eight pencils, and seven reaches, and prints every cell that breaks a law.
- Two instruments stand against the counter and share no code with it: the torus reduction, which counts the roots of one real equation in one variable and multiplies by `a`, and a flood, which rasters the curves as walls and counts the faces they cut the plane into.

## THE REDUCTION

- Over `b` orbits the centre turns once for `u` in `[0, 2 pi)`, and the pencil at seat `p` draws `z(u) = A e^(i b u) + r p e^(-i M u)` inside and `A e^(i b u) + r p e^(i M u)` outside, with `A` the centre path radius and `M = a - b` inside, `M = a + b` outside, so `A = r M / b` in both cases.
- Two pencils meet where `z_p(u) = z_q(v)`; with `s = u + v` and `d = u - v` and `e^(i b u) - e^(i b v) = 2i e^(i b s/2) sin(b d/2)` the equation becomes `2i A e^(i a s/2) sin(b d/2) = -r(p e^(i M d/2) - q e^(-i M d/2))` inside, where `b + M = a`, and the same with `e^(-i a s/2)` outside, where `b - M = -a`.
- Taking moduli kills `s` and leaves one equation in `d` alone, `2 A abs(sin(b d/2)) = r abs(p e^(i M d/2) - q e^(-i M d/2))`; for every solution `d` the phase then fixes `a s / 2` modulo `2 pi`, giving exactly `a` values of `s` in `[0, 4 pi)`.
- So with `x = d/2` and `N` the number of `x` in `[0, 2 pi)` solving that equation, two distinct curves meet at `a N / 2` parameter pairs and one curve meets itself at `a (N - 2) / 4`, the two dropped solutions `x = 0` and `x = pi` being the diagonal `u = v`.
- Two thresholds fall out of it: `A b > r abs(p) M` is `abs(p) < 1`, the loop threshold, and `A > r abs(p)` is `abs(p) b < M`, the seat inside the centre path. With `0 < abs(p)` for the centre seat they make one window, `0 < abs(p) < min(1, M/b)`, which is three hypotheses and not one: the loop threshold binds whenever `a > 2b` and outside every track, the seat threshold whenever `a < 2b` inside.
- The seat at the wheel's centre, `abs(p) = 0`, is outside every statement below: it draws the centre circle `b` times over, so its own crossing count is `0` and not `a(b - 1)`, and its `2ab` parameter pairs with another curve sit at `2a` points, each hit once a lap.

## THE PROOF

- **One curve crosses itself `a(b - 1)` times when `0 < abs(p)` and both thresholds hold.** Set `C = r abs(p)`. For `p = q` the reduction reads `A abs(sin(b x)) = C abs(sin(M x))`, that is `A sin(b x) = -+ C sin(M x)`, the zeros of the imaginary parts of `F = A e^(i b x) - C e^(i M x)` and `G = A e^(i b x) + C e^(i M x)`.
- `arg F` turns at the rate `(A^2 b + C^2 M - A C (b + M) cos((b - M) x)) / abs(F)^2`, which is positive for every `x` if and only if `(A - C)(A b - C M) > 0`, and the same rate with the sign of the cosine reversed governs `G`. That product is positive on two components; below both thresholds is the one where `A > C`, and it is the component that matters, since on the other, `abs(p) > 1` with the seat outside the centre path, the argument climbs by `2 pi M` and each branch takes `2M` zeros, not `2b`.
- Below both thresholds `arg F` climbs by `2 pi b` over a period, so `Im F` has exactly `2b` zeros, and `Im G` likewise; with `C > 0` the two share a zero only where `sin(b x) = 0 = sin(M x)`, which since `gcd(b, M) = gcd(a, b) = 1` means `x = 0` or `x = pi`. At `C = 0` the two branches are the same function and share all `2b`, which is the centre seat and its count of nothing.
- Hence `N = 4b - 2` and the self count is `a (N - 2) / 4 = a(b - 1)`.
- **Two distinct curves of one radius cross `2ab` times when `0 < abs(p) = abs(q)` and both thresholds hold.** The right side is `2 r abs(p) abs(sin(M x + f))` with `2f` the angle between the seats, the turning rate is unchanged because a phase only shifts the cosine, each branch again has exactly `2b` zeros, and the two share none unless `sin(b x) = 0 = sin(M x + f)`, so `N = 4b` and the count is `a N / 2 = 2ab`.
- **Two distinct curves of different radii, both with `0 < abs(p)`, cross `2ab` times under a sufficient bound.** The equation squares to `4 A^2 sin^2(b x) = r^2(abs(p) - abs(q))^2 + 4 r^2 abs(p) abs(q) sin^2(M x + f)`, whose right side `rho` never vanishes and obeys `abs(rho') <= 2 r sqrt(abs(p) abs(q)) M`, so `2A sin(b x) - rho` can only turn where `abs(cos(b x)) <= k = r sqrt(abs(p) abs(q)) M / (A b)`; if `2 A sqrt(1 - k^2) > r(abs(p) + abs(q))` it is positive there, hence strictly monotone on both flanks of each of the `b` humps, giving two zeros a hump, `N = 4b`, and again `2ab`.
- That bound is sufficient and far from necessary: inside `7/3` the two design at `abs(p) = 0.950` and `0.672` reads `1.604` against `1.622` and so fails it, while the count is 42, which is `2ab`.
- **A whole design.** With `k` distinct curves, every seat obeying the hypotheses above, the roulette carries `2ab C(k, 2) + k a(b - 1)` crossings, and as many distinct nodes when no three branches meet, which the multiplicity bucket is there to check.

## THE REGIONS

- At a generic reach the roulette is a connected 4-regular plane graph on its nodes, so Euler cuts the plane into `points + 2` regions, the unbounded one among them; connected because every two curves cross, `2ab > 0`, and 4-regular because every node is a transversal double point.
- With a crowded node the graph is not 4-regular and the count is `branches - points + 2`, which is the same number where nothing is crowded.
- The flood confirms both at 1600 and at 2400 pixels: inside `3/1` one seat 2, `5/2` 7, `7/3` 16 at `abs(p) = 0.5` and at `0.9`, `3/1` two seats 8, `5/2` two seats 32, and the alignment cell inside `7/3` at seat `0.527`, where 1288 crossings sit at 1148 nodes carrying 2352 branches, 1206.
- At `k = 1` and `b = 1` there is no node at all and the 2 regions are Jordan's, not Euler's.
- Past the seat threshold the flood tracks the counted crossings and never the law: inside `7/4` two seats at `abs(p) = 0.900` and `0.636` count 84 crossings and flood 86 regions, where the law would say 98 and 100.

## THE CREST

- Past the seat threshold the curve is not the one the proof describes, and the count falls in steps, each step the loss of one hump of `R(x) = abs(sin(b x)) / abs(sin(M x))` under the level `C/A`.
- The steps are exactly the local maxima of `R`, so the self law holds up to the crest, the least local maximum of `R`, which is never below 1 because `R` reaches 1 at the midpoint of any two consecutive zeros of `sin(b x)`, `abs(sin(b x)) = 1` there and `abs(sin(M x)) <= 1`.
- Below the crest the count is `a(b - 1)`; at the top of the band, as `C/A` climbs to `b/M`, it is `a(a - b)`; between them it is a multiple of `a` that never rises.
- The pair law is not so patient: one seat past the seat threshold is enough to lose it, inside `7/4` at seats `0.900` and `0.636` reading 42 against `2ab = 56` while both self counts hold at 21 and the crest, `1.333`, is not yet reached.
- At exactly `C/A = 1` the curve runs through the centre and `a` branches meet there, so the counts printed at that one reach, 25 at `7/5` and 31 at `7/6`, are neither the law nor a multiple of `a`; no census cell sits on it.
- Inside a track with `a >= 2b`, and outside every track, the seat threshold is past the loop threshold and the crest never bites.

## RUN

- `cargo run --release --manifest-path research/spirograph/lab/rs/roulette-nodes/Cargo.toml`
- About two minutes; prints only, writes nothing.
- The verbs `scratch`, `law`, `gap`, `torus`, `centre` and `regions` run one section each.

## WITNESSES

- The law: 5553 cells and 4455175 crossings over both tracks, `a <= 20`, `b <= 10`, four designs and seven reaches, with 143 aligned cells printed, 576 near misses, 26 cells needing a further doubling and none left unsettled.
- Two cells disagree, both by a few crossings on a pair at a near tangency, both reading the law at their base sample count and losing crossings under refinement: inside `20/3` carpet at seat `0.250`, four pairs at 116 against 120, settled 3664 against 3680; outside `10/9` carpet at seat `0.400`, two pairs at 178 against 180, settled 5676 against 5680. The torus reads 3680 and 5680, at two million roots and at eight million.
- The self law by ratio, one pencil inside at `abs(p) = 0.5`: `3/1` and `4/1` and `5/1` and `7/1` at 0, `5/2` at 5, `7/2` at 7, `7/3` at 14, `8/3` at 16, `11/4` at 33, against `a(b - 1)`.
- The pair law: the carpet's eight fills inside `a/1` at tile reach `0.30`, seat `0.20`, give `168, 224, 280, 336, 392` for `a = 3..7`, against `2ab C(8, 2) = 56a`, with no self crossing.
- The whole law: the carpet's eight fills inside `7/3` give `1176 + 112 = 1288` at every tile reach from `0.2` to `1.2`, seat `0.13` to `0.80`, against `2ab C(8, 2) + 8 a(b - 1)`.
- Two curves inside `7/4` at tile reach `0.90`, seats `0.600` and `0.424`, give 56 crossings and 21 self crossings each, against `2ab = 56` and `a(b - 1) = 21`.
- The alignment: the carpet inside `7/3` at tile reach `0.79`, seats `0.527` and `0.373`, carries 28 crowded nodes each holding 6 crossings, which is four branches through one point, and none at `0.78` or `0.80`; the census cell nearest it, seat `0.950`, crowds 28 nodes at three branches.
- The centre seat: inside `5/2` a seat at `0.5` beside a seat at the wheel's centre reads 20 crossings against `2ab = 20`, but they sit at `2a = 10` points and the flood reads 17 regions against the law's 27; the centre curve's own polyline self count is meaningless, 5999 at 6000 samples and 2999 at 3000.
- The crest: 213 cells over 31 ratios, every count a multiple of `a`, no rung rising, every cell reading `a(b - 1)` below the crest and less above it, and all 31 ladders ending at `a(a - b)`.
- The crest bites: inside `7/6` the self count is 35 at seat `0.158`, 21 at `0.175` and 7 at `0.9`, while `a(b - 1) = 35`; inside `7/5` it is 28 at `0.42` and 14 at `0.44`; the crest is 1 at `4/3` and `6/5` and `1.089` at `5/4`, and where `a = 2b - 1` the two ends agree and the ladder is flat.
- The torus against the polyline: twelve named cells, eleven agreeing, the twelfth the outside `1/10` carpet cell at seat `0.100` where the polyline reads 172 against the torus and the law at 156.
- The crate tests `two_curves_on_one_orbit_cross_twice_the_ratio_whatever_their_seats`, `a_curve_below_the_threshold_crosses_itself_a_times_b_less_one`, `an_alignment_reach_crowds_the_nodes_and_the_points_fall_short` and `the_coincident_pencils_collapse_to_one_pencil_a_curve`.
