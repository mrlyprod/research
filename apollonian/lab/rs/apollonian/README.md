# Apollonian

- Grows an integral Apollonian circle packing from its root quadruple in exact integers, checks Descartes on every quadruple it makes, identifies the strip packing's line-tangent circles as the Ford circles, and prints the curvature census.
- A circle is the triple `(k, k x, k y)`: curvature, and the centre scaled by it. A line is `k = 0` with `(k x, k y)` its outward unit normal. Every entry is an integer for every packing here.
- The one move is the reflection `v' = 2(v_1 + v_2 + v_3) - v` applied to all three coordinates at once, which swaps the two circles tangent to the other three; it needs no square root, so the whole run is integer arithmetic.
- Six exact integer identities hold on every quadruple, all of them the bilinear form `B(u, v) = (sum u_i)(sum v_i) - 2 sum u_i v_i` on the four columns: `B(k, k) = B(k, kx) = B(k, ky) = B(kx, ky) = 0` and `B(kx, kx) = B(ky, ky) = -4`. The first is Descartes; the rest are the position half of the same theorem, and the reflection preserves all six because it lies in the form's orthogonal group.
- The two roots are the strip packing `(0, 0, 2, 2)`, two lines a unit apart with the unit-diameter circles at `x = 0` and `x = 1` between them, and the bounded packing `(-1, 2, 2, 3)`, whose outer circle is the unit circle.
- The tree is the standard one: from the root every swap is taken, from any later quadruple every swap but the one just undone, and a branch stops when the new curvature passes `T`. Each node is one new circle, so the census is the node count and the root quadruple is not in it.

## THE STRIP

- The strip packing is periodic in `x`, so the census is taken on one period: the two swaps at the root that replace a line are kept and the two that translate are dropped, which leaves exactly the circles whose centre has `0 < x < 1`. The two root circles at `x = 0` and `x = 1` are the only ones of the whole packing on the boundary and are never counted, so `strayed` tests the strict inequality, `k x <= 0 || k x >= k`.
- A circle of the packing is tangent to the line `y = 0` exactly when `k y = 1`, and to `y = 1` exactly when `k y = k - 1`, so the line-tangent census is read off the integers with no geometry.
- Both line tests are for circles of positive curvature: the line `y = 1` is `(0, 0, 1)` and would pass `k y = 1`, the line `y = 0` is `(0, 0, -1)` and would pass `k y = k - 1`. The tree never hands a line to either test, since a line is only ever replaced.
- `ford` walks the Stern-Brocot tree of Farey pairs in exact integers and checks at every mediant that the reflection returns `(k, k x, k y) = (2 r^2, 2 p r, 1)` for the mediant `p/r` of the parents `a/b` and `c/d`, that the quadruple `(0, 2 b^2, 2 d^2, 2 r^2)` satisfies Descartes, that the other root of that quadratic is `2 (b - d)^2`, and that the two parents are tangent, `(a d - b c)^2 = 1`.
- The same walk sums `floor(Q/b)` over the nodes of `(0, 1)`, adds the node `0/1` at brightness `Q`, and compares the total with `Q(Q + 1)/2`, the Farey stack's brightness on `[0, 1)`.

## RUN

- `cargo run --release --manifest-path research/apollonian/lab/rs/apollonian/Cargo.toml`
- The verbs `ford`, `strip`, `census` and `design` run one section each; about fourteen seconds in total, `census` twelve of them. Prints only, writes nothing.
- `N(T)` in every table is the node count, so the root quadruple is excluded: four circles in the bounded column, one in the strip column.

## WITNESSES

- `ford`: at `b <= 4000` the walk makes 4863601 mediants against `sum_{b <= 4000} phi(b) - 1 = 4863601`, with 0 broken quadruples, 0 coordinate misses and 0 non-tangent parents; brightness 8002000 against `4000 * 4001 / 2`. Same at `b <= 50`, `200` and `1000`: 773, 12231, 304191 nodes and 1275, 20100, 500500 brightness.
- `strip`: 20770674 circles to `T = 2097152`, tree depth 1023, 0 broken quadruples, 0 duplicate circles on the `T <= 2048` control (2448 circles, 2448 distinct), and 0 circles outside `0 < k x < k`. Curvatures fall in the residues `0, 2, 8, 18` mod 24, at 4144636, 6223160, 6241134 and 4161744.
- Of those, 318963 have `k y = 1` and every one of them is a Ford circle: `k/2` a perfect square `b^2`, `k x` a multiple `2 a b` with `gcd(a, b) = 1` and `0 <= a <= b`, 0 off-Ford. The same run reads 318963 tangent to the far line.
- `strip` on the decade grid: `N(T)` is `2, 48, 950, 19298, 390478, 7899138` at `T = 10^1 .. 10^6`, exponent by ratio `1.3802, 1.2965, 1.3078, 1.3061, 1.3060`.
- The line-tangent count at curvature at most `2 Q^2` is `sum_{b <= Q} phi(b) - 1` at `Q = 32`, `181` and `1024`: 323, 10059 and 318963, and the top-tangent count is the same by the strip's reflection symmetry.
- `strip` on the octave grid, exponent as the ratio `log(N(T_2)/N(T_1)) / log(T_2/T_1)`: `1.2925, 1.2722, 1.2716, 1.2925, 1.3073, 1.3011, 1.3064, 1.3050, 1.3056` over the octave grid to `T = 2097152`.
- `census`: the bounded packing `(-1, 2, 2, 3)` to `T = 10^7` gives 555198593 circles, depth 3162, 0 broken quadruples and 67163 distinct circles against 67163 counted on the `T <= 10^4` control. `N(T)` is `5, 165, 3325, 67163, 1359167, 27463391, 555198593` at `T = 10^1 .. 10^7`.
- The same census's exponent by ratio is `1.5185, 1.3043, 1.3053, 1.3061, 1.3055, 1.3057`, and the flat reading `log N(T)/log T` is `1.2492` at `T = 10^7`, still climbing.
- The bounded packing's curvatures fall in the eight residues `2, 3, 6, 11, 14, 15, 18, 23` mod 24.
- `design`: `base^delta` at the literature `delta` misses every integer for `2 <= base <= 100`, the nearest being `52^delta = 174.005426001`, then `68`, `89`, `49`, `23`, `20`; the worst gap is `0.488109816` at base 47.
- The same verb reads the miss as a dimension: the nearest design dimension in the window is `log 351/log 89 = 1.305694144`, off `delta` by `0.000007416`, then `log 247/log 68 = 1.305694579` at `0.000007851`.
