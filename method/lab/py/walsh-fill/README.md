# walsh-fill

- Computes the Walsh coefficients `hat f(S) = 2^-dim sum_(c in F) (-1)^(|S cap c|)` of every base-2 design at `dim 1..4` by one Hadamard product, and their level sums `W_j = sum_(|S| = j) hat f(S)`, kept as the integers `2^dim W_j`.
- Corner `i` is the binary digits of `i`, most significant first, as `mrly_math::bang::factory::code_to_corners` reads a code.
- `section_expansion` counts filled cells on literal grids, every nonempty code at `dim 1..3` and 2000 seeded codes at `dim 4`, sides `1..9`, plus a literal corner histogram against all 65536 codes at `dim 4`, and checks `2^dim fill(N) = sum_j 2^dim W_j N^(dim-j)` at odd `N` and `w N^dim` at even `N`.
- `section_profile` checks `2^dim W_j = sum_w a_w K_j(w)`, `K_j` the Krawtchouk polynomial, on every code at `dim 1..4`, and `K^2 = 2^dim I`.
- `section_mirror` checks `fill_F(-N) = (-1)^dim fill_F'(N)` on every code at `dim 1..4`, tests `fill_F(-N) = +- void_F(N)` against the profile criterion `a_j + a_(dim-j) = C(dim, j)`, the count formula and the self-dual designs.
- `section_roots` maps every root `r != 1/2` of `P_F(n)` to `mu = 1/(2r - 1)` and evaluates `R(mu) = sum_j W_j mu^j` there, on the 77 profiles at `dim 1..3`, and checks `R > 0` on a grid of `(-1, 1)`.
- `section_drift` checks `dim/2 - mean = W_1/(2 W_0)` on every nonempty code at `dim 1..4`, and `2^dim W_1 = U`, the bichromatic edge count, exactly on the down-sets.
- `section_stability` checks `2^dim sum_(c in F) fill_(F+c)(N) = sum_j 2^(2 dim) W^j N^(dim-j)`, `W^j` the squared-coefficient weight, and `sum_c fill_(F+c)(N) = w N^dim` over all flips, by literal counts at odd sides `1..9`, every code at `dim 1..3`.
- `section_entropy` prints `-log2` of the fill ratio at sides `2, 3, 9, ..., 729` for `dim 1 code 1`, `dim 1 code 2`, `dim 2 code 7`, `dim 3 code 23`.
- `section_threshold` checks `2^dim W_0 = S(dim, t)` and `2^dim W_1 = (t+1) C(dim, t+1)` for the rule "at most `t` odd" at `dim 1..8` by a Hadamard transform of each rule, prints the majority ratios, and groups `S(dim, t)/2^dim` over `0 <= t < dim <= 2048` by exact value.
- `section_flat` finds the codes with `W_j = 0` for `0 < j < dim` at `dim 1..4`, checks them against the level-parity criterion, and intersects them with the 896 bent codes at `dim 4`.

## RUN

```
uv run python research/method/lab/py/walsh-fill/walsh.py
```

- From the repo root; one core, about 4 seconds: expansion 0.1s, mirror 0.6s, drift 0.6s, threshold 2.1s, the rest under 0.2s each. A verb name runs one section.

## WITNESSES

- `method.md` The coin, Theorem: 339045 expansion checks at odd sides and 271236 even-side checks, 0 mismatches; `W` of `dim 2 code 7`, `dim 2 code 11`, `dim 3 code 23`, `dim 3 code 232`.
- `method.md` The coin, `N -> -N` swaps fill and void exactly at a balanced profile: swapping designs 2, 4, 40, 2800 at `dim 1..4`, equal to the profile criterion and the count formula, self-dual 2, 4, 16, 256 inside; `dim 3 code 27`, corners `000, 001, 011, 100`, the least non-self-dual swap, `dim 2 code 3` swaps, `dim 3 code 1` does not.
- `method.md` The coin, The roots are zeros of the biased mean: 77 profiles, worst residual `2.7e-15`.
- `method.md` The coin, The drift is half the slope of the log biased mean at the fair coin: drift identity on 65808 codes; `2^dim W_1 = U` on the nonempty down-sets, 2, 5, 19, 167 at `dim 1..4`, and on no other code.
- `method.md` The coin, Noise stability is the fill summed over the flips by filled corners: 2730 literal checks; the four polynomials of the `dim 3 code 23` flip orbit.
- `method.md` The coin, A level costs `log2(2^dim/w)` bits at infinite side: side-3 costs `0.584962501`, `1.584962501`, `0.169925001`, `0.432959407` at `dim 1 code 1`, `dim 1 code 2`, `dim 2 code 7`, `dim 3 code 23`, and `N` times the gap at side 729 against `(W_1/W_0)/ln 2`.
- `method.md` The coin, Threshold rules and Shared limit ratios: majority `1/2, 3/4, 1/2, 11/16, 1/2, 21/32, 1/2, 163/256` at `dim 1..8`; `5/16` and `11/16` at `dim 4`; beyond `1/2`, held by the 1024 odd majorities, 26 shared values to `dim 2048`, each held by exactly two rules.
- `method.md` The coin, Lower-order-free designs: 4, 8, 4, 140 at `dim 1..4`, none of the 896 bent codes at `dim 4`.
