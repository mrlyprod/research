# node-error

- The error of the spun stacks' node count, recounted to norm bound 8388607 in Rust and read over dyadic windows, Gaussian and Eisenstein.
- It reads the claim of `lab/py/totient-constant` and `lab/py/eisenstein-stack` that the deviation of the count from `c N^2` is of size `N log N`: that reading is Refuted here, and the measured scale is `N^(5/4)`.

## DEFINITIONS

- `K` is `Q(i)` or `Q(sqrt -3)`, class number one, `w` units, `rho = 2 pi/(w sqrt |D_K|)`, `chi` the character of `K`; ideals and associate classes are the same monoid.
- `f(n)` is the sum of `Phi(a)` over the ideals of norm `n`, `S(N) = sum_{n <= N} f(n)` the lit nodes at norm bound `N`, `c = rho/(2 zeta(2) L(2, chi))` and `E(N) = S(N) - c N^2`.
- `a(n) = sum_{d | n} chi(d)` counts the ideals of norm `n`, `A(N) = sum_{n <= N} a(n)`, `P(N) = A(N) - rho N` is the ideal-count error (the circle problem's), and `T(N) = sum_{n <= N} n a(n)` is the sum of the ideal norms.
- `E_1(N) = T(N) - (rho/2) N^2` is the `b = 1` term of the Mobius layer `S(N) = sum_b mu_K(b) T(N/N(b))`, and `E_1(N) = N P(N) - int_0^N P` exactly, by Abel summation.
- For real `x > 0`, `E(x) = sum_b mu_K(b) E_1(x/N(b))` exactly over all ideals `b`, the terms with `N(b) > x` being `-(rho/2) (x/N(b))^2` because `T` is 0 below 1; `layers` holds the coefficients `mu * (mu chi)` of that sum by norm.
- No lattice and no factoring: `f = (n a(n)) * m` with `m = mu * (mu chi)`, the coefficients of `zeta_K(s-1)/zeta_K(s) = (zeta(s-1)/zeta(s)) (L(s-1, chi)/L(s, chi))` read from `1/zeta_K = (1/zeta)(1/L)`; `a`, `m`, `f` come from three divisor sieves.
- `c` is assembled from `rho` and `L(2, chi)` summed to the cut, whose tail is under `1e-13`; it prints `0.260634696495` and `0.235217881630`, the values of `lab/py/totient-constant`.
- Windows are dyadic, `N` in `[2^k, 2^(k+1))`, every integer `N` in each; a statistic is taken over the integers of the window after scaling pointwise by the named power of `N`, and `lo` means the scale uses the window's left edge `2^k`.
- `E` is read at integer `N`; at real `x` in `[N, N+1)` the count is `S(N)`, so `E(x)` differs from `E(N)` by under `2 c (N + 1)`, below every scale here.

## WHAT IT PRINTS

- Per field, one row per window `2^6 .. 2^22`: the count of integers, `rms(E)/N^(5/4)`, the same over `lo^(5/4)`, the mean of `E/N^(5/4)`, `max |E|/N^(5/4)`, `max |E|/N^(4/3)` with its `N`, `rms(E)/(N ln N)` and `max |E|/(N ln N)`.
- Per field, the decomposition of `E` by the `b = 1` term: `rms(E_1)`, `rms(N P)`, `rms(int_0^N P)`, `rms(E - E_1)` over `N^(5/4)` and the correlation of `E` with `E_1`.
- A fit of `ln rms(E)` on `ln N` as an exponent, over all windows, windows `2^14` up, and the last four.
- The cut `N = 8388607`, the statistics on the sliver within `0.02` of it, `E/N^(5/4)` at the cut, and the record `|E|/N^(5/4)` with where it sits.

## RUN

- `cargo run --release --manifest-path research/stack/lab/rs/node-error/Cargo.toml` from the repo root.
- One core, about 4 s, about 600 MB, standard library only, writes nothing.
- `cargo test --release --manifest-path research/stack/lab/rs/node-error/Cargo.toml` pins `S(N)` at the nine norm bounds 50 to 102400 to the counts of the class route in `lab/py/totient-constant` (672 to 2732257376 Gaussian, 630 to 2466558234 Eisenstein), `A(2000)` to a brute lattice-point count of both norm forms, and `E(n)` to the layer sum `sum_m layers(m) E_1(n/m)` within `0.01` at `n = 100, 250, 400` in both fields.

## WITNESSES

- The scale is `N^(5/4)`: from window `2^12` to `2^22`, `rms(E)/N^(5/4)` stays in `0.2708 .. 0.2750` Gaussian and `0.2507 .. 0.2535` Eisenstein, `0.2729` and `0.2512` in the top window; scaled by the left edge it reads `0.4684` and `0.4311`, the `0.47` and `0.43` of the audit; the fitted exponent is `1.2497` and `1.2500` over `2^14 .. 2^22`.
- The scale is not `N log N`: `rms(E)/(N ln N)` climbs from `0.1963` to `0.8712` Gaussian and `0.1812` to `0.8018` Eisenstein between windows `2^6` and `2^22`, and `max |E|/(N ln N)` from `0.4526` to `4.4817` and `0.4377` to `3.6467`, against the `0.44` of the norm bounds to 3200 that read as bounded.
- The mean-square upper half is Proved (`O(N^(5/4))` for the rms over `[N, 2N)`, the last rows of `claims/spun-stack.md`): the measured `rms(E)/N^(5/4)` is flat, at most `0.2750` Gaussian and `0.2535` Eisenstein over `2^12 .. 2^22`, so the proved order is the measured one.
- The proved scale `N^(4/3)` holds with room: `max |E|/N^(4/3)` never exceeds `0.5425` Gaussian and `0.4395` Eisenstein over all windows, and is `0.4002` and `0.3124` in the top window.
- The `b = 1` term carries the scale but not alone: `rms(E_1)/N^(5/4)` is `0.3984` and `0.3218` in the top window, `rms(N P)/N^(5/4)` the same to three digits, `rms(int_0^N P)/N^(5/4)` falls to `0.0050` and `0.0034`, and `rms(E - E_1)/N^(5/4)` is `0.2059` and `0.1239`, correlation of `E` with `E_1` `0.8775` and `0.9359`, so the other layers are smaller but of the same order, and `N P(N)` with `P(N)` of size `N^(1/4)` is of the order `N^(5/4)`.
- The cut: the top window is complete, every integer to `8388607` counted; the sliver within `0.02` of the cut reads `rms(E)/N^(5/4)` `0.2693` and `0.2517`, `max |E|/N^(4/3)` `0.3022` and `0.2739`, and `E/N^(5/4)` `-0.2341` and `+0.3981` at the cut itself; the record `|E|/N^(5/4)` is `1.4517` at `N = 5193840` and `1.2164` at `N = 4159518`, at `0.62` and `0.50` of the cut, so nothing extreme sits at it.
