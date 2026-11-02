---
title: Design masks
lead: A Toeplitz attention mask whose offsets form a digit design reaches every offset in two hops exactly when its digit sumset covers the digits; the sparsest ties the strided mask in density and loses to it on tiny Shakespeare, which reads the nearest offsets.
figure: research-masks
slug: masks
---

## Offset masks

A Toeplitz mask lets query `i` see key `j` when `i - j` lies in an offset set `D`; two layers reach `D + D`. The digit design `D(F)` holds the integers with all base-`b` digits in `F`, `0` in `F`. Generator: [the lab](../lab/py/design-masks/); rows: [the claims](../claims/offset-masks.md).

**For `T >= b`, `D(F)` cut to `[0, T)` reaches every offset below `T` in two hops exactly when `F + F` holds every digit `0..b-1`; otherwise the least offset missed is the least digit outside `F + F`. Proved.** If `F + F` holds every digit, split each digit of `n < T` as `n_i = f_i + g_i` in `F`; then `x = sum f_i b^i` and `y = sum g_i b^i` lie in `D(F)`, sum to `n` and are at most `n`. Otherwise let `m` be the least digit outside `F + F`: each `n < m` is `f + g` with `f, g` in `F`, single digits of `D(F)`, while `m = x + y` in `D(F)` forces `x, y <= m < b`, so `x, y` in `F`, impossible. For `F = {0, 1}` at base `3` this is the classical fact that every integer is a sum of two with ternary digits `0, 1`; nothing new is claimed.

**Full 2-hop coverage at side `81` needs at least `13` offsets. Proved.** `k` offsets make at most `k(k+1)/2` sums, and `78 < 81`.

**Census, Verified.** `8` of `112` distinct digit designs reach every offset. The sparsest, `d3-01`, holds `16` offsets and `976` of `3321` causal pairs; strided `l = 7` holds `18` and `975`. No random set of `16` offsets reaches full coverage in `400` draws. An integer program beats both: `16` offsets, `869` pairs, solver bound `629`. `general-2`, the 2d design that beat its window on copy, is the Toeplitz mask of `0` and the odd offsets.

**Tiny Shakespeare, Verified.** Validation nats, `600` steps, `3` paired seeds: `window-8` `1.886`, `strided-7` `1.928`, `d5-012` `1.951`, dense `1.976`, `d3-01` `1.978`, Cantor `d3-02` `2.215`. Strided beats the three designs by `0.043` per seed over `6` seeds, t `11.1`, df `5`. The designs beat their random-offset twins only through twins missing offset `1`. Against its window, missing `1` costs `0.29` to `0.59` nats, missing `2` costs `0.04` to `0.09`, holding both at most `0.03`. Every window beats dense, Refuting the predicted tie within `0.02`.

**What it shows.** A two-layer character model at `600` steps reads its nearest offsets, so this ranks learning speed, not capacity; full 2-hop coverage bought nothing measurable. The next test is a mask that holds `0..2` and spends its far offsets like the integer-program set.
