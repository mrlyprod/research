# Farey Discrepancy

Three verbs on one object: the node set of the scale stack, and the Franel-Landau discrepancy read on it. With `rho_1 < ... < rho_m` the nodes of a set in `(0,1]` ascending and `delta_j = rho_j - j/m`, the meter is `S2 = sum delta_j^2` and `S1 = sum |delta_j|`.

## VERBS

- `stack` - the lit nodes of the scale stack are the Farey nodes: the literal reduced lines of every scale `n <= Q`, the sorted window list and the mediant walk are one set, sized by `sum_{k<=Q} phi(k)`, at `Q = 10, 30, 60, 125`. Each node `a/b` is drawn by `floor(Q/b)` of those scales, and scale `n` lights `phi(n)` nodes no earlier scale lit.
- `meter` - the Franel-Landau meter on `F_Q`, seven rungs `Q = 125 .. 8000`: `S2*Q`, `S1/sqrt(Q)` and the local exponents `d ln S / d ln Q`, with the sorted window route cross-checked against the mediant walk at the first four rungs. The top rung is `Q = 8000` with 19455782 nodes.
- `design` - the same meter restricted to a digit design `S_F`, the whole numbers whose every base digit lies in a digit set. Two conventions, `a` and `b` both in `S_F` (strict) and `b` alone in `S_F` (denominator), against the unrestricted control, at base 3 `{0,1}` to `Q = 3^11` and base 10 without the digit 9 to `Q = 10^5`. Per lane it prints `S2` and `S1` with their consecutive-rung ratios and local exponents `e_2` and `e_1`, never a fit; the transplanted normalisers `S2 Q^(e-alpha)` and `S1/Q^(alpha/2)`; the constants `S1/card` and `S2/card`; and the widest gap between consecutive nodes with its left endpoint.

## HOW IT COUNTS

- One Stern-Brocot mediant walk per rung enumerates `F_Q` in ascending order; the walk itself carries no node list and the three lanes filter it as it runs, so each keeps its own rank and the restricted lanes ride the control for free. The verb's memory is the `O(Q)` sieves it builds per rung, a few megabytes at the top.
- The count is checked against a sieve that never enumerates a fraction: `sum_{b in S_F} phi(b)` for the denominator convention, and `sum_{b in S_F} sum_{d | b} mu(d) #{multiples of d in S_F up to b}` for the strict one. The `chk` column is that comparison, `PASS` at every rung of every table.
- Rungs are powers of the base so the design's set is self-similar at each one.
- The cut is the control's cost, one walk step per node of `F_Q`: 9538759028 steps at `Q = 3^11`, 3039650754 at `Q = 10^5`.

## RUN

- `cargo run --release --manifest-path research/farey/lab/rs/farey-discrepancy/Cargo.toml` from the repo root runs all three verbs; append `-- stack`, `-- meter` or `-- design` for one.
- `stack` and `meter` take under a second at 0.7 GB peak; `design` takes a minute at a few megabytes. Prints only, writes nothing.
- `cargo test --release --manifest-path research/farey/lab/rs/farey-discrepancy/Cargo.toml` checks both restricted counts against brute-force enumeration at base 3, base 10 without 9 and the full set.

## WITNESSES

- `farey.md` WHERE THE LINES LAND - brightness `floor(Q/b)` on every lit node, and the node table at `N = 30`.
- `farey.md` PRIMES ARE THE MAXIMALLY NOVEL SCALES - the new nodes of scale `n` number `phi(n)` for `n = 2..30`.
- `farey.md` FRANEL AND LANDAU, 1924 - `card F_Q = sum phi(k)` at `Q = 10, 30, 60, 125`.
- `farey.md` THE METER READS WHAT RH PREDICTS - nodes, `S2*Q`, `S1/sqrt(Q)` and the local exponents at the seven rungs.
- `farey.md` THE METER ON A DIGIT DESIGN - the two exponent tables, the count checks, the constants `S1/card` and `S2/card`, the widest gap between consecutive nodes, and the `S2*Q` and `S1/Q^(alpha/2)` readings of the denominator lane.
- `farey.md` FAREY ORDER IS THE STACK, NOT THE DESIGN - brightness by literal stacking at `Q = 30` on all 278 lit fractions and up to `Q = 125`.
