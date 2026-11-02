# four-baselines

- The digit designs of [progressions](../../../notes/progressions.md) against four baselines counted exactly at the same `N = b^L <= 10^12`, `b <= 40`: the greedy base-`3` set, the best Salem-Spencer digit class, the best Behrend shell, and `r_3(N)` where known.
- A design is `F^L` for `F` a largest set progression-free mod `b` containing `0`, found by a depth-first search mod `b`; every such set is kept, their counts asserted equal to the lift column of the census, and a sample of each base is checked by the carry criterion and by brute force at level `3`.
- The baselines are digit DPs over the prefix of `N - 1`: greedy by ternary digits, Behrend by the square sum with `d` to `40` in base `2d - 1`, the digit class by the digit multiset with at most `40000` types per `(d, L)`; `r_3` is the count of the A065825 terms, exact to `N = 209`.
- A fixed digit sum is not a baseline: `2, 6, 10` in base `5` with digits below `3` share digit sum `2`, and a control asserts it.
- `spans` gives the interval of `N` on which a level of a design, a set inside `[0, S]`, beats greedy, Behrend and the `r_3` floor; `literature` sets two sizes from the comment of A065825 against greedy.

## CONTROLS

- Brute-force `r_3(n)` at `n <= 23` and the A065825 count against A003002; greedy against the binary-to-ternary inverse and enumeration; Behrend and digit-class DPs against enumeration with their maximisers progression-free; a memoised recount of Behrend at four large `(d, N)`; truncated design counts against enumeration over every maximiser.

## RUN

- `uv run python -I research/progressions/lab/py/four-baselines/four-baselines.py`, `70` seconds.
- Prints the controls, the designs, the table of `332` values of `N`, the winners, the neighbours, the spans, the literature points, the closest losses and the caps; writes the table to `data/research/progressions/lab/py/four-baselines/four-baselines.txt`.

## WITNESSES

- One winner in the table: `N = 37^7 = 94931877133`, design `10000000` against greedy `8912896`, Behrend `3540559`, digit class `2018016`; `21` ties at `N = 3^L`.
- `172` further rows where a cut design of another base is above the four, base `37` at `150`, the lead `1.5259` at `N = 3^17`.
- The span of the level-`7` set of base `37`: `N` in `[92294880546, 98807184868]`.
