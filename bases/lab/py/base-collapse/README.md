# base-collapse

- Computes the block digit set, the exact dimension and the brute-force check for a multiplicatively dependent family of bases, all in exact integer and rational arithmetic.
- A cell is a list of bases `b_i = r^(e_i)` sharing one root `r`, each with a digit set containing `0`; the study finds `r`, the exponents `e_i` and `M = lcm(e_i)`, and collapses the family to one design in base `B = r^M`.
- Cells: `G` bases `4, 8` on `{0,1}` and `{0,1,2,3}`; `H` bases `4, 16` on `{0,1}` and `{0,1,4,5}`; `N` bases `9, 27` on `{0,1,2}` and `{0,...,8}`; `T` bases `4, 8, 16` on `{0,1}`, `{0,1,2,3}` and `{0,...,7}`; `I` bases `4, 16` on `{0,1,2}` and the full digit set, the cell whose block count is not a power of the root.
- Height for every brute-force check is `10^13`.

## VERBS

- `blocks` builds the block digit set `A` twice and asserts the two agree: once by sieving all `r^M` base-`r` words against the per-group constraint, once by testing every integer below `r^M` for membership in each original design. It also runs the sharpness cell `b_1 = 2` on `{1}` with `b_2 = 4` on the full set, where `0` is missing from a digit set and the collapse over-counts. Runtime under `0.05` s.
- `dim` prints the exact dimension `log_r(card A) / M`, the parts `log_r(card A_i) / e_i`, the raw sum `sum_i dim A_i - (m - 1)`, the naive budget `max(0, raw sum)` and the gap. The dimension prints as a rational when `card A` is a power of `r` and as the exponent `log_r(card A) / M` with its float otherwise; a cell with a part that is not a power of `r`, such as `I`, gets the dimension and no budget. Runtime under `0.05` s.
- `check` runs three checks per cell to `10^13`: every element of `F(B, A)` below the height passes a digit test in each original base; the joint count from enumerating the lowest-dimension original design and filtering it by the others equals the count of `F(B, A)`, which with the first check gives set equality; and the joint count below `B^level` equals `(card A)^level` at every level with `B^level <= 10^13`. Runtime `16` s, one core, peak resident `252` MB run alone and `253` MB on the three-verb run.
- The enumerated sides are `2^22 - 1` elements at cells `G`, `H` and `T` and `3^14 - 1` at cell `N`; nothing is stored but the frontier.

## RUN

```
uv run python research/bases/lab/py/base-collapse/collapse.py
uv run python research/bases/lab/py/base-collapse/collapse.py blocks dim
```

## WITNESSES

- `bases.md:88` dependence is an equivalence relation and one class is the integer powers of its least member
- `bases.md:92` the collapse theorem, its block digit set and its exact dimension formula
- `bases.md:101-107` the four proof steps, from `e_i | M` to the count law and the open set condition
- `bases.md:105` the sharpness cell `A = {0, 1, 3}` and the elements `4` and `5` in the collapse and not in the joint set
- `bases.md:107` cell `I`: `card A = 9` at `r = 2` and `M = 4`, dimension `log_2(9) / 4 = log_2(3) / 2`, irrational
- `bases.md:115-120` the table rows: `A = {0, 1, 16, 17}` and `1/3` against `1/6`; `A = {0, 1, 4, 5}` and `1/2` against `0`; nine blocks and `1/3` against `1/6`; sixteen blocks and `1/3` against a raw sum of `-1/12` read at `0`
- `bases.md:124` the three checks to `10^13` and the counts `16384`, `1048576`, `6561`, `4096` at level 7, 10, 4, 3
