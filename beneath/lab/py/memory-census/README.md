# Memory Census

- The memory dial: a design accepts a digit word when every window of `k` consecutive digits is allowed, and `k = 1` is the memoryless design of the same code.
- The census of every width-`k` rule at base 2 for `(dim, k)` in `{(1,1), (1,2), (1,3), (1,4), (2,1), (2,2)}`: class counts under two groups, Perron roots as exact algebraic numbers, and the memory number `kappa`.

## THE DIAL

- The alphabet is the `2^dim` digit vectors `d in {0,1}^dim`, each read as the corner integer `c = sum_i d[i] 2^i` with axis `i` on bit `i`.
- A width-`k` rule is a set `W` of windows `(d_1, ..., d_k)`; the window integer is `w = sum_(j=1..k) c_j 2^(dim (k - j))`, the first digit most significant.
- The rule code has bit `w` set iff window `w` is allowed, so a rule at `(dim, k)` is a code in `[0, 2^(2^(k dim)))`: `4, 16, 256, 65536` rules at `dim 1` and `k = 1..4`, `16` and `65536` at `dim 2` and `k = 1, 2`.
- A word `d_1 d_2 ... d_level` is read coarsest digit first and is accepted iff every window of `k` consecutive digits is allowed; every word of length `level < k` is accepted.
- `N_W(level)` is the number of accepted words of length `level`, and `#W` the number of allowed windows.
- At `k = 1` the code is a bang code and the accepted words of length `level` are the cells of the design of that code at that level, so the memory dial at `k = 1` is the plain design. **Verified** (memory-census, the `(dim,1)` rows of `census.csv`).
- `mrly_math::bang::universe::corners(dim)` emits the corner vector row first and `corner_index` folds it most significant first, so a crate design's corner integer at `dim 2` is `c = x + 2y`, bit `0` the column: the same map this census uses, and no code label moves. **Verified** (the `demos` test `width_one_is_the_plane_design_cell_for_cell` pins codes `11` and `13`, which the axis swap exchanges and whose designs differ).

## THE TRANSFER MATRIX

- The states are the `2^(dim(k-1))` windows of width `k - 1`; `A_W[s][t] = 1` iff `s` and `t` overlap in `k - 2` digits and the `k`-window they form is allowed. At `k = 1` there is one state and `A_W = [#W]`.
- `N_W(level) = 1^T A_W^(level-k+1) 1` for `level >= k - 1`. **Proved** (a word of length `level` is a path of `level - k + 1` steps: its width-`(k-1)` windows in order are the states, each consecutive pair forms one width-`k` window, and every such path is an accepted word).
- `rho(W)` is the Perron root, the spectral radius of `A_W`; it is a root of the characteristic polynomial and dominates every eigenvalue in modulus. **Proved** (Perron-Frobenius for a nonnegative matrix).
- `lim_level N_W(level)^(1/level) = rho(W)` whenever `rho(W) > 0`. **Proved** (`max_(s,t) (A^n)[s][t] <= 1^T A^n 1 <= S^2 max_(s,t) (A^n)[s][t]` with `S` the number of states, and the spectral radius is the limit of the `n`-th root of the matrix norm).
- `log_2 rho(W)` is the growth exponent of `N_W(level)`, which is what this study measures. That it is also the Hausdorff dimension per axis of the accepted set is **Conjecture** here, not proved and not cited.

## THE MEMORY NUMBER

- `kappa(W) = log_2(#W) / k - log_2 rho(W)`, defined when `rho(W) > 0`; a rule with `rho(W) = 0` accepts finitely many words and is called dead, and the empty rule `W = {}` is dead with `#W = 0`, so no dead rule carries a `kappa`. The crate returns `0.0` on the empty rule as a convention, not as a value of the formula.
- `kappa(W) >= 0` for every rule. **Proved** (an accepted word of length `mk` splits into `m` disjoint allowed windows, so `N_W(mk) <= #W^m`; take `m`-th roots and the limit to get `rho^k <= #W`).
- `kappa(W) = 0` when `W = F^k` is a product with `F` a non-empty set of digit vectors, the memoryless case. **Proved** (every word over `F` is accepted, so `N_W(level) = #F^level`, `rho = #F >= 1` and `#W = #F^k`). At `F = {}` the rule is dead and `kappa` is not defined.
- So `kappa` measures how much of the window budget the rule wastes: zero means the windows compose freely, and the further above zero, the more the rule forbids by interference rather than by banning digits.
- In the whole census `kappa(W) = 0` holds on exactly the non-empty product classes and on nothing else: `2` classes at every `(1,k)` and `5` at every `(2,k)`, with `0` counterexamples over `19563` live classes. **Verified** (memory-census, `census.csv` columns `classes_kappa0` and `kappa0_nonproduct`; the test is exact, `kappa = 0` iff the minimal polynomial of `rho` divides `x^k - #W`).
- For `k >= 2` the largest `kappa` in the census is attained at `rho = 1`, by the largest rule of zero entropy, and the maximum is a tie over `1, 3, 4` classes at `(1,2), (1,3), (1,4)` and `3` at `(2,2)`. **Verified** (memory-census, `census.csv` columns `kappa_max`, `kappa_max_code`, `kappa_max_ties`).
- At `k = 1` every live rule is a product, `kappa` is identically `0`, and the maximum is attained on every live class at once, the full rule at `rho = 2^dim` included. **Verified** (memory-census, the `(dim,1)` rows).
- `kappa_max = 0` at `k = 1`, `log_2(3)/2 = 0.792481` at `(1,2)` on code `11`, `log_2(6)/3 = 0.861654` at `(1,3)` on code `175`, `log_2(13)/4 = 0.925110` at `(1,4)` on code `49071`, and `log_2(10)/2 = 1.660964` at `(2,2)` on code `36079`, each the least code of its tie. **Verified** (memory-census, `census.csv`).
- The window budget of zero entropy, the largest `#W` a live class with `rho = 1` carries, is `1, 3, 6, 13` at `dim 1` and `k = 1..4` and `1, 10` at `dim 2`: a rule may allow that many windows and still accept subexponentially many words, and code `11` at `(1,2)`, which allows `00, 01, 11`, accepts every `0^a 1^b` at run lengths `a, b >= 0`, `N_W(level) = level + 1`. **Verified** (memory-census, `census.csv` column `rho1_windows`).

## THE GROUP

- `G_(dim,k)` is `B_dim`, the signed permutations of the `dim` axes, acting by the *same* element on every digit of a window, together with window reversal `(d_1, ..., d_k) -> (d_k, ..., d_1)`.
- Every element preserves `N_W(level)` for every `level`. **Proved** (a diagonal `B_dim` element relabels the digit alphabet by a bijection, which relabels the states and conjugates `A_W` by a permutation matrix; reversal sends a word to its reverse, a bijection of accepted words, and transposes `A_W`).
- Reversal therefore fixes the characteristic polynomial as well, since `A` and `A^T` share it.
- A width-`k` rule at `dim` axes is the same object as a subset of the `k dim`-cube: the same censuses of the `k dim`-cube counted with a smaller group. `#G_(dim,k) = 2^(dim+1) dim!` for `k >= 2` against `#B_(k dim) = 2^(k dim) (k dim)!`; at `k = 1` reversal is trivial and `#G_(dim,1) = 2^dim dim!`.
- The distinct window permutations are `2, 4, 4, 4` at `dim 1` and `k = 1..4` (reversal is trivial at `k = 1`) and `8, 16` at `dim 2` and `k = 1, 2`. **Verified** (memory-census, `census.csv` columns `group_G`, `group_B`).
- `G_(1,4) < G_(2,2) < B_4` as permutation groups of the `4`-cube: flipping all four bits is diagonal `B_2` flipping both axes in both blocks, and the `(1,4)` reversal is the `(2,2)` block swap composed with the diagonal axis swap. **Proved** (the generators are exhibited), and the class counts nest as `16960 > 4660 > 402`. **Verified** (memory-census; `402` is A000616 at `4`).

## THE CENSUS

- One row per `(dim, k)` in `census.csv`, one row per distinct Perron minimal polynomial in `classes.csv`.
- Classes under `G_(dim,k)`: `3, 9, 88, 16960` at `dim 1` and `k = 1..4`; `6, 4660` at `dim 2` and `k = 1, 2`. **Verified** (memory-census, orbit walk, agreeing with an independent Burnside average on every row).
- Classes under diagonal `B_dim` alone, no reversal: `3, 10, 136, 32896` at `dim 1`; `6, 8548` at `dim 2`. **Verified** (memory-census, same two verbs).
- At `k = 1` the two groups agree and the counts are A000616 at `dim`: `3` at `dim 1` and `6` at `dim 2`, the design census reproduced. **Verified** (memory-census, the `(dim,1)` rows).
- At `dim 1` and no reversal the count is `2^(2^k - 1) + 2^(2^(k-1) - 1)`. **Proved** (Burnside on the group of order 2: the flip acts on the `2^k` windows as `w -> 2^k - 1 - w`, all `2^(k-1)` two-cycles, so `(2^(2^k) + 2^(2^(k-1)))/2`).
- Against A000616 at `k dim`, the ratio of classes under `G_(dim,k)` is `1, 3/2, 4, 42.19...` at `dim 1` and `k = 1..4` and `1, 11.59...` at `dim 2`: the groups agree at `k = 1` and the smaller group splits cube classes from `k = 2` on. **Verified** (memory-census, `census.csv` column `a000616`).
- Burnside extends the class counts past the orbit walk at no cost: under `G_(1,k)` for `k = 1..8` they are `3, 9, 88, 16960, 1074036736, 4611686053860868096, 85070591730234617055658644612208132096, 28948022309329048855892746252171977006958709724020498949042189405102555529216`. **Verified** (memory-census, Burnside extension).
- The `G_(1,k)` count has the closed form `a(2m) = 2^(2^(2m)-2) + 2^(2^(2m-1)-2) + 2^(2^(2m-1)+2^(m-1)-1)` for `m >= 1` and `a(2m+1) = 2^(2^(2m+1)-2) + 2^(2^(2m)-1) + 2^(2^(2m)+2^m-2)` for `m >= 0`. **Proved** (Burnside over the order-4 group: the digit flip fixes no window, reversal fixes the `2^ceil(k/2)` palindromes, and flip-reversal fixes the `2^(k/2)` antipalindromes at even `k` and none at odd `k`; memory-census, verb `burnside`, where the cycle index and the closed form agree at `k = 1..11`).
- Under `G_(2,k)` for `k = 1..4`: `6, 4660, 1152921592116822016, 7237005577332262213973186563042994284449319951280334537682043176672785596416`. **Verified** (memory-census, Burnside extension).
- Dead classes, those with `rho = 0`: `1, 2, 13, 2093` at `dim 1` and `1, 53` at `dim 2`. **Verified** (memory-census, `census.csv` column `dead_classes`).

## THE PERRON ROOTS

- Every distinct transfer matrix is batched into one `gp -q` script: characteristic polynomial from an exact integer Faddeev-LeVerrier in Python, then `factor` over `Q` and `polrootsreal` on each factor in PARI; the factor carrying the largest real root is the minimal polynomial of `rho`.
- Distinct characteristic polynomials: `3, 6, 23, 431` at `dim 1` and `k = 1..4`; `5, 333` at `dim 2`. Distinct minimal polynomials of `rho`: `3, 4, 10, 177` at `dim 1`; `5, 185` at `dim 2`. **Verified** (memory-census, `census.csv`).
- At `dim 1` and `k >= 2` every transfer matrix has determinant in `{-1, 0, 1}`, so the constant term of every characteristic polynomial is `0`, `1` or `-1`. **Proved** (rows `s` and `s + 2^(k-2)` are both supported on the columns `2s` and `2s+1` taken modulo `2^(k-1)`, and those column pairs partition the columns as `s` runs over `0..2^(k-2)-1`, so the matrix is a row permutation of a block diagonal matrix with `2^(k-2)` blocks of size `2 x 2` over `{0,1}`; memory-census, verb `lemmas`, where the determinant and the signed block product agree and land in `{-1, 0, 1}` over all `16, 256, 65536` rules at `k = 2, 3, 4`).
- Distinct minimal polynomials of `rho` are distinct `rho`, so the minimal polynomial count is a count of growth rates. **Proved** (every conjugate of `rho(W)` is a root of the characteristic polynomial of `A_W` and so an eigenvalue of `A_W`, hence at most `rho(W)` in modulus, so two conjugate Perron roots are equal in modulus and, both being nonnegative, equal; memory-census, verb `lemmas`, where no conjugate exceeds `rho` on any of the `463` characteristic polynomials at `dim 1` and `k = 1..4` and the `3, 4, 10, 177` minimal polynomials carry `3, 4, 10, 177` distinct `rho`).
- The golden ratio `x^2 - x - 1`, `rho = 1.618033988749`, is the Perron root of code `7` at `(1,2)`, the rule forbidding the window `11`. **Verified** (memory-census, `classes.csv`).
- At `(1,3)` the four named roots land on the four expected codes, each the least code of its class: `x^3 - x^2 - x - 1` tribonacci `1.839286755214` on code `127`, `x^2 - x - 1` golden on code `55`, `x^3 - x^2 - 1` supergolden `1.465571231876` on code `23`, `x^3 - x - 1` plastic `1.324717957244` on code `54`. **Verified** (memory-census, `classes.csv`).
- Their class counts at `(1,3)` are `1, 7, 5, 4` and their rule counts `2, 16, 14, 12`; at `(1,4)` the same four polynomials carry `7, 193, 526, 696` classes. **Verified** (memory-census, `classes.csv`).
- At `(2,2)` the same four occur on codes `327, 19, 323, 326` with `121, 588, 54, 48` classes, so the named roots are not a `dim 1` accident. **Verified** (memory-census, `classes.csv`).
- The two largest non-integer roots at `(1,4)` are `x^4 - x^3 - x^2 - x - 1` at `1.927561975482` on code `32767` and `x^4 - 2x^3 + x^2 - 2x + 1` at `1.883203505913` on code `64511`; at `(2,2)` they are `x^2 - 3x - 3` at `3.791287847477` on code `32767` and `x^2 - 4x + 1` at `3.732050807568` on code `49151`. **Verified** (memory-census, `classes.csv`).
- Every Perron root at width `k` occurs again at width `k + 1`: `0` missing on all four tested steps. **Proved** (`W' = {(d_1..d_(k+1)) : (d_1..d_k) in W and (d_2..d_(k+1)) in W}` accepts the same words of length `k + 1` and above, which is all `rho` needs; at `level = k` exactly `W'` has no window and accepts every word) and **Verified** (memory-census, the `nest` lines).
- A root is *strict* when every conjugate is smaller in modulus; the decision is numeric, PARI's complex roots against a `1e-20` gap, not exact, and the `strict` column is empty on the dead row `x`. The non-strict roots are `12` of `177` at `(1,4)` and `5` of `185` at `(2,2)`, and every one of them is `p(x^m)` for some `m >= 2` with `p` the minimal polynomial of a strict root already in the census. **Verified** (memory-census, `census.csv` columns `weak_perron_polys`, `weak_are_radicals`).
- That identification is a census observation, not a theorem here: the peripheral spectrum of a nonnegative matrix is `rho` times roots of unity, and the step from there to a minimal polynomial in `x^m` needs the conjugate set of `rho` stable under `x -> zeta x`, which this study does not argue, and these transfer matrices are reducible in general.
- So every root of a live rule is either `1` or a radical of a strict Perron number in `(1, 2^dim]`; the dead rules carry `rho = 0`, on `1, 2, 13, 2093` classes at `dim 1` and `1, 53` at `dim 2`, and `0` is a radical of nothing. Witnesses: `x^2 - 2`, `x^3 - 2` and `x^4 - 2` at `(1,4)`, and `x^4 - x^2 - 1` and `x^6 - x^3 - 1`, the golden ratio's square and cube roots. **Verified** (memory-census, `classes.csv`).

## OEIS

- `3, 10, 136, 32896, 2147516416`, the `dim 1` count with no reversal, greps three hits in the local dump: A055708, A056006, A191363. All three are lists of integers with a `sigma` property and the agreement is a coincidence of the closed form `2^(m-1)(2^m + 1)` at `m = 2^(k-1)`; none is this census. **Verified** (grep `,3,10,136,32896` in the local dump, then the name of each hit read at source).
- `3, 9, 88, 16960` (classes under `G_(1,k)`), `6, 4660` (under `G_(2,k)`), `3, 4, 10, 177` (distinct Perron minimal polynomials at `dim 1`) and `3, 6, 23, 431` (distinct characteristic polynomials at `dim 1`) each grep to zero hits. **Verified** (grep `,3,9,88,16960`, `,6,4660,`, `,3,4,10,177,`, `,3,6,23,431,` in the local dump).
- Each grep is a comma-prefixed consecutive run against the local dump, so a zero hit is the absence of that run from the dump, not the absence of the sequence from OEIS.

## RUN

- `python3 research/beneath/lab/py/memory-census/memory.py` from the repository root; the standard library and `gp` only, no third-party package.
- Writes `census.csv` and `classes.csv` beside this file.
- Runtimes on a laptop: orbit walk `0.04s` at `(1,4)` and `0.14s` at `(2,2)`, characteristic polynomials `3.61s` at `(1,4)`, PARI `0.17s` at `(1,4)` and `0.06s` at `(2,2)`, Burnside extension `0.01s`, whole study `4.3s`.
- `python3 research/beneath/lab/py/memory-census/memory.py burnside` prints the class count under `G_(1,k)` at `k = 1..11` from the cycle index of the order-4 group on the `2^k` windows, checks the closed form against every term, and prints the orbit count of the same group on the windows themselves, A005418 at `k`; `0.07s`.
- `python3 research/beneath/lab/py/memory-census/memory.py lemmas` checks the determinant lemma over every rule at `k = 2, 3, 4` against a Bareiss determinant and against the signed block product, and the Perron lemma over every characteristic polynomial at `k = 1..4`; `13.4s`.

## WITNESSES

- [beneath](../../../notes/memory-dial.md), The census - the six-row table of rules, classes under both groups, dead classes and minimal polynomials, and the `19563` live classes: `census.csv` and the `live classes` line of the run.
- [beneath](../../../notes/memory-dial.md), The coupling - `kappa = 0` on exactly the non-empty product classes, `2` at every `(1,k)` and `5` at every `(2,k)`, with `0` counterexamples: `census.csv` columns `classes_kappa0` and `kappa0_nonproduct`.
- [beneath](../../../notes/memory-dial.md), The coupling - the maxima `log_2(3)/2`, `log_2(6)/3`, `log_2(13)/4` and `log_2(10)/2` on the least codes `11`, `175`, `49071` and `36079`, with ties of `1, 3, 4, 3`: `census.csv` columns `kappa_max`, `kappa_max_code`, `kappa_max_ties`.
- [beneath](../../../notes/memory-dial.md), The famous constants are one notch in - the four named roots on the least codes `7`, `23`, `54`, `127` at `dim 1` and `19`, `323`, `326`, `327` at `(2,2)`: `classes.csv`.
- [beneath](../../../notes/memory-dial.md), What the dial buys - the distinct minimal polynomials of `rho`, `3, 4, 10, 177` at `dim 1` and `5, 185` at `dim 2`, against `3, 6, 23, 431` and `5, 333` characteristic polynomials: `census.csv`.
- [beneath](../../../notes/memory-dial.md), The memory dial - the closed form for the classes under `G_(1,k)` at `k = 1..11`, the determinant lemma at `dim 1` and `k = 2, 3, 4`, and the Perron lemma at `k = 1..4`: the verbs `burnside` and `lemmas`.

## COLUMNS

- `census.csv` one row per `(dim, k)`: rules, classes under both groups, group orders, polynomial counts, the `kappa` extremes with the least code attaining the maximum, the size of that tie, the windows it allows and the zero-entropy window budget, the weak-root counts, the Burnside cross-check, dead classes and per-verb runtimes.
- `classes.csv` one row per distinct minimal polynomial of `rho`: degree, `rho` to 12 digits, class and rule counts, the least code carrying it, its window count, the `kappa` range over its classes, strictness (empty on the dead row `x`) and the name when the polynomial proves one.
