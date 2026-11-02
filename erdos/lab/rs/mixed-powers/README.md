# mixed-powers

- Computes, for a finite base set `D` of integers `>= 3` and a level `k >= 1`, whether every large integer is a sum of distinct terms of the multiset `M_k(D)` of powers `d^j`, `d` in `D`, `j >= k`, one copy per base: Erdos problem 124 at one level.
- The certificate is the surplus certificate of `research/erdos/notes/problem-124.md`: the sums `P_n` of the `n` least terms are held as one bit array, each term folded in by one shift-or pass, and after each term `T = 1 + max{x <= S_n/2 : x not in P_n}` is read by a word scan down from `S_n/2`; when `T <= a_(n+1)` and the step `a_(m+1) <= S_m - 2T + 1` holds until the surplus `(sigma - 1) a_(m+1) >= C_k + 2T - 1` takes over, the level is complete and `F = T - 1` is the largest integer that is no sum. The step check runs over the terms up to `10^30` in exact integer arithmetic, `sigma` and `C_k` over the common denominator of the `d - 1`.
- The census enumerates the base sets in `[3, R]` minimal for `sigma > 1` and `gcd = 1`, every other such set containing one, and certifies each at `k = 1, 2, ...` until a cell passes the bit cap.
- The wide certificate keeps only `P_n` meet `[0, W]` as bits and the holes of `P_n` in `(W, S_n/2]` as a sorted list capped at `2^25`, exact for `S_n < 2^63`, by the update rule proved on the note, so a level closes on `W` bits where the full array needs `S_n`; at the close it sorts every non-sum by its route, below the least term, window or pair. A capped cell returns `Cap` either because the list passes `2^25` or because the candidates of one step past `S_n` do, without saying which.
- The route verb reads, per minimal base set, the least `sum(M_1(D) below N)/N` over terms `N` in `[10^E, 10^30]`, truncated to six places, the quantity that bounds two-part partitions below `sigma = 2`.
- The 3-adic verbs hold Erdos problem 406 in the relaxation of Lagarias and Abram and Lagarias: `C(1, M_1, ..., M_r)` is the set of 3-adic integers `lambda` with every digit of `lambda` and of each `M_i lambda` in `{0, 1}`, read by the carry automaton whose state is the vector of carries of the `M_i lambda`; the dimension is `log_3` of the largest spectral radius over its reachable strongly connected components, bracketed by the least and largest ratio of a power iterate, and the prefix test decides `C = {0}` from the `M_i mod 3^D`, the powers `2^m` reduced modulo `3^D`, so the exponents matter only modulo `2 3^(D - 1)`, with `D <= 79`, exact in `u128`.

## VERBS

- `census R K BITS THREADS`: every minimal base set in `[3, R]` at `k = 1..K` with a `2^BITS`-bit array per cell; prints `F`, the count of positive non-sums up to `F`, the index `n` and term `a_n` of the certificate, `S_n`, the term where the surplus takes over, and the depth histogram. Runtime `2.9` s at `census 10 5 32 4`, `8.7` s at `census 12 5 32 4`, on `4` threads, `512` MB per thread.
- `deep R K BITS THREADS`: the census through the wide certificate at `W = 2^BITS`, each cell also printing its non-sums by route and the route of `F`, and per level the cells whose `F` comes by each route. Runtime `21.3` s at `deep 10 6 31 4`, `2.8` s at `deep 12 3 31 4`, on `4` threads.
- `wcell D K BITS`: one cell through the wide certificate. Runtime `1.8` s at `wcell 3,4,5 5 33`, `0.9` s at `wcell 3,4,6 4 33`, `3.6` to `5.4` s at `W = 2^33` and `13.1` to `14.7` s at `W = 2^34` on the `k = 5` cells of the note, `2` GB at `W = 2^34`.
- `wcontrol BITS`: every cell of the minimal sets in `[3, 10]` at `k <= 4` that the full certificate settles within `2^30` bits, against the wide certificate at `W = 2^BITS`. Runtime `6.8` s at `wcontrol 16`: `109` agree, none differ. It samples only `k <= 4`, bases up to `10` and `S_n <= 2^30`, never the cap path or the scale of the `k = 5` cells.
- `split D K X COUNT`: the non-sums in `[X, X + COUNT)` by a meet in the middle over the terms up to `X + COUNT`, independent of the bit array. Runtime `0.3` s at `split 3,4,5 5 31333640736 64`.
- `chainless BITS`: the multiset `4, 5, 6, 7, 12, 18, 31, ...` of the note, `a_(n+2) = S_n - y_n - y_(n+1)`, up to `2^BITS`: windows, the least of `2 (S_n - a_(n+1)) - a_(n+1)`, the values of `a_(n+3) - a_(n+2) - a_(n+1)`, the chain of non-sums `a_(n+1) + y_n` and every non-sum by route. Runtime under `0.1` s at `chainless 27`.
- `forced BITS`: the multiset of `chainless` with `a_(n+2) = S_n + 2` or `S_n + 3` forced at `n >= 6`, `n = 0 mod 6, 8, 10`: open windows, chained pairs, the least of `2 (S_n - a_(n+1)) - a_(n+1)`, the non-sums below the last term and those in the last `12` term intervals. Runtime under `0.1` s at `forced 26`.
- `windowed BITS`: the multiset `V = 4, 5, 6, 7, 13, 19, 47, ...` of the note, `a_(n+2) = S_n - u_n - u_(n+1)` with `u_n` set from the chain holes `x_n = a_(n+1) + u_n`, up to `2^BITS`: sortedness, the windows and their least width over the term, chained pairs, windows holding a non-sum, the least of `12 (S_n - a_(n+1)) - a_(n+1)`, the largest `a_(n+1) - 3 a_n`, the values of `a_(n+3) - a_(n+2) - a_(n+1)` at `n = 3 mod 5` from `3`, the chain of non-sums `x_n` and every non-sum by route. Runtime under `0.1` s at `windowed 27`.
- `jumps R K BITS THREADS`: every minimal base set in `[3, R]` at `k = 1..K` through the wide certificate at `W = 2^BITS`, logging `S_n` and `T_n` at every step: checks the jump law of the note on every step, and prints per cell the rises of `T`, the birth step `b` of `F`, `F - a_(b+1)`, whether `T_b - 1 > R_b`, `R_b = S_b - a_(b+1)`, the route of `F`, the peak of `T_n/a_(n+1)` after the last base enters, and per capped cell the run of rises open at the cap; then per level the counts, the largest and least peaks, the gap from the last rise to the certificate, the cells run and the violations of the jump law and the route lemma; it stops a base set at its first capped level. Runtime `16.6` s at `jumps 10 6 31 4` on `4` threads, `3.5` GB peak resident, `1.8` GB on `2` threads; the holes past the window ride as `u64` in two buffers swapped each step, `2^28` bytes each, so a capped cell holds `W/8` bytes of bits and `2^29` of holes.
- `trace D K BITS`: one cell step by step: `a_n`, its base, `S_n`, `a_(n+1)`, `R_n = S_n - a_(n+1)`, `T_n`, `T_n/a_(n+1)`, and `rises` on the rows `n` with `T_(n+1) > T_n`. Runtime under `0.1` s at `trace 3,4,5 3 26`.
- `cell D K BITS`: one base set, comma separated, at one level. Runtime `1.0` s at `cell 3,4,5 4 34`, `1.3` s at `cell 3,6,9,10,12 2 34`, `2` GB.
- `control`: nine cells against a plain knapsack over every term up to `4 F`, same `F` and same count of non-sums, and a `gcd = 3` set that never certifies below `2^24`. Runtime `0.9` s.
- `set D K B`: the same base set read as a set of powers, ties counted once, by a plain knapsack to `B`; prints the largest non-sum below `B`. Runtime under `0.01` s at `set 3,5,6,9 1 20000`, which reads `649` against the multiset's `F = 22`.
- `windows D K E F`: over the terms `a_(n+1)` in `[10^E, 10^F]` of `M_K(D)`, the count of windows `a_(n+2) > S_n` and of pairs `m < n` of windows where the upper window `(S_m, a_(m+2))` meets the lower window `(S_n - a_(n+1), a_(n+2) - a_(n+1))`. Runtime under `0.01` s at `windows 3,4,5 1 6 30`: `27` of `123`, no pair.
- `graham T P Q BITS`: the set `floor(T (P/Q)^n)`, `n >= 1`, up to `2^BITS`, exact in integers; the count of windows and of windows holding a non-sum, and the largest non-sum up to `2^(BITS-1)`. Runtime under `0.1` s at `graham 2 5 3 26`: `31` windows, all holding a non-sum, largest `23559582`.
- `refute BITS`: the set `2 F_m - 1`, `m >= 2`, up to `2^BITS`, which refutes the surplus lemma L of the note: the least of `2 (S_n - a_(n+1)) - a_(n+1)`, the windows and those holding a non-sum, and the non-sums below `10^6`. Runtime under `0.1` s at `refute 27`: `36` windows, all holding a non-sum.
- `ternary N HI THREADS`: `g_3(n)` of Erdos problem 817 for `n = 1..N`, exhaustive: for each largest element in turn, a search from the top down over sets whose `3^n` ternary sums are distinct, one bit array of sums, pruned by `g_3(r)` of the run before and by `sum a_i^2 >= (9^n - 1)/8`. Runtime `5.9` s at `ternary 6 200 8`.
- `probe N TOP THREADS`: the same search at one largest element. Runtime `36` s at `probe 7 380 8`.
- `band LO HI W THREADS`: the least largest element in `[LO, HI]` of an admissible 7-set whose six largest elements lie within `W` of it, exhaustive inside that band. Runtime `2.8` s at `band 419 475 70 8`, which finds `{302, 409, 447, 459, 465, 466, 474}`, and `114` s at `band 419 473 100 8`, which finds none.
- `offsets LO HI`: the offset families `{0, 1, 3, 8, 22, 60}`, `{0, 2, 6, 9, 23, 61}`, `{0, 2, 5, 7, 21, 60}` below a largest element plus one free offset. Runtime `1.3` s.
- `sweep N LO HI THREADS`: the least largest element in `[LO, HI]` of an admissible `N`-set, `N <= 7`, exhaustive: the difference set `D(H)` of the chosen elements as one signed bit array, the last two elements word-parallel, pruned by `g_3(r)`, the variance bound and the concentration cut of the note, which a binary search applies at every node owing `r >= 3` elements; prints per top the nodes by depth, the cut tests and the time. Runtime `342` s at `sweep 7 169 473 8`, which finds none, against `36` s for `probe 7 380 8` alone.
- `tally N LO HI THREADS PLAIN`: every admissible `N`-set with largest element in `[LO, HI]` through the cut search, each checked from the definition and against the cut on all its prefixes; with `PLAIN = 1` also through the search without the cut, comparing the lists. Runtime `0.3` s at `tally 6 168 176 8 1`, `12.4` s at `tally 7 474 474 8 0`.
- `four N HI THREADS`: `g_4(n)` of Erdos problem 817 for `n = 1..N`, the least `N` with an `n`-set in `[1, N]` whose subset sums hold no non-trivial 4-term progression, exhaustive: for each largest element in turn, a search from the top down where each node keeps the candidates `y` with `S U (S + y)` free of 4-term progressions, `S` the sums so far, children testing only their parent's candidates, pruned by `g_4(r)` of the run before; prints every extremal set, each rechecked from the definition. Runtime `2.0` s at `four 6 200 8`, `63` s at `four 7 130 8` on `8` threads.
- `fourtop G TOP ALL`: the same search at one largest element, `G` the comma-separated `0, g_4(1), ..., g_4(n - 1)`; `ALL = 0` stops at the first set. Runtime `283` s at `fourtop 0,1,3,5,14,40,79 180 0`.
- `fourplain N HI`: the least largest element of an `N`-set in `[1, HI]` with no 4-term progression in its subset sums, by a plain search from the bottom up that recomputes the subset sums and tests every pair as the first two terms at each node, no pruning but heredity. Runtime `38.4` s at `fourplain 6 100`.
- `route R E`: the least ratio per minimal base set in `[3, R]` over terms in `[10^E, 10^30]`. Runtime under `0.01` s.
- `tight K BITS HOLES`: `{3, 4, 7}`, where `sigma = 1`, at `k = 1..K` through the wide certificate at `W = 2^BITS` with a list of `2^HOLES` middle holes, default `25`: at each step with `T <= a_(n+1)` it sets `L = C_k + 2T - 1` and `B = max(8, ceil(ln(3L)/(ln 4 - 1/5)) + 1)`, checks the step on every term up to `4^B`, and past `4^B` leaves the step to Ellison's bound on `abs(2^x - 3^y)`, the tight certificate of the note; prints `F`, the count of positive non-sums up to `F`, the index and term of the close, `S_n`, `L`, `B`, the terms checked, the least slack, the margin `(B + 1)(ln 4 - 1/5) - ln(3L)` and the steps tried, or the term where the cap stops. Runtime `0.4` s at `tight 4 31`, `12` s and `4.9` GB at `tight 5 34 26`, which caps; `10.2` s at `tight 4 20`.
- `ones R SIZE`: the base sets of at most `SIZE` bases in `[3, R]` with `sigma = 1`, `gcd = 1` and no two bases powers of one integer, by a search over one common denominator; prints the count, those with two smooth bases `2^alpha 3^beta` and the pair the rule takes, the two least smooth bases, and those without, and `rho`, the least `abs(2^x - 3^y)/2^x` over `x` in Ellison's exceptional set and `y >= 1`. Runtime under `0.1` s at `ones 47 6`.
- `pair D K BITS HOLES`: one base set with `sigma = 1` and two smooth bases at `k = 1..K` through the wide certificate, the paired certificate of `research/erdos/notes/sigma-one.md`: at each step with `T <= a_(n+1)` it sets `L = C_k + 2T - 1` and `H = max(ceil((f - 1) L 2^19/7153), ((f - 1) L)^(1/theta))`, rounded up, `f` the larger base of the pair, and checks the step on every term up to `H`; prints the pair, `C_k`, `F`, the count of positive non-sums up to `F`, the index and term of the close, `S_n`, `L`, `H`, the terms checked, the least slack and the steps tried, or the cap, or a horizon past `10^30`. Runtime `2.8` s at `pair 3,4,10,19 5 34 26`, `6.0` s at `pair 3,5,6,21 5 35 26`.
- `pairs R SIZE K BITS HOLES THREADS`: `pair` on every set of `ones R SIZE` with two smooth bases, each up to its first open level, then per level the cells closed and the largest `F` and `H`. Runtime `89.2` s and `7.9` GB at `pairs 47 6 5 35 26 2`, `35.2` s at `pairs 47 6 5 33 25 3`, `4.1` s at `pairs 47 6 4 24 25 4`.
- `pcontrol R SIZE K`: every cell of `pairs` at `k = 1..K` and `W = 2^31` against a plain knapsack to `4 F + 4`, same `F` and same count of non-sums. Runtime `0.9` s at `pcontrol 47 6 2`.
- `near E F BITS`: the neighbour pairs of `e^a` among the powers of `f`, `a` up to `A(e, f) = ceil((81 ln 2 + 2520 ln e ln f)/ln e) + 1`, the bound that Laurent's Corollary 2 puts on `a` when `0 < abs(e^a - f^b) < 2^80`, in exact integer arithmetic on `64`-bit limbs; prints those with gap below `2^80`, lists those below `2^BITS`, and the largest smaller power and the largest ratio of smaller power to gap. Runtime under `0.01` s at `near 3 5 24`, `0.03` s at `near 31 46 0`.
- `lpair D K BITS HOLES`: `pair` with the two-log rule of `research/erdos/notes/sigma-one.md`, for any base set with `sigma = 1` and no two bases powers of one integer: at each step with `T <= a_(n+1)` it takes, over every pair of bases `e < f`, the largest smaller power `H` of a neighbour pair of `near` with gap below `(f - 1) L`, and checks the step on every term up to the least `H`. Runtime `0.3` s at `lpair 3,5,7,13 4 31 25`.
- `logs R SIZE K BITS HOLES THREADS ALL`: `lpair` on every set of `ones R SIZE` without two smooth bases, or on all with `ALL = 1`, each up to its first open level, then per level the cells closed and the largest `F` and `H`. Runtime `2.8` s and `1.6` GB at `logs 47 6 5 31 25 4 0`, `6.9` s and `2.2` GB at `logs 47 6 5 33 25 2 0`, `5.6` s at `logs 47 6 4 31 25 4 1`.
- `lcontrol R SIZE K ALL`: `pcontrol` for the cells of `logs`. Runtime `1.3` s at `lcontrol 47 6 2 0`, `9.0` s at `lcontrol 47 6 2 1`.
- `cantor E1,E2,... CAP`: the set `C(1, 2^E1, 2^E2, ...)` through its carry automaton, at most `4` exponents below `126`: the ternary digits of each multiplier, the states, the cyclic components and those with more edges than nodes, the radius bracket, the dimension, two first-return cycles at a state of the largest component with a prefix reaching it, checked by running the carries directly, the dimension of the Cantor set they generate, and the admissible prefixes of length `1..28`. Runtime under `0.01` s at `cantor 8` and `cantor 2,8`.
- `carry M1,M2,... CAP`: the same for integer multipliers. Runtime under `0.01` s at `carry 94`.
- `tuples R TOP THREADS CAP`: every `R`-subset of the even exponents in `[2, TOP]` through the automaton, capped at `CAP` states: the positive ones with radius and dimension, the countable nonzero ones, the count of sets equal to `{0}`, the capped ones, the most states. Runtime `128.8` s at `tuples 1 40 8 16777216`, which caps at `40`, `16.2` s at `tuples 2 40 8 16777216`, under `0.1` s at `tuples 3 124 8 16777216`.
- `deaths R TOP D THREADS`: every `R`-subset of the even exponents in `[2, TOP]` through the prefix test modulo `3^D`: the count with `C = {0}`, the sets that pass at length `D`, the deepest death and the histogram of death lengths. Runtime `19.1` s at `deaths 3 3000 79 8`, `0.7` s at `deaths 3 1000 79 8`, under `0.01` s at `deaths 2 40 79 8`.
- `classes R D THREADS`: every ordered `R`-tuple of exponent classes modulo `2 3^(D - 1)` through the prefix test at length `D`: the count that dies, the count that survives against the first-moment bound `2^((R + 1)(D - 1))`, and the lower density of `R`-tuples with `C = {0}`, truncated. Runtime `10.0` s at `classes 3 7 8`, `0.42` s at `classes 3 6 8`, `0.22` s at `classes 2 8 8`.

## RUN

```
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- census 10 5 32 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- census 12 5 32 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- cell 3,4,5 4 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- cell 3,6,9,10,12 2 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- control
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- set 3,5,6,9 1 20000
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- route 10 12
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- deep 10 6 31 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- deep 12 3 31 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,4,8,9 5 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,5,7,9 5 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,4,9,10 5 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,5,6,9 5 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,6,8,9,10 5 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,4,6 4 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,4,6 5 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,5,8,9 5 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,4,5 6 34
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,6,8,9,12 3 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcell 3,6,9,10,12 3 33
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- wcontrol 16
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- split 3,4,5 5 31333640736 64
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- chainless 27
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- forced 26
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- windowed 27
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- jumps 10 6 31 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- trace 3,4,5 3 26
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- windows 3,4,5 1 6 30
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- graham 2 5 3 26
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- refute 27
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- ternary 6 200 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- band 419 475 70 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- band 419 473 100 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- sweep 7 169 473 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 474 474 8 0
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 6 168 176 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 169 348 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 349 445 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 446 458 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 459 467 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 468 473 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tally 7 474 474 8 1
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tight 4 31
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tight 4 20
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tight 5 34 26
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- ones 47 6
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- pair 3,4,10,19 5 34 26
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- pairs 47 6 5 35 26 2
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- pairs 47 6 5 33 25 3
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- pairs 47 6 4 24 25 4
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- pcontrol 47 6 2
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- split 3,4,10,19 5 271729742861 64
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- set 3,4,7 2 20000000
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- split 3,4,7 3 166025260 64
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- split 3,4,7 4 14961204826 64
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- cantor 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- cantor 2,8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- cantor 2,18
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tuples 1 40 8 16777216
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tuples 2 40 8 16777216
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- tuples 3 124 8 16777216
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- deaths 2 40 79 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- deaths 3 3000 79 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- classes 3 7 8
cargo run --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml -- classes 2 8 8
```

- `cargo test --release --manifest-path research/erdos/lab/rs/mixed-powers/Cargo.toml`, `18` tests, under `1` s after the build.

## READS

- `30` minimal base sets in `[3, 10]`, all certified at `k = 1, 2, 3`, `21` at `k = 4`, `5` at `k = 5`; `103` in `[3, 12]`, all at `k = 1, 2` with `{3, 6, 9, 10, 12}` at `k = 2` from the cell verb, `F = 1473914231`.
- `F({3, 4, 5}, k)` at `k = 1..4`: `79, 77613, 4330731, 1075364603`, with `11, 1128, 45704, 1785062` positive non-sums up to `F`; `F({3, 4, 6}, k)` at `k = 1..3`: `986, 242113, 58941162`.
- The one base set in `[3, 12]` with `sigma = 1` and `gcd = 1` is `{3, 4, 7}`, outside the certificate.
- Through the wide certificate: all `30` minimal sets in `[3, 10]` at `k = 4`, `28` at `k = 5`, `{3, 5, 6, 7}` and `{4, 5, 6, 7, 9}` at `k = 6`; `F({3, 4, 5}, 5) = 31333640736` with `100335626` non-sums; `F({3, 4, 6}, 4) = 14059625243`; `101` of the `103` in `[3, 12]` at `k = 3`.
- From `k = 4` on every certified cell makes `F` by the pair route.
- `V` up to `2^27`: `37` terms, `7` windows at `n = 5, 10, ..., 35`, none chained, all holding a non-sum, least width `0.255` of the term; `33` chain values, all non-sums.
- The jump law holds on every step of the `173` cells `jumps` reaches at `W = 2^31`; `T_n/a_(n+1)` after the last base enters peaks at `1.9430` and passes `1.1727` in every certified cell with `k >= 2`; all `28` capped cells sit in one open run of rises.
- Route ratios over `[10^12, 10^30]` for the `30` minimal base sets in `[3, 10]`: between `1.105953` and `1.537469`.
- `F({3, 4, 7}, k)` at `k = 1..4` by the tight certificate: `581, 3982888, 166025260, 14961204826`, with `37, 5207, 391070, 23496512` positive non-sums up to `F`, `B = 8, 16, 19, 23`; `k = 5` caps at `W = 2^34` with `2^26` holes at term `40`.
- Smooth-pair sets, `sigma = 1`, at most six bases in `[3, 47]`: `52` of `62`, all closed at `k = 1..4`, `51` at `k = 5`, `{3, 4, 7}` capping at term `40` at `W = 2^35`; `F({3, 4, 10, 19}, k) = 251, 66091, 4138541, 1052160335, 271729742861` at `k = 1..5`; largest `H` `186717655760522`.
- `g_3(n) = 1, 3, 8, 22, 60, 168, 474` at `n = 1..7` for Erdos problem 817; `474` has one extremal set and `168` two.
- `C(1, 2^a, 2^b, 2^c) = {0}` for all `561375500` triples of even exponents up to `3000`, every prefix test dead by length `50`, deepest at `(2, 8, 332)`; up to `124` every triple automaton has at most `35` states.
- Pairs of even exponents up to `40`: `177` empty, `3` countable, `10` of positive dimension, the largest `0.228391` at `(2, 8)`.
- `dim C(1, 2^8) = 0.306871`, radius `1.400924990`; `dim C(1, 2^m)` in `[0.259662, 0.262985]` at `m = 16..38`.
- Classes of triples dead by length `D`: `0.691358, 0.818777, 0.894353, 0.937762, 0.963267` at `D = 3..7`.

## WITNESSES

- `research/erdos/claims/erdos.md`, every row citing this study:
- `route 10 12`: the two-part obstruction at twelve digits.
- `census 10 5 32 4`: cofinite at `k = 1, 2, 3` for every base set in `[3, 10]`.
- `census 12 5 32 4` and `cell 3,6,9,10,12 2 34`: cofinite at `k = 1, 2` for every base set in `[3, 12]`.
- `census 10 5 32 4`, `cell 3,4,5 4 34` and `control`: the largest non-sums of `{3, 4, 5}` and `{3, 4, 6}` and the knapsack agreement.
- `census 10 5 32 4` and `cell 3,4,5 4 34`: the Conjecture row on `{3, 4, 5}` at every level.
- `cell 3,5,6,9 1 20` and `set 3,5,6,9 1 20000`: sets and multisets differ on tie cells.
- `graham 2 5 3 26`: the seed of the incompleteness of `floor(2 (5/3)^n)`.
- `windows`: the power multisets run no window chain.
- `census 10 5 32 4`, `cell 3,4,5 4 34` and `cell 3,4,5 5 34`: the certificate for `{3, 4, 5}` at `k = 5` passes `2^34` bits.
- `refute 27`: the Refuted row on the surplus lemma L.
- `chainless 27` and its crate test: the Refuted row on the chain-excluded lemma.
- `forced 26`: the six forced-window variants of the note, none with a non-sum in its last `12` term intervals.
- `windowed 27` and its crate test: the Refuted row on the windowed lemma.
- `jumps 10 6 31 4`: the jump law on every step of the `173` cells it reaches, the route of `F` by its birth step, the peaks of `T_n/a_(n+1)` and the open runs at the cut.
- `deep 10 6 31 4`, the `wcell` runs and `wcontrol 16`: cofinite at `k = 1..4` for every base set in `[3, 10]`, `28` of `30` minimal sets at `k = 5`, and the route of every `F`.
- `deep 10 6 31 4` and `split 3,4,5 5 31333640736 64`: `F({3, 4, 5}, 5)`.
- `deep 12 3 31 4`: `101` of `103` minimal sets in `[3, 12]` at `k = 3`.
- `ternary 6 200 8`: `g_3(n)` for Erdos problem 817 at `n = 1..6`.
- `band 419 475 70 8` and the crate test on the set: `g_3(7) <= 474`.
- `band 419 473 100 8`: no set beating `474` with its six largest elements within `100` of the top.
- `sweep 7 169 473 8`: `g_3(7) = 474`.
- `tally 7 LO HI 8 1` over `[169, 348]`, `[349, 445]`, `[446, 458]`, `[459, 467]`, `[468, 473]`, `[474, 474]`: the search without the cut agrees, about `27` minutes in all.
- `tally 7 474 474 8 0` and `tally 6 168 176 8 1`: the extremal sets at `n = 6, 7`, and the cut removing none of the sets it was tried on, with the crate tests on the cut.
- `four 6 200 8`, `fourplain N 100` for `N = 1..6` and the crate test on the four-term values: `g_4(n) = 1, 3, 5, 14, 40, 79` at `n = 1..6`.
- `four 7 130 8` and the crate test on the four-term values: `131 <= g_4(7) <= 237`.
- `tight 4 31`, `tight 4 20`, `set 3,4,7 1 2000000`, `set 3,4,7 2 20000000`, the `split` runs at `3,4,7` and the crate test on the sigma-one set: `{3, 4, 7}` cofinite at `k = 1..4` and its largest non-sums.
- The certificate, exactness, monotone, wide certificate, two-part, window, two-route, one-inequality, surplus-does-not-heal, chain, birth, jump and route-of-`F` rows are proofs on `research/erdos/notes/problem-124.md`; the sigma-one deficit, tight certificate and sigma-one exactness rows, the smooth-pair reduction and the two-log reduction are proofs on `research/erdos/notes/sigma-one.md`.
- `ones 47 6`: the `62` base sets of at most six bases in `[3, 47]` with `sigma = 1`, `gcd = 1` and no two bases powers of one integer, the `52` with two smooth bases, and `rho = 7153/2^19`, in the smooth-pair reduction row.
- `pairs 47 6 5 35 26 2`, `pairs 47 6 5 33 25 3`, `pairs 47 6 4 24 25 4`, `pcontrol 47 6 2`, the `split` runs at `3,4,10,19`, `3,5,6,21` and `3,6,8,15,21,29`, and the crate test on the smooth pair: the `52` smooth-pair sets cofinite at `k <= 4`, `51` at `k = 5`, and the largest non-sums of `{3, 4, 10, 19}`.
- `near 3 5 24`, `near 3 7 0`, `near 5 7 0`, `near 31 46 0` and the crate test on the close powers: the neighbour pairs below `2^80` and `abs(3^a - 5^b) >= min(2^80, min(3^a, 5^b)/47)`.
- `logs 47 6 5 31 25 4 0`, `logs 47 6 5 33 25 2 0`, `logs 47 6 4 24 25 4 0`, `logs 47 6 4 31 25 4 1`, `lcontrol 47 6 2 0`, `lcontrol 47 6 2 1`, the `split` runs at `3,5,7,13`, `lpair 3,5,7,13 6 34 25` and the crate test on the two-log rule: the `10` sets without two smooth bases cofinite at `k = 1..5`, and the two rules agreeing on the `52`.
- The concentration-cut row is a proof on `research/erdos/notes/problems-817-406.md`.
- `deaths 3 3000 79 8`, `tuples 3 124 8 16777216` and the crate test on the two deciders: no four exponents within a span of `3000` for a nonzero 3-adic `lambda`.
- `tuples 2 40 8 16777216`, `deaths 2 40 79 8` and `cantor 2,18`: the pairs up to `40`.
- `cantor 8` and the crate test on the prefix counts: the Table 5.2 entry at `2^8`.
- `classes 3 3 8` to `classes 3 7 8` and the crate test on the first moment: the class densities.
- `tuples 1 40 8 16777216`: the Conjecture row on `dim C(1, 2^m)`.
- The prefix-test and almost-every-set-is-empty rows are proofs on `research/erdos/notes/problems-817-406.md`.
