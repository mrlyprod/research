# gasket-witness-weights

- Reparametrises the gasket residual `R(level)` by the witness `z` instead of by the multiplier pair `(s,t)`: every off-diagonal collinear pair of `G_level` is `(sz, tz)` with `gcd(s,t) = 1` and `z` a positive integer vector, uniquely.
- Checks the free-digit automaton `B(s,t)` is a constrained tensor square of a single-coordinate carry automaton `C(s,t)`: `T = S (x) S - U (x) U - V (x) V + W (x) W`, so `B` needs no four-tuple state graph.
- Checks the box construction `P_level(s,t) = #{z >= (1,1) : z_1 + z_2 <= (3^level-1)/(2 max(s,t)), sz, tz in G_level}` against `B(s,t)` and against the tensor iteration at multipliers where a forward build of `B` blows up.
- Checks the witness weight `w = z_1 + z_2` is never below 4, that the weight layers scale exactly as `R_{3w}(level) = R_w(level-1)`, and that the weight-four layer is counted by the no-adjacent-ones set `F_level` with `|F_level| = Fib(level+1) - 1`.
- Checks every multiplier pair above `(3^level-1)/10` contributes exactly 4 ordered collinear pairs.
- Checks the golden ceiling `M_level(z) <= Fib(level+1) - 1` on every coprime direction in a box, at every level to 40, with `z = (1,3)` the only attainer, and on six families chosen to favour a breach at every level to 45.
- Carries the second moment `E(level)` and the residual `R(level)` to level 17 by grouping all `3^level` points by primitive direction, and checks both `R/3^level` and `R/phi^(2 level)` fall at every level through 17.
- Rebuilds the ray automaton in the direction coordinate: a multiplier word is a word over the increments `{0, z_2, -z_1}` summing to zero, so `M_level(z) + 1` counts the closed paths of a carry automaton whose states live in `[-z_1/2, z_2/2]`.
- Checks the branch structure that turns the golden ceiling into a theorem: out-degree at most two, branch states in one residue class mod 3, and the two successors of a branch state differing by `q/3` where `q` is the unique member of `{z_1, z_2, w}` divisible by 3.
- Checks the proved cases over the box and over the de-duplicated union of the six adversarial families, names the directions left over and runs them to level 60, and verifies nine explicit Fibonacci certificates in exact integer arithmetic.
- Checks the renewal criterion `Sum_{j>=2} f_j Fib(level+1-j) <= Fib(level-1)` on the first-return counts of every occupied direction of the box.
- Solves the golden potential in exact `Q(sqrt5)` arithmetic: `u(0) = 1` and `phi u(c) = sum of u over the successors of c` for every live `c != 0`, read off as `U(z) = sum of u over the successors of the start state other than itself`.
- Checks the criterion `U(z) <= phi^-2`, which implies `M_level(z) <= Fib(level+1) - 1` at every level, on the box, the six adversarial families and a stress list of high `v_3` directions, and confirms every solution against both inequalities of the criterion rather than only against the linear system.
- Checks the occupancy residue law: writing `q = 3^k q1` for the coordinate divisible by 3 and `p` for the other, a direction carrying mass at any level has `q1 = p mod 3`, so the congruence alone empties a large share of the directions with `3 | z1 z2`.
- Checks the short first returns: no first return has length between 2 and `v_3(q)`; `f_2` is nonzero only at `{1,3}` and `f_3` only at `{1,9}`, `{1,12}`, `{3,10}`, `{4,9}`, each equal to 1.
- Checks the degree potential: `pi(c) = 1` where a live state has out-degree two, `phi^-1` where it has out-degree one, `pi(0) = 1`, is a super-solution of the golden criterion, and sweeping it under the same operator gives a decreasing chain of exact upper bounds on `U`; the sweep runs over consecutive depths, so the depth reported is the least one that works.
- Checks the burst identity that blocks every valuation-graded potential: `U` is `phi^-(k-1)` times the sum of `u` over the `2^(k-1)` burst-floor states, all of valuation 0, while `u(p) = phi^-1` at the valuation-0 state `p`.
- Checks the golden partition bound at `v_3(q) = 1`: `U <= phi^-1 (1 - phi^-max(t,2))` with `t = v_3(q1 - p)`, hence `U <= phi^-2` on the arithmetic class `t <= 2`.
- Checks the two automaton-free restatements of the criterion: `Sum_level (M_level + 1) phi^-level <= phi^4`, and `Sum_m phi^-l(m) <= phi` over the multipliers `m` of the direction, where `l(m)` is the number of base-3 digits of `(z1 + z2) m`.
- Checks the route that dies: the state maximum `G(level) = max_c N(c,level)` does not obey `G(level) <= G(level-1) + G(level-2)`.

## RUN

- `uv run --with numpy python research/coprime/lab/py/gasket-witness-weights/witness.py`
- About 165 seconds on one core, peak under 3GB at level 17; prints one line per law and exits nonzero on the first failure.
- numpy is used only for the level 15 to 17 array sort; every other law is stdlib.

## WITNESSES

- 473 coprime pairs below 40: the tensor square reproduces the `B(s,t)` return counts at every level to 9, zero mismatches.
- Carry states against reachable `B` states: `(365,1094)` 729 against 26931, `(41,122)` 81 against 835, `(25,52)` 38 against 393, `(31,40)` 35 against 354.
- 812 coprime pairs: the box construction agrees with `B(s,t)` at level 9, zero mismatches; it agrees with the tensor iteration at `(365,1094)`, `(41,122)`, `(122,123)`, `(1,2460)` and `(2431,2458)` at level 9 and level 12.
- `R(level) = 20, 88, 432, 1624, 5512, 15896, 46064, 124928, 335704, 863848` for level 4..13.
- `R_{3w}(level) = R_w(level-1)` on all 1869 weight layers divisible by 3 across level 5..13, zero failures.
- No witness of weight below 4 at any level to 13.
- `R_4(level) = 12, 36, 108, 336, 988, 2596, 6672, 17480, 45720` for level 4..12, each equal to twice the number of ordered coprime non-3-power-ratio pairs drawn from `F_level`, and `|F_level| = Fib(level+1) - 1` at each of those levels.
- Pairs above `(3^level-1)/10` number 18, 57, 163, 402, 1019, 2702, 7060, 18607 at level 6..13; every one contributes exactly 4.
- `M_level(1,3) = Fib(level+1) - 1` for `level <= 40`; over 13158 coprime directions with `z_1 <= 120` and `z_1 <= z_2 <= 240` there is no `level <= 40` and no `z` with `M_level(z) > Fib(level+1) - 1`, and `(1,3)` is the sole attainer at level 40.
- Six adversarial families at every `level <= 45`, coprime counts 16940, 253, 2998, 1499, 1199, 1199, total 24088, zero breaches: binary base-3 pairs below `3^8`, no-adjacent-ones pairs below `3^7`, `(1,t)` with `t < 3000`, consecutive below 1500, `(s,3s-1)` and `(s,3s+1)` with `s < 1200`. The shelf script `gasket-ray-machine/scripts/verify.py` ships the same six with the binary family at `3^7`, 11369 directions, to stay inside its time budget.
- `E(level) = 4003372, 11679626, 34050692, 99800950, 292848756` and `R(level) = 863848, 2211960, 5549452, 14100688, 35354824` for level 13..17.
- `R(level)/3^level` peaks at `0.8401158` at level 8 and falls at every level to `0.2737709` at level 17; `R(level)/phi^(2 level)` peaks at `3.2378233` at level 12 and falls at every level to `2.7724831`; the level ratio `R(level+1)/R(level)` reads `2.5073119` at level 17, below `phi^2 = 2.6180339`.
- The direction automaton reproduces the gasket-digit automaton on all 947 coprime pairs below 40 at every level to 20, zero mismatches.
- Of the 13158 coprime `z` with `z_1 <= 120` and `z_1 <= z_2 <= 240`, exactly 218 have a live automaton beyond the start state; every one has out-degree at most two, one branch class, and carries inside `[-z_1/2, z_2/2]`, the largest live set 37 states. The branch argument settles 206 of them, 107 of which have `v_3(q) = 1`. The twelve left are `(1,9)`, `(1,27)`, `(1,81)`, `(1,90)`, `(4,117)`, `(9,73)`, `(9,82)`, `(9,235)`, `(10,81)`, `(13,108)`, `(27,217)`, `(27,226)`.
- The first three of those twelve are shift rays, closed by `Fib(p+2) Fib(q+2) = Fib(p+q+3) - Fib(p+1) Fib(q+1)`; the other nine carry Fibonacci certificates of denominators 18, 40, 381, 18, 2013, 18, 40, 2013, 34.
- The six adversarial families overlap: their 24088 coprime members are 23435 distinct directions, of which 717 lie in the box and 22718 are new. Of the 22718: 20945 have zero mass, 1693 fall to the branch argument, 3 are shift rays, and 77 are left to the enumeration. Those 77 hold at every level to 60, zero breaches, worst ratio to the ceiling below `0.1516`.
- The renewal criterion holds on all 218 occupied directions of the box at every level to 46; `f_1 = 1` everywhere and `f_2 = 1` only at `(1,3)`.
- `M_level(z) <= D_level(w)` on all 829 coprime directions with `z_1 <= 30`, `z_1 <= z_2 <= 60` at level 12; the bound is far weaker than the ceiling, `D_24(w) = 4196351, 1683971, 613817, 228519` at `w = 4, 10, 28, 82` against `Fib(25) - 1 = 75024`.
- The state maximum at `(1,9)` runs `1, 1, 1, 2, 4, 6, 9`, so `G(4) = 4 > G(3) + G(2) = 3`; the recursion fails for 8 directions of the box.
- The golden potential is nonnegative and satisfies both criterion inequalities on all 218 occupied directions of the box; 214 have `U <= phi^-2` and the four failures are the shift rays `(1,3)`, `(1,9)`, `(1,27)`, `(1,81)`, each with `U = phi - 1` exactly.
- `U` takes 57 distinct values over the box, in order `phi - 1` on the four shift rays, `2 - phi` at `(1,12)`, `(3,10)`, `(4,9)`, then `(9 phi - 14)/2` at `(1,90)`, `(9,82)`, `(10,81)`.
- Over 36037 directions - the box, the six families and a stress list of high `v_3(q)` directions - 1995 are occupied and 1987 satisfy `U <= phi^-2`; the 8 failures are exactly the shift rays `(1,3^j)` for `j = 1..8`, no direction has `U` in the open interval `(2 - phi, phi - 1)`, the maximum among the passing directions is `2 - phi` attained only at `(1,12)`, `(3,10)`, `(4,9)`, and the largest live set is 256 states.
- Of the 13158 coprime box directions, 6566 have `3` dividing neither coordinate and carry no mass at any level; 6374 more have `3` dividing a coordinate and still carry none; 218 are occupied.
- Adversarial cross-check on all 218: `f_1 = 1`, the exact partial sum `Sum_{j=2}^{46} f_j phi^-j` from the independent first-return enumeration is dominated by `phi^-1 U` from the linear solve, and `U <= phi^-1` at every direction, with equality only on the shift rays.
- Of the 36037 directions, 20193 have `3 | z1 z2` and 15556 of those match `q1 = p mod 3`; all 1995 occupied directions match, zero violations, so the congruence empties 4637 directions with no automaton built.
- No occupied direction has a first return of length between 2 and `v_3(q)`, zero failures over 1995.
- The burst identity `U = phi^-(k-1) Sum_m u(q1 m)` over the `2^(k-1)` burst-floor states, and `u(p) = phi^-1` whenever `2p <= q`, hold exactly on all 1995.
- The degree potential is a super-solution on 1968 of the 1995 occupied directions, 1902 of which have no double branching and 66 of which do, so the tool reaches strictly past the branch case.
- Sweeping it settles `U <= phi^-2` on 1966 directions - least depth 1 on 1804, 3 on 101, 4 on 44, 5 on 11, 6 on 6 - and leaves 29: the 8 shift rays `(1,3^j)` and `(1,756)`, `(1,2196)`, `(1,2214)`, `(1,2268)`, `(1,2430)`, `(9,2188)`, `(10,2187)`, `(13,1080)`, `(13,3267)`, `(27,730)`, `(27,2188)`, `(28,729)`, `(28,2187)`, `(40,1053)`, `(81,2188)`, `(82,2187)`, `(91,2214)`, `(121,3159)`, `(243,2188)`, `(244,2187)`, `(819,2539)`.
- 757 occupied directions have `v_3(q) = 1`; the bound `U <= phi^-1 (1 - phi^-max(t,2))` holds at every one, is attained exactly at `(1,12)` and `(3,10)`, and 512 of them fall in the class `t <= 2` where it gives `U <= phi^-2` with zero failures.
- The two automaton-free restatements agree with the linear solve on 111 directions with `z1 < 40`, `z2 < 120`, checked in exact `Q(sqrt5)` arithmetic.
