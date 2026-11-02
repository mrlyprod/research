# band-return-times

- The band automaton of a direction `(z1, z2)` with `z1 + z2 = w` and `3 | z1`: states are the integers `j` in `[-(z2-1)//2, (z1-1)//2]`, the moves are `j -> (j + a)/3` over the increments `a` in `{0, z1, -z2}` whose quotient is integral, and the band is invariant.
- Out-degree is 2 on `j = 0 mod 3`, 1 on `j = z2 mod 3` and 0 on the third class, so the mean out-degree is `1 + (n0 - n2)/N` for the counts of those two classes among the `N` band states, and never 1 as an identity.
- `N` consecutive integers split the three classes mod 3 to within one, so `|mean - 1| <= 1/N` always and the mean is exactly 1 whenever `n0 = n2`, which `3 | N` gives and does not exhaust: `3 | N` is sufficient and not necessary.
- The census the verb prints over `w` in `{13, 40, 100, 101, 121, 257, 364, 1093}`: 591 coprime directions, 465 of them exactly critical, 60 of those 465 with `3` not dividing `N`.
- The two smallest witnesses, also printed by the verb: `(z1, z2) = (3, 1)` has states `{0, 1}`, degrees `2, 1` and mean `3/2`, the smallest non-critical direction; `(3, 2)` has degrees `2, 0` and mean 1 at `N = 2` with `3` not dividing `N`.
- A walk leaves 0 by the forced increment `z1` and a return is a walk back to 0; a return of length `n` spells a multiplier `m` with `m z1` and `m z2` binary in base 3 on disjoint supports, so `m w` is binary of length `n`.
- The first return time `d(z)` is the base-3 length of the shortest binary lift `m w`, and for `w = R_k = (3^k - 1)/2` the depth `ceil(d(z) / k)` is the block depth of the pincer.
- At the horizon `n = bk`, for a block depth `b`, the slot profile of the column transfer is uniform, so the transfer is one fixed matrix in `k` per residue `j` and the return count `L(k, bk)` obeys a constant-coefficient linear recurrence in `k` whose dominant root is the block rate `lam_b`.
- Every column sum of every block deviates from `2^b/3` by one of two values fixed by the parity of `b`, `-1/3` or `+2/3` at even `b` and `-2/3` or `+1/3` at odd `b`, so `lam_b` lies in `2^b/3 + [-1/3, 2/3]` at even `b` and in `2^b/3 + [-2/3, 1/3]` at odd `b`.
- The lift half of the same pincer is the submask count `N_K(m)`, the submasks of a binary `K` divisible by a divisor `m` of it; its solution set is closed under complement in `K`, under disjoint union and under nested difference, so it is always even and counts the distinct unions of pairwise disjoint irreducible submasks.
- The decomposition into irreducibles is not unique, so `N_K(m) <= #packings <= 2^(iota)` for the irreducible count `iota`, the right equality holding exactly when the irreducibles are pairwise disjoint; a power-of-two count does not force disjointness and the two ladders differ.
- A column transfer with a state set free of `k` is a linear representation of `N` as a series over the column word, so its state count is at least the Hankel rank; that rank reads `3, 7, 14, 31, 62, 126, 253` against the full `3, 7, 15, 31, 63, 127, 255` at word length `1..7` on each side, so no such transfer with fewer than 253 states exists at `k <= 15` and its unboundedness past there is observed and not proved.
- The reversed reading of that matrix is its transpose at equal Hankel side lengths `p = q` and is no second check; the irreducible count `iota` has the same full-rank Hankel, so the bounded object is not the computable one.
- The headline is two-sided at every `b`, `|3 lam_b / 2^b - 1| <= 2^(1 - b)`: the parity refines which edge is which, not the rate; the printed excess `3 lam_b - 2^b` is positive at every `b` the verb reaches, so the upper edge is the live one.

## WHAT IT COMPUTES

- `automaton`: the state count, the edge count, the out-degree profile and the mean out-degree of the band automaton over the coprime directions of the given weights, with the band invariance and the identity `mean = 1 + (n0 - n2) / N` asserted, the count of exactly critical directions and of those among them with `3` not dividing `N`, and the two smallest witnesses `(3, 1)` and `(3, 2)`.
- `returns`: `L(k, n)`, the number of primitive returns of the weight `R_k` inside horizon `n`, by a transfer matrix on the column carry; the per-block growth `lam_b` of `L(k, bk)`; the excess `rho` over the free model `2^(n-1) / R_k`; and the lengths at which no return exists.
- `hist`: the first-return histogram of every coprime direction of weight `R_k`, its bulk against the multiplier model `Sum_{max T = t} 2^k / m_T`, the depth histogram and the survival of the deep tail, and the deepest return against `sqrt(R_k)`.
- `critical`: the deepest and the median first return over a ladder of weights prime to 3, against the scale `sqrt(w)` a critical walk on a band of `w/2` states predicts, with the leave-one-out range of the `d_max` slope, the median exponent under cuts on the occupied count `Z`, the two-predictor fit on `log w` and `log Z`, and the least-squares diagnostics.
- `model`: `D(k, N) = Sum over the primitive lifts of length at most N of 2^L / m` at lift length `L`, the equidistribution count of return pairs, summed by the same transfer matrix and bracketed by the length of the lift.
- `ladder`: the fixed matrix at `n = bk`, its head lengths `r0(j)`, the exact characteristic polynomial of `L(k, bk)` in `k`, the minimal polynomial of `lam_b` with its radical form where the degree allows, `lam_b` as a certified interval by exact bisection, the parity column-sum bracket that holds it, and `max_j rho(B(b, j))` against `lam_b`.
- `lift`: the depth-2 lift census, `M_k = Sum_T N_K(m)` by meet in the middle with the Hankel rank of the sequence itself, the irreducible supply `Sum_T iota_T` and the largest one, the share of `T` whose count is trivial, whose count is a power of two and whose irreducibles are pairwise disjoint, the residue classes `m_T mod 9` of the ladder with the two-line argument that class 4 is empty, and the Hankel rank of `N` and of `iota` over the column word.
- `check`: the transfer against a brute enumeration of binary multiples, the block-ladder identities, and the breadth-first return time against the shortest lift found by brute force.

## HOW TO RUN

- `uv run python research/coprime/lab/py/band-return-times/returns.py check`
- `uv run python research/coprime/lab/py/band-return-times/returns.py automaton --weights 13 40 100 101 121 257 364 1093`
- `uv run python research/coprime/lab/py/band-return-times/returns.py returns --kmin 3 --bmax 9 --kbig 160`
- `uv run python research/coprime/lab/py/band-return-times/returns.py hist --kmax 14`
- `uv run python research/coprime/lab/py/band-return-times/returns.py critical --wmin 2000 --wmax 400000 --step 1.15`
- `uv run python research/coprime/lab/py/band-return-times/returns.py model --kmin 8 --kmax 16 --bmax 8`
- `uv run python research/coprime/lab/py/band-return-times/returns.py ladder --bmax 14`
- `uv run python research/coprime/lab/py/band-return-times/returns.py lift --kmax 12 --kbig 17 --nmax 7`
- `lift` costs `2^(1.5k)` on the `M_k` ladder and `2^(2k-1)` on the irreducible census, 170 s at the defaults, 80 s of it the Hankel rank at word length 7.
- `hist` costs `w^1.5`, 43 s at `k = 14`; `ladder` costs a second a step to `b = 14` and a minute to `b = 20`; every other verb is under 30 s.

## WITNESSES

- the-window-occupancy-and-the-band.md, THE WINDOW AT DIMENSION ONE: the lift half as a union count over irreducible submasks, and the return count of the critical band automaton, the block ladder at every depth, the block rate `lam_b` as an exact algebraic number, the return-time support and the depth reading of the deep tail.
