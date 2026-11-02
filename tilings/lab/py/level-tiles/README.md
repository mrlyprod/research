# level-tiles

- Decides, for a base `b` and a digit set `F`, which levels `A_n = F + b F + ... + b^(n-1) F` tile the integers and which are spectral, from the zero set of the digit transform `hat F` at roots of unity.
- The zero set `Z(F) = {m : Phi_m | F(x)}` is computed exactly by integer division by `Phi_m` for every `m` with `phi(m) <= deg F`; the level zero set comes from the index lemma, never from a factorisation of the level.
- T1 and T2 of Coven-Meyerowitz are read on the level exponents; T1 at every level is decided by the closed condition of `notes/tiling-levels.md`, T2 at every level by the levels up to the horizon `1 + G(2 omega - 1)`.
- Independent checks: PARI `gp -q` factors the masks; an exhaustive exact-cover search looks for a complement in `Z/N`; the Coven-Meyerowitz complement is built and `A + B` is checked to be all residues; a clique search decides rational spectra in `Z/L`, `L` the lcm of the level zero set, and PARI certifies when no other spectrum can exist.
- Exact integer arithmetic throughout; floats only in the final re-check of a found spectrum.

## RUN

```
uv run python research/tilings/lab/py/level-tiles/level_tiles.py
uv run python research/tilings/lab/py/level-tiles/level_tiles.py census measure
```

- No argument runs every verb, about 90 seconds on one core. Needs PARI (`gp`).
- `hand` (1 s): the level `{0, 2, 6, 8}` of `{0, 2}` at base 3, its mask `Phi_4^2 Phi_12` from PARI and from the index lemma, T1 read `4` against `2`, and no complement in `Z/N` for `4 | N <= 128`.
- `index` (13 s): 752 masks, every digit set with `0` at bases 2..8, levels 1..3 and level 4 at bases <= 4, factored by PARI; the cyclotomic part equals the index lemma with multiplicity on all of them.
- `census` (8 s): every digit set with `0` and at least two digits at bases 2..16 (65519 sets); per base the counts tiling `Z/b`, tiling at every level, passing T1 at every level, the counterexamples to the guess raw and after dividing by the gcd, the T2 horizon and the depth histogram; asserts on every set the T1 depth formula against the T1 depth read level by level to level 16, the tiling depth against the least of the T1 depth, the T2 depth and 8 (the T1 depth at prime-power sizes), and the prime-power-base law.
- `search` (35 to 50 s): 494 levels (bases <= 8, levels <= 2) against an exhaustive complement search over `|A| | N <= 512`, 492 decided and 2 past the step cap of 2000000; 803 tile verdicts (bases <= 12, levels <= 3) certified by the Coven-Meyerowitz complement.
- `witness` (11 s): the five named digit sets level by level with PARI factorisation, exponents, verdict and complement; every multiple `vF`, `v = 1..12`, of `{0, 1, 8, 9}` and `{0, 3, 8, 11}` searched for a complement in `Z/12`, with the exponents of `2` in `S_vF`; three T1 depths against the tiling depth read to level 10; then the T2 search over every base to 30 and every size with two primes, at most 200000 digit sets a cell.
- `spectral` (4 s): 6866 levels with `|A| <= 36` at bases <= 12, levels <= 3, against the tiling verdict; the clique search decides rational spectra, and PARI marks the 77 digit sets with a root on the unit circle that is not a root of unity, whose 117 levels without a rational spectrum stay undecided.
- `measure` (8 s): 576422 product-form digit sets `{0..N-1} + m {0..L-1}` at base `p <= 24`, `N, L <= 12`, `N <= m <= p^2`, the every-level verdict against the spectral condition of Liu-Wang-Zheng Theorem 1.3; and consecutive digits at bases to 64 against `N | p`.

## WITNESSES

- `notes/tiling-levels.md` "The index lemma": 752 of 752 masks (`index`).
- "The hand check": the Cantor level (`hand`).
- "T1 at every level": the asserts on 65519 sets, the depth histogram, and the depths `2` of `{0, 1, 4, 5}` at base 6 and `3` of `{0, 1, 8, 9}` at base 10 (`census`, `witness`).
- "The guess is false": `{0, 2}` at base 6, `{0, 1, 8, 9}` and `{0, 3, 8, 11}` at base 12, the only two counterexamples after the gcd at bases <= 16, both mirror-symmetric (`census`); no multiple of either tiles `Z/12` (`witness`).
- "T2 at every level": `{0, 1, 2, 6, 7, 8}` at base 18, the counts 0, 16, 0, 4, 258 at bases 12, 18, 20, 24, 30, and the horizon checked six levels past itself (`witness`, `census`).
- `notes/tilings.md` "The census": the totals 65519, 588, 609, 641 and 32 (`census`); 492 of 492 and 803 of 803 (`search`).
- `notes/limit-measure.md` "Spectral levels": 6749 decided levels, spectral iff tiling on all, 117 undecided (`spectral`).
- "The limit measure": 576422 of 576422 and 2016 of 2016 (`measure`); the T1 depth `8` of `{0, 1} + 256{0, 1}` at base 18 (`witness`).
