# eisenstein-visibility

- Recomputes the base-3 page constants at 60-digit working precision, printed to 43 digits: pi by Machin, trigamma `psi1(x) = sum 1/(n+x)^2` by Euler-Maclaurin in exact rationals (100 explicit terms, Bernoulli tail through `B_20`), then `L(2, chi_-3) = (psi1(1/3) - psi1(2/3))/9`, Catalan `G = (psi1(1/4) - psi1(3/4))/16`, `zeta(2) = pi^2/6`, `zeta_K(2) = zeta(2) L`, `zeta_Qi(2) = zeta(2) G` and both reciprocals.
- Cross-checks by the reflection identity `psi1(1/3) + psi1(2/3) = 4 pi^2/3`, prints the implied `Cl2(pi/3)`, and runs the two-term relation search as continued fractions of `L sqrt(3)/pi^2` and `L/pi^2` with `G/pi^2` and `zeta(2)/pi^2` as controls, reporting each last convergent denominator.
- Repeats the same constants by independent float routes: the paired Dirichlet series for the `L`-function, the Eisenstein lattice sum `(1/6) sum 1/N(z)^2` and the Gaussian `(1/4) sum 1/(a^2+b^2)^2` truncated at norm `4e6` with their analytic tails, and the `gcd(a,b) = 1` census of `[1,3000]^2`.
- Builds `Z[omega]` from `omega^2 = -1 - omega` (norm `a^2 - a b + b^2`, nearest-integer division by the conjugate, Euclidean gcd), sieves 400000 seeded random pairs for the coprime fraction, checks `base^level - base^(level-1)` against `phi(base^level)` for base 2, 3, 5 and the base 6 counterexample, and draws the hexagonal visibility figure by per-pixel coverage blending.

## RUN

`uv run python research/bases/lab/py/eisenstein-visibility/eisenstein_visibility.py` then `uv run python research/bases/lab/py/eisenstein-visibility/makefig.py`

About 3 s and 0.1 s. Domain is the source domain: prec 60, 200000 Dirichlet pairs, lattice norm `<= 4000000` at radius 4400, census `[1,3000]^2`, sieve 400000 pairs in `[-500,500]^2`, figure `880x760` over `[-30,30]^2`. `makefig.py` writes `bases-fig.png` to `data/research/bases/lab/py/eisenstein-visibility/`.

## WITNESSES

- bases.md:16 `6/pi^2 = 0.607927...`; bases.md:17 the `3000 x 3000` grid gives `0.608042`
- bases.md:51 forty digits asked, 43 printed: bases.md:59 `1.644934066848226`, bases.md:60 `0.781302412896486`, bases.md:61 `1.285190955484149`, bases.md:62 `0.778094489175179`, bases.md:63 `0.915965594177219`, bases.md:64 `0.663700804613853`
- bases.md:53-55 lattice sum `1.2851908043 + 1.51e-07 = 1.2851909555` against `zeta(2) L = 1.2851909555`, ten digits; the direct Dirichlet sum for bases.md:60 is `0.781302412896`
- bases.md:67-68 Eisenstein sieve `311237/400000 = 0.77809` on 400000 pairs, against the predicted `0.778094`
- bases.md:75-76 CF of `L sqrt(3)/pi^2` is `[0, 7, 3, 2, 2, 3, 2, 23, 2, 7, 1, 12, 1, 17, 1, 29, 1, 1, 2, 1, 6, 12]`, ordinary quotients, last convergent denominator `1805284405980`
- bases.md:77 control `zeta(2)/pi^2` has CF `[0, 5, 1, 10^59]`, convergent `1/6`
- bases.md:85 `Cl2(pi/3) = 1.0149416064...`
- bases.md:21 the figure `bases-fig`; the run also reports `364/612` in-frame points visible `= 0.5948`

## NOTE

- bases.md:76 states a bound, not a value: the printed last convergent denominator `1805284405980` is what clears `10^11`.
- bases.md:173 claims byte-for-byte redrawing; the regenerated figure is not: it differs from the pinned one by at most 12 levels on 1.1 percent of its pixels.
- The recorded transcript labels the reflection residual `~1e-40` where the value is `-1.1035169449471404E-42`; bases.md:174 records that already, so the label is dropped here and the identity printed plainly.
