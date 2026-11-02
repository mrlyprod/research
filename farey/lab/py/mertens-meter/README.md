# mertens-meter

- Weights the Farey stack by `mu(n)`, the node `a/b` carrying `sum_{k <= N/b} mu(kb)`, and counts how often that equals `M(floor(N/b))`.
- Checks the global readout `sum_{n <= x} M(floor(x/n)) = 1` at every `x` through 20000.
- Builds `M(x)/sqrt(x)` for `x = 1..50000` from a linear Mobius sieve, resamples it uniformly in `log x` on 8192 points, Hann-windows it and takes the real FFT power spectrum.
- Reads the axis as `gamma = 2 pi f`, keeps local maxima above three times the band median over `8 < gamma < 55`, and matches each of the first eight nontrivial zeta zeros to its nearest peak.

## RUN

`uv run python research/farey/lab/py/mertens-meter/mertens_meter.py`

Under a second. Domain is the source domain: `N_mertens = 50000`, `N_stack = 200`, `N_readout = 20000`, 8192 log samples.

## WITNESSES

- farey.md:128 to farey.md:135, the eight rows in order: `14.1347 | 13.94 | 0.20`, `21.0220 | 20.90 | 0.12`, `25.0109 | 24.97 | 0.04`, `30.4249 | 30.19 | 0.23`, `32.9351 | 32.52 | 0.42`, `37.5862 | 37.74 | 0.16`, `40.9187 | 40.64 | 0.27`, `43.3271 | 42.97 | 0.36`
- research/claims/:53 the readout `sum_{n <= N} M(floor(N/n)) = 1` holds at every `N` through 20000, zero breaches
- research/claims/:346 detected `13.94, 20.90, 24.97, 30.19, 32.52, 37.74, 40.64, 42.97`, errors `0.04` to `0.42`

## NOTE

- The bin width is `0.5806`, so every error above sits inside one bin: this is a rendering, not a measurement.
- farey.md:117 calls the node brightness `M(floor(N/b))` and farey.md:118 marks that **Proved**; the generator sums `mu(kb)`, and the two agree at only 64 of 200 denominators. Repeated at farey.md:233, research/claims/:53 and research/claims/:346.
