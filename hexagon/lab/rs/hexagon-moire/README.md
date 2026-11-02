# Hexagon Moire

- The stack of diagonal cube slices, one hexagon per odd side n, for the carpet, net, tree and void designs.
- Cut ink of every layer as an exact rational against the four closed forms, and the void centre cell.
- Exact area-weighted Pearson correlations of layer pairs on the full hexagon: the doubling sign law, its Richardson limit, adjacent and gcd-echo pairs, the other families.
- The pair limit in closed form: the cut cell law s_y = s_x + s_z + w, the phase map n = a m + c, and the exact rational correlation of every family against every phase map, with the finite-layer gap beside it.
- The ideal-frame stack at 3601 samples per axis: quarter-line bands, the one-sided plateau at X = 1/4, the void star.
- The twisted triangle averages and the coarse crosshair law at 200003 samples.
- The ghost star of the rendered carpet stack on a 1200 by 2399 raster, and its log-corrected decay in the ideal and the lattice frame.
- The cell frame of the ghost star: the arm's rational ink law to n = 2001, the decay ladders at L = 0 mod 4, L = 2 mod 4 and odd L out to L = 6400, and the band half-widths with their held-out predictions.
- The width family in closed form: the derived tail weight E(K) against its 8-periodic run, the per-layer identity in exact rationals at 14 distinct half-widths broken out by n mod 8, and the summed ladders with the constant, both 1/L^2 branches and the sliding window's own L residue trap.
- The constants from the `mrly_num` series and the partial character sums at N = 53 and 55.
- The one-layer law: the six ink-law coefficients of each family, the summed identity against the counted hexagons in exact rationals, and the approach in both classes of the layer count.
- The flat carpet stack: pair-term, carpet and void ink at generic points of 20 primitive lines against the ray law, the crosshair period sums, the solid diagonals, the gcd sums by their Jordan grouping against the double sum, the odd-harmonic Euler sums, the four exact variance constants at N = 999999 against their limits, and the odd-denominator fraction census.
- The eternal points of the xor pair: a rational grid of the cut searched at every odd n to 55 and to 2001, with the permutations of (3/4, 3/4, 0) beside.
- Every domain is the source's own; nothing was shrunk. About a minute.

## RUN

- `cargo run --release --manifest-path research/hexagon/lab/rs/hexagon-moire/Cargo.toml` from the repo root

## WITNESSES

- hexagon.md the cut ink of one layer: the ink laws match 28/28 layers in exact rationals for all four designs at every odd n <= 55, and the void centre cell is ink in 28 of 28 against 14 for the other three.
- hexagon.md layer pairs and the doubling law: the sign law holds on all 18 pairs from (3,5) to (601,1201); the full-hexagon blind pairs read (5,9) = -0.14179450 and (5,7) = -0.08542646.
- hexagon.md layer pairs and the doubling law, the magnitude: Richardson r_inf from -0.11715991 to -0.11711630, branches 0.1171270 and 0.1171274, against 253/2160 = 0.11712963 and 19/162 = 0.11728395.
- hexagon.md the other pair families: tree and void doubling at (201,401) +0.00049938 and +0.00107260; tree adjacent -0.07016649, void adjacent +0.00786627; echo (67,201) carpet +0.21473346, tree +0.14984772, void +0.07721358.
- hexagon.md the quarter line and the crosshairs: plateau 0.12210, 0.12340, 0.12412, 0.12446 at N = 151, 301, 601, 1201; the X = 1/4 bands +0.066465 and -0.090881 at N = 55.
- hexagon.md the quarter line and the crosshairs, the coarse model: A = 1/3 excess -0.083726 against -0.083333 at N = 5555; A = 1/7 -0.058065 at N = 55 against -0.036036 at N = 5555.
- hexagon.md the void star: arms 0.513667, 0.486492, 0.491844, 0.491525 and the Z arm 0.399416 over a mean 0.278670.
- hexagon.md the twist that kills the ray family: the twisted average falls like 1/N at every x, the A = C excess -0.003580 at N = 55 and -0.000052 at N = 5555.
- hexagon.md the ghost star: star minus background -0.09385, -0.03102, -0.01756 at 5, 28, 56 layers; excess times L runs -0.7779 to -1.2519 in the ideal frame and -1.0212 to -1.3645 in the lattice frame from L = 28 to L = 400, ln L slopes -0.1807 and -0.1911 in the ideal frame and -0.1252 and -0.1233 in the lattice frame over the windows 100 to 200 and 200 to 400.
- hexagon.md the ghost star, the cell frame: the arm ink law 1/2 + chi_8(n)/(2n) matches 1001/1001 layers at odd n <= 2001 and the three arms agree in 111/111 at odd n <= 221; excess times L + (ln L)/4 runs -0.2939128437 at L = 28 to -0.2937605886 at L = 6400 against the closed form -0.2937605857, ln L slope -0.250000.
- hexagon.md the ghost star, the parity split: the residual times L^2 runs -0.11975831 to -0.11979188 at L = 100 to 6400 with ratio 4.00 per doubling against -23/192 = -0.11979167; the L = 2 mod 4 ladder reads +0.13020825 at L = 3202 against +25/192 = +0.13020833; the odd ladder misses by +1/8 and what is left times L reads -0.250708, -0.250244, +0.249934, +0.249984, -0.250004 at L = 27, 99, 401, 1601, 6399.
- hexagon.md the ghost star, the width family: the twelve swept half-widths read -1/4, -1/6, -1/10, -3/28, -5/36, -3/22, -3/26, -9/68, -5/42, -13/100, -17/132 and -33/260 at L = 800, each the derived rational -(K + b)/(4(2K + 1)) and each found again by a blind best-rational search.
- hexagon.md the ghost star, the width family checked: W = 0 to 12 recomputed to L = 3200 land within 5.1e-8 with every gap falling by four per doubling; the held-out W = 14, 18, 28, 36 land within 4.0e-8 of -7/60, -5/38, -7/58 and -9/74 at L = 1600; W = 1 and W = 3 read -0.249999 and -0.166666, the same point set and the same rational as W = 0 and W = 2.
- hexagon.md the constants: pi/4 + pi^2/32 = 1.0938233009 with 1.0826656664 at N = 55; (pi + pi^2)/16 = 0.8131998159 against pi^2 ln2/(7 zeta3) = 0.8130217042; G/8 = 0.1144956993 with 0.1144757884 at N = 55 and G/8 - 1/8 with -0.0104828892 at N = 53.
- hexagon.md the constants: L(2, chi_-3) = 0.7813024129 printed from its series and matched by nothing above.
- hexagon.md the one-layer law: M(mean - A - c eps) = -B S - d s1 + e s3 - f s2 holds against the counted hexagons in exact rationals at every layer count to N = 55, all four families, 7 of 7 in each class n = 1, 3, 5, 7 mod 8.
- hexagon.md the one-layer law, the limits: carpet G/8 and G/8 - 1/8, net -G/8 and 1/8 - G/8, tree pi/48 + pi^2/48 + G/6 and void (pi + pi^2)/16 at both classes of the layer count.
- hexagon.md the one-layer law, the approach: at M = 3200 and 3201 the gap times M^2 reads -0.01562500 and +0.01562500 against -+1/64, and the gap times M reads -0.18750000 and -0.06250000 for the void against -3/16 and -1/16, -0.06250651 and -0.02082683 for the tree against -1/16 - 1/(48M) and -1/48 + 1/(48M), and -0.31249999 for the carpet split against -5/16.
- hexagon.md the width family: excess_W(n) = kappa chi + (m + q chi_8(n))/n + chi/(8 n^2) matches the counted band in exact rationals 1354 of 1354, at 14 distinct half-widths from W = 0 to 64 and every odd n from K to 201, with the four classes n = 1, 3, 5, 7 mod 8 counted apart; the nine odd rows repeat their even neighbours and are not counted twice.
- hexagon.md the width family summed: C_W = m(ln2 + gamma/2) + q L(1, chi_8) - G/8 + Delta_W holds to 1e-9 at L = 1600 for seven widths, W = 0 returning -0.2937605857 and Delta_W = 2/7 at W = 6.
- hexagon.md the width family, the two 1/L^2 branches: residual times L^2 reads -0.1197915 and +0.1302082 at W = 0 and +0.0027574 and +0.0174632 at W = 16, against -q/4 + 1/64 + m/48 and +q/4 + 1/64 + m/48.
- hexagon.md the width family, the L residue trap: the sliding window from L/2 to L reads -0.24999980 at L = 1600 and -0.43078703 at L = 1602, against m/2 = -1/4 and the limit m/2 + kappa/ln 2 = -0.43033688.
- hexagon.md the width family, the derived tail weight: E(K) = #(|j| <= K, j = 3,4,5 mod 8) + floor((K+2)/4) - K agrees with the 8-periodic run 0 -1 -1 0 1 2 2 1 at 201 of 201 values K = 0..200.
- hexagon.md the exact doubling constant: the limit integral returns covariance 253/9216 and variance 15/64 at all four branches, so r = 253/2160, with (301,601) reading -0.11745304 and (601,1201) -0.11729091, gap times m at -0.0974 and -0.0969.
- hexagon.md every pair limit: carpet and net read -+253/2160, -11/135 and +29/135, the tree 0, -61/864 and +4/27, the void 0, +7/864 and +2/27; the tree and void doubling covariances are exactly 0.
- hexagon.md every pair limit, the finite-layer gap: (249,251) reads -0.08150224 against -11/135 and (99,297) reads +0.21476417 against +29/135, gaps 2.1e-5 and 5.1e-5 with gap times m at -0.0052 and -0.0050.
- hexagon.md every pair limit, the second residue class: (103,205) reads +0.11914004 and (203,405) reads +0.11814528 against +253/2160, gap times m at +0.207 and +0.206, about twice the m = 1 mod 4 class's.
- flat-carpet-stack.md the ray law: 20 lines at 3 generic points and 500000 odd layers read within 3.7e-5 of -sgn(ab)(-1)^c/|ab| on the odd-odd lines and 2.1e-5 of 0 on the rest, carpet excess within 1.7e-5 of -pair/4, void excess within 1.9e-5 of +pair/2; 181 of 181 reduced a/b with b <= 24 match (-1)^a over b for odd b and 0 for even b; u + v = 1 is solid void ink at 300000 of 300000 layers.
- flat-carpet-stack.md the gcd sums: Jordan grouping against the gcd double sum at N = 6001 to 2.8e-11; sum of h_c/c^2 over odd c 1.3810359531 and over all n 1.9069358482 against 7 zeta(3)/16 + pi^2 ln2/8 and 7 zeta(3)/8 + pi^2 ln2/8; T2 = 1.7102721160 = pi^2 ln2/4; T4 = 1.1122336970; S2/N 0.8130215039 and S4/N 0.5536131179 at N = 999999 against 0.8130217042 and 0.5536124372.
- flat-carpet-stack.md the variance constants at N = 999999: tree 0.4064818975, void 0.2768055212, carpet 0.2724275136, net 0.2724571445 against 0.4065108521, 0.2768062186, 0.2724569807 and 0.2724569807; reduced fractions with odd denominator b <= 999999 number 202642133233 against (2/pi^2) B^2 = 202641962000.1.
- diagonal-slice-stack.md the eternal points: 34 grid points of constant parity to n = 55 and 5 to n = 2001, the permutations of (1/4, 1/4, 1), (0, 3/4, 3/4) and (1/4, 3/4, 1/2); (3/4, 3/4, 0) alternates 0 1 0 1 and (3/4, 0, 3/4) lies outside the layer.
