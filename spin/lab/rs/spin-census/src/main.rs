mod design;
mod eisenstein;
mod gaussian;
mod mass;
mod orbit;
mod powder;
mod shadow;
mod shape;

use design::{bit_cells, canonical, carry, plane, square_group, BASE};
use mass::{distance, horizon, ripple, scaling_error, shells, Ripple};

const SUBJECTS: [u128; 7] = [79, 95, 127, 239, 255, 495, 511];
const BINS: usize = 24;
const LOW: f64 = 27.0;

fn dimension(code: u128) -> f64 {
    (code.count_ones() as f64).ln() / (BASE as f64).ln()
}

fn list(values: &[f64], places: usize) -> String {
    values
        .iter()
        .map(|value| format!("{value:.places$}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn corner_bit(table: &[(usize, usize)]) -> usize {
    table
        .iter()
        .position(|cell| *cell == (0, 0))
        .expect("a corner digit")
}

fn centre_bit(table: &[(usize, usize)]) -> usize {
    table
        .iter()
        .position(|cell| *cell == (1, 1))
        .expect("a centre digit")
}

fn read(
    code: u128,
    level: usize,
    digit: (usize, usize),
    table: &[(usize, usize)],
) -> (Ripple, f64, u64) {
    let grid = plane(code, BASE, level);
    let side = grid.shape[0] as f64;
    let bulk = shells(&grid, digit);
    let far = horizon(code, digit, table) * side;
    let error = scaling_error(&bulk, LOW, far, code.count_ones() as f64);
    let curve = ripple(&bulk, LOW, far, dimension(code), BINS);
    (curve, error, bulk.total())
}

fn anchors(table: &[(usize, usize)]) {
    println!("ANCHORS");
    println!("bit to cell map: {table:?}");
    let group = square_group();
    let mut seen = vec![false; 512];
    let mut classes = 0;
    for code in 0..512u128 {
        if seen[code as usize] {
            continue;
        }
        classes += 1;
        for member in design::orbit(&group, code, table) {
            seen[member as usize] = true;
        }
    }
    println!(
        "square group order {} classes over 512 codes {classes}",
        group.len()
    );
    let grid = plane(127, BASE, 4);
    let bulk = shells(&grid, (1, 1));
    let side = grid.shape[0];
    let rings = mrly_math::spin::profile(&design::floats(&grid), side, 6000).unwrap();
    let crate_mass = mrly_math::spin::mass(&rings, side);
    println!(
        "code 127 level 4: shell total {} crate profile mass {:.2} gap {:.2e}",
        bulk.total(),
        crate_mass,
        (crate_mass - bulk.total() as f64).abs() / bulk.total() as f64
    );
    let step = mrly_math::spin::reach(side) / 5999.0;
    for radius in [9.0f64, 27.0, 40.5] {
        let cut = (radius / step) as usize;
        let partial: f64 = (1..=cut)
            .map(|k| {
                let (a, b) = (k as f64 * step, (k - 1) as f64 * step);
                std::f64::consts::PI * (a + b) * (rings[k] as f64 + rings[k - 1] as f64) / 2.0
                    * step
            })
            .sum();
        println!(
            "  radius {radius:>5.1} shell count {:>6} crate profile integral {:>9.2} ratio {:.4}",
            bulk.at(radius),
            partial,
            partial / bulk.at(radius) as f64
        );
    }
}

fn holes(table: &[(usize, usize)]) {
    println!();
    println!("THE CENTRE HOLE, EXACT NEAREST FILLED CELL TO THE RASTER CENTRE");
    let centre = centre_bit(table);
    for code in SUBJECTS {
        let grid = plane(code, BASE, 5);
        let side = grid.shape[0];
        let four = mass::nearest_cell(&grid, (1, 1));
        let sixth = (BASE as u64).pow(4).pow(2);
        println!(
            "  code {code:>3} fill {} centre digit {} four times the squared distance {four} against (side/3)^2 = {sixth} at least {} distance {:.6} side/6 {:.6}",
            code.count_ones(),
            if code >> centre & 1 == 1 { "filled" } else { "empty " },
            four >= sixth,
            (four as f64).sqrt() / 2.0,
            side as f64 / 6.0
        );
    }
}

fn spin_dimension(table: &[(usize, usize)], level: usize) {
    println!();
    println!(
        "SPIN DIMENSION AT LEVEL {level}, FIXED POINT OF THE CORNER DIGIT AND OF THE CENTRE DIGIT"
    );
    let corner = corner_bit(table);
    let centre = centre_bit(table);
    for code in SUBJECTS {
        for (name, bit, digit) in [
            ("corner", corner, (0usize, 0usize)),
            ("centre", centre, (1, 1)),
        ] {
            if code >> bit & 1 == 0 {
                println!("  code {code:>3} {name}: digit empty, no fixed point");
                continue;
            }
            let (curve, error, total) = read(code, level, digit, table);
            println!(
                "  code {code:>3} fill {} {name}: D exact {:.6} slope {:.6} gap {:.2e} periods {} scaling error {:.2e} ripple swing {:.5} drift {:.5} mass {total}",
                code.count_ones(),
                dimension(code),
                curve.slope,
                (curve.slope - dimension(code)).abs(),
                curve.periods,
                error,
                curve.swing,
                curve.drift
            );
        }
    }
}

fn acid(table: &[(usize, usize)], level: usize) {
    println!();
    println!("THE EQUAL DIMENSION PAIRS AT LEVEL {level}, CORNER FIXED POINT");
    let corner = corner_bit(table);
    let mut kept: Vec<(u128, Ripple)> = Vec::new();
    for code in SUBJECTS {
        if code >> corner & 1 == 0 {
            continue;
        }
        let (curve, _, _) = read(code, level, (0, 0), table);
        println!("  code {code:>3} ripple {}", list(&curve.curve, 4));
        kept.push((code, curve));
    }
    for left in 0..kept.len() {
        for right in left + 1..kept.len() {
            if kept[left].0.count_ones() != kept[right].0.count_ones() {
                continue;
            }
            let gap = distance(&kept[left].1.curve, &kept[right].1.curve);
            let bar = kept[left].1.drift.max(kept[right].1.drift);
            println!(
                "  {} against {}: fill {} ripple gap {:.5} drift bar {:.5} ratio {:.1}",
                kept[left].0,
                kept[right].0,
                kept[left].0.count_ones(),
                gap,
                bar,
                gap / bar
            );
        }
    }
}

fn ripple_census(table: &[(usize, usize)], level: usize) {
    println!();
    println!("RIPPLE CENSUS AT LEVEL {level}: EVERY CODE WITH A FILLED CORNER DIGIT");
    let corner = corner_bit(table);
    let group = square_group();
    let transpose = group
        .iter()
        .find(|map| map[1] == BASE && map[BASE] == 1)
        .expect("the transpose")
        .clone();
    let flip: Vec<u128> = (0..512u128)
        .map(|code| carry(&transpose, code, table))
        .collect();
    let mut curves: Vec<(u128, Ripple)> = Vec::new();
    for code in 1..512u128 {
        if code >> corner & 1 == 0 {
            continue;
        }
        let (curve, _, _) = read(code, level, (0, 0), table);
        curves.push((code, curve));
    }
    let index = |code: u128| curves.iter().position(|row| row.0 == code);
    let mut mirror: f64 = 0.0;
    for (code, curve) in &curves {
        if let Some(at) = index(flip[*code as usize]) {
            mirror = mirror.max(distance(&curve.curve, &curves[at].1.curve));
        }
    }
    println!(
        "  codes read {} transpose control, worst ripple gap {:.2e}",
        curves.len(),
        mirror
    );
    let stamp = |code: u128| code.min(flip[code as usize]);
    let mut collisions: Vec<(f64, f64, u128, u128)> = Vec::new();
    let mut closest = (f64::INFINITY, 0u128, 0u128);
    for left in 0..curves.len() {
        for right in left + 1..curves.len() {
            let (a, b) = (curves[left].0, curves[right].0);
            if a.count_ones() != b.count_ones() || stamp(a) == stamp(b) {
                continue;
            }
            if a != stamp(a) || b != stamp(b) {
                continue;
            }
            let gap = distance(&curves[left].1.curve, &curves[right].1.curve);
            let bar = curves[left].1.drift.max(curves[right].1.drift);
            if gap < closest.0 {
                closest = (gap, a, b);
            }
            if gap < bar {
                collisions.push((gap, bar, a, b));
            }
        }
    }
    collisions.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite"));
    println!(
        "  transpose class pairs of equal fill whose ripples sit inside their own drift bar: {}",
        collisions.len()
    );
    for (gap, bar, left, right) in &collisions {
        println!(
            "    classes {left:>3} and {right:>3} fill {} gap {gap:.5} bar {bar:.5} ratio {:.2}",
            left.count_ones(),
            gap / bar
        );
    }
    let swing = |code: u128| curves[index(code).expect("a read code")].1.swing;
    println!(
        "  closest distinct class pair: {} and {} fill {} gap {:.5} with ripple swings {:.5} and {:.5}",
        closest.1,
        closest.2,
        closest.1.count_ones(),
        closest.0,
        swing(closest.1),
        swing(closest.2)
    );
}

fn powder_rings(level: usize, pad: usize, only: &[u128]) {
    println!();
    println!(
        "POWDER RINGS AT LEVEL {level}, PAD {pad}, RING AVERAGED POWER AGAINST THE FREQUENCY INDEX"
    );
    let side = (BASE as f64).powi(level as i32);
    let low = 3.0 * pad as f64 / side;
    let high = pad as f64 / 8.0;
    println!("  band {low:.1} to {high:.1} in frequency index, bins 240, phases {BINS}, slide window 3 periods by a quarter period");
    for code in SUBJECTS {
        if !only.is_empty() && !only.contains(&code) {
            continue;
        }
        let grid = plane(code, BASE, level);
        let read = powder::powder(&grid, pad, low, high, 240, BINS);
        let dim = dimension(code);
        println!(
            "  code {code:>3} fill {} D {:.6} band slope {:.5} against -D {:.5} gap {:.4} slide {:.5} to {:.5} spread {:.4} porod -3 gap {:.4} log period swing {:.4}",
            code.count_ones(),
            dim,
            read.slope,
            -dim,
            (read.slope + dim).abs(),
            read.low,
            read.high,
            read.high - read.low,
            (read.slope + 3.0).abs(),
            read.swing
        );
    }
}

fn spin_spectrum(table: &[(usize, usize)]) {
    println!();
    println!("SPIN SPECTRUM, P_m OVER ALL 512 BASE THREE CODES");
    let group = square_group();
    let stamp: Vec<u128> = (0..512u128)
        .map(|code| canonical(&group, code, table))
        .collect();
    let first = orbit::census(1, 1024, 12);
    let second = orbit::census(2, 768, 12);
    let mut pairs: Vec<(u128, u128, f64, f64)> = Vec::new();
    for left in 1..512u128 {
        for right in left + 1..512u128 {
            if left.count_ones() != right.count_ones() {
                continue;
            }
            if stamp[left as usize] == stamp[right as usize] {
                continue;
            }
            let one = orbit::gap(&first[left as usize], &first[right as usize]);
            let two = orbit::gap(&second[left as usize], &second[right as usize]);
            if one < 1e-9 && two < 1e-9 {
                pairs.push((left, right, one, two));
            }
        }
    }
    println!(
        "  pairs outside one square class agreeing at levels 1 and 2 to 1e-9: {}",
        pairs.len()
    );
    let mut classes: Vec<(u128, u128)> = pairs
        .iter()
        .map(|(a, b, _, _)| (stamp[*a as usize], stamp[*b as usize]))
        .collect();
    classes.sort_unstable();
    classes.dedup();
    println!("  distinct class pairs among them: {}", classes.len());
    for (left, right) in classes.iter().take(8) {
        let third = orbit::gap(
            &orbit::spectrum(*left, 3, 1024, 24),
            &orbit::spectrum(*right, 3, 1024, 24),
        );
        println!(
            "    classes {left} and {right} fill {} level 3 gap {:.2e} verdict {}",
            left.count_ones(),
            third,
            if third < 1e-9 {
                "isospectral"
            } else {
                "separated at level 3"
            }
        );
    }
    let mut buckets: Vec<Vec<u128>> = Vec::new();
    for code in 1..512u128 {
        let mut placed = false;
        for bucket in buckets.iter_mut() {
            let head = bucket[0] as usize;
            if first[head].len() == first[code as usize].len()
                && orbit::agree(&first[head], &first[code as usize], 1e-9)
                && orbit::agree(&second[head], &second[code as usize], 1e-9)
            {
                bucket.push(code);
                placed = true;
                break;
            }
        }
        if !placed {
            buckets.push(vec![code]);
        }
    }
    println!(
        "  distinct spin spectra over the 511 nonempty codes at levels 1 and 2 together: {}",
        buckets.len()
    );
    let big = buckets.iter().map(|bucket| bucket.len()).max().unwrap_or(0);
    println!("  largest spectral bucket holds {big} codes");
}

fn sponge_shadow(top: usize) {
    println!();
    println!("THE SPONGE SHADOW: LATTICE LINES IN DIRECTION (a,b,c) MEETING THE LEVEL L SPONGE, AND THE SAME FOR THE SOLID CUBE");
    let digits = shadow::digits();
    let cube = shadow::cube_digits();
    println!(
        "  sponge digits {} cube digits {}",
        digits.len(),
        cube.len()
    );
    let mut full = Vec::new();
    for view in shadow::views(3) {
        let counts: Vec<usize> = (1..=top)
            .map(|level| shadow::shadow(level, view, &digits))
            .collect();
        let solid: Vec<usize> = (1..=top - 1)
            .map(|level| shadow::shadow(level, view, &cube))
            .collect();
        let share: Vec<f64> = solid
            .iter()
            .zip(&counts)
            .map(|(all, hit)| *hit as f64 / *all as f64)
            .collect();
        let opaque = share.iter().all(|value| (value - 1.0).abs() < 1e-12);
        if opaque {
            full.push(view);
        }
        println!(
            "  view {view:?} sponge {counts:?} cube {solid:?} share {} {}",
            list(&share, 5),
            if opaque { "opaque" } else { "see through" }
        );
    }
    println!(
        "  views whose lattice lines the sponge blocks completely to level {}: {full:?}",
        top - 1
    );
    let deep = shadow::shadow(top, [1, 1, 1], &cube);
    println!(
        "  the space diagonal at level {top}: sponge {} cube {deep} equal {}",
        shadow::shadow(top, [1, 1, 1], &digits),
        shadow::shadow(top, [1, 1, 1], &digits) == deep
    );
}

fn gaussian_farey(top: usize) {
    println!();
    println!("THE GAUSSIAN FAREY: RADII NEW AT SCALE n");
    let cap = 2 * top * top;
    let least = gaussian::least_factors(cap);
    let direct = gaussian::union_counts(top, false);
    let boxed = gaussian::union_counts(top, true);
    let rule: Vec<usize> = (1..=top)
        .map(|scale| gaussian::new_disc(scale, &least))
        .collect();
    let mismatch = direct
        .iter()
        .zip(&rule)
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(index, _)| index + 1)
        .collect::<Vec<_>>();
    println!(
        "  disc reading, new radii at n = 1..20: {:?}",
        &direct[..20.min(direct.len())]
    );
    println!(
        "  the square free rule reproduces the disc reading at every n up to {top}: {}",
        mismatch.is_empty()
    );
    if !mismatch.is_empty() {
        println!(
            "  first mismatches at n = {:?}",
            &mismatch[..8.min(mismatch.len())]
        );
    }
    println!(
        "  box reading, new radii at n = 1..20: {:?}",
        &boxed[..20.min(boxed.len())]
    );
    let split: Vec<usize> = (1..=top)
        .filter(|n| direct[n - 1] != boxed[n - 1])
        .collect();
    println!(
        "  first scale where the box and the disc disagree: {:?}",
        split.first()
    );
    let primitives: Vec<usize> = (1..=20)
        .map(|scale| {
            (1..=2 * scale * scale)
                .filter(|norm| gaussian::primitive_norm(*norm, &least))
                .count()
        })
        .collect();
    println!("  norms below 2n^2 with a primitive representation, n = 1..20: {primitives:?}");
    let prefix = gaussian::two_square_prefix(cap, &least);
    let sieve: Vec<usize> = (1..=top)
        .map(|scale| gaussian::mobius_count(scale, &prefix, &least))
        .collect();
    let broken: Vec<usize> = (1..=top)
        .filter(|n| sieve[n - 1] != direct[n - 1])
        .collect();
    println!("  the Mobius identity new(n) = sum_d mu(d) B(2n^2/d^2) over d | rad(n) holds to n = {top}: {}", broken.is_empty());
    if !broken.is_empty() {
        println!(
            "  it first fails at n = {:?}",
            &broken[..4.min(broken.len())]
        );
    }
    let reach = 192usize;
    let wide = gaussian::least_factors(2 * reach * reach);
    let long = gaussian::two_square_prefix(2 * reach * reach, &wide);
    println!(
        "  the radical six family, new(n)/B(2n^2) climbing to the Jordan factor {:.5}:",
        gaussian::jordan(6, &wide)
    );
    for scale in [6usize, 12, 24, 48, 96, 192] {
        let count = gaussian::mobius_count(scale, &long, &wide);
        let all = long[2 * scale * scale];
        println!(
            "    n {scale:>3} new {count:>6} all {all:>6} ratio {:.5}",
            count as f64 / all as f64
        );
    }
}

fn gaussian_limit(top: usize) {
    println!();
    println!("THE GAUSSIAN FAREY LIMIT: new(n) sqrt(ln n) / n^2 AGAINST sqrt 2 K J_2(rad n) FOR n UP TO {top}");
    let cap = 2 * top * top;
    let least = gaussian::least_factors(cap);
    let prefix = gaussian::two_square_prefix(cap, &least);
    println!(
        "  K = {:.16} (A064533), every k up to 2 n^2 classified, cut n <= {top}, no radius sampled",
        gaussian::LANDAU_RAMANUJAN
    );
    let landau = prefix[cap] as f64 * (cap as f64).ln().sqrt() / cap as f64;
    println!("  Landau factor at the cut: B({cap}) sqrt(ln {cap}) / {cap} = {landau:.6} against K");
    let radicals = [2usize, 3, 5, 6, 10, 30, 210];
    let mut by_class = Vec::new();
    for radical in radicals {
        let members: Vec<usize> = (2..=top)
            .filter(|n| gaussian::radical(*n, &least) == radical)
            .collect();
        println!(
            "  radical {radical:>3}: {} scales up to {top}, target sqrt2 K J_2 = {:.6}",
            members.len(),
            gaussian::landau_target(radical, &least)
        );
        let mut ladder = Vec::new();
        for shift in [0u32, 2, 4, 6] {
            let bound = top >> shift;
            if let Some(n) = members.iter().rev().find(|n| **n <= bound) {
                if !ladder.contains(n) {
                    ladder.push(*n);
                }
            }
        }
        for n in &ladder {
            let count = gaussian::mobius_count(*n, &prefix, &least);
            let value = gaussian::scaled(*n, count);
            let target = gaussian::landau_target(*n, &least);
            println!(
                "    n {n:>5} new {count:>9} scaled {value:.6} ratio {:.5}",
                value / target
            );
        }
        by_class.push((radical, members));
    }
    let band: Vec<(usize, f64)> = (top / 2 + 1..=top)
        .map(|n| {
            let count = gaussian::mobius_count(n, &prefix, &least);
            (
                n,
                gaussian::scaled(n, count) / gaussian::landau_target(n, &least),
            )
        })
        .collect();
    let low = band
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("a band");
    let high = band
        .iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .expect("a band");
    println!(
        "  every n in ({}, {top}]: ratio least {:.5} at n = {} (radical {}), greatest {:.5} at n = {} (radical {})",
        top / 2,
        low.1,
        low.0,
        gaussian::radical(low.0, &least),
        high.1,
        high.0,
        gaussian::radical(high.0, &least)
    );
    println!(
        "  classes listed: {:?}",
        by_class
            .iter()
            .map(|(radical, members)| (*radical, members.len()))
            .collect::<Vec<_>>()
    );
}

fn eisenstein_farey(top: usize) {
    println!();
    println!("THE EISENSTEIN FAREY: RADII NEW AT SCALE n IN THE DISC OF RADIUS SQRT 3");
    let cap = eisenstein::WINDOW * top * top;
    let least = gaussian::least_factors(cap);
    let direct = eisenstein::union_counts(top);
    let rule: Vec<usize> = (1..=top)
        .map(|scale| eisenstein::new_rule(scale, &least))
        .collect();
    let prefix = eisenstein::loeschian_prefix(cap, &least);
    let sieve: Vec<usize> = (1..=top)
        .map(|scale| eisenstein::mobius_count(scale, &prefix, &least, None))
        .collect();
    let without: Vec<usize> = (1..=top)
        .map(|scale| eisenstein::mobius_count(scale, &prefix, &least, Some(3)))
        .collect();
    let first = |left: &[usize], right: &[usize]| -> Option<usize> {
        (1..=top).find(|n| left[n - 1] != right[n - 1])
    };
    println!(
        "  lattice reading over every point of Z[omega], new radii at n = 1..30: {:?}",
        &direct[..30.min(direct.len())]
    );
    println!(
        "  the square free rule reproduces the lattice reading at every n up to {top}: {}",
        first(&direct, &rule).is_none()
    );
    println!(
        "  the Mobius identity new(n) = sum_d mu(d) L(3n^2/d^2) over d | rad(n) holds to n = {top}: {}",
        first(&direct, &sieve).is_none()
    );
    println!(
        "  the same sum with the prime 3 left out of the radical first fails at n = {:?}",
        first(&direct, &without)
    );
    println!(
        "  scales n <= {top} divisible by 3 on which the full sum holds: {} of {}",
        (3..=top)
            .step_by(3)
            .filter(|n| sieve[n - 1] == direct[n - 1])
            .count(),
        top / 3
    );
}

fn eisenstein_limit(top: usize) {
    println!();
    println!("THE EISENSTEIN FAREY LIMIT: new(n) sqrt(ln n) / n^2 AGAINST (3 / sqrt 2) K_3 J_2(rad n) FOR n UP TO {top}");
    let cap = eisenstein::WINDOW * top * top;
    let least = gaussian::least_factors(cap);
    let prefix = eisenstein::loeschian_prefix(cap, &least);
    println!(
        "  K_3 = {:.16}, every k up to 3 n^2 classified, cut n <= {top}, no radius sampled",
        eisenstein::HEX_CONSTANT
    );
    let partial = eisenstein::constant_from_primes(cap, &least);
    println!(
        "  K_3 = sqrt(B / (2 sqrt 3)), B = prod over p = 2 mod 3 of (1 - p^-2)^-1 taken over p <= {cap}: {partial:.12}, short of K_3 by {:.2e}, every omitted factor exceeding 1",
        eisenstein::HEX_CONSTANT - partial
    );
    let landau = prefix[cap] as f64 * (cap as f64).ln().sqrt() / cap as f64;
    println!(
        "  Landau factor at the cut: L({cap}) sqrt(ln {cap}) / {cap} = {landau:.6} against K_3"
    );
    let radicals = [2usize, 3, 5, 6, 7, 10, 30, 210];
    for radical in radicals {
        let members: Vec<usize> = (2..=top)
            .filter(|n| gaussian::radical(*n, &least) == radical)
            .collect();
        println!(
            "  radical {radical:>3}: {} scales up to {top}, target (3/sqrt2) K_3 J_2 = {:.6}",
            members.len(),
            eisenstein::target(radical, &least)
        );
        let mut ladder = Vec::new();
        for shift in [0u32, 2, 4, 6] {
            let bound = top >> shift;
            if let Some(n) = members.iter().rev().find(|n| **n <= bound) {
                if !ladder.contains(n) {
                    ladder.push(*n);
                }
            }
        }
        for n in &ladder {
            let count = eisenstein::mobius_count(*n, &prefix, &least, None);
            let value = gaussian::scaled(*n, count);
            println!(
                "    n {n:>5} new {count:>9} scaled {value:.6} ratio {:.5}",
                value / eisenstein::target(*n, &least)
            );
        }
    }
    let band: Vec<(usize, f64)> = (top / 2 + 1..=top)
        .map(|n| {
            let count = eisenstein::mobius_count(n, &prefix, &least, None);
            (
                n,
                gaussian::scaled(n, count) / eisenstein::target(n, &least),
            )
        })
        .collect();
    let low = band
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("a band");
    let high = band
        .iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .expect("a band");
    println!(
        "  every n in ({}, {top}]: ratio least {:.5} at n = {} (radical {}), greatest {:.5} at n = {} (radical {})",
        top / 2,
        low.1,
        low.0,
        gaussian::radical(low.0, &least),
        high.1,
        high.0,
        gaussian::radical(high.0, &least)
    );
    println!(
        "  the radical 3 family, new(n) / L(3n^2) climbing to the Jordan factor 8/9 = 0.88889:"
    );
    let mut scale = 3usize;
    while scale <= top {
        let count = eisenstein::mobius_count(scale, &prefix, &least, None);
        let all = prefix[eisenstein::WINDOW * scale * scale];
        println!(
            "    n {scale:>5} new {count:>9} all {all:>9} ratio {:.5}",
            count as f64 / f64::from(all)
        );
        scale *= 3;
    }
}

fn main() {
    let table = bit_cells();
    if std::env::args().any(|a| a == "shape") {
        shape::control();
        shape::base_five();
        shape::cube_three();
        return;
    }
    if std::env::args().any(|a| a == "limit") {
        gaussian_limit(4096);
        return;
    }
    if std::env::args().any(|a| a == "hex") {
        eisenstein_farey(64);
        eisenstein_limit(4096);
        return;
    }
    anchors(&table);
    holes(&table);
    spin_dimension(&table, 6);
    acid(&table, 7);
    ripple_census(&table, 6);
    ripple_census(&table, 7);
    powder_rings(7, 4096, &[]);
    powder_rings(7, 8192, &[127, 255, 495]);
    spin_spectrum(&table);
    sponge_shadow(5);
    gaussian_farey(64);
}
