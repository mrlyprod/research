use mrly_num::factor::mobius_sieve;
use mrly_math::memory::{allowed_windows, kappa, perron, Rule};
use std::env;
use std::time::Instant;

// THE WINDOW PROFILE

fn induced_tables() -> (Vec<u8>, Vec<u8>) {
    let mut to2 = vec![0u8; 256];
    let mut to1 = vec![0u8; 256];
    for p in 0..256usize {
        let mut a = 0u8;
        let mut b = 0u8;
        for w in 0..8usize {
            if (p >> w) & 1 == 0 {
                continue;
            }
            a |= 1 << (w >> 1);
            a |= 1 << (w & 3);
            b |= 1 << ((w >> 2) & 1);
            b |= 1 << ((w >> 1) & 1);
            b |= 1 << (w & 1);
        }
        to2[p] = a;
        to1[p] = b;
    }
    (to2, to1)
}

fn small_profiles(n: usize) -> (u8, u8, u8) {
    match n {
        1 => (0, 0, 2),
        2 => (0, 4, 3),
        3 => (0, 8, 2),
        _ => unreachable!(),
    }
}

fn digit_profiles(n: usize) -> ([u8; 3], [u8; 3]) {
    let bits = usize::BITS as usize - n.leading_zeros() as usize;
    let digits: Vec<usize> = (0..bits).rev().map(|i| (n >> i) & 1).collect();
    let mut plain = [0u8; 3];
    let mut zeroed = [0u8; 3];
    let mut padded = vec![0usize];
    padded.extend_from_slice(&digits);
    for k in 1..=3usize {
        for run in digits.windows(k) {
            let mut w = 0usize;
            for &d in run {
                w = (w << 1) | d;
            }
            plain[k - 1] |= 1 << w;
        }
        for run in padded.windows(k) {
            let mut w = 0usize;
            for &d in run {
                w = (w << 1) | d;
            }
            zeroed[k - 1] |= 1 << w;
        }
    }
    (plain, zeroed)
}

// THE ZETA TRANSFORM

fn subset_sums_u64(source: &[u64], bits: u32) -> Vec<u64> {
    let mut out = source.to_vec();
    for b in 0..bits {
        for w in 0..out.len() {
            if (w >> b) & 1 == 1 {
                out[w] += out[w ^ (1 << b)];
            }
        }
    }
    out
}

fn subset_sums_i64(source: &[i64], bits: u32) -> Vec<i64> {
    let mut out = source.to_vec();
    for b in 0..bits {
        for w in 0..out.len() {
            if (w >> b) & 1 == 1 {
                out[w] += out[w ^ (1 << b)];
            }
        }
    }
    out
}

// THE PHASE GRID

fn phase_grid(depth: u32, cap: usize) -> Vec<(String, usize)> {
    let steps = [
        1.0f64,
        1.189_207_115_002_721,
        std::f64::consts::SQRT_2,
        1.681_792_830_507_429,
    ];
    let mut out: Vec<(String, usize)> = Vec::new();
    for level in 8..=depth {
        for (j, step) in steps.iter().enumerate() {
            let x = ((1u64 << level) as f64 * step).floor() as usize;
            if x > cap {
                continue;
            }
            let label = format!("{}.{:02}", level, j * 25);
            if out.last().map(|(_, v)| *v) == Some(x) {
                continue;
            }
            out.push((label, x));
        }
    }
    out
}

// THE CLASS REPRESENTATIVES

fn window_maps(width: usize) -> Vec<Vec<usize>> {
    let windows = 1usize << width;
    let flip: Vec<usize> = (0..windows).map(|w| windows - 1 - w).collect();
    let rev: Vec<usize> = (0..windows)
        .map(|w| {
            let mut out = 0usize;
            for j in 0..width {
                out |= ((w >> j) & 1) << (width - 1 - j);
            }
            out
        })
        .collect();
    let identity: Vec<usize> = (0..windows).collect();
    let mut maps = vec![identity.clone(), flip.clone(), rev.clone()];
    maps.push((0..windows).map(|w| flip[rev[w]]).collect());
    maps.sort();
    maps.dedup();
    maps
}

fn representatives(width: usize) -> Vec<usize> {
    let windows = 1usize << width;
    let maps = window_maps(width);
    let codes = 1usize << windows;
    let mut rep = vec![usize::MAX; codes];
    for code in 0..codes {
        if rep[code] != usize::MAX {
            continue;
        }
        for map in &maps {
            let mut image = 0usize;
            for (w, &target) in map.iter().enumerate().take(windows) {
                if (code >> w) & 1 == 1 {
                    image |= 1 << target;
                }
            }
            rep[image] = code;
        }
    }
    rep
}

// THE SWEEP

struct Reading {
    mass: Vec<Vec<u64>>,
    meter: Vec<Vec<i64>>,
    peak: Vec<Vec<i64>>,
}

struct Sweep {
    phases: Vec<(String, usize)>,
    read: [Reading; 3],
    mertens: Vec<i64>,
}

fn sweep(mu: &[i8], profile: &[u8], phases: &[(String, usize)], anchors: &[usize]) -> Sweep {
    let (to2, to1) = induced_tables();
    let mut count = [vec![0u64; 4], vec![0u64; 16], vec![0u64; 256]];
    let mut mass = [vec![0i64; 4], vec![0i64; 16], vec![0i64; 256]];
    let mut current = [vec![0i64; 4], vec![0i64; 16], vec![0i64; 256]];
    let mut peak = [vec![0i64; 4], vec![0i64; 16], vec![0i64; 256]];
    let mut read = [
        Reading {
            mass: Vec::new(),
            meter: Vec::new(),
            peak: Vec::new(),
        },
        Reading {
            mass: Vec::new(),
            meter: Vec::new(),
            peak: Vec::new(),
        },
        Reading {
            mass: Vec::new(),
            meter: Vec::new(),
            peak: Vec::new(),
        },
    ];
    let mut mertens = Vec::new();
    let mut next_phase = 0usize;
    let mut next_anchor = 0usize;
    let cap = phases.last().map(|(_, x)| *x).unwrap_or(0);
    for n in 1..=cap {
        let (p1, p2, p3) = if n >= 4 {
            let p = profile[n] as usize;
            (to1[p] as usize, to2[p] as usize, p)
        } else {
            let (a, b, c) = small_profiles(n);
            (c as usize, b as usize, a as usize)
        };
        count[0][p1] += 1;
        count[1][p2] += 1;
        count[2][p3] += 1;
        let m = mu[n] as i64;
        if m != 0 {
            mass[0][p1] += m;
            mass[1][p2] += m;
            mass[2][p3] += m;
            for (slot, (seed, bits)) in [(p1, 3usize), (p2, 15), (p3, 255)].iter().enumerate() {
                let free = bits & !seed;
                let mut sub = free;
                loop {
                    let rule = seed | sub;
                    let value = current[slot][rule] + m;
                    current[slot][rule] = value;
                    let size = value.abs();
                    if size > peak[slot][rule] {
                        peak[slot][rule] = size;
                    }
                    if sub == 0 {
                        break;
                    }
                    sub = (sub - 1) & free;
                }
            }
        }
        if next_anchor < anchors.len() && n == anchors[next_anchor] {
            mertens.push(current[0][3]);
            next_anchor += 1;
        }
        if next_phase < phases.len() && n == phases[next_phase].1 {
            for slot in 0..3 {
                let bits = [2u32, 4, 8][slot];
                let a = subset_sums_u64(&count[slot], bits);
                let z = subset_sums_i64(&mass[slot], bits);
                assert_eq!(z, current[slot], "zeta meter against the running meter");
                read[slot].mass.push(a);
                read[slot].meter.push(z);
                read[slot].peak.push(peak[slot].clone());
            }
            next_phase += 1;
        }
    }
    Sweep {
        phases: phases.to_vec(),
        read,
        mertens,
    }
}

// THE MEMORYLESS CONTROLS

fn design_meter(mu: &[i8], levels: usize, digits: &[u64]) -> (i64, i64, u64) {
    let mut powers = vec![1u64; levels + 1];
    for i in 1..=levels {
        powers[i] = powers[i - 1] * 3;
    }
    let lead: Vec<u64> = digits.iter().copied().filter(|&d| d != 0).collect();
    let base = digits.len();
    let mut meter = 0i64;
    let mut peak = 0i64;
    let mut mass = 0u64;
    for length in 1..=levels {
        let tails = base.pow((length - 1) as u32);
        for head in &lead {
            for tail in 0..tails {
                let mut value = head * powers[length - 1];
                let mut rest = tail;
                for &power in powers.iter().take(length - 1) {
                    value += digits[rest % base] * power;
                    rest /= base;
                }
                meter += mu[value as usize] as i64;
                if meter.abs() > peak {
                    peak = meter.abs();
                }
                mass += 1;
            }
        }
    }
    (meter, peak, mass)
}

// THE DIRECT RECOUNT

fn recount(limit: usize) -> ([Vec<u64>; 3], [Vec<u64>; 3]) {
    let mut plain = [vec![0u64; 4], vec![0u64; 16], vec![0u64; 256]];
    let mut zeroed = [vec![0u64; 4], vec![0u64; 16], vec![0u64; 256]];
    for n in 1..=limit {
        let (p, z) = digit_profiles(n);
        for k in 1..=3usize {
            let bits = (1usize << (1 << k)) - 1;
            for target in [
                (p[k - 1] as usize, &mut plain[k - 1]),
                (z[k - 1] as usize, &mut zeroed[k - 1]),
            ] {
                let (seed, store) = target;
                let free = bits & !seed;
                let mut sub = free;
                loop {
                    store[seed | sub] += 1;
                    if sub == 0 {
                        break;
                    }
                    sub = (sub - 1) & free;
                }
            }
        }
    }
    (plain, zeroed)
}

fn accepts_agrees(limit: usize) {
    for n in 1..=limit {
        let bits = usize::BITS as usize - n.leading_zeros() as usize;
        let word: Vec<usize> = (0..bits).rev().map(|i| (n >> i) & 1).collect();
        let (p, _) = digit_profiles(n);
        for k in 1..=3usize {
            let codes = 1usize << (1 << k);
            for code in 0..codes {
                let rule = Rule::new(1, k, code as u64).unwrap();
                let want = rule.accepts(&word);
                let have = (p[k - 1] as usize) & !code == 0;
                assert_eq!(
                    want, have,
                    "profile membership against Rule::accepts at n {n} k {k} code {code}"
                );
            }
        }
    }
}

// THE REPORT

fn ratio(value: i64, mass: u64) -> String {
    if mass == 0 {
        return "-".to_string();
    }
    format!("{:.6}", value as f64 / (mass as f64).sqrt())
}

struct Facts {
    windows: usize,
    rho: f64,
    coupling: f64,
    closed: bool,
    zero_mass: u64,
    mass: u64,
    meter: i64,
    peak: i64,
}

fn print_rows(tag: &str, width: usize, code: usize, sweep: &Sweep) {
    let slot = width - 1;
    for (index, (label, _)) in sweep.phases.iter().enumerate() {
        let a = sweep.read[slot].mass[index][code];
        let m = sweep.read[slot].meter[index][code];
        let p = sweep.read[slot].peak[index][code];
        println!(
            "{tag} k={width} code={code} x={label} A={a} M={m} Mmax={p} r={} rmax={}",
            ratio(m, a),
            ratio(p, a)
        );
    }
}

type Control = (&'static str, &'static [u64], &'static [(usize, i64, i64)]);

fn main() {
    let started = Instant::now();
    let depth: u32 = env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(30);
    let cap = 1usize << depth;
    let mu = mobius_sieve(cap);
    let sieved = started.elapsed().as_secs_f64();
    let mut profile = vec![0u8; cap + 1];
    for n in 4..=cap {
        profile[n] = profile[n >> 1] | (1u8 << (n & 7));
    }
    let anchors: Vec<usize> = (1..=8)
        .map(|e| 10usize.pow(e))
        .filter(|&x| x <= cap)
        .collect();
    let phases = phase_grid(depth, cap);
    let walked = Instant::now();
    let sweep = sweep(&mu, &profile, &phases, &anchors);
    let walk = walked.elapsed().as_secs_f64();
    drop(profile);

    let mertens_want = [-1i64, 1, 2, -23, -48, 212, 1037, 1928];
    assert_eq!(
        sweep.mertens,
        mertens_want[..anchors.len()].to_vec(),
        "the full line against A084237"
    );
    println!(
        "control mertens anchors 10^1..10^{} = {:?} (A084237)",
        anchors.len(),
        sweep.mertens
    );

    let designs: [Control; 3] = [
        (
            "0 1",
            &[0, 1],
            &[(14, 11, 105), (16, 149, 173), (18, -30, 312)],
        ),
        ("1 2", &[1, 2], &[(18, -1461, 1582)]),
        ("0 2", &[0, 2], &[]),
    ];
    for (name, digits, pins) in designs {
        for levels in [14usize, 16, 18, 20] {
            let top = digits.iter().copied().max().unwrap() as usize
                * (3usize.pow(levels as u32) - 1)
                / 2;
            if top > cap {
                println!("control design base3 digits {name} L={levels} needs {top} which is past 2^{depth}, not pinned here");
                continue;
            }
            let (m, p, a) = design_meter(&mu, levels, digits);
            let pinned = pins.iter().find(|(l, _, _)| *l == levels);
            println!(
                "control design base3 digits {name} L={levels} A={a} M={m} Mmax={p} pinned={}",
                if pinned.is_some() { "yes" } else { "no" }
            );
            if let Some((_, want_m, want_p)) = pinned {
                assert_eq!(
                    (m, p),
                    (*want_m, *want_p),
                    "base 3 design against lab/rs/mobius-designs"
                );
            }
        }
    }

    let limit = (1usize << 20).min(cap);
    let (plain, zeroed) = recount(limit);
    accepts_agrees(1 << 12);
    let cut = phases.iter().position(|(_, x)| *x == limit).unwrap();
    for (slot, want) in plain.iter().enumerate() {
        assert_eq!(
            sweep.read[slot].mass[cut], *want,
            "the profile recurrence against a direct digit recount at 2^20"
        );
    }
    println!("control recount 2^20 agrees with the profile recurrence on all 276 rules");
    println!("control Rule::accepts agrees with profile containment on all 276 rules below 2^12");

    let last = sweep.phases.len() - 1;
    let mut facts: Vec<Vec<Facts>> = Vec::new();
    for width in 1..=3usize {
        let slot = width - 1;
        let codes = 1usize << (1 << width);
        let mut row = Vec::new();
        for code in 0..codes {
            let rule = Rule::new(1, width, code as u64).unwrap();
            row.push(Facts {
                windows: allowed_windows(&rule),
                rho: perron(&rule),
                coupling: kappa(&rule),
                closed: plain[slot][code] == zeroed[slot][code],
                zero_mass: zeroed[slot][code],
                mass: sweep.read[slot].mass[last][code],
                meter: sweep.read[slot].meter[last][code],
                peak: sweep.read[slot].peak[last][code],
            });
        }
        facts.push(row);
    }

    let reps = representatives(3);
    let class_reps: Vec<usize> = (0..256).filter(|&c| reps[c] == c).collect();
    println!(
        "classes k=3 representatives={} (G_(1,3) orbits)",
        class_reps.len()
    );

    let full: Vec<f64> = (0..sweep.phases.len())
        .map(|i| sweep.read[0].peak[i][3] as f64 / (sweep.read[0].mass[i][3] as f64).sqrt())
        .collect();
    let band_low = full.iter().cloned().fold(f64::INFINITY, f64::min);
    let band_high = full.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let strays = |floor: u64| -> Vec<(usize, usize, u64)> {
        (1..=3usize)
            .flat_map(|w| (0..(1usize << (1 << w))).map(move |c| (w, c)))
            .filter(|&(w, c)| facts[w - 1][c].mass >= floor)
            .filter(|&(width, code)| {
                let slot = width - 1;
                (0..sweep.phases.len()).all(|index| {
                    let a = sweep.read[slot].mass[index][code];
                    if a < floor {
                        return true;
                    }
                    let r = sweep.read[slot].peak[index][code] as f64 / (a as f64).sqrt();
                    r < band_low || r > band_high
                })
            })
            .map(|(w, c)| (w, c, facts[w - 1][c].mass))
            .collect()
    };

    let last_label = sweep.phases[last].0.clone();
    let score =
        |w: usize, c: usize| facts[w - 1][c].peak as f64 / (facts[w - 1][c].mass as f64).sqrt();
    let sweepwide = |w: usize, c: usize| {
        (0..sweep.phases.len())
            .filter(|&i| sweep.read[w - 1].mass[i][c] > 0)
            .map(|i| {
                (
                    sweep.read[w - 1].peak[i][c] as f64
                        / (sweep.read[w - 1].mass[i][c] as f64).sqrt(),
                    i,
                )
            })
            .fold(
                (0.0f64, 0usize),
                |best, next| {
                    if next.0 > best.0 {
                        next
                    } else {
                        best
                    }
                },
            )
    };
    let against = |w: usize, c: usize| {
        (0..sweep.phases.len())
            .filter(|&i| sweep.read[w - 1].mass[i][c] > 0)
            .map(|i| {
                let mine = sweep.read[w - 1].peak[i][c] as f64
                    / (sweep.read[w - 1].mass[i][c] as f64).sqrt();
                let line =
                    sweep.read[0].peak[i][3] as f64 / (sweep.read[0].mass[i][3] as f64).sqrt();
                (mine / line, i)
            })
            .fold(
                (0.0f64, 0usize),
                |best, next| {
                    if next.0 > best.0 {
                        next
                    } else {
                        best
                    }
                },
            )
    };
    let rank = |floor: u64| {
        let mut out: Vec<(usize, usize)> = (1..=3usize)
            .flat_map(|w| (0..(1usize << (1 << w))).map(move |c| (w, c)))
            .filter(|&(w, c)| facts[w - 1][c].mass >= floor)
            .collect();
        out.sort_by(|&(wa, a), &(wb, b)| {
            score(wb, b)
                .partial_cmp(&score(wa, a))
                .unwrap()
                .then((wa, a).cmp(&(wb, b)))
        });
        out
    };
    let mut named: Vec<usize> = vec![23, 54, 126, 127, 255];
    named.extend(
        strays(10_000)
            .iter()
            .filter(|(w, _, _)| *w == 3)
            .map(|(_, c, _)| *c),
    );
    for floor in [64u64, 10_000] {
        let list = rank(floor);
        let high: Vec<(usize, usize)> = list.iter().take(10).copied().collect();
        let low: Vec<(usize, usize)> = list.iter().rev().take(10).copied().collect();
        named.extend(
            high.iter()
                .chain(low.iter())
                .filter(|(w, _)| *w == 3)
                .map(|(_, c)| *c),
        );
        println!(
            "top all widths rmax at phase {last_label} over {} rules of at least {floor} elements: {}",
            list.len(),
            high.iter()
                .map(|&(w, c)| format!("k{w}code{c}:{:.6}", score(w, c)))
                .collect::<Vec<_>>()
                .join(",")
        );
        println!(
            "bottom all widths rmax at phase {last_label} over {} rules of at least {floor} elements: {}",
            list.len(),
            low.iter()
                .map(|&(w, c)| format!("k{w}code{c}:{:.6}", score(w, c)))
                .collect::<Vec<_>>()
                .join(",")
        );
        let span_low = list
            .iter()
            .map(|&(w, c)| score(w, c))
            .fold(f64::INFINITY, f64::min);
        let span_high = list
            .iter()
            .map(|&(w, c)| score(w, c))
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "span all widths floor={floor} rules={} rmax at phase {last_label} runs [{span_low:.6}, {span_high:.6}]",
            list.len()
        );
        let mut wide: Vec<(f64, usize, usize, usize)> = list
            .iter()
            .map(|&(w, c)| {
                let (v, i) = sweepwide(w, c);
                (v, w, c, i)
            })
            .collect();
        wide.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let wide_low = wide.iter().map(|r| r.0).fold(f64::INFINITY, f64::min);
        println!(
            "span sweepwide floor={floor} rules={} rmax over all {} phases runs [{wide_low:.6}, {:.6}], the top on k{} code {} at phase {}",
            wide.len(),
            sweep.phases.len(),
            wide[0].0,
            wide[0].1,
            wide[0].2,
            sweep.phases[wide[0].3].0
        );
        println!(
            "wide top ten sweepwide rmax floor={floor}: {}",
            wide.iter()
                .take(10)
                .map(|r| format!("k{}code{}:{:.6}@{}", r.1, r.2, r.0, sweep.phases[r.3].0))
                .collect::<Vec<_>>()
                .join(",")
        );
        let tail = sweep.phases.len() - sweep.phases.len() / 4;
        let late = wide.iter().filter(|r| r.3 >= tail).count();
        println!(
            "latepeak floor={floor} rules whose sweepwide rmax falls in the last quarter of the phases, from {}: {late} of {}",
            sweep.phases[tail].0,
            wide.len()
        );

        let mut factors: Vec<(f64, usize, usize, usize)> = list
            .iter()
            .map(|&(w, c)| {
                let (v, i) = against(w, c);
                (v, w, c, i)
            })
            .collect();
        factors.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        println!(
            "factor same phase against the full line, floor={floor}, largest over all phases: {}",
            factors
                .iter()
                .take(10)
                .map(|r| format!("k{}code{}:{:.6}@{}", r.1, r.2, r.0, sweep.phases[r.3].0))
                .collect::<Vec<_>>()
                .join(",")
        );
        let at_last: Vec<f64> = list
            .iter()
            .map(|&(w, c)| score(w, c) / score(1, 3))
            .collect();
        println!(
            "factor same phase at {last_label} floor={floor} runs [{:.6}, {:.6}] over {} rules",
            at_last.iter().cloned().fold(f64::INFINITY, f64::min),
            at_last.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            at_last.len()
        );
    }

    for width in 1..=3usize {
        let slot = width - 1;
        let codes = 1usize << (1 << width);
        for code in 0..codes {
            let f = &facts[slot][code];
            let mut track = Vec::new();
            for (index, (label, _)) in sweep.phases.iter().enumerate() {
                if !label.ends_with(".00") {
                    continue;
                }
                let level: u32 = label[..label.len() - 3].parse().unwrap();
                if !level.is_multiple_of(4) && level != depth {
                    continue;
                }
                let a = sweep.read[slot].mass[index][code];
                let p = sweep.read[slot].peak[index][code];
                track.push(format!("{level}:{}", ratio(p, a)));
            }
            println!(
                "rule k={width} code={code} rep={} W={} rho={:.9} kappa={:.6} zeroclosed={} zeroA={} A={} M={} Mmax={} r={} rmax={} track={}",
                if width == 3 { reps[code] } else { code },
                f.windows,
                f.rho,
                f.coupling,
                if f.closed { "yes" } else { "no" },
                f.zero_mass,
                f.mass,
                f.meter,
                f.peak,
                ratio(f.meter, f.mass),
                ratio(f.peak, f.mass),
                track.join(",")
            );
        }
    }

    for code in 0..4usize {
        print_rows("row", 1, code, &sweep);
    }
    for code in 0..16usize {
        print_rows("row", 2, code, &sweep);
    }
    named.sort();
    named.dedup();
    for code in &named {
        print_rows("row", 3, *code, &sweep);
    }

    let fibbinary: Vec<u64> = vec![
        1, 2, 4, 5, 8, 9, 10, 16, 17, 18, 20, 21, 32, 33, 34, 36, 37, 40, 41, 42,
    ];
    let golden: Vec<u64> = (1..=64u64)
        .filter(|&n| Rule::new(1, 2, 7).unwrap().accepts(&binary_word(n)))
        .collect();
    assert_eq!(
        golden[..fibbinary.len()],
        fibbinary[..],
        "code 7 at k=2 against A003714"
    );
    println!(
        "control code 7 at k=2 opens {:?} which is A003714 without its zero",
        &golden[..12]
    );
    let mersenne: Vec<u64> = (1..=1024u64)
        .filter(|&n| Rule::new(1, 2, 11).unwrap().accepts(&binary_word(n)))
        .collect();
    assert_eq!(
        mersenne,
        vec![1, 3, 7, 15, 31, 63, 127, 255, 511, 1023],
        "code 11 at k=2 is A000225"
    );
    assert_eq!(
        facts[1][11].mass, depth as u64,
        "code 11 at k=2 holds one element per level"
    );
    println!(
        "control code 11 at k=2 opens the Mersenne numbers A000225, A=2^m-1 count {}",
        facts[1][11].mass
    );

    let mut climbing: Vec<(usize, usize)> = Vec::new();
    let mut climbing_meter: Vec<(usize, usize)> = Vec::new();
    for width in 1..=3usize {
        let slot = width - 1;
        let codes = 1usize << (1 << width);
        for (code, fact) in facts[slot].iter().enumerate().take(codes) {
            if fact.mass < 1000 {
                continue;
            }
            let mut up = true;
            for index in (last - 7)..=last {
                let a = sweep.read[slot].mass[index][code] as f64;
                let p = sweep.read[slot].peak[index][code] as f64;
                let b = sweep.read[slot].mass[index - 1][code] as f64;
                let q = sweep.read[slot].peak[index - 1][code] as f64;
                if p / a.sqrt() <= q / b.sqrt() {
                    up = false;
                    break;
                }
            }
            if up {
                climbing.push((width, code));
            }
            let mut rises = true;
            for index in (last - 7)..=last {
                let a = sweep.read[slot].mass[index][code] as f64;
                let p = sweep.read[slot].meter[index][code].abs() as f64;
                let b = sweep.read[slot].mass[index - 1][code] as f64;
                let q = sweep.read[slot].meter[index - 1][code].abs() as f64;
                if p / a.sqrt() <= q / b.sqrt() {
                    rises = false;
                    break;
                }
            }
            if rises {
                climbing_meter.push((width, code));
            }
        }
    }
    println!(
        "climbing rmax rises at every one of the last eight phases on {} rules: {:?}",
        climbing.len(),
        climbing
    );
    println!(
        "climbing abs r rises at every one of the last eight phases on {} rules: {:?}",
        climbing_meter.len(),
        climbing_meter
    );

    let line = &facts[0][3];
    println!(
        "band fullline k=1 code=3 A={} M={} Mmax={} r={} rmax={}",
        line.mass,
        line.meter,
        line.peak,
        ratio(line.meter, line.mass),
        ratio(line.peak, line.mass)
    );
    for width in 2..=3usize {
        let slot = width - 1;
        let codes = 1usize << (1 << width);
        let edges = [0.0f64, 1e-9, 0.1, 0.2, 0.3, 0.5, 0.8, 2.0];
        for floor in [64u64, 10_000] {
            for pair in edges.windows(2) {
                let (lo, hi) = (pair[0], pair[1]);
                let members: Vec<usize> = (0..codes)
                    .filter(|&c| {
                        let f = &facts[slot][c];
                        f.mass >= floor && f.rho > 0.0 && f.coupling >= lo && f.coupling < hi
                    })
                    .collect();
                if members.is_empty() {
                    continue;
                }
                let scores: Vec<f64> = members
                    .iter()
                    .map(|&c| facts[slot][c].peak as f64 / (facts[slot][c].mass as f64).sqrt())
                    .collect();
                let mean = scores.iter().sum::<f64>() / scores.len() as f64;
                let least = scores.iter().cloned().fold(f64::INFINITY, f64::min);
                let most = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let argmax = members[scores.iter().position(|&s| s == most).unwrap()];
                let argmin = members[scores.iter().position(|&s| s == least).unwrap()];
                println!(
                "band k={width} floor={floor} kappa=[{lo},{hi}) rules={} rmax_mean={:.6} rmax_min={:.6} on code {argmin} rmax_max={:.6} on code {argmax}",
                members.len(),
                mean,
                least,
                most
            );
            }
        }
        let mut by_rho: Vec<(f64, usize)> = (0..codes)
            .filter(|&c| facts[slot][c].mass >= 64 && facts[slot][c].rho > 0.0)
            .map(|c| ((facts[slot][c].rho * 1e4).round() / 1e4, c))
            .collect();
        by_rho.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let mut seen: Vec<f64> = Vec::new();
        for (r, c) in &by_rho {
            if seen.contains(r) {
                continue;
            }
            seen.push(*r);
            let members: Vec<usize> = by_rho
                .iter()
                .filter(|(q, _)| q == r)
                .map(|(_, c)| *c)
                .collect();
            let scores: Vec<f64> = members
                .iter()
                .map(|&m| facts[slot][m].peak as f64 / (facts[slot][m].mass as f64).sqrt())
                .collect();
            let mean = scores.iter().sum::<f64>() / scores.len() as f64;
            println!(
                "rho k={width} rho={r:.4} rules={} least={c} rmax_mean={:.6} rmax_min={:.6} rmax_max={:.6}",
                members.len(),
                mean,
                scores.iter().cloned().fold(f64::INFINITY, f64::min),
                scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            );
        }
    }

    let mut running = 0i64;
    let mut ceiling = 0i64;
    let mut early = Vec::new();
    for (n, &sign) in mu.iter().enumerate().take(401).skip(1) {
        running += sign as i64;
        if running.abs() > ceiling {
            ceiling = running.abs();
        }
        if [1usize, 5, 13, 31, 200, 256].contains(&n) {
            early.push(format!("{n}:{:.6}", ceiling as f64 / (n as f64).sqrt()));
        }
    }
    println!(
        "gridstart the full line rmax below the grid reads {}",
        early.join(",")
    );

    let mut only_left = 0u64;
    let mut only_right = 0u64;
    let mut both = 0u64;
    let mut first_gap = 0usize;
    for n in 1..=(1usize << 20) {
        let (p, _) = digit_profiles(n);
        let left = (p[1] as usize) & !14usize == 0;
        let right = (p[2] as usize) & !126usize == 0;
        match (left, right) {
            (true, true) => both += 1,
            (true, false) => only_left += 1,
            (false, true) => {
                only_right += 1;
                if first_gap == 0 {
                    first_gap = n;
                }
            }
            _ => {}
        }
    }
    println!(
        "pair k=2 code 14 against k=3 code 126 below 2^20: shared={both} only14={only_left} only126={only_right} symmetric={} least in 126 not 14 is {first_gap}",
        only_left + only_right
    );

    for width in 1..=3usize {
        let codes = 1usize << (1 << width);
        let mut agree = Vec::new();
        for code in 0..codes {
            let rule = Rule::new(1, width, code as u64).unwrap();
            let mut same = true;
            'words: for length in 1..=14usize {
                for value in 0..(1usize << length) {
                    let word: Vec<usize> = (0..length).rev().map(|i| (value >> i) & 1).collect();
                    let trimmed: Vec<usize> =
                        word.iter().copied().skip_while(|&d| d == 0).collect();
                    if rule.accepts(&word) != rule.accepts(&trimmed) {
                        same = false;
                        break 'words;
                    }
                }
            }
            if same {
                agree.push(code);
            }
        }
        println!(
            "reading k={width} the word language agrees with the integer set on {} of {codes} codes over every word to length 14: {:?}",
            agree.len(),
            agree
        );
    }

    let edges = [0.0f64, 1e-9, 0.1, 0.2, 0.3, 0.5, 0.8, 2.0];
    for floor in [64u64, 10_000] {
        for only in [false, true] {
            for pair in edges.windows(2) {
                let (lo, hi) = (pair[0], pair[1]);
                let members: Vec<(usize, usize)> = (1..=3usize)
                    .flat_map(|w| (0..(1usize << (1 << w))).map(move |c| (w, c)))
                    .filter(|&(w, c)| {
                        let f = &facts[w - 1][c];
                        (!only || w == 3)
                            && f.mass >= floor
                            && f.rho > 0.0
                            && f.coupling >= lo
                            && f.coupling < hi
                    })
                    .collect();
                if members.is_empty() {
                    continue;
                }
                let scores: Vec<f64> = members.iter().map(|&(w, c)| score(w, c)).collect();
                let mean = scores.iter().sum::<f64>() / scores.len() as f64;
                let least = scores.iter().cloned().fold(f64::INFINITY, f64::min);
                let most = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let argmax = members[scores.iter().position(|&s| s == most).unwrap()];
                let argmin = members[scores.iter().position(|&s| s == least).unwrap()];
                println!(
                    "kappaband widths={} floor={floor} kappa=[{lo},{hi}) rules={} rmax_mean={mean:.6} at phase {last_label} rmax_min={least:.6} on k{} code {} rmax_max={most:.6} on k{} code {} members={}",
                    if only { "3" } else { "1,2,3" },
                    members.len(),
                    argmin.0,
                    argmin.1,
                    argmax.0,
                    argmax.1,
                    members
                        .iter()
                        .map(|&(w, c)| format!("k{w}c{c}"))
                        .collect::<Vec<_>>()
                        .join(" ")
                );
            }
        }
    }

    let top_at = full.iter().position(|&v| v == band_high).unwrap();
    let low_at = full.iter().position(|&v| v == band_low).unwrap();
    println!(
        "band control the full line rmax runs [{band_low:.6}, {band_high:.6}] over every phase, the top at phase {} and the floor at phase {}",
        sweep.phases[top_at].0, sweep.phases[low_at].0
    );
    for floor in [64u64, 10_000] {
        let census = (1..=3usize)
            .flat_map(|w| (0..(1usize << (1 << w))).map(move |c| (w, c)))
            .filter(|&(w, c)| facts[w - 1][c].mass >= floor)
            .count();
        let outside = strays(floor);
        println!(
            "band outside floor={floor} census={census} rules never inside the control band at any phase where they hold {floor} elements: {} {:?}",
            outside.len(),
            outside
        );
    }

    let closed_counts: Vec<usize> = (0..3)
        .map(|slot| {
            let codes = 1usize << (1 << (slot + 1));
            (0..codes).filter(|&c| facts[slot][c].closed).count()
        })
        .collect();
    println!("zeroclosed counts k=1,2,3 = {closed_counts:?} of 4, 16, 256 tested below 2^20");
    for width in 1..=3usize {
        let slot = width - 1;
        let codes = 1usize << (1 << width);
        for (code, fact) in facts[slot].iter().enumerate().take(codes) {
            let want = match width {
                1 => code == 0 || code & 1 == 1,
                2 => (code >> 1) & 1 == 1,
                _ => (code >> 2) & 1 == 1 && (code >> 3) & 1 == 1,
            };
            assert_eq!(
                fact.closed, want,
                "the zero-closed criterion at k {width} code {code}"
            );
        }
    }
    println!("zeroclosed criterion holds code for code: k=1 the codes allowing the digit 0 with the empty code, k=2 the codes allowing 01, k=3 the codes allowing both 010 and 011");
    println!(
        "run depth={depth} N={cap} sieve={sieved:.1}s sweep={walk:.1}s total={:.1}s phases={}",
        started.elapsed().as_secs_f64(),
        sweep.phases.len()
    );
}

fn binary_word(n: u64) -> Vec<usize> {
    let bits = 64 - n.leading_zeros() as usize;
    (0..bits).rev().map(|i| ((n >> i) & 1) as usize).collect()
}
