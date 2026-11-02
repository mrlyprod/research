use coprime_terms::brute;
use coprime_terms::design::Design;
use coprime_terms::engine::{caps, count_one, methods, terms_with, Mode};
fn parse(args: &[String], index: usize, fallback: u32) -> u32 {
    args.get(index)
        .map(|s| s.parse().unwrap())
        .unwrap_or(fallback)
}

pub fn run(args: &[String]) {
    let args: Vec<String> = args.to_vec();
    let name = args.get(1).map(|s| s.as_str()).unwrap_or("menger");
    let brute_top = parse(&args, 2, 6);
    let sample_level = parse(&args, 3, 14);
    let threads: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    let design = Design::named(name).expect("unknown design");
    let fill = design.fill();
    println!("design {} D {} k {}", design.name, design.dimension, fill);

    let engine = terms_with(&design, brute_top, threads, Mode::Auto);
    for level in engine.iter() {
        let direct = brute::count(&design, level.level, threads);
        println!(
            "brute level {} engine {} enumerated {} {}",
            level.level,
            level.value,
            direct,
            if level.value == direct as i128 {
                "MATCH"
            } else {
                "FAIL"
            }
        );
    }

    for level in 3..=8u32 {
        let span = 3u64.pow(level);
        let mut bad = 0;
        let mut seen = 0;
        for modulus in 1..span {
            if modulus % 3 == 0 {
                continue;
            }
            let all = methods(&design, level, modulus);
            seen += 1;
            if all.iter().any(|v| *v != all[0]) {
                bad += 1;
                if bad == 1 {
                    println!("  first split at modulus {} {:?}", modulus, all);
                }
            }
        }
        println!(
            "methods level {} moduli {} splits {} {}",
            level,
            seen,
            bad,
            if bad == 0 { "MATCH" } else { "FAIL" }
        );
    }

    let pinned: [(Mode, u32); 6] = [
        (Mode::Direct, 10),
        (Mode::Convolve, 10),
        (Mode::Zeta, 11),
        (Mode::Bitset, 8),
        (Mode::Rows, 10),
        (Mode::Cube, 10),
    ];
    let reference = terms_with(&design, 11, threads, Mode::Auto);
    for (mode, top) in pinned.iter() {
        let run = terms_with(&design, *top, threads, *mode);
        let same = run
            .iter()
            .zip(reference.iter())
            .all(|(a, b)| a.value == b.value);
        println!(
            "pinned {:?} levels 1..{} {}",
            mode,
            top,
            if same { "MATCH" } else { "FAIL" }
        );
    }

    let span = 3u64.pow(sample_level);
    let bounds = caps(sample_level, design.dimension);
    println!(
        "sample level {} caps residue {} bitset {} rows {}",
        sample_level, bounds.residue, bounds.bitset, bounds.rows
    );
    let mut probes: Vec<u64> = Vec::new();
    for shift in 0..14u32 {
        let base = span >> shift;
        for step in 0..3u64 {
            let candidate = base.saturating_sub(step * 7).max(1);
            if !candidate.is_multiple_of(3) && !probes.contains(&candidate) {
                probes.push(candidate);
            }
        }
    }
    let wide = 3u64.pow(sample_level.saturating_sub(11)) + 1;
    for extra in [
        1u64, 2, 5, 7, 11, 13, 17, 41, 61, 101, 157, 437, 439, 440, 443, 446, wide,
    ] {
        if !probes.contains(&extra) {
            probes.push(extra);
        }
    }
    probes.sort_unstable();
    let mut bad = 0;
    for modulus in probes.iter() {
        let reach = span / modulus + 1;
        let middle = reach <= 20_000 || *modulus == wide;
        let mut seen: Vec<(String, u128)> = Vec::new();
        if design.dimension == 3 {
            seen.push((
                "convolve".to_string(),
                count_one(&design, sample_level, *modulus, Mode::Convolve),
            ));
            seen.push((
                "cube".to_string(),
                count_one(&design, sample_level, *modulus, Mode::Cube),
            ));
            if middle {
                seen.push((
                    "zeta".to_string(),
                    count_one(&design, sample_level, *modulus, Mode::Zeta),
                ));
                seen.push((
                    "rows".to_string(),
                    count_one(&design, sample_level, *modulus, Mode::Rows),
                ));
            }
            if reach <= 4_000 {
                seen.push((
                    "direct".to_string(),
                    count_one(&design, sample_level, *modulus, Mode::Direct),
                ));
                seen.push((
                    "bitset".to_string(),
                    count_one(&design, sample_level, *modulus, Mode::Bitset),
                ));
            }
        } else {
            seen.push((
                "zeta".to_string(),
                count_one(&design, sample_level, *modulus, Mode::Zeta),
            ));
            if reach <= 40_000 {
                seen.push((
                    "direct".to_string(),
                    count_one(&design, sample_level, *modulus, Mode::Direct),
                ));
            }
        }
        if modulus
            .checked_pow(design.dimension as u32)
            .map(|v| v <= 8_000_000)
            .unwrap_or(false)
        {
            seen.push((
                "residue".to_string(),
                count_one(&design, sample_level, *modulus, Mode::Residue),
            ));
        }
        let agree = seen.iter().all(|(_, v)| *v == seen[0].1);
        if !agree {
            bad += 1;
            println!("  modulus {} {:?}", modulus, seen);
        }
    }
    println!(
        "probes level {} count {} splits {} {}",
        sample_level,
        probes.len(),
        bad,
        if bad == 0 { "MATCH" } else { "FAIL" }
    );

    if design.dimension == 3 {
        let deep = 19u32;
        let shelf = 3u64.pow(deep - 11) + 1;
        let mut ring_bad = 0;
        for (modulus, pinned, with_rows) in [
            (shelf, 58164940132507u128, true),
            (440, 61610728675376604u128, false),
        ] {
            let mut seen: Vec<(String, u128)> = vec![
                (
                    "cube".to_string(),
                    count_one(&design, deep, modulus, Mode::Cube),
                ),
                ("pinned".to_string(), pinned),
            ];
            if with_rows {
                seen.push((
                    "rows".to_string(),
                    count_one(&design, deep, modulus, Mode::Rows),
                ));
            }
            if !seen.iter().all(|(_, v)| *v == seen[0].1) {
                ring_bad += 1;
                println!("  ring modulus {} {:?}", modulus, seen);
            }
        }
        println!(
            "ring level {} moduli 2 splits {} {}",
            deep,
            ring_bad,
            if ring_bad == 0 { "MATCH" } else { "FAIL" }
        );
    }

    let mut peel = 0;
    let mut peel_bad = 0;
    if design.zero_filled() {
        for modulus in [1u64, 2, 5, 7, 11, 23, 47, 101, 1009, 20011] {
            if modulus >= 3u64.pow(sample_level - 1) {
                continue;
            }
            let mode = if design.dimension == 3 {
                Mode::Convolve
            } else {
                Mode::Zeta
            };
            let high = count_one(&design, sample_level, 3 * modulus, mode);
            let low = count_one(&design, sample_level - 1, modulus, mode);
            peel += 1;
            if high != low {
                peel_bad += 1;
                println!("  peel modulus {} {} vs {}", modulus, high, low);
            }
        }
    }
    println!(
        "peel checks {} splits {} {}",
        peel,
        peel_bad,
        if peel_bad == 0 { "MATCH" } else { "FAIL" }
    );

    let whole = count_one(&design, sample_level, 1, Mode::Residue);
    let expect = (fill as u128).pow(sample_level);
    println!(
        "fill level {} residue {} k^n {} {}",
        sample_level,
        whole,
        expect,
        if whole == expect { "MATCH" } else { "FAIL" }
    );
}
