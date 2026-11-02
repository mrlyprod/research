mod blocks;
mod census;
mod design;
mod mass;
mod pairs;
mod radii;
mod terms;

use design::Design;

fn run_terms(design: &Design, top: u32, threads: usize) {
    let stored = terms::stored(design);
    let delta = design.density();
    println!(
        "{} k {} delta {:.6} live domain n <= {}",
        design.name, design.fill, delta, top
    );
    let mut last_gap = 0.0;
    let mut agree = true;
    for term in terms::terms(design, top, threads) {
        let ratio = term.value as f64 / (design.fill as f64).powi(term.level as i32);
        let gap = ratio - delta;
        let kept = stored
            .iter()
            .find(|(n, _)| *n == term.level)
            .map(|(_, v)| *v);
        let tag = match kept {
            Some(v) if v == term.value => "stored",
            Some(_) => "DIFFERS",
            None => "fresh",
        };
        agree &= tag != "DIFFERS";
        println!(
            "n {:2} A {:>26} A/k^n {:.6} gap {:+.2e} halving {:.2} {:.1}s {}",
            term.level,
            term.value,
            ratio,
            gap,
            if last_gap != 0.0 { gap / last_gap } else { 0.0 },
            term.seconds,
            tag
        );
        last_gap = gap;
    }
    for (n, v) in stored.iter().filter(|(n, _)| *n > top) {
        let gap = *v as f64 / (design.fill as f64).powi(*n as i32) - delta;
        println!("n {:2} A {:>26} gap {:+.2e} stored only", n, v, gap);
    }
    let brute_top = match design.dimension {
        3 => 5,
        _ => 8,
    };
    let brute: Vec<u64> = (1..=brute_top).map(|n| terms::brute(design, n)).collect();
    let matches = brute.iter().enumerate().all(|(i, b)| {
        stored
            .iter()
            .any(|(n, v)| *n as usize == i + 1 && *v == *b as i128)
    });
    println!(
        "{} live terms agree with stored {}, brute force to n = {} agrees {}",
        design.name, agree, brute_top, matches
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);
    match args.first().map(String::as_str) {
        Some("terms") => {
            let design = Design::named(&args[1]).expect("carpet, menger or vicsek");
            run_terms(design, args[2].parse().expect("top level"), threads);
        }
        Some("radii") => {
            radii::report(args[1].parse().expect("height"), threads);
        }
        Some("blocks") => match args.len() {
            3 => blocks::census(
                args[1].parse().expect("height"),
                args[2].parse().expect("length"),
                None,
            ),
            5 if args[3] == "probe" => blocks::census(
                args[1].parse().expect("height"),
                args[2].parse().expect("length"),
                Some(args[4].bytes().map(|b| b - b'0').collect()),
            ),
            n => blocks::single(
                args[1].parse().expect("a"),
                args[2].parse().expect("b"),
                args[3].parse().expect("length"),
                (n == 6 && args[4] == "probe").then(|| args[5].bytes().map(|b| b - b'0').collect()),
            ),
        },
        Some("mass") => {
            mass::report(
                args[1].parse().expect("a"),
                args[2].parse().expect("b"),
                args[3].parse().expect("top level"),
            );
        }
        Some("pairs") => {
            pairs::report(args[1].parse().expect("height"), threads);
        }
        Some("ray") => {
            radii::single(args[1].parse().expect("a"), args[2].parse().expect("b"));
        }
        Some("census") => {
            for level in args[1..].iter().map(|a| a.parse().expect("level")) {
                census::report(level, level == 13);
            }
        }
        _ => {
            run_terms(&design::CARPET, 18, threads);
            run_terms(&design::MENGER, 16, threads);
            run_terms(&design::VICSEK, 18, threads);
            for level in 13..=16 {
                census::report(level, level == 13);
            }
        }
    }
}
