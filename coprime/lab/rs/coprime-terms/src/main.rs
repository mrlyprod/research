mod check;

use coprime_terms::design::Design;
use coprime_terms::engine::{terms_each, Mode};

fn terms(args: &[String]) {
    let name = args.get(1).map(|s| s.as_str()).unwrap_or("menger");
    let top: u32 = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(8);
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    let design = Design::named(name).expect("unknown design");
    println!(
        "{} q 3 D {} k {} threads {}",
        design.name,
        design.dimension,
        design.fill(),
        threads
    );
    println!("n value seconds");
    println!("0 0 0.000");
    terms_each(&design, top, threads, Mode::Auto, &mut |level| {
        println!("{} {} {:.3}", level.level, level.value, level.seconds);
        use std::io::Write;
        std::io::stdout().flush().ok();
    });
}

fn profile(args: &[String]) {
    let name = args.get(1).map(|s| s.as_str()).unwrap_or("menger");
    let level: u32 = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(14);
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(8);
    let design = Design::named(name).expect("unknown design");
    coprime_terms::engine::profile(&design, level, threads);
}

fn main() {
    let all: Vec<String> = std::env::args().collect();
    let rest = &all[1..];
    match rest.first().map(|s| s.as_str()) {
        Some("check") => check::run(rest),
        Some("profile") => profile(rest),
        _ => terms(rest),
    }
}
