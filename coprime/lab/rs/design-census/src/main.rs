mod boundary;
mod classes;
mod coprime;
mod fills;
mod quasi;
mod tables;

use mrly_math::bang::Code;
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("limits") => boundary::limits(),
        Some("dim4") => boundary::dim4(),
        Some("terms") => boundary::terms(
            args[1].parse().expect("base"),
            args[2].parse().expect("dimension"),
            Code::from(args[3].parse::<u128>().expect("code")),
        ),
        _ => {
            let here = Path::new(env!("CARGO_MANIFEST_DIR"));
            coprime::report();
            println!();
            fills::orbits_report();
            println!();
            fills::report(&here.join("sequences.csv"));
            println!();
            classes::report(&here.join("counts.csv"));
        }
    }
}
