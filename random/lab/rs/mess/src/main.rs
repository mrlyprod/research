mod census;
mod frontier;
mod json;
mod laws;
mod mess;
mod sources;
mod typical;

use mess::{ceiling, symmetric_table, table};
use mrly_life::elementary::output;

fn rule_thirty() {
    println!("RULE 30");
    let values: Vec<u8> = (0..8usize)
        .map(|x| output(30, (x >> 2 & 1) as u8, (x >> 1 & 1) as u8, (x & 1) as u8))
        .collect();
    let t = table(&values, 3);
    println!(
        "l xor (c or r): degree={} ones={} nl={} autocorrelation_max={} differential={} ceiling(3)={}",
        t.degree,
        t.ones,
        t.nonlinearity,
        t.autocorrelation_max,
        t.differential,
        ceiling(3)
    );
    let majority = table(&symmetric_table(&[0, 0, 1, 1], 3), 3);
    println!(
        "majority on 3: degree={} ones={} nl={} differential={}",
        majority.degree, majority.ones, majority.nonlinearity, majority.differential
    );
}

fn frontier_report() -> Vec<frontier::Front> {
    let mut fronts = Vec::new();
    println!("FRONTIER over all 2^(k+1) symmetric functions");
    println!("k ceiling best_nl best_count full_degree_nl full_degree_count balanced_nl balanced_count functions_balanced");
    for k in 1..=frontier::KFRONT {
        let front = frontier::run(k);
        let (best_nl, best_count) = frontier::best(&front);
        assert_eq!(
            best_nl as u128,
            ceiling(k),
            "k={k} the ceiling is not reached"
        );
        let full = &front.rows[k];
        let balanced_walsh = front
            .rows
            .iter()
            .map(|r| r.balanced_best_walsh)
            .min()
            .unwrap();
        let balanced_count: u64 = front
            .rows
            .iter()
            .filter(|r| r.balanced_best_walsh == balanced_walsh)
            .map(|r| r.balanced_best_count)
            .sum();
        let balanced_total: u64 = front.rows.iter().map(|r| r.balanced_count).sum();
        let degrees: Vec<String> = front
            .rows
            .iter()
            .map(|r| {
                if r.count == 0 {
                    "-".to_string()
                } else {
                    format!(
                        "{}x{}",
                        frontier::nonlinearity(k, r.best_walsh),
                        r.best_count
                    )
                }
            })
            .collect();
        println!(
            "{k} {} {best_nl} {best_count} {} {} {} {balanced_count} {balanced_total}",
            ceiling(k),
            if full.count == 0 {
                0
            } else {
                frontier::nonlinearity(k, full.best_walsh)
            },
            full.best_count,
            if balanced_total == 0 {
                0
            } else {
                frontier::nonlinearity(k, balanced_walsh)
            },
        );
        println!("  by degree 0..{k} (nl x count): {}", degrees.join(" "));
        fronts.push(front);
    }
    fronts
}

fn typical_report() {
    println!(
        "TYPICAL symmetric function, {} draws, seed {:#x}",
        typical::DRAWS,
        typical::SEED
    );
    for (label, tree) in [("coin g", false), ("Source::Random(0..draws)", true)] {
        println!("{label}");
        for k in typical::KS {
            let t = typical::run(k, tree);
            println!(
                "k={} median nl/ceiling={:.4} min={:.4} max={:.4} full degree share={:.3} balanced share={:.3} median degree={}",
                t.k, t.median_ratio, t.min_ratio, t.max_ratio, t.full_degree, t.balanced, t.median_degree
            );
        }
    }
}

fn main() {
    let sources = sources::fixed();
    let entries = census::build(&sources);
    census::report_classes(&entries);
    println!();
    println!(
        "brute-force table agrees on {} cells of the classes (k <= 14)",
        census::brute(&entries)
    );
    let (cells, tight, total) = census::thin_bound(&entries);
    println!("NL <= min(wt, 2^k - wt) and NL <= ceiling on {cells} cells of the classes; NL = min(wt, 2^k - wt) on {tight} of {total} non-affine class cells with k >= 8");
    laws::report(&entries);
    println!();
    census::report_designs(&entries);
    println!();
    census::report_trend(&entries);
    println!();
    census::report_tight(&entries);
    println!();
    census::report_density(&entries);
    println!();
    census::best_by_k(&entries);
    println!();
    census::summary(&entries);
    println!();
    census::report_dropped(&entries);
    println!();
    census::report_degrees(&entries);
    println!();
    rule_thirty();
    println!();
    typical_report();
    println!();
    json::write(&entries);
    println!();
    let fronts = frontier_report();
    json::write_frontier(&fronts);
}
