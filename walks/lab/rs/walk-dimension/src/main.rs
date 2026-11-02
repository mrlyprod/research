mod design;
mod figure;
mod gasket;
mod spectral;
mod walkers;

use design::{components, plane, sponge, Graph, BASE};
use mrly_core::{Rng, Tensor};
use std::path::Path;

const SUBJECTS: [u128; 7] = [79, 95, 127, 239, 255, 495, 511];

struct Subject {
    name: String,
    fill: usize,
    base: usize,
    df: f64,
    levels: (usize, usize),
    grid: fn(usize) -> Tensor,
    walk_level: usize,
}

fn base3(code: u128) -> Subject {
    let name = format!("bang dim 2, base 3, code {code}");
    let fill = code.count_ones() as usize;
    Subject {
        name,
        fill,
        base: BASE,
        df: (fill as f64).ln() / (BASE as f64).ln(),
        levels: (4, 5),
        grid: match code {
            79 => |level| plane(79, BASE, level),
            95 => |level| plane(95, BASE, level),
            127 => |level| plane(127, BASE, level),
            239 => |level| plane(239, BASE, level),
            255 => |level| plane(255, BASE, level),
            495 => |level| plane(495, BASE, level),
            _ => |level| plane(511, BASE, level),
        },
        walk_level: 5,
    }
}

fn subjects() -> Vec<Subject> {
    let mut out: Vec<Subject> = SUBJECTS.iter().map(|code| base3(*code)).collect();
    out.push(Subject {
        name: "bang dim 2, code 7".to_string(),
        fill: 3,
        base: 2,
        df: 3f64.ln() / 2f64.ln(),
        levels: (6, 7),
        grid: |level| plane(7, 2, level),
        walk_level: 8,
    });
    out.push(Subject {
        name: "bang dim 3, base 3, code 23".to_string(),
        fill: 20,
        base: BASE,
        df: 20f64.ln() / 3f64.ln(),
        levels: (3, 4),
        grid: sponge,
        walk_level: 4,
    });
    out
}

fn giant_graph(grid: &Tensor) -> Graph {
    Graph::of(&components(grid).giant)
}

fn pair(subject: &Subject, coarse: usize, fine: usize) -> (Vec<f64>, usize) {
    let low = giant_graph(&(subject.grid)(coarse));
    let high = giant_graph(&(subject.grid)(fine));
    let nodes = high.nodes();
    (
        spectral::exponents(
            &spectral::low(&low),
            &spectral::low(&high),
            subject.base as f64,
        ),
        nodes,
    )
}

fn list(values: &[f64], places: usize) -> String {
    values
        .iter()
        .map(|value| format!("{value:.places$}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn anchors() {
    println!("ANCHORS, SPECTRAL");
    let low = spectral::low(&gasket::build(7).graph);
    let high = spectral::low(&gasket::build(8).graph);
    let tau: Vec<f64> = (1..=4).map(|mode| low[mode] / high[mode]).collect();
    println!(
        "gasket tau levels 7 to 8, modes 1 to 4: {} (exact 5)",
        list(&tau, 5)
    );
    println!(
        "gasket d_w from tau: {} (exact {:.6})",
        list(
            &tau.iter().map(|t| t.ln() / 2f64.ln()).collect::<Vec<_>>(),
            5
        ),
        5f64.ln() / 2f64.ln()
    );
    for (code, name, levels) in [(511u128, "solid", (4usize, 5usize)), (7, "path", (5, 6))] {
        let mut spectra = Vec::new();
        for level in [levels.0, levels.1] {
            let graph = Graph::of(&plane(code, BASE, level));
            let side = BASE.pow(level as u32) as f64;
            let values = spectral::low(&graph);
            let exact = 2.0 - 2.0 * (std::f64::consts::PI / side).cos();
            println!(
                "{name} level {level}: nodes {} lambda_2 {:.12} closed form {:.12} gap {:.1e}",
                graph.nodes(),
                values[1],
                exact,
                (values[1] - exact).abs()
            );
            spectra.push(values);
        }
        let est = spectral::exponents(&spectra[0], &spectra[1], BASE as f64);
        println!(
            "{name} d_w levels {} to {} from lambda_2: {:.5}",
            levels.0, levels.1, est[0]
        );
    }
}

fn census() -> Vec<u128> {
    println!();
    println!("CLASSES AND SPANNING");
    let group = design::group();
    let classes = design::classes(&group);
    let orbit_sum: usize = classes.iter().map(|class| class.1).sum();
    println!(
        "group order {} classes {} orbit sum {orbit_sum}",
        group.len(),
        classes.len()
    );
    let canonical = |code: u128| *design::orbit(&group, code).iter().next().expect("an orbit");
    let mut spanning: Vec<u128> = Vec::new();
    let mut spanning_all: Vec<u128> = Vec::new();
    let mut single: Vec<u128> = Vec::new();
    let mut strays: Vec<(u128, usize, f64)> = Vec::new();
    for code in 1..512u128 {
        let at5 = components(&plane(code, BASE, 5));
        if at5.spanning {
            spanning.push(code);
            if (1..5).all(|level| components(&plane(code, BASE, level)).spanning) {
                spanning_all.push(code);
            }
            if at5.count == 1 {
                single.push(code);
            } else {
                strays.push((code, at5.count, at5.share));
            }
        }
    }
    let reps_of = |codes: &[u128]| {
        let mut reps: Vec<u128> = codes.iter().map(|code| canonical(*code)).collect();
        reps.sort_unstable();
        reps.dedup();
        reps
    };
    let reps = reps_of(&spanning);
    println!(
        "codes whose giant component touches all four walls at level 5: {} of 511, in {} classes with reps {:?}",
        spanning.len(),
        reps.len(),
        reps
    );
    println!(
        "of those, spanning at every level 1 to 5: {}",
        spanning_all.len()
    );
    let strict = reps_of(&single);
    println!(
        "codes that are one component touching all four walls at level 5: {} of 511, in {} classes with reps {:?}",
        single.len(),
        strict.len(),
        strict
    );
    println!("spanning codes carrying stray components at level 5, with component count and giant share:");
    for (code, count, share) in &strays {
        println!(
            "  code {code} class {} components {count} share {share:.4}",
            canonical(*code)
        );
    }
    for rep in &reps {
        let members: Vec<u128> = spanning
            .iter()
            .copied()
            .filter(|code| canonical(*code) == *rep)
            .collect();
        let size = design::orbit(&group, *rep).len();
        println!(
            "  class {rep} fill {} orbit {size} spanning members {:?}",
            rep.count_ones(),
            members
        );
    }
    println!("representatives: giant share at levels 4 and 5, components at 5, spanning at 5");
    for (rep, size) in &classes {
        if *rep == 0 {
            continue;
        }
        let at4 = components(&plane(*rep, BASE, 4));
        let at5 = components(&plane(*rep, BASE, 5));
        println!(
            "  rep {rep:>3} fill {} orbit {size:>2} share {:.4} {:.4} components {:>5} spanning {}",
            rep.count_ones(),
            at4.share,
            at5.share,
            at5.count,
            at5.spanning
        );
    }
    spanning
}

fn spectral_census(subjects: &[Subject]) -> Vec<f64> {
    println!();
    println!(
        "SPECTRAL CENSUS, d_w FROM LAMBDA_2 AT THE UPPER LEVEL PAIR, DRIFT FROM THE PAIR BELOW"
    );
    let mut out = Vec::new();
    for subject in subjects {
        let (coarse, fine) = subject.levels;
        let (below, _) = pair(subject, coarse - 1, coarse);
        let (step, nodes) = pair(subject, coarse, fine);
        let dw = step[0];
        let kmin = step.iter().copied().fold(f64::INFINITY, f64::min);
        let kmax = step.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        println!(
            "{:>16} fill {:>2} d_f {:.4} nodes {nodes:>6} levels {coarse} to {fine} modes [{}] below [{}] d_w {dw:.4} spread {:.4} drift {:.4} d_s {:.4}",
            subject.name,
            subject.fill,
            subject.df,
            list(&step, 4),
            list(&below, 4),
            kmax - kmin,
            (dw - below[0]).abs(),
            2.0 * subject.df / dw
        );
        out.push(dw);
    }
    out
}

fn anomaly() {
    println!();
    println!("CODE 127, THE SLOW MODE AT LEVEL 4");
    let grid = plane(127, BASE, 4);
    let graph = giant_graph(&grid);
    let mode = spectral::slow_mode(&graph);
    let block = 27usize;
    let mut inside = 0.0;
    let mut total = 0.0;
    let (mut row_sum, mut col_sum) = (0.0, 0.0);
    for (node, flat) in graph.cells.iter().enumerate() {
        let at = design::coords(*flat, &graph.shape);
        let weight = mode[node] * mode[node];
        total += weight;
        row_sum += weight * at[0] as f64;
        col_sum += weight * at[1] as f64;
        if at[0] >= 2 * block && at[1] < block {
            inside += weight;
        }
    }
    println!(
        "weight of mode 1 inside the block under tile cell (2,0), rows 54 to 80 cols 0 to 26: {:.4}; centroid row {:.1} col {:.1} of 81",
        inside / total,
        row_sum / total,
        col_sum / total
    );
}

fn walker_census(subjects: &[Subject]) -> Vec<f64> {
    println!();
    println!(
        "ANCHORS, WALKERS, {} BLIND ANTS, SEED {}",
        walkers::WALKERS,
        walkers::SEED
    );
    let mut rng = Rng::new(walkers::SEED);
    let solid = plane(511, BASE, 5);
    let side = BASE.pow(5) as f64;
    let trace = walkers::grid_walk(&solid, &mut rng, walkers::WALKERS, 5.0);
    let (dw, drift) = walkers::fit(&trace, (side / 6.0).powi(2));
    println!("solid level 5 d_w {dw:.4} drift {drift:.4} (exact 2)");
    let path = plane(7, BASE, 6);
    let side = BASE.pow(6) as f64;
    let trace = walkers::grid_walk(&path, &mut rng, walkers::WALKERS, 6.0);
    let (dw, drift) = walkers::fit(&trace, (side / 6.0).powi(2));
    println!("path level 6 d_w {dw:.4} drift {drift:.4} (exact 2)");
    let gasket = gasket::build(8);
    let trace = walkers::gasket_walk(&gasket, &mut rng, walkers::WALKERS);
    let (dw, drift) = walkers::fit(&trace, (256.0f64 / 6.0).powi(2));
    let exact = 5f64.ln() / 2f64.ln();
    println!(
        "gasket level 8 d_w {dw:.4} drift {drift:.4} (exact {exact:.4}, off by {:.2}%)",
        (dw - exact).abs() / exact * 100.0
    );
    println!();
    println!("WALKER CENSUS, d_w FROM THE DISPLACEMENT FIT, DRIFT BETWEEN THE WINDOW HALVES");
    let mut out = Vec::new();
    for subject in subjects {
        let giant = components(&(subject.grid)(subject.walk_level)).giant;
        let nodes = giant.size() - giant.count(0);
        let narrow = giant.shape.iter().min().copied().unwrap_or(0) as f64;
        let trace = walkers::grid_walk(&giant, &mut rng, walkers::WALKERS, 5.0);
        let (dw, drift) = walkers::fit(&trace, (narrow / 6.0).powi(2));
        println!(
            "{:>16} level {} nodes {nodes:>6} d_w {dw:.4} drift {drift:.4} d_s {:.4}",
            subject.name,
            subject.walk_level,
            2.0 * subject.df / dw
        );
        out.push(dw);
    }
    out
}

fn main() {
    let subjects = subjects();
    anchors();
    census();
    let spectral = spectral_census(&subjects);
    anomaly();
    let walker = walker_census(&subjects);
    println!();
    println!("CENSUS TABLE: design, fill, d_f, d_w spectral, d_w walkers, d_s from the walkers");
    for (index, subject) in subjects.iter().enumerate() {
        println!(
            "{:>16} {:>2} {:.4} {:.3} {:.3} {:.2}",
            subject.name,
            subject.fill,
            subject.df,
            spectral[index],
            walker[index],
            2.0 * subject.df / walker[index]
        );
    }
    println!(
        "d_w at or above 2 on every subject: {}",
        spectral.iter().chain(&walker).all(|dw| *dw >= 2.0 - 0.01)
    );
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = here.ancestors().nth(5).expect("the root");
    let folder = root
        .join("data")
        .join(here.strip_prefix(root).expect("under the root"));
    std::fs::create_dir_all(&folder).expect("the data folder");
    let figure = folder.join("walks-fig.png");
    figure::write(
        &figure,
        &figure::Series {
            spectral,
            walker,
            fractal: subjects.iter().map(|subject| subject.df).collect(),
        },
    );
    println!(
        "wrote {}",
        figure.strip_prefix(root).expect("under the root").display()
    );
}
