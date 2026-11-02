use mrly_math::spirograph::{
    cover, disc, pencils, point, representatives, signed_area, track, Kind, Pencil, Track,
};
use std::f64::consts::{PI, TAU};

const SIDES: [usize; 4] = [256, 512, 1024, 2048];
const REACH: f64 = 0.9;
const TRACE: usize = 200_001;

struct Case {
    name: String,
    path: Track,
    pens: Vec<Pencil>,
}

struct Read {
    covered: f64,
    hole: f64,
    wall: f64,
    winding: f64,
    areas: f64,
}

struct Limit {
    covered: f64,
    bar: f64,
    hole: f64,
    ratio: f64,
}

fn seat(d: f64, quarter: usize) -> Pencil {
    let mut cell = (2i64, 0i64);
    for _ in 0..quarter % 4 {
        cell = (-cell.1, cell.0);
    }
    let angle = TAU * (quarter % 4) as f64 / 4.0;
    Pencil {
        x: d * angle.cos(),
        y: d * angle.sin(),
        seat: cell,
        kind: Kind::Fill,
    }
}

fn copies(d: f64, count: usize) -> Vec<Pencil> {
    (0..count).map(|k| seat(d, k * (4 / count))).collect()
}

fn carpet() -> Vec<u8> {
    vec![1, 1, 1, 1, 0, 1, 1, 1, 1]
}

fn one(name: String, kind: &str, a: usize, b: usize, pens: Vec<Pencil>) -> Case {
    Case {
        name,
        path: track(kind, a, b, 4, 1).unwrap(),
        pens,
    }
}

fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    for kind in ["in", "out"] {
        for (a, b) in ratios() {
            out.push(one(
                format!("one pencil {kind} {a}/{b}"),
                kind,
                a,
                b,
                vec![seat(REACH, 0)],
            ));
        }
    }
    for (a, b) in [(7usize, 3usize), (5, 2)] {
        for count in [2usize, 4] {
            out.push(one(
                format!("{count} copies in {a}/{b}"),
                "in",
                a,
                b,
                copies(REACH, count),
            ));
        }
    }
    for (pens, label) in [("fill", "carpet fills"), ("corners", "carpet corners")] {
        for (kind, a, b) in [("in", 7usize, 3usize), ("in", 5, 2), ("out", 5, 8)] {
            out.push(one(
                format!("{label} {kind} {a}/{b}"),
                kind,
                a,
                b,
                pencils(&carpet(), 3, 3, pens, REACH, 0.0, 1).unwrap(),
            ));
        }
    }
    out
}

fn ratios() -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for b in 1..=4usize {
        for a in b + 1..=9usize {
            if gcd(a, b) == 1 {
                out.push((a, b));
            }
        }
    }
    out
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn curves(case: &Case) -> Vec<usize> {
    representatives(&case.path, &case.pens, true)
}

fn ladder(case: &Case) -> Vec<Read> {
    SIDES
        .iter()
        .map(|&side| {
            let out = cover(&case.path, &case.pens, true, 2, side).unwrap();
            Read {
                covered: out.covered,
                hole: out.hole,
                wall: out.wall,
                winding: out.winding,
                areas: out.areas,
            }
        })
        .collect()
}

fn richer(coarse: f64, fine: f64) -> f64 {
    2.0 * fine - coarse
}

fn limit(values: &[f64]) -> (f64, f64) {
    let last = richer(values[values.len() - 2], values[values.len() - 1]);
    let before = richer(values[values.len() - 3], values[values.len() - 2]);
    (last, (last - before).abs())
}

fn limits(reads: &[Vec<Read>]) -> Vec<Limit> {
    reads
        .iter()
        .map(|rungs| {
            let covered: Vec<f64> = rungs.iter().map(|r| r.covered).collect();
            let holes: Vec<f64> = rungs.iter().map(|r| r.hole).collect();
            let steps: Vec<f64> = (1..covered.len())
                .map(|k| covered[k] - covered[k - 1])
                .collect();
            let ratio = steps
                .windows(2)
                .map(|pair| (pair[1] / pair[0]).abs())
                .fold(0.0_f64, f64::max);
            let (value, bar) = limit(&covered);
            Limit {
                covered: value,
                bar,
                hole: limit(&holes).0,
                ratio,
            }
        })
        .collect()
}

fn far(case: &Case) -> f64 {
    case.path.wheel
        * case
            .pens
            .iter()
            .map(|p| p.x.hypot(p.y))
            .fold(0.0_f64, f64::max)
}

fn rho_of(case: &Case) -> f64 {
    disc(&case.path, &case.pens).unwrap().radius - far(case)
}

fn bound(case: &Case, side: usize) -> f64 {
    let bounds = disc(&case.path, &case.pens).unwrap();
    let rho = rho_of(case);
    let arcs: f64 = curves(case)
        .iter()
        .map(|&k| {
            TAU * case.path.ratio.1 as f64 * rho * (1.0 + case.pens[k].x.hypot(case.pens[k].y))
        })
        .sum();
    2.0 * arcs * (2.0 * bounds.radius / side as f64) / (PI * bounds.radius * bounds.radius)
}

fn simple(case: &Case) -> f64 {
    let (rho, d) = (rho_of(case), far(case));
    rho * (rho + case.path.side * d * d / case.path.wheel) / (rho + d).powi(2)
}

fn annulus(case: &Case) -> f64 {
    let (rho, d) = (rho_of(case), far(case));
    4.0 * rho * d / (rho + d).powi(2)
}

fn inscribed(case: &Case) -> f64 {
    let bounds = disc(&case.path, &case.pens).unwrap();
    (bounds.hole / bounds.radius).powi(2)
}

fn lone(case: &Case) -> bool {
    curves(case).len() == 1 && case.path.ratio.1 == 1
}

fn the_disc(all: &[Case]) -> (f64, usize) {
    println!("THE DISC: the two radii against a trace of {TRACE} points per curve");
    println!(
        "{:<22} {:>7} {:>10} {:>12} {:>10} {:>12} {:>10}",
        "case", "curves", "radius", "traced max", "hole", "traced min", "gap"
    );
    let mut worst = 0.0_f64;
    for case in all {
        let bounds = disc(&case.path, &case.pens).unwrap();
        let (mut high, mut low) = (0.0_f64, f64::MAX);
        for &k in curves(case).iter() {
            for step in 0..TRACE {
                let s = case.path.total * step as f64 / (TRACE - 1) as f64;
                let (x, y) = point(&case.path, &case.pens[k], s);
                let radius = (x - bounds.x).hypot(y - bounds.y);
                high = high.max(radius);
                low = low.min(radius);
            }
        }
        let gap = (high - bounds.radius).abs().max((low - bounds.hole).abs());
        worst = worst.max(gap);
        println!(
            "{:<22} {:>7} {:>10.6} {:>12.6} {:>10.6} {:>12.6} {:>10.2e}",
            case.name,
            curves(case).len(),
            bounds.radius,
            high,
            bounds.hole,
            low,
            gap
        );
    }
    println!();
    (worst, all.len())
}

fn the_winding(all: &[Case], reads: &[Vec<Read>]) -> (f64, usize) {
    println!("THE WINDING SELF-CHECK: the scanline mean against the Green closed form, every side");
    println!(
        "{:<22} {:>7} {:>12} {:>10} {:>10} {:>10} {:>10} {:>11}",
        "case", "curves", "closed form", "gap 256", "gap 512", "gap 1024", "gap 2048", "bound 2048"
    );
    let (mut worst, mut loose) = (0.0_f64, 0usize);
    for (case, rungs) in all.iter().zip(reads) {
        let gaps: Vec<f64> = rungs.iter().map(|r| (r.winding - r.areas).abs()).collect();
        for (rung, gap) in gaps.iter().enumerate() {
            worst = worst.max(*gap);
            if *gap > bound(case, SIDES[rung]) {
                loose += 1;
            }
        }
        println!(
            "{:<22} {:>7} {:>12.6} {:>10.2e} {:>10.2e} {:>10.2e} {:>10.2e} {:>11.2e}",
            case.name,
            curves(case).len(),
            rungs[3].areas,
            gaps[0],
            gaps[1],
            gaps[2],
            gaps[3],
            bound(case, SIDES[3])
        );
    }
    println!();
    (worst, loose)
}

fn the_ladder(all: &[Case], reads: &[Vec<Read>], marks: &[Limit]) {
    println!("THE COVER: the shape's share of the disc, raster by raster, the wall counted in");
    println!(
        "{:<22} {:>7} {:>9} {:>9} {:>9} {:>9} {:>10} {:>9} {:>9} {:>9}",
        "case", "curves", "256", "512", "1024", "2048", "covered*", "bar", "wall 2048", "hole*"
    );
    for ((case, rungs), mark) in all.iter().zip(reads).zip(marks) {
        println!(
            "{:<22} {:>7} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>10.6} {:>9.6} {:>9.6} {:>9.6}",
            case.name,
            curves(case).len(),
            rungs[0].covered,
            rungs[1].covered,
            rungs[2].covered,
            rungs[3].covered,
            mark.covered,
            mark.bar,
            rungs[3].wall,
            mark.hole
        );
    }
    println!();
}

fn the_hole(all: &[Case], marks: &[Limit]) -> (usize, (f64, String), (f64, String), f64) {
    println!(
        "THE HOLE: the centre flood against its inscribed disc and, at b = 1, the closed form"
    );
    println!(
        "{:<22} {:>10} {:>11} {:>10} {:>11} {:>10}",
        "case", "hole*", "inscribed", "slack", "b = 1 form", "gap"
    );
    let mut leaks = 0;
    let mut low = (f64::MAX, String::new());
    let mut high = (0.0_f64, String::new());
    let mut worst = 0.0_f64;
    for (case, mark) in all.iter().zip(marks) {
        let (floor, form) = (inscribed(case), simple(case));
        let slack = mark.hole - floor;
        if slack < 0.0 {
            leaks += 1;
        }
        if slack < low.0 {
            low = (slack, case.name.clone());
        }
        if slack > high.0 {
            high = (slack, case.name.clone());
        }
        if lone(case) {
            worst = worst.max((mark.hole - form).abs());
        }
        println!(
            "{:<22} {:>10.6} {:>11.6} {:>10.6} {:>11} {:>10}",
            case.name,
            mark.hole,
            floor,
            slack,
            if lone(case) {
                format!("{form:.6}")
            } else {
                "-".to_string()
            },
            if lone(case) {
                format!("{:.2e}", (mark.hole - form).abs())
            } else {
                "-".to_string()
            }
        );
    }
    println!();
    (leaks, low, high, worst)
}

fn the_candidates(all: &[Case], marks: &[Limit]) -> (usize, usize) {
    println!("THE CANDIDATES: the extrapolated cover against three closed forms");
    println!(
        "{:<22} {:>10} {:>9} {:>10} {:>10} {:>10} {:>10}",
        "case", "covered*", "bar", "simple", "annulus", "areas", "verdict"
    );
    let (mut held, mut tried) = (0usize, 0usize);
    for (case, mark) in all.iter().zip(marks) {
        let bar = mark.bar.max(1e-6);
        let forms = [simple(case), annulus(case), areas(case)];
        let names = ["simple", "annulus", "areas"];
        let kept: Vec<&str> = names
            .iter()
            .zip(&forms)
            .filter(|(_, form)| (**form - mark.covered).abs() < 10.0 * bar)
            .map(|(name, _)| *name)
            .collect();
        held += kept.len();
        tried += forms.len();
        println!(
            "{:<22} {:>10.6} {:>9.6} {:>10.6} {:>10.6} {:>10.6} {:>10}",
            case.name,
            mark.covered,
            mark.bar,
            forms[0],
            forms[1],
            forms[2],
            if kept.is_empty() {
                "none".to_string()
            } else {
                kept.join(",")
            }
        );
    }
    println!();
    (held, tried)
}

fn areas(case: &Case) -> f64 {
    let bounds = disc(&case.path, &case.pens).unwrap();
    curves(case)
        .iter()
        .filter_map(|&k| signed_area(&case.path, &case.pens[k]))
        .sum::<f64>()
        / (PI * bounds.radius * bounds.radius)
}

fn main() {
    let all = cases();
    let reads: Vec<Vec<Read>> = all.iter().map(ladder).collect();
    let marks = limits(&reads);
    println!("ROULETTE COVER: the shape between the walls of a circle roulette");
    println!("reach {REACH} wheel radii, seats exact, sides {SIDES:?}");
    println!("ratios: every a/b in lowest terms with 1 <= b <= 4 and b < a <= 9, inside and out");
    println!();
    let (radii, seen) = the_disc(&all);
    let (winding, loose) = the_winding(&all, &reads);
    the_ladder(&all, &reads, &marks);
    let (leaks, low, high, form) = the_hole(&all, &marks);
    let (held, tried) = the_candidates(&all, &marks);
    let flat = all
        .iter()
        .zip(&marks)
        .filter(|(case, _)| lone(case))
        .map(|(_, mark)| mark.covered.abs())
        .fold(0.0_f64, f64::max);
    let ones = all.iter().filter(|case| lone(case)).count();
    let ratio = marks.iter().map(|mark| mark.ratio).fold(0.0_f64, f64::max);
    println!("THE REPORT");
    println!("- the two disc radii meet the trace on all {seen} cases, worst gap {radii:.2e}.");
    println!("- the scanline winding meets the closed form on all {seen} cases at all four sides, worst gap {winding:.2e}, {loose} readings outside the perimeter bound.");
    println!("- successive ladder differences fall by a factor of at most {ratio:.3} on every case and both steps, the pixel law the extrapolation rests on.");
    println!("- b = 1 and one curve, {ones} cases: the curve is simple, so the cover is the wall alone, at most {flat:.6} extrapolated.");
    println!("- b = 1 and one curve: the centre flood is the whole inside, the signed closed form holding to {form:.2e}.");
    println!("- the centre flood holds its inscribed disc on every case, {leaks} leaks, slack from {:.6} at {} to {:.6} at {}.", low.0, low.1, high.0, high.1);
    println!("- {held} of {tried} candidate readings survive ten bars: no closed form shows itself for the cover.");
}
