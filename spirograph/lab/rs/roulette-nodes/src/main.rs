use mrly_math::roulette::{nodes, spread, Nodes};
use mrly_math::spirograph::{pencils, trace, track, Pencil, Track};
use mrly_num::factor::gcd;
use std::f64::consts::TAU;

const TOL: f64 = 4e-4;
const GRAIN: usize = 400;
const FLOOR: usize = 6000;
const CEIL: usize = 24000;
const ROOTS: usize = 2_000_000;
const DEEP: usize = 5;
const GRIT: usize = 12;

// THE DESIGNS

struct Design {
    name: &'static str,
    types: Vec<u8>,
    width: usize,
    height: usize,
}

impl Design {
    fn far(&self) -> f64 {
        let (w, h) = (self.width as i64, self.height as i64);
        let mut far: f64 = 0.0;
        for i in 0..self.height {
            for j in 0..self.width {
                if self.types[i * self.width + j] == 0 {
                    continue;
                }
                let (i, j) = (i as i64, j as i64);
                let seat = ((2 * j + 1 - w) as f64 / 2.0, (h - 2 * i - 1) as f64 / 2.0);
                far = far.max(seat.0.hypot(seat.1));
            }
        }
        far
    }

    fn reach(&self, target: f64) -> f64 {
        target * (self.width as f64 / 2.0).hypot(self.height as f64 / 2.0) / self.far()
    }

    fn fits(&self, target: f64) -> bool {
        (0.05..=2.0).contains(&self.reach(target))
    }

    fn seats(&self, target: f64) -> Vec<Pencil> {
        pencils(
            &self.types,
            self.width,
            self.height,
            "fill",
            self.reach(target),
            0.0,
            1,
        )
        .unwrap()
    }
}

fn designs() -> Vec<Design> {
    vec![
        Design {
            name: "pin",
            types: vec![1, 0, 0],
            width: 3,
            height: 1,
        },
        Design {
            name: "two",
            types: vec![1, 0, 0, 0, 0, 0, 0, 1, 0],
            width: 3,
            height: 3,
        },
        Design {
            name: "plus",
            types: vec![0, 1, 0, 1, 0, 1, 0, 1, 0],
            width: 3,
            height: 3,
        },
        Design {
            name: "carpet",
            types: vec![1, 1, 1, 1, 0, 1, 1, 1, 1],
            width: 3,
            height: 3,
        },
    ]
}

// THE RUN

struct Run {
    curves: usize,
    least: f64,
    most: f64,
    count: Nodes,
}

fn beat(kind: &str, a: usize, b: usize) -> usize {
    if kind == "in" {
        a - b
    } else {
        a + b
    }
}

fn wide(kind: &str, a: usize, b: usize) -> f64 {
    beat(kind, a, b) as f64 / b as f64
}

fn bound(kind: &str, a: usize, b: usize) -> f64 {
    wide(kind, a, b).min(1.0)
}

fn grain(kind: &str, a: usize, b: usize) -> usize {
    (GRAIN * beat(kind, a, b).max(b)).clamp(FLOOR, CEIL)
}

fn lay(kind: &str, a: usize, b: usize) -> Track {
    track(kind, a, b, 4, 1).unwrap()
}

fn sweep(kind: &str, a: usize, b: usize, design: &Design, target: f64, samples: usize) -> Run {
    grip(kind, a, b, design, target, samples, TOL)
}

fn grip(
    kind: &str,
    a: usize,
    b: usize,
    design: &Design,
    target: f64,
    samples: usize,
    tol: f64,
) -> Run {
    let path = lay(kind, a, b);
    let seats = spread(&path, &design.seats(target), true);
    let reaches: Vec<f64> = seats.iter().map(|p| p.x.hypot(p.y)).collect();
    Run {
        curves: seats.len(),
        least: reaches.iter().cloned().fold(f64::MAX, f64::min),
        most: reaches.iter().cloned().fold(0.0, f64::max),
        count: nodes(&path, &seats, samples, tol).unwrap(),
    }
}

fn law(a: usize, b: usize, k: usize) -> usize {
    a * b * k * (k - 1) + k * a * (b - 1)
}

fn ratios(kind: &str, top: usize, deep: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for b in 1..=deep {
        for a in 1..=top {
            if gcd(a as u128, b as u128) == 1 && (kind == "out" || a > b) {
                out.push((a, b));
            }
        }
    }
    out
}

// THE ROOTS

fn roots(wide: f64, b: usize, m: usize, cp: f64, cq: f64, phase: f64, steps: usize) -> usize {
    let gap = (cp - cq).abs();
    let mean = 4.0 * cp * cq;
    let mut out = 0;
    for sign in [1.0, -1.0] {
        let f = |k: usize| {
            let x = TAU * (k as f64 + 0.5) / steps as f64;
            let ripple = (m as f64 * x + phase).sin();
            2.0 * wide * (b as f64 * x).sin() - sign * (gap * gap + mean * ripple * ripple).sqrt()
        };
        let mut prev = f(steps - 1);
        for k in 0..steps {
            let now = f(k);
            if (prev < 0.0) != (now < 0.0) {
                out += 1;
            }
            prev = now;
        }
    }
    out
}

fn torus(
    kind: &str,
    a: usize,
    b: usize,
    p: (f64, f64),
    q: (f64, f64),
    steps: usize,
) -> (usize, usize) {
    let (cp, cq) = (p.0.hypot(p.1), q.0.hypot(q.1));
    let phase = (p.1.atan2(p.0) - q.1.atan2(q.0)) / 2.0;
    let count = roots(wide(kind, a, b), b, beat(kind, a, b), cp, cq, phase, steps);
    (count, if cp == cq && phase == 0.0 { 4 } else { 0 })
}

fn main() {
    let want = std::env::args().nth(1);
    let pick = |name: &str| match want.as_deref() {
        Some(verb) => verb == name,
        None => true,
    };
    if pick("scratch") {
        scratch();
    }
    if pick("law") {
        census();
    }
    if pick("gap") {
        gap();
    }
    if pick("torus") {
        reduction();
    }
    if pick("centre") {
        centre();
    }
    if pick("regions") {
        regions();
    }
}

// THE SCRATCH

fn at_reach(design: &Design, reach: f64) -> f64 {
    let circum = (design.width as f64 / 2.0).hypot(design.height as f64 / 2.0);
    reach * design.far() / circum
}

fn scratch() {
    let all = designs();
    let named = |name: &str| all.iter().position(|d| d.name == name).unwrap();
    let (pin, two, carpet) = (
        &all[named("pin")],
        &all[named("two")],
        &all[named("carpet")],
    );
    println!("THE SCRATCH");
    let mut row = Vec::new();
    for a in [3, 4, 5, 6, 7] {
        let run = sweep("in", a, 1, carpet, at_reach(carpet, 0.3), FLOOR);
        row.push(format!(
            "{a}/1 k{} {} {}",
            run.curves,
            run.count.paired(),
            run.count.selved()
        ));
    }
    println!(
        "- carpet fills inside at reach 0.30, ratio a/1: {}",
        row.join("  ")
    );
    println!("  against 2ab C(k,2) = 56a and k a(b - 1) = 0");
    let mut row = Vec::new();
    for reach in [0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.2] {
        let target = at_reach(carpet, reach);
        let run = sweep("in", 7, 3, carpet, target, grain("in", 7, 3));
        row.push(format!(
            "{reach:.1}|p|{:.2}:{}+{}{}",
            run.most,
            run.count.paired(),
            run.count.selved(),
            if run.count.crowded > 0 { "*" } else { "" }
        ));
    }
    println!(
        "- carpet fills inside 7/3 by reach, pairs + selves: {}",
        row.join("  ")
    );
    println!("  against 2ab C(8,2) + 8 a(b - 1) = 1176 + 112 = 1288 below the threshold");
    for target in [0.5, 0.9] {
        let mut row = Vec::new();
        for (a, b) in [
            (3, 1),
            (4, 1),
            (5, 1),
            (7, 1),
            (5, 2),
            (7, 2),
            (7, 3),
            (8, 3),
            (11, 4),
        ] {
            let run = sweep("in", a, b, pin, target, grain("in", a, b));
            row.push(format!("{a}/{b}:{}", run.count.selved()));
        }
        println!("- one pencil inside at |p| = {target}: {}", row.join("  "));
    }
    for target in [0.5, 0.9] {
        let mut row = Vec::new();
        for (a, b) in [(3, 1), (4, 1), (5, 1), (7, 1), (5, 2), (7, 3)] {
            let run = sweep("out", a, b, pin, target, grain("out", a, b));
            row.push(format!("{a}/{b}:{}", run.count.selved()));
        }
        println!("- one pencil outside at |p| = {target}: {}", row.join("  "));
    }
    let run = sweep("in", 7, 4, two, at_reach(two, 0.9), grain("in", 7, 4));
    println!(
        "- two carpet curves inside 7/4 at reach 0.90: k {} pairs {} selves {:?} against 56 and 21",
        run.curves,
        run.count.paired(),
        run.count.selves
    );
    let mut row = Vec::new();
    for reach in [0.76, 0.78, 0.79, 0.80, 0.81, 0.84] {
        let run = sweep(
            "in",
            7,
            3,
            carpet,
            at_reach(carpet, reach),
            grain("in", 7, 3),
        );
        row.push(format!(
            "{reach:.2}:{}+{} most {} crowded {}",
            run.count.paired(),
            run.count.selved(),
            run.count.most,
            run.count.crowded
        ));
    }
    println!(
        "- carpet fills inside 7/3 near the alignment: {}",
        row.join("  ")
    );
}

// THE LAW

fn faults(run: &Run, a: usize, b: usize) -> Option<String> {
    let (want, pair) = (a * (b - 1), 2 * a * b);
    let stray = run.count.selves.iter().any(|&count| count != want)
        || run.count.pairs.iter().any(|&count| count != pair);
    if !stray {
        return None;
    }
    Some(format!(
        "selves {:?} against {want}, pairs {:?} against {pair}",
        run.count.selves, run.count.pairs
    ))
}

fn branches(most: usize) -> usize {
    let mut n = 1;
    while n * (n - 1) / 2 < most {
        n += 1;
    }
    n
}

fn settle(
    kind: &str,
    a: usize,
    b: usize,
    design: &Design,
    target: f64,
    samples: usize,
) -> (Run, usize, bool) {
    let mut steps = samples;
    let mut run = sweep(kind, a, b, design, target, steps);
    let mut seen = vec![run.count.total()];
    for turn in 1..=DEEP {
        steps = 2 * steps - 1;
        run = sweep(kind, a, b, design, target, steps);
        seen.push(run.count.total());
        let last = seen.len() - 1;
        if last >= 2 && seen[last] == seen[last - 1] && seen[last - 1] == seen[last - 2] {
            return (run, turn, true);
        }
    }
    (run, DEEP, false)
}

fn census() {
    println!("THE LAW");
    let all = designs();
    let targets = [0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95];
    let (mut cells, mut counted, mut shown) = (0, 0, 0);
    let (mut broken, mut crowded, mut shaky, mut near, mut wobbly) = (0, 0, 0, 0, 0);
    let mut blurred = 0;
    for kind in ["in", "out"] {
        for (a, b) in ratios(kind, 20, 10) {
            let edge = bound(kind, a, b);
            for design in &all {
                for target in targets {
                    if target >= edge || !design.fits(target) {
                        continue;
                    }
                    let samples = grain(kind, a, b);
                    let (run, turns, steady) = settle(kind, a, b, design, target, samples);
                    cells += 1;
                    counted += run.count.total();
                    shaky += usize::from(turns > 2);
                    wobbly += usize::from(!steady);
                    let head = format!(
                        "{kind} {a}/{b} {} |p| {:.3} k {}",
                        design.name, run.most, run.curves
                    );
                    if !steady {
                        println!(
                            "- wobbly {head}: {} nodes over {DEEP} doublings, never three counts alike",
                            run.count.total()
                        );
                    }
                    let mut aligned = false;
                    if run.count.crowded > 0 {
                        let tight = grip(kind, a, b, design, target, samples, TOL / 10.0);
                        aligned = tight.count.crowded > 0;
                        if aligned {
                            crowded += 1;
                            if shown < 8 {
                                shown += 1;
                                println!(
                                    "- crowded {head}: {} of {} nodes take {} branches",
                                    tight.count.crowded,
                                    tight.count.total(),
                                    branches(tight.count.most)
                                );
                            }
                        } else {
                            near += 1;
                        }
                    }
                    let fault = faults(&run, a, b);
                    if fault.is_some() || run.count.total() != law(a, b, run.curves) {
                        let fault = fault.unwrap_or_else(|| "the total misses".to_string());
                        if aligned {
                            blurred += 1;
                            println!(
                                "- blurred {head}: {fault}, {} nodes take {} branches",
                                run.count.crowded,
                                branches(run.count.most)
                            );
                        } else {
                            broken += 1;
                            println!(
                                "- broken {head}: {fault}, settled at {} against {} at the base {} samples and the law {}",
                                run.count.total(),
                                grip(kind, a, b, design, target, samples, TOL).count.total(),
                                samples,
                                law(a, b, run.curves)
                            );
                        }
                    }
                }
            }
        }
    }
    println!("- cells {cells}, nodes {counted}, broken {broken}, crowded {crowded} of which {blurred} disagree and are excluded, near {near}");
    println!("- every count is three sample counts alike, the last two doubled: {shaky} cells needed a further doubling, {wobbly} never settled");
}

// THE GAP

fn crest(b: usize, m: usize) -> f64 {
    let steps = 400_000;
    let ridge = |k: usize| {
        let x = std::f64::consts::PI * k as f64 / steps as f64;
        (b as f64 * x).sin().abs() / (m as f64 * x).sin().abs()
    };
    let (mut back, mut here) = (ridge(1), ridge(2));
    let mut least = f64::MAX;
    for k in 3..steps {
        let next = ridge(k);
        if here > back && here >= next {
            least = least.min(here);
        }
        (back, here) = (here, next);
    }
    least
}

fn gap() {
    println!("THE GAP");
    let all = designs();
    let named = |name: &str| all.iter().position(|d| d.name == name).unwrap();
    let (pin, two) = (&all[named("pin")], &all[named("two")]);
    let levels = [0.5, 0.9, 1.02, 1.1, 1.3, 1.6, 2.0, 2.6, 3.4, 4.5, 6.0];
    let (mut cells, mut stray, mut broken, mut bumps, mut tails) = (0, 0, 0, 0, 0);
    let mut ladders = Vec::new();
    for (a, b) in ratios("in", 20, 10) {
        let (m, edge) = (a - b, wide("in", a, b));
        if edge >= 1.0 {
            continue;
        }
        let top = crest(b, m).min(b as f64 / m as f64);
        let mut rungs: Vec<(f64, usize)> = Vec::new();
        for level in levels {
            let target = level * edge;
            if target >= 1.0 || !pin.fits(target) {
                continue;
            }
            let (run, ..) = settle("in", a, b, pin, target, grain("in", a, b));
            let count = run.count.selved();
            cells += 1;
            stray += usize::from(count % a != 0);
            if (level < top) != (count == a * (b - 1)) {
                broken += 1;
            }
            if let Some(&(_, last)) = rungs.last() {
                bumps += usize::from(count > last);
            }
            rungs.push((level, count));
        }
        if let Some(&(level, count)) = rungs.last() {
            tails += usize::from(level > 1.0 && count == a * m);
        }
        let row: Vec<String> = rungs
            .iter()
            .map(|(level, count)| format!("{level}:{}a", count / a))
            .collect();
        ladders.push(format!(
            "{a}/{b} b - 1 = {} m = {m} crest {top:.3}: {}",
            b - 1,
            row.join(" ")
        ));
    }
    println!("- cells {cells}, counts not a multiple of a {stray}, cells the crest misreads {broken}, rungs that rise {bumps}");
    println!("- the ladder by C/A, the seat over the centre path, one pencil inside:");
    for line in ladders.iter().take(6) {
        println!("  {line}");
    }
    println!("  ladders {}, ending at a(a - b) {tails}", ladders.len());
    let mut strays = Vec::new();
    for (a, b) in ratios("in", 20, 10) {
        let edge = wide("in", a, b);
        if edge >= 1.0 {
            continue;
        }
        for level in [1.1, 1.6, 2.4] {
            let target = level * edge;
            if target >= 1.0 || !two.fits(target) {
                continue;
            }
            let (run, ..) = settle("in", a, b, two, target, grain("in", a, b));
            if run.least <= edge || run.curves < 2 {
                continue;
            }
            let seen: Vec<String> = run.count.pairs.iter().map(|c| c.to_string()).collect();
            if run.count.pairs.iter().any(|&c| c != 2 * a * b) {
                strays.push(format!(
                    "{a}/{b} |p| {:.2} to {:.2}: {} against 2ab {}",
                    run.least,
                    run.most,
                    seen.join(" "),
                    2 * a * b
                ));
            }
        }
    }
    let seats = two.seats(0.9);
    let path = lay("in", 7, 4);
    let kept = spread(&path, &seats, true);
    let hit = nodes(&path, &kept, grain("in", 7, 4), TOL).unwrap();
    println!(
        "- one seat past the edge is enough: in 7/4 the two design at |p| {:.3} and {:.3}, the edge at 0.750 and the crest at {:.3}, reads pair {} against 2ab 56 while both selves hold at {:?} against a(b - 1) 21",
        kept[0].x.hypot(kept[0].y),
        kept[1].x.hypot(kept[1].y),
        crest(4, 3).min(4.0 / 3.0),
        hit.paired(),
        hit.selves
    );
    println!("- two curves past the edge do not read 2ab, and read no one number a ratio:");
    for line in strays.iter().take(5) {
        println!("  {line}");
    }
    println!("  cells with a stray pair count {}", strays.len());
    println!("- at C/A = 1 exactly the curve runs through the centre, a branches meeting there, so the reach below is neither the law nor a multiple of a and no census cell sits on it:");
    for (a, b) in [(7, 5), (7, 6), (9, 5)] {
        let edge = wide("in", a, b);
        let mut row = Vec::new();
        for step in 0..9 {
            let target = edge * (0.8 + 0.05 * step as f64);
            let (run, ..) = settle("in", a, b, pin, target, grain("in", a, b));
            row.push(format!("{:.3}:{}", run.most, run.count.selved()));
        }
        println!(
            "- across the edge at {a}/{b}, (a - b)/b = {edge:.3}, a(b - 1) = {}: {}",
            a * (b - 1),
            row.join(" ")
        );
    }
}

// THE TORUS

fn reduction() {
    println!("THE TORUS");
    let all = designs();
    for (kind, a, b, name, target) in [
        ("in", 7, 3, "pin", 0.5),
        ("in", 7, 3, "pin", 0.9),
        ("in", 11, 4, "two", 0.5),
        ("in", 7, 5, "two", 0.3),
        ("in", 7, 5, "pin", 0.9),
        ("in", 7, 6, "pin", 0.9),
        ("out", 5, 2, "two", 0.9),
        ("out", 7, 3, "carpet", 0.4),
        ("out", 1, 10, "carpet", 0.1),
        ("in", 7, 4, "two", 0.6),
        ("in", 20, 3, "carpet", 0.25),
        ("out", 10, 9, "carpet", 0.4),
    ] {
        let design = all.iter().find(|d| d.name == name).unwrap();
        let path = lay(kind, a, b);
        let seats = spread(&path, &design.seats(target), true);
        let count = nodes(&path, &seats, grain(kind, a, b), TOL).unwrap();
        let mut reads = [0, 0];
        for (turn, steps) in [ROOTS, 4 * ROOTS].into_iter().enumerate() {
            for (i, one) in seats.iter().enumerate() {
                for (j, other) in seats.iter().enumerate().skip(i) {
                    let (r, drop) = torus(kind, a, b, (one.x, one.y), (other.x, other.y), steps);
                    reads[turn] += if i == j {
                        a * (r - drop) / 4
                    } else {
                        a * r / 2
                    };
                }
            }
        }
        println!(
            "- {kind} {a}/{b} {name} |p| {target} k {}: polyline {} against torus {} at {ROOTS} roots and {} at {}, and the law {}",
            seats.len(),
            count.total(),
            reads[0],
            reads[1],
            4 * ROOTS,
            law(a, b, seats.len())
        );
    }
}

// THE CENTRE

fn core() -> Design {
    Design {
        name: "core",
        types: vec![1, 0, 0, 0, 1, 0, 0, 0, 0],
        width: 3,
        height: 3,
    }
}

fn centre() {
    println!("THE CENTRE");
    let hub = core();
    for (a, b) in [(3, 1), (5, 2), (7, 3)] {
        let path = lay("in", a, b);
        let seats = spread(&path, &hub.seats(0.5), true);
        let samples = grain("in", a, b);
        let count = nodes(&path, &seats, samples, TOL).unwrap();
        let coarse = nodes(&path, &seats, samples / 2, TOL).unwrap();
        println!(
            "- in {a}/{b}, a seat at |p| 0.5 beside a seat at the wheel's centre: the pair reads {} crossings against 2ab {}, and the flood reads {} regions at 1600 and {} at 2400, against 2a + a(b - 1) + 2 = {} and the law's 2ab + a(b - 1) + 2 = {}",
            count.paired(),
            2 * a * b,
            flood(&path, &seats, 24000, 1600),
            flood(&path, &seats, 24000, 2400),
            2 * a + a * (b - 1) + 2,
            2 * a * b + a * (b - 1) + 2
        );
        println!(
            "  the centre curve is one circle traced b times, so its own polyline self count is meaningless: {} at {samples} samples against {} at {}, and the seat off centre keeps its {} against a(b - 1) {}",
            count.selves[1],
            coarse.selves[1],
            samples / 2,
            count.selves[0],
            a * (b - 1)
        );
    }
}

// THE REGIONS

fn flood(path: &Track, seats: &[Pencil], samples: usize, side: usize) -> usize {
    let trail = trace(path, seats, samples).unwrap();
    let (mut low, mut high) = ([f64::MAX; 2], [f64::MIN; 2]);
    for pair in trail.chunks_exact(2) {
        for k in 0..2 {
            low[k] = low[k].min(f64::from(pair[k]));
            high[k] = high[k].max(f64::from(pair[k]));
        }
    }
    let span = (high[0] - low[0]).max(high[1] - low[1]);
    let place = |p: [f64; 2]| {
        let step = |v: f64, k: usize| {
            (((v - low[k]) / span * (side - 5) as f64) as usize + 2).min(side - 1)
        };
        (step(p[0], 0), step(p[1], 1))
    };
    let mut wall = vec![false; side * side];
    for k in 0..seats.len() {
        for i in 0..samples - 1 {
            let base = (k * samples + i) * 2;
            let one = [f64::from(trail[base]), f64::from(trail[base + 1])];
            let two = [f64::from(trail[base + 2]), f64::from(trail[base + 3])];
            let (from, to) = (place(one), place(two));
            let far = from.0.abs_diff(to.0).max(from.1.abs_diff(to.1));
            for step in 0..=2 * far {
                let share = step as f64 / (2 * far).max(1) as f64;
                let point = [
                    one[0] + (two[0] - one[0]) * share,
                    one[1] + (two[1] - one[1]) * share,
                ];
                let (x, y) = place(point);
                wall[x * side + y] = true;
            }
        }
    }
    let mut seen = vec![false; side * side];
    let mut faces = 0;
    for start in 0..side * side {
        if wall[start] || seen[start] {
            continue;
        }
        let mut stack = vec![start];
        seen[start] = true;
        let mut size = 0;
        while let Some(cell) = stack.pop() {
            size += 1;
            let (x, y) = (cell / side, cell % side);
            for (dx, dy) in [(1, 0), (side - 1, 0), (0, 1), (0, side - 1)] {
                let (nx, ny) = ((x + dx) % side, (y + dy) % side);
                let next = nx * side + ny;
                if !wall[next] && !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        faces += usize::from(size >= GRIT);
    }
    faces
}

fn regions() {
    println!("THE REGIONS");
    let all = designs();
    let named = |name: &str| all.iter().position(|d| d.name == name).unwrap();
    let (pin, two, carpet) = (
        &all[named("pin")],
        &all[named("two")],
        &all[named("carpet")],
    );
    for (kind, a, b, design, target, note) in [
        ("in", 3, 1, pin, 0.5, "one seat, no crossing, Jordan"),
        ("in", 5, 2, pin, 0.5, "one seat"),
        ("in", 7, 3, pin, 0.5, "one seat"),
        ("in", 7, 3, pin, 0.9, "one seat"),
        ("in", 3, 1, two, 0.5, "two seats"),
        ("in", 5, 2, two, 0.5, "two seats"),
        ("in", 7, 4, two, 0.9, "two seats past the seat threshold"),
        (
            "in",
            7,
            3,
            carpet,
            0.527,
            "eight seats at an alignment reach",
        ),
    ] {
        let path = lay(kind, a, b);
        let seats = spread(&path, &design.seats(target), true);
        let count = nodes(&path, &seats, grain(kind, a, b), TOL).unwrap();
        println!(
            "- {kind} {a}/{b} {} at |p| {target}, {note}: crossings {} points {} branches {} so E - V + 2 = {}, flood {} at 1600 and {} at 2400, the law {}",
            design.name,
            count.total(),
            count.points,
            count.branches,
            count.branches - count.points + 2,
            flood(&path, &seats, 24000, 1600),
            flood(&path, &seats, 24000, 2400),
            law(a, b, seats.len())
        );
    }
}
