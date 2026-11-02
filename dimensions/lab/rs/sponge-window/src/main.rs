use mrly_math::three::sponge::{deep, dimension, distance, exact, profile, COVER, EDGE};
use std::time::Instant;

const PLUS: f64 = 7.0 / 27.0;
const SIDE: usize = 40;
const RUN: usize = 800;
const TOP: usize = 4;
const COARSE: usize = 72;
const CENTRE: usize = 120;
const STEPS: usize = 3000;
const DEPTH: i32 = 40;

// CROSSING

struct Hole {
    level: usize,
    start: f64,
    side: f64,
    rows: Vec<f64>,
}

fn crossing(line: f64, top: usize) -> Vec<Hole> {
    let (mut local, mut start, mut rows) = (3.0 * line, 0.0, vec![0.0]);
    let mut out = Vec::new();
    for level in 1..=top {
        let side = 3f64.powi(-(level as i32) - 1);
        let digit = (3.0 * local).floor().clamp(0.0, 2.0);
        local = 3.0 * local - digit;
        if digit == 1.0 {
            let here = rows.iter().map(|q| q + side).collect();
            out.push(Hole {
                level,
                start: start + side,
                side,
                rows: here,
            });
        }
        let choices: &[f64] = if digit == 1.0 {
            &[0.0, 2.0]
        } else {
            &[0.0, 1.0, 2.0]
        };
        rows = rows
            .iter()
            .flat_map(|q| choices.iter().map(move |b| q + b * side))
            .collect();
        start += digit * side;
    }
    out
}

fn integrand(a: f64, radius: f64, hole: &Hole, line: f64) -> f64 {
    let (half, reach) = (hole.side / 2.0, line - hole.start);
    let fall = (radius * radius - (hole.start + a).powi(2)).max(0.0).sqrt();
    let width = reach.min(hole.side - fall) - a.max(fall);
    if fall >= half || width <= 0.0 {
        return 0.0;
    }
    4.0 * (half - fall) * width
}

fn simpson(
    f: &dyn Fn(f64) -> f64,
    (a, fa): (f64, f64),
    (b, fb): (f64, f64),
    fm: f64,
    whole: f64,
    depth: usize,
) -> f64 {
    let m = (a + b) / 2.0;
    let (lm, rm) = ((a + m) / 2.0, (m + b) / 2.0);
    let (flm, frm) = (f(lm), f(rm));
    let left = (m - a) * (fa + 4.0 * flm + fm) / 6.0;
    let right = (b - m) * (fm + 4.0 * frm + fb) / 6.0;
    if depth == 0 || (left + right - whole).abs() <= 1e-22 {
        return left + right + (left + right - whole) / 15.0;
    }
    simpson(f, (a, fa), (m, fm), flm, left, depth - 1)
        + simpson(f, (m, fm), (b, fb), frm, right, depth - 1)
}

fn quadrature(radius: f64, hole: &Hole, line: f64) -> f64 {
    let f = |a: f64| integrand(a, radius, hole, line);
    let reach = line - hole.start;
    let low = ((radius * radius - hole.side * hole.side / 4.0)
        .max(0.0)
        .sqrt()
        - hole.start)
        .clamp(0.0, reach);
    let cuts = 64;
    let mut total = 0.0;
    for k in 0..cuts {
        let a = low + (reach - low) * k as f64 / cuts as f64;
        let b = low + (reach - low) * (k + 1) as f64 / cuts as f64;
        let (fa, fb, fm) = (f(a), f(b), f((a + b) / 2.0));
        total += simpson(
            &f,
            (a, fa),
            (b, fb),
            fm,
            (b - a) * (fa + 4.0 * fm + fb) / 6.0,
            40,
        );
    }
    total
}

// RASTER

fn box_low(radius: f64, hole: &Hole, line: f64) -> f64 {
    (radius * radius - hole.side * hole.side / 4.0)
        .max(0.0)
        .sqrt()
        .max(hole.start)
        .min(line)
}

fn rows_raster(radius: f64, hole: &Hole, line: f64) -> [f64; 3] {
    let low = box_low(radius, hole, line);
    let (step, pace) = ((line - low) / SIDE as f64, hole.side / RUN as f64);
    let half = (2.0 * step * step + pace * pace).sqrt() / 2.0;
    let (mut mid, mut inner, mut outer) = (0.0, 0.0, 0.0);
    for row in &hole.rows {
        for a in 0..SIDE {
            for b in 0..SIDE {
                for c in 0..RUN {
                    let u = 1.0 / 3.0 + low + (a as f64 + 0.5) * step;
                    let v = 1.0 / 3.0 + low + (b as f64 + 0.5) * step;
                    let d = distance([row + (c as f64 + 0.5) * pace, u, v]);
                    mid += f64::from(u8::from(d > radius));
                    inner += f64::from(u8::from(d > radius + half));
                    outer += f64::from(u8::from(d > radius - half));
                }
            }
        }
    }
    let cell = step * step * pace;
    [mid * cell, inner * cell, outer * cell]
}

fn outside(radius: f64, holes: &[Hole], line: f64) -> f64 {
    let cell = 1.0 / (3.0 * COARSE as f64);
    let mut worst: f64 = 0.0;
    for a in 0..COARSE {
        for b in 0..COARSE {
            for c in 0..COARSE {
                let along = (a as f64 + 0.5) * cell;
                let (y, z) = ((b as f64 + 0.5) * cell, (c as f64 + 0.5) * cell);
                let (u, v) = (y.min(1.0 / 3.0 - y), z.min(1.0 / 3.0 - z));
                if u > line || v > line {
                    continue;
                }
                let inside = holes.iter().any(|h| {
                    let low = box_low(radius, h, line);
                    u >= low && v >= low && h.rows.iter().any(|q| along > *q && along < q + h.side)
                });
                if !inside {
                    worst = worst.max(distance([along, 1.0 / 3.0 + y, 1.0 / 3.0 + z]));
                }
            }
        }
    }
    worst
}

fn centre_raster(radius: f64) -> [f64; 3] {
    let lap = (radius * radius - EDGE * EDGE).max(0.0).sqrt();
    let step = (EDGE - lap) / CENTRE as f64;
    let half = step * 3f64.sqrt() / 2.0;
    let (mut mid, mut inner, mut outer) = (0.0, 0.0, 0.0);
    for a in 0..CENTRE {
        for b in 0..CENTRE {
            for c in 0..CENTRE {
                let at = [a, b, c].map(|k| 1.0 / 3.0 + lap + (k as f64 + 0.5) * step);
                let d = distance(at);
                mid += f64::from(u8::from(d > radius));
                inner += f64::from(u8::from(d > radius + half));
                outer += f64::from(u8::from(d > radius - half));
            }
        }
    }
    let cell = 8.0 * step.powi(3);
    [mid * cell, inner * cell, outer * cell]
}

// PERIOD

fn periodic(radius: f64) -> f64 {
    let (mut total, mut weight, mut reach) = (20.0 / 27.0, 1.0, radius);
    for _ in 0..=DEPTH {
        total += weight * exact(reach);
        weight *= 27.0 / 20.0;
        reach /= 3.0;
    }
    radius.powf(dimension() - 3.0) * total
}

fn golden(f: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64, sign: f64) -> f64 {
    let ratio = (5f64.sqrt() - 1.0) / 2.0;
    for _ in 0..200 {
        let (x, y) = (b - ratio * (b - a), a + ratio * (b - a));
        if sign * f(x) > sign * f(y) {
            b = y;
        } else {
            a = x;
        }
    }
    (a + b) / 2.0
}

fn main() {
    let clock = Instant::now();
    println!("THRESHOLDS");
    for m in 1..=5 {
        let half = 3f64.powi(-m - 1) / 2.0;
        println!(
            "delta_{m} = sqrt(1/36 + s_{m}^2/4) = {:.12}",
            (1.0 / 36.0 + half * half).sqrt()
        );
    }
    println!(
        "arms swallowed at sqrt(10)/18 = {:.12}",
        10f64.sqrt() / 18.0
    );
    println!("plus swallowed at sqrt(2)/6 = {:.12}", COVER);
    println!();
    println!(
        "LEVELS: closed deep, lab quadrature by level, row raster by level (midpoint, bracket)"
    );
    let radii = [
        0.22,
        0.2,
        0.18,
        0.174,
        0.17,
        0.168,
        0.1675,
        0.167,
        0.16675,
        0.1667,
        0.16667,
        EDGE,
        0.148,
        1.0 / 8.0,
        1.0 / 12.0,
    ];
    let mut worst_ratio: f64 = 0.0;
    let mut worst_centre: f64 = 0.0;
    let mut worst_gap: f64 = 0.0;
    for &radius in &radii {
        let line = radius.min(EDGE);
        let holes = crossing(line, 12);
        let parts: Vec<f64> = holes
            .iter()
            .map(|h| h.rows.len() as f64 * quadrature(radius, h, line))
            .collect();
        let sum: f64 = parts.iter().sum();
        worst_gap = worst_gap.max((deep(radius) - sum).abs());
        println!(
            "radius {radius:.6}  T = {:.12}  deep = {:.6e}  lab sum = {:.6e}  gap = {:.1e}",
            exact(radius),
            deep(radius),
            sum,
            (deep(radius) - sum).abs()
        );
        for (h, part) in holes.iter().zip(&parts).filter(|(h, _)| h.level <= TOP) {
            if *part == 0.0 {
                println!(
                    "  level {} x {}: closed 0, raster skipped",
                    h.level,
                    h.rows.len()
                );
                continue;
            }
            let [mid, inner, outer] = rows_raster(radius, h, line);
            if radius >= EDGE {
                worst_ratio = worst_ratio.max((mid / part - 1.0).abs());
            }
            assert!(
                inner <= *part && *part <= outer,
                "bracket misses level {} at {radius}",
                h.level
            );
            println!(
                "  level {} x {}: lab {:.6e}  raster {:.6e}  in [{:.6e}, {:.6e}]",
                h.level,
                h.rows.len(),
                part,
                mid,
                inner,
                outer
            );
        }
        if radius >= EDGE {
            let far = outside(radius, &holes, line);
            assert!(
                far <= radius,
                "an uncovered point outside the boxes at {radius}"
            );
            let [mid, inner, outer] = centre_raster(radius);
            let closed = PLUS - exact(radius) - 24.0 * deep(radius);
            assert!(
                inner <= closed && closed <= outer,
                "centre bracket misses at {radius}"
            );
            worst_centre = worst_centre.max((mid / closed - 1.0).abs());
            println!("  arm cells outside the boxes: largest distance {far:.9} <= radius");
            println!("  centre deficit: closed {closed:.9e}  raster {mid:.9e}  in [{inner:.6e}, {outer:.6e}]");
        }
    }
    println!("largest relative gap on [1/6, sqrt(2)/6], row raster midpoint against the level integral: {worst_ratio:.2e}");
    println!(
        "largest gap at all {} radii, crate deep against the lab quadrature: {worst_gap:.1e}",
        radii.len()
    );
    println!("largest relative gap on [1/6, sqrt(2)/6], centre raster midpoint against the closed deficit: {worst_centre:.2e}");
    println!();
    println!("PERIOD: p from exact T on u = ln(1/eps) over one period from ln(1/(sqrt(2)/6))");
    let start = (1.0 / COVER).ln();
    let span = 3f64.ln();
    let at = |u: f64| periodic((-u).exp());
    let grid = |n: usize| -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| start + span * (i as f64 + 0.5) / n as f64)
            .map(|u| (u, at(u)))
            .collect()
    };
    let walk = grid(STEPS);
    let step = span / STEPS as f64;
    let (top_u, _) =
        walk.iter().copied().fold(
            (0.0, f64::NEG_INFINITY),
            |b, x| if x.1 > b.1 { x } else { b },
        );
    let (bottom_u, _) =
        walk.iter()
            .copied()
            .fold((0.0, f64::INFINITY), |b, x| if x.1 < b.1 { x } else { b });
    let crest = golden(&at, top_u - step, top_u + step, 1.0);
    let trough = golden(&at, bottom_u - step, bottom_u + step, -1.0);
    let (high, low) = (at(crest), at(trough));
    println!("maximum p = {high:.9} at eps = {:.9}", (-crest).exp());
    println!("minimum p = {low:.9} at eps = {:.9}", (-trough).exp());
    println!(
        "swing (max - min)/min = {:.6} %",
        100.0 * (high - low) / low
    );
    let mean = |n: usize| grid(n).iter().map(|x| x.1).sum::<f64>() / n as f64;
    println!(
        "logarithmic mean of p over the period: {:.9} ({} steps), {:.9} ({} steps)",
        mean(STEPS),
        STEPS,
        mean(2 * STEPS),
        2 * STEPS
    );
    let window: Vec<f64> = walk
        .iter()
        .filter(|x| (-x.0).exp() > EDGE)
        .map(|x| x.1)
        .collect();
    let window_high = window.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let window_low = window.iter().copied().fold(f64::INFINITY, f64::min);
    println!("on the phases (1/6, sqrt(2)/6]: p runs over [{window_low:.9}, {window_high:.9}], {} of {STEPS} steps", window.len());
    let falling = window.windows(2).all(|pair| pair[1] > pair[0]);
    println!("p falls as eps grows, at every step across (1/6, sqrt(2)/6]: {falling}");
    for eps in [1.0 / 12.0, 1.0 / 8.0, EDGE, 0.2, COVER] {
        let edge = profile(eps).map_or("none".to_string(), |v| format!("{v:.9}"));
        println!(
            "p({eps:.9}) = {:.9}   upper edge (Deep = 0) {edge}",
            periodic(eps)
        );
    }
    println!();
    println!("runtime {:.1} s", clock.elapsed().as_secs_f64());
}
