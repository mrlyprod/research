use crate::lattice::{cell, column, row, Rule};
use crate::sums::{mean, odds};

#[derive(Clone, Copy)]
pub enum Line {
    X,
    Z,
    W,
    D,
    K,
}

impl Line {
    pub fn name(self) -> &'static str {
        match self {
            Line::X => "X",
            Line::Z => "Z",
            Line::W => "X+Z",
            Line::D => "X-Z",
            Line::K => "2X+Z",
        }
    }

    fn at(self, x: f64, z: f64) -> f64 {
        match self {
            Line::X => x,
            Line::Z => z,
            Line::W => x + z,
            Line::D => x - z,
            Line::K => 2.0 * x + z,
        }
    }
}

pub struct Frame {
    pub side: usize,
    pub layers: usize,
    pub values: Vec<f64>,
    pub whole: Vec<bool>,
    pub points: Vec<f64>,
}

pub fn points(side: usize) -> Vec<f64> {
    (0..side)
        .map(|step| (step as f64 + 0.5) / side as f64)
        .collect()
}

pub fn stack(limit: usize, rule: &Rule, side: usize) -> Frame {
    let points = points(side);
    let mut values = vec![0.0f64; side * side];
    let mut whole = vec![true; side * side];
    let mut layers = 0;
    let mut line = vec![0.0f64; side];
    let mut reach = vec![true; side];
    for number in odds(limit) {
        let n = number as i64;
        let columns: Vec<i64> = points.iter().map(|p| column(*p, n)).collect();
        let rows: Vec<i64> = points.iter().map(|p| row(*p, n)).collect();
        let mut start = 0;
        while start < side {
            let x = columns[start];
            let mut end = start;
            while end < side && columns[end] == x {
                end += 1;
            }
            for (col, z) in rows.iter().enumerate() {
                let found = cell(rule, n, x, *z);
                line[col] = f64::from(found.unwrap_or(false));
                reach[col] = found.is_some();
            }
            let clipped = reach.iter().any(|ok| !ok);
            for slot in start..end {
                let base = slot * side;
                for (value, add) in values[base..base + side].iter_mut().zip(&line) {
                    *value += add;
                }
                if clipped {
                    for (keep, ok) in whole[base..base + side].iter_mut().zip(&reach) {
                        *keep &= ok;
                    }
                }
            }
            start = end;
        }
        layers += 1;
    }
    for value in values.iter_mut() {
        *value /= layers as f64;
    }
    Frame {
        side,
        layers,
        values,
        whole,
        points,
    }
}

impl Frame {
    pub fn window(&self, line: Line, keep: impl Fn(f64) -> bool) -> f64 {
        let mut kept = Vec::new();
        for row in 0..self.side {
            for col in 0..self.side {
                let at = row * self.side + col;
                if self.whole[at] && keep(line.at(self.points[row], self.points[col])) {
                    kept.push(self.values[at]);
                }
            }
        }
        mean(&kept)
    }

    pub fn hexmean(&self) -> f64 {
        let kept: Vec<f64> = self
            .values
            .iter()
            .zip(&self.whole)
            .filter(|(_, keep)| **keep)
            .map(|(value, _)| *value)
            .collect();
        mean(&kept)
    }

    pub fn band(&self, line: Line, position: f64, width: f64) -> [f64; 3] {
        let left = self.window(line, |at| at >= position - width && at < position);
        let right = self.window(line, |at| at > position && at <= position + width);
        let centre = self.window(line, |at| (at - position).abs() <= width);
        [left, right, centre]
    }
}

pub fn plateau(limit: usize, rule: &Rule, side: usize) -> (f64, f64) {
    let frame = stack(limit, rule, side);
    let background = frame.hexmean();
    let mut sums = [vec![0.0f64; side], vec![0.0f64; side]];
    let mut whole = [vec![true; side], vec![true; side]];
    for number in odds(limit) {
        let n = number as i64;
        for (slot, x) in [n - 1, n].into_iter().enumerate() {
            for (at, point) in frame.points.iter().enumerate() {
                match cell(rule, n, x, row(*point, n)) {
                    Some(fill) => sums[slot][at] += f64::from(fill),
                    None => whole[slot][at] = false,
                }
            }
        }
    }
    let strip = |slot: usize| {
        let kept: Vec<f64> = sums[slot]
            .iter()
            .zip(&whole[slot])
            .filter(|(_, keep)| **keep)
            .map(|(value, _)| value / frame.layers as f64)
            .collect();
        mean(&kept) - background
    };
    (strip(0), strip(1))
}

pub fn quarter(rule: &Rule) {
    println!("carpet ideal frame at N = 55, 3601 samples per axis, one-sided bands of width 0.018");
    let frame = stack(55, rule, 3601);
    let background = frame.hexmean();
    println!("  hexagon mean {background:.6}");
    let spec = [
        (Line::X, 0.25),
        (Line::X, 0.75),
        (Line::X, 1.0 / 3.0),
        (Line::X, 0.2),
        (Line::X, 0.5),
        (Line::Z, 0.25),
        (Line::D, 0.25),
        (Line::W, 1.25),
    ];
    for (line, position) in spec {
        let [left, right, _] = frame.band(line, position, 0.018);
        println!(
            "  {} = {position:.4}: left {:+.6}  right {:+.6}",
            line.name(),
            left - background,
            right - background
        );
    }
    println!("quarter line X = 1/4, one-sided strips straight at the line, 3601 samples");
    for limit in [151usize, 301, 601, 1201] {
        let (left, right) = plateau(limit, rule, 3601);
        println!(
            "  N = {limit:4}: left {left:+.5}  right {right:+.5}  mean magnitude {:.5}  (limit 1/8)",
            (left.abs() + right.abs()) / 2.0
        );
    }
}

pub fn void_star(rule: &Rule) {
    println!(
        "void ideal frame at N = 55, 3601 samples per axis, centred bands of half-width 0.004"
    );
    let frame = stack(55, rule, 3601);
    let background = frame.hexmean();
    println!("  hexagon mean {background:.6}  (limit 1/4)");
    for (line, position) in [
        (Line::X, 0.5),
        (Line::Z, 0.5),
        (Line::W, 1.0),
        (Line::D, 0.0),
        (Line::K, 1.5),
    ] {
        let [_, _, centre] = frame.band(line, position, 0.004);
        println!(
            "  {} = {position:.2}: ink {centre:.6}  ratio to mean {:.3}",
            line.name(),
            centre / background
        );
    }
}

pub fn fading(rule: &Rule) {
    println!("carpet main diagonal X = Z, band half-width 0.01, 2801 samples per axis");
    let mut scaled = Vec::new();
    for count in [28usize, 100, 200, 400] {
        let frame = stack(2 * count - 1, rule, 2801);
        let excess = frame.window(Line::D, |at| at.abs() <= 0.01) - frame.hexmean();
        println!(
            "  L = {count:3}: excess {excess:+.6}  excess * L {:+.4}",
            excess * count as f64
        );
        scaled.push(excess * count as f64);
    }
    println!(
        "  slope of excess * L against ln L: {:+.4} (100 to 200)  {:+.4} (200 to 400)",
        (scaled[2] - scaled[1]) / 2f64.ln(),
        (scaled[3] - scaled[2]) / 2f64.ln()
    );
}
