use census::save;
use mrly_core::error::Result;
use mrly_core::{json, Json};
use mrly_math::design::elements;
use mrly_num::zeta::{Complex, Line};

const NAME: &str = "research-zeta";
const PEEL: u32 = 8;
const DEPTH: usize = 10;
const CUT: f64 = 16.0;
const LMAX: usize = 16;
const NODES: usize = 34;
const HEIGHT: f64 = 40.0;
const WIDE: f64 = 1.20;
const REACH: f64 = 0.70;
const DISC: f64 = 0.45;
const RUN: f64 = 0.03;
const RISE: f64 = 0.1;
const ONE: Complex = Complex::new(1.0, 0.0);

struct Ladder {
    q: f64,
    lq: f64,
    k: f64,
    alpha: f64,
    gamma: Vec<f64>,
    qpow: Vec<f64>,
    head: Vec<f64>,
    blogs: Vec<f64>,
    drops: Vec<Vec<f64>>,
    logs: Vec<f64>,
    fall: Vec<f64>,
}

fn power(log: f64, w: Complex) -> Complex {
    Complex::exp(Complex::new(-w.re * log, -w.im * log))
}

impl Ladder {
    fn new(q: u64, digits: &[u64]) -> Ladder {
        let all = elements(q, digits, DEPTH);
        let split = q.pow(PEEL - 1);
        let stop = q.pow(PEEL);
        let head: Vec<f64> = all
            .iter()
            .filter(|&&n| n < split)
            .map(|&n| n as f64)
            .collect();
        let blogs: Vec<f64> = all
            .iter()
            .filter(|&&n| n >= split && n < stop)
            .map(|&n| (n as f64).ln())
            .collect();
        let logs: Vec<f64> = all
            .iter()
            .filter(|&&n| n >= split)
            .map(|&n| (n as f64).ln())
            .collect();
        let fall: Vec<f64> = all
            .iter()
            .filter(|&&n| n >= split)
            .map(|&n| 1.0 / n as f64)
            .collect();
        let drops = (0..NODES)
            .map(|node| {
                blogs
                    .iter()
                    .map(|&log| (-(node as f64) * log).exp())
                    .collect()
            })
            .collect();
        let k = digits.len() as f64;
        let mut gamma = vec![k];
        for l in 1..=LMAX {
            gamma.push(digits.iter().map(|&a| (a as f64).powi(l as i32)).sum());
        }
        let qpow = (0..=LMAX).map(|l| (q as f64).powi(-(l as i32))).collect();
        Ladder {
            q: q as f64,
            lq: (q as f64).ln(),
            k,
            alpha: k.ln() / (q as f64).ln(),
            gamma,
            qpow,
            head,
            blogs,
            drops,
            logs,
            fall,
        }
    }

    fn cofactor(&self, s: Complex) -> Complex {
        let top = ((CUT - s.re).ceil().max(1.0) as usize).min(NODES);
        let mut g = vec![Complex::default(); NODES + 1];
        let mut pw: Vec<Complex> = self
            .logs
            .iter()
            .map(|&l| power(l, s + top as f64))
            .collect();
        let mut acc = Complex::default();
        for v in &pw {
            acc = acc + *v;
        }
        g[top] = acc;
        for slot in g.iter_mut().take(NODES + 1).skip(top + 1) {
            let mut acc = Complex::default();
            for (p, &drop) in pw.iter_mut().zip(&self.fall) {
                *p = *p * drop;
                acc = acc + *p;
            }
            *slot = acc;
        }
        let phase: Vec<Complex> = self.blogs.iter().map(|&log| power(log, s)).collect();
        let mut numerator = Complex::default();
        for node in (0..top).rev() {
            let w = s + node as f64;
            let qw = power(self.lq, w);
            let mut num = Complex::default();
            for (p, &drop) in phase.iter().zip(&self.drops[node]) {
                num = num + *p * drop;
            }
            let mut binom = ONE;
            for l in 1..=LMAX {
                binom = binom * (-(w + (l as f64 - 1.0))) * (1.0 / l as f64);
                if node + l > NODES {
                    break;
                }
                num = num + binom * qw * (self.qpow[l] * self.gamma[l]) * g[node + l];
            }
            g[node] = num / (ONE - qw * self.k);
            if node == 0 {
                numerator = num;
            }
        }
        let mut d = Complex::default();
        for &n in &self.head {
            d = d + power(n.ln(), s);
        }
        (ONE - power(self.lq, s) * self.k) * d + numerator
    }

    fn centres(&self) -> Vec<(Complex, bool)> {
        let mut out = Vec::new();
        for t in ordinates(self.lq, HEIGHT + 1.0) {
            let zero = Complex::new(self.alpha, t);
            let r0 = self.cofactor(zero) * (1.0 / self.lq);
            out.push((zero, Complex::abs(r0) > 1e-12));
            let one = Complex::new(self.alpha - 1.0, t);
            let r1 = one * (self.gamma[1] / (self.k * (self.q - 1.0))) * r0;
            out.push((one, Complex::abs(r1) > 1e-12));
        }
        out
    }
}

fn ordinates(lq: f64, reach: f64) -> Vec<f64> {
    let step = std::f64::consts::TAU / lq;
    (0..)
        .map(|j| j as f64 * step)
        .take_while(|t| *t < reach)
        .collect()
}

fn refine(l: &Ladder, seed: Complex) -> Option<Complex> {
    let h = 1e-6;
    let mut s = seed;
    for _ in 0..40 {
        let v = l.cofactor(s);
        if Complex::abs(v) < 1e-13 {
            return Some(s);
        }
        let a = l.cofactor(s + Complex::new(h, 0.0));
        let b = l.cofactor(s - Complex::new(h, 0.0));
        let d = (a - b) * (1.0 / (2.0 * h));
        if Complex::abs(d) < 1e-30 {
            return None;
        }
        let step = v / d;
        if Complex::abs(step) > 1.0 {
            return None;
        }
        s = s - step;
    }
    if Complex::abs(l.cofactor(s)) < 1e-10 {
        Some(s)
    } else {
        None
    }
}

fn hunt(l: &Ladder, lo: f64, hi: f64) -> Vec<Complex> {
    let cols = ((hi - lo) / RUN).round() as usize;
    let rows = (HEIGHT / RISE).round() as usize;
    let at = |i: usize, j: usize| {
        Complex::new(
            lo + (hi - lo) * i as f64 / cols as f64,
            0.02 + (HEIGHT - 0.04) * j as f64 / rows as f64,
        )
    };
    let lanes = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let grid = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..lanes)
            .map(|lane| {
                let at = &at;
                scope.spawn(move || {
                    (lane..=cols)
                        .step_by(lanes)
                        .map(|i| {
                            let column: Vec<f64> = (0..=rows)
                                .map(|j| Complex::abs(l.cofactor(at(i, j))))
                                .collect();
                            (i, column)
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let mut grid = vec![0.0f64; (cols + 1) * (rows + 1)];
        for handle in handles {
            for (i, column) in handle.join().expect("a lane of the hunt") {
                grid[i * (rows + 1)..(i + 1) * (rows + 1)].copy_from_slice(&column);
            }
        }
        grid
    });
    let mut found: Vec<Complex> = Vec::new();
    for i in 1..cols {
        for j in 1..rows {
            let here = grid[i * (rows + 1) + j];
            let mut least = true;
            for di in 0..3 {
                for dj in 0..3 {
                    if (di, dj) != (1, 1) && grid[(i + di - 1) * (rows + 1) + j + dj - 1] <= here {
                        least = false;
                    }
                }
            }
            if !least {
                continue;
            }
            if let Some(root) = refine(l, at(i, j)) {
                if root.re < lo || root.re > hi || root.im < 0.02 || root.im > HEIGHT {
                    continue;
                }
                if found.iter().all(|z| Complex::abs(*z - root) > 1e-4) {
                    found.push(root);
                }
            }
        }
    }
    found.sort_by(|a, b| a.im.partial_cmp(&b.im).unwrap());
    found
}

fn split(l: &Ladder, zeros: &[Complex]) -> Panel {
    let centres = l.centres();
    let mut teeth = Vec::new();
    let mut hollow = Vec::new();
    let mut family = Vec::new();
    for z in zeros {
        let mut best = f64::INFINITY;
        let mut live = false;
        for (c, alive) in &centres {
            let d = Complex::abs(*z - *c);
            if d < best {
                best = d;
                live = *alive;
            }
        }
        if best >= DISC {
            family.push(*z);
        } else if live {
            teeth.push(*z);
        } else {
            hollow.push(*z);
        }
    }
    Panel {
        alpha: l.alpha,
        lq: l.lq,
        teeth,
        hollow,
        family,
    }
}

struct Panel {
    alpha: f64,
    lq: f64,
    teeth: Vec<Complex>,
    hollow: Vec<Complex>,
    family: Vec<Complex>,
}

fn points(zeros: &[Complex]) -> Json {
    json!(zeros.iter().map(|z| [z.re, z.im]).collect::<Vec<_>>())
}

fn record(p: &Panel) -> Json {
    json!({
        "columns": [p.alpha, p.alpha - 1.0],
        "ordinates": ordinates(p.lq, HEIGHT),
        "family": points(&p.family),
        "teeth": points(&p.teeth),
        "hollow": points(&p.hollow),
    })
}

fn panels() -> (Panel, Panel) {
    let design = Ladder::new(3, &[0, 1]);
    let full = Ladder::new(2, &[0, 1]);
    let known = Line::new().zeros(6);
    assert!(Complex::abs(full.cofactor(Complex::new(0.5, known[0]))) < 1e-8);

    let dz = hunt(&design, design.alpha - 0.92, design.alpha + 3.02);
    let fz = hunt(&full, full.alpha - 0.92, full.alpha + 3.02);
    let f = split(&full, &fz);
    let d = split(&design, &dz);

    assert_eq!(fz.len(), 10);
    assert_eq!(f.teeth.len(), 0);
    assert_eq!(f.hollow.len(), 4);
    assert_eq!(f.family.len(), 6);
    for z in &f.hollow {
        assert!((z.re - full.alpha).abs() < 1e-6);
    }
    for (z, g) in f.family.iter().zip(&known) {
        assert!((z.re - 0.5).abs() < 1e-6 && (z.im - g).abs() < 1e-6);
    }

    assert_eq!(dz.len(), 13);
    assert_eq!(d.hollow.len(), 0);
    assert_eq!(d.teeth.len() + d.family.len(), 13);
    assert!(d.teeth.iter().all(|z| (z.re - design.alpha).abs() > 1e-3));
    assert!(d.family.iter().any(|z| (z.re - 0.391038600).abs() < 1e-6));
    assert!(d
        .family
        .iter()
        .any(|z| (z.re + 0.273079611).abs() < 1e-6 && (z.im - 39.262315320).abs() < 1e-6));
    let least = d.family.iter().map(|z| z.re).fold(f64::INFINITY, f64::min);
    let most = d
        .family
        .iter()
        .map(|z| z.re)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((least + 0.273079611).abs() < 1e-6 && (most - 0.391038600).abs() < 1e-6);
    assert!(dz
        .iter()
        .all(|z| z.re > design.alpha - WIDE && z.re < design.alpha + REACH));
    assert!(fz
        .iter()
        .all(|z| z.re > full.alpha - WIDE && z.re < full.alpha + REACH));
    (d, f)
}

fn main() -> Result<()> {
    let (d, f) = panels();
    save(NAME, &json!({"design": record(&d), "full": record(&f)}))
}
