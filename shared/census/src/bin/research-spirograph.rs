use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::bang::Code;
use mrly_math::roulette;
use mrly_math::spirograph::{distinct, nodes as law, pencils, trace, track, Pencil, Track};
use mrly_math::two;
use std::collections::HashSet;

const NAME: &str = "research-spirograph";
const CODE: u128 = 495;
const SIDE: usize = 3;
const LEVEL: usize = 1;
const BASE: usize = 3;
const RING: usize = 7;
const WHEEL: usize = 3;
const REACH: f64 = 0.9;
const SAMPLES: usize = 8000;
const CURVES: usize = 8;
const CROSSINGS: usize = 1288;
const TOL: f64 = 4e-4;
const GRID: usize = 256;

fn main() -> Result<()> {
    let tile = two::create(Code::from(CODE), SIDE, LEVEL, 0, BASE)?;
    let types = tile.types().bytes()?.to_vec();
    let path = track("in", RING, WHEEL, 4, 1)?;
    let pens = pencils(&types, tile.width(), tile.height(), "fill", REACH, 0.0, 1)?;
    assert_eq!(distinct(&path, &pens, true), CURVES);
    assert_eq!(law(&path, &pens, true), Some(CROSSINGS as u64));
    assert_eq!(
        roulette::nodes(&path, &pens, SAMPLES, TOL)?.total(),
        CROSSINGS
    );
    let marks = crossings(&path, &pens, SAMPLES)?;
    assert_eq!(marks.len(), CROSSINGS);
    let nodes: Vec<_> = marks.iter().map(|&(x, y)| json!([x, y])).collect();
    save(NAME, &json!({"nodes": nodes}))
}

fn crossings(path: &Track, pens: &[Pencil], samples: usize) -> Result<Vec<(f64, f64)>> {
    let points = trace(path, pens, samples)?;
    let steps = samples - 1;
    let at = |k: usize, i: usize| {
        let base = 2 * (k * samples + i);
        (f64::from(points[base]), f64::from(points[base + 1]))
    };
    let (mut low, mut high) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for pair in points.chunks_exact(2) {
        let (x, y) = (f64::from(pair[0]), f64::from(pair[1]));
        low = (low.0.min(x), low.1.min(y));
        high = (high.0.max(x), high.1.max(y));
    }
    let size = (high.0 - low.0, high.1 - low.1);
    let cell = |p: (f64, f64)| {
        let column = ((p.0 - low.0) / size.0 * GRID as f64) as usize;
        let row = ((p.1 - low.1) / size.1 * GRID as f64) as usize;
        (column.min(GRID - 1), row.min(GRID - 1))
    };
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); GRID * GRID];
    for k in 0..pens.len() {
        for i in 0..steps {
            let (a, b) = (at(k, i), at(k, i + 1));
            let (x0, y0) = cell((a.0.min(b.0), a.1.min(b.1)));
            let (x1, y1) = cell((a.0.max(b.0), a.1.max(b.1)));
            for x in x0..=x1 {
                for y in y0..=y1 {
                    buckets[x * GRID + y].push((k * steps + i) as u32);
                }
            }
        }
    }
    let mut seen: HashSet<(u32, u32)> = HashSet::new();
    let mut out = Vec::new();
    for bucket in &buckets {
        for (u, &left) in bucket.iter().enumerate() {
            for &right in &bucket[u + 1..] {
                let (p, q) = (left as usize, right as usize);
                let ((kp, ip), (kq, iq)) = ((p / steps, p % steps), (q / steps, q % steps));
                let gap = ip.abs_diff(iq);
                if kp == kq && (gap == 1 || gap == steps - 1) {
                    continue;
                }
                let (a, b) = (at(kp, ip), at(kp, ip + 1));
                let (c, d) = (at(kq, iq), at(kq, iq + 1));
                if side(a, b, c) * side(a, b, d) < 0
                    && side(c, d, a) * side(c, d, b) < 0
                    && seen.insert((left.min(right), left.max(right)))
                {
                    out.push(meet(a, b, c, d));
                }
            }
        }
    }
    Ok(out)
}

fn side(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> i32 {
    let turn = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    match turn.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => 1,
        Some(std::cmp::Ordering::Less) => -1,
        _ => 0,
    }
}

fn meet(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> (f64, f64) {
    let run = (b.0 - a.0, b.1 - a.1);
    let other = (d.0 - c.0, d.1 - c.1);
    let denominator = run.0 * other.1 - run.1 * other.0;
    let step = ((c.0 - a.0) * other.1 - (c.1 - a.1) * other.0) / denominator;
    (a.0 + step * run.0, a.1 + step * run.1)
}
