use std::collections::HashSet;
use std::env;

// THE OBJECTS

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct Circle {
    k: i64,
    x: i64,
    y: i64,
}

type Quad = [Circle; 4];

fn form(u: [i64; 4], v: [i64; 4]) -> i128 {
    let su: i128 = u.iter().map(|&a| a as i128).sum();
    let sv: i128 = v.iter().map(|&a| a as i128).sum();
    let dot: i128 = (0..4).map(|i| u[i] as i128 * v[i] as i128).sum();
    su * sv - 2 * dot
}

fn curvatures(q: &Quad) -> [i64; 4] {
    [q[0].k, q[1].k, q[2].k, q[3].k]
}

fn abscissae(q: &Quad) -> [i64; 4] {
    [q[0].x, q[1].x, q[2].x, q[3].x]
}

fn ordinates(q: &Quad) -> [i64; 4] {
    [q[0].y, q[1].y, q[2].y, q[3].y]
}

fn descartes(q: &Quad) -> bool {
    let (k, x, y) = (curvatures(q), abscissae(q), ordinates(q));
    form(k, k) == 0 && form(k, x) == 0 && form(k, y) == 0
}

fn frame(q: &Quad) -> bool {
    let (x, y) = (abscissae(q), ordinates(q));
    form(x, x) == -4 && form(y, y) == -4 && form(x, y) == 0
}

fn reflect(q: &Quad, i: usize) -> Circle {
    let (mut k, mut x, mut y) = (0i64, 0i64, 0i64);
    for (j, c) in q.iter().enumerate() {
        if j != i {
            k += c.k;
            x += c.x;
            y += c.y;
        }
    }
    Circle {
        k: 2 * k - q[i].k,
        x: 2 * x - q[i].x,
        y: 2 * y - q[i].y,
    }
}

fn swap(q: &Quad, i: usize) -> Quad {
    let mut out = *q;
    out[i] = reflect(q, i);
    out
}

fn strip_root() -> Quad {
    [
        Circle { k: 0, x: 0, y: -1 },
        Circle { k: 0, x: 0, y: 1 },
        Circle { k: 2, x: 0, y: 1 },
        Circle { k: 2, x: 2, y: 1 },
    ]
}

fn bounded_root() -> Quad {
    [
        Circle { k: -1, x: 0, y: 0 },
        Circle { k: 2, x: 1, y: 0 },
        Circle { k: 2, x: -1, y: 0 },
        Circle { k: 3, x: 0, y: 2 },
    ]
}

// THE CENSUS

struct Run {
    circles: u64,
    quads: u64,
    broken: u64,
    backward: u64,
    strayed: u64,
    bins: Vec<u64>,
    res: [u64; 24],
    depth: usize,
    bottom: u64,
    top: u64,
    offford: u64,
}

fn grow(
    root: Quad,
    first: &[usize],
    steps: &[i64],
    strip: bool,
    watch: Option<&mut HashSet<Circle>>,
) -> Run {
    let cap = *steps.last().unwrap();
    let mut run = Run {
        circles: 0,
        quads: 1,
        broken: 0,
        backward: 0,
        strayed: 0,
        bins: vec![0; steps.len()],
        res: [0; 24],
        depth: 0,
        bottom: 0,
        top: 0,
        offford: 0,
    };
    let mut seen = watch;
    if !descartes(&root) || !frame(&root) {
        run.broken += 1;
    }
    let mut stack: Vec<(Quad, usize, usize)> = Vec::new();
    for &i in first {
        let c = reflect(&root, i);
        if c.k <= cap {
            stack.push((swap(&root, i), i, 1));
        }
    }
    while let Some((q, last, depth)) = stack.pop() {
        run.quads += 1;
        run.depth = run.depth.max(depth);
        if !descartes(&q) || !frame(&q) {
            run.broken += 1;
        }
        let c = q[last];
        run.circles += 1;
        run.res[c.k.rem_euclid(24) as usize] += 1;
        let at = steps.iter().position(|&s| c.k <= s).unwrap();
        run.bins[at] += 1;
        if strip {
            if c.y == 1 {
                run.bottom += 1;
                if !is_ford(c.k, c.x) {
                    run.offford += 1;
                }
            }
            if c.y == c.k - 1 {
                run.top += 1;
            }
            if c.x <= 0 || c.x >= c.k {
                run.strayed += 1;
            }
        }
        if let Some(set) = seen.as_mut() {
            if !set.insert(c) {
                run.backward += 1;
            }
        }
        for j in 0..4 {
            if j == last {
                continue;
            }
            let n = reflect(&q, j);
            if n.k > cap {
                continue;
            }
            if n.k <= q[j].k {
                run.backward += 1;
                continue;
            }
            stack.push((swap(&q, j), j, depth + 1));
        }
    }
    run
}

fn totals(bins: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(bins.len());
    let mut acc = 0u64;
    for &b in bins {
        acc += b;
        out.push(acc);
    }
    out
}

fn is_ford(k: i64, x: i64) -> bool {
    if k % 2 != 0 {
        return false;
    }
    let s = k / 2;
    let b = (s as f64).sqrt().round() as i64;
    if b * b != s || b <= 0 {
        return false;
    }
    if x % (2 * b) != 0 {
        return false;
    }
    let a = x / (2 * b);
    a >= 0 && a <= b && gcd(a, b) == 1
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

// THE FORD IDENTIFICATION

struct Ford {
    nodes: u64,
    broken: u64,
    missed: u64,
    untangent: u64,
    bright: u128,
    deepest: i64,
}

fn ford(top: i64) -> Ford {
    let mut f = Ford {
        nodes: 0,
        broken: 0,
        missed: 0,
        untangent: 0,
        bright: 0,
        deepest: 0,
    };
    let root = strip_root();
    let base = [root[0], root[2], root[3], root[1]];
    let mut stack = vec![(base, 0i64, 1i64, 1i64, 1i64)];
    while let Some((q, a, b, c, d)) = stack.pop() {
        let (p, r) = (a + c, b + d);
        if r > top {
            continue;
        }
        if (a * d - b * c) * (a * d - b * c) != 1 {
            f.untangent += 1;
        }
        let m = reflect(&q, 3);
        let child = [q[0], q[1], m, q[2]];
        if !descartes(&child) || !frame(&child) {
            f.broken += 1;
        }
        if m.k != 2 * r * r || m.x != 2 * p * r || m.y != 1 {
            f.missed += 1;
        }
        let lhs = (2 * b * b + 2 * d * d + 2 * r * r) as i128;
        let rhs: i128 = 2
            * ((2 * b * b) as i128 * (2 * b * b) as i128
                + (2 * d * d) as i128 * (2 * d * d) as i128
                + (2 * r * r) as i128 * (2 * r * r) as i128);
        if lhs * lhs != rhs {
            f.broken += 1;
        }
        if 2 * (2 * b * b + 2 * d * d) - 2 * (b - d) * (b - d) != 2 * r * r {
            f.broken += 1;
        }
        f.nodes += 1;
        f.deepest = f.deepest.max(r);
        f.bright += (top / r) as u128;
        stack.push(([q[0], q[1], m, q[2]], a, b, p, r));
        stack.push(([q[0], m, q[2], q[1]], p, r, c, d));
    }
    f.bright += top as u128;
    f
}

fn totients(top: usize) -> Vec<u64> {
    let mut phi: Vec<u64> = (0..=top as u64).collect();
    for i in 2..=top {
        if phi[i] == i as u64 {
            let mut j = i;
            while j <= top {
                phi[j] -= phi[j] / i as u64;
                j += i;
            }
        }
    }
    phi
}

// THE VERBS

const DELTA: f64 = 1.3056867280498771846459862068510;

fn show(name: &str, run: &Run, steps: &[i64]) {
    let tot = totals(&run.bins);
    println!(
        "{name}: {} circles, {} quadruples, depth {}, broken {}, backward {}, strayed {}",
        run.circles, run.quads, run.depth, run.broken, run.backward, run.strayed
    );
    println!("  T  N(T)  log N / log T  ratio");
    let mut last: Option<(i64, u64)> = None;
    for (i, &t) in steps.iter().enumerate() {
        let n = tot[i];
        if n == 0 {
            continue;
        }
        let e = (n as f64).ln() / (t as f64).ln();
        let r = match last {
            Some((pt, pn)) if pn > 0 => {
                format!(
                    "{:.4}",
                    ((n as f64 / pn as f64).ln()) / ((t as f64 / pt as f64).ln())
                )
            }
            _ => "-".to_string(),
        };
        println!("  {t}  {n}  {e:.4}  {r}");
        last = Some((t, n));
    }
    let live: Vec<String> = (0..24)
        .filter(|&r| run.res[r] > 0)
        .map(|r| format!("{r}:{}", run.res[r]))
        .collect();
    println!("  mod 24 {}", live.join(" "));
}

fn main() {
    let verb = env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let all = verb == "all";

    if all || verb == "ford" {
        for top in [50i64, 200, 1000, 4000] {
            let f = ford(top);
            let phi = totients(top as usize);
            let want: u64 = phi[1..=top as usize].iter().sum::<u64>() - 1;
            let bright = (top as u128) * (top as u128 + 1) / 2;
            println!(
                "ford b<={top}: nodes {} want {} broken {} missed {} untangent {} bright {} want {}",
                f.nodes, want, f.broken, f.missed, f.untangent, f.bright, bright
            );
            println!(
                "  curvature at b={top} is {} deepest {}",
                2 * top * top,
                f.deepest
            );
        }
    }

    if all || verb == "strip" {
        let steps = [
            2i64, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152,
        ];
        let mut set = HashSet::new();
        let small = grow(strip_root(), &[0, 1], &steps[..6], true, Some(&mut set));
        println!(
            "strip control T<=2048: {} circles, {} distinct, bottom {} top {}",
            small.circles,
            set.len(),
            small.bottom,
            small.top
        );
        let run = grow(strip_root(), &[0, 1], &steps, true, None);
        show("strip octaves", &run, &steps);
        let decades = [10i64, 100, 1000, 10000, 100000, 1000000];
        let dec = grow(strip_root(), &[0, 1], &decades, true, None);
        show("strip decades", &dec, &decades);
        println!(
            "  bottom-tangent {} top-tangent {} off-Ford {}",
            run.bottom, run.top, run.offford
        );
        for q in [32i64, 181, 1024] {
            let f = grow(strip_root(), &[0, 1], &[2 * q * q], true, None);
            let phi = totients(q as usize);
            let want: u64 = phi[1..=q as usize].iter().sum::<u64>() - 1;
            println!("  line-tangent below 2*{q}^2: {} want {}", f.bottom, want);
        }
    }

    if all || verb == "census" {
        let steps = [10i64, 100, 1000, 10000, 100000, 1000000, 10000000];
        let mut set = HashSet::new();
        let small = grow(
            bounded_root(),
            &[0, 1, 2, 3],
            &steps[..4],
            false,
            Some(&mut set),
        );
        println!(
            "bounded control T<=10000: {} circles, {} distinct",
            small.circles,
            set.len()
        );
        let run = grow(bounded_root(), &[0, 1, 2, 3], &steps, false, None);
        show("bounded (-1,2,2,3)", &run, &steps);
    }

    if all || verb == "design" {
        let mut best: Vec<(f64, i64, f64)> = Vec::new();
        for q in 2i64..=100 {
            let v = (q as f64).powf(DELTA);
            let n = v.round();
            best.push(((v - n).abs(), q, v));
        }
        best.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        println!("design: q^delta against the integers, delta = {DELTA}");
        for &(g, q, v) in best.iter().take(6) {
            let d = v.round().ln() / (q as f64).ln();
            println!(
                "  q {q}  q^delta {v:.9}  gap {g:.9}  log N / log q {d:.9}  off {:.9}",
                (d - DELTA).abs()
            );
        }
        let w = best.last().unwrap();
        println!("  worst gap {:.9} at q {}", w.0, w.1);
        let mut near: Vec<(f64, i64, f64)> = best
            .iter()
            .map(|&(_, q, v)| {
                let d = v.round().ln() / (q as f64).ln();
                ((d - DELTA).abs(), q, d)
            })
            .collect();
        near.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let n = near[0];
        println!(
            "  nearest design dimension {:.9} at q {} off {:.9}",
            n.2, n.1, n.0
        );
    }
}
