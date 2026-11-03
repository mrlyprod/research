use mrly_math::spirograph::{distinct, pencils, point, trace, track, Kind, Pencil};
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy)]
struct Wheel {
    a: usize,
    b: usize,
    side: f64,
    amp: f64,
}

fn wheel(a: usize, b: usize, inside: bool) -> Wheel {
    let side = if inside { -1.0 } else { 1.0 };
    Wheel {
        a,
        b,
        side,
        amp: (a as f64 / b as f64 + side).abs(),
    }
}

fn lean(w: &Wheel, rho: f64, x: f64) -> f64 {
    let k = (rho - w.amp) / (rho + w.amp);
    if x >= PI - 1e-13 {
        return PI / 2.0 + if k > 0.0 { PI / 2.0 } else { -PI / 2.0 };
    }
    x / 2.0 + (k * (x / 2.0).tan()).atan()
}

fn sweep(w: &Wheel, rho: f64, x: f64) -> f64 {
    w.b as f64 * w.side * x / w.a as f64 + lean(w, rho, x)
}

fn ring(w: &Wheel, rho: f64, x: f64) -> f64 {
    (w.amp * w.amp + rho * rho + 2.0 * w.amp * rho * x.cos()).sqrt()
}

fn mark_at(w: &Wheel, rho: f64, r: f64) -> f64 {
    let c = ((r * r - w.amp * w.amp - rho * rho) / (2.0 * w.amp * rho)).clamp(-1.0, 1.0);
    4.0 * w.a as f64 * sweep(w, rho, c.acos()) / PI
}

#[derive(Clone, Copy)]
struct Seat {
    half: (f64, f64),
    corner: bool,
}

fn carpet() -> Vec<Seat> {
    let mut out = Vec::new();
    for i in 0..3i64 {
        for j in 0..3i64 {
            if i == 1 && j == 1 {
                continue;
            }
            let (hx, hy) = ((j - 1) as f64, (1 - i) as f64);
            out.push(Seat {
                half: (hx, hy),
                corner: hx != 0.0 && hy != 0.0,
            });
        }
    }
    out
}

fn reach_unit() -> f64 {
    1.0 / (1.5_f64).hypot(1.5)
}

fn seat_arm(seat: &Seat, t: f64) -> f64 {
    let u = t * reach_unit();
    (seat.half.0 * u).hypot(seat.half.1 * u)
}

fn seat_oct(w: &Wheel, seat: &Seat) -> i64 {
    let dial = match (seat.half.0 as i64, seat.half.1 as i64) {
        (1, 0) => 0,
        (1, 1) => 1,
        (0, 1) => 2,
        (-1, 1) => 3,
        (-1, 0) => 4,
        (-1, -1) => 5,
        (0, -1) => 6,
        _ => 7,
    };
    (-(w.b as i64) * (w.side as i64) * dial).rem_euclid(8)
}

fn band(w: &Wheel, rho: f64) -> (f64, f64) {
    ((w.amp - rho).abs(), w.amp + rho)
}

#[derive(Clone, Copy)]
struct Curve {
    arm: f64,
    oct: i64,
    corner: bool,
}

fn curves(w: &Wheel, seats: &[Seat], t: f64) -> (Vec<Curve>, Vec<usize>) {
    let mut out: Vec<Curve> = Vec::new();
    let mut owner = Vec::new();
    for seat in seats {
        let arm = seat_arm(seat, t);
        let oct = seat_oct(w, seat);
        let at = out
            .iter()
            .position(|c: &Curve| (c.arm - arm).abs() < 1e-12 && c.oct == oct);
        match at {
            Some(k) => owner.push(k),
            None => {
                owner.push(out.len());
                out.push(Curve {
                    arm,
                    oct,
                    corner: seat.corner,
                });
            }
        }
    }
    (out, owner)
}

fn meeting(cs: &[Curve], m: i64, n: i64, flat: bool) -> (usize, usize, i64, Vec<usize>) {
    let mut best = (0usize, 0usize, 0i64, Vec::new());
    for slot in 0..8i64 {
        let mut ends = 0usize;
        let mut seen: Vec<usize> = Vec::new();
        for (k, c) in cs.iter().enumerate() {
            let mark = if c.corner { m } else { n };
            let signs: &[i64] = if flat && !c.corner { &[1] } else { &[1, -1] };
            for &sign in signs {
                if (c.oct + sign * mark).rem_euclid(8) == slot {
                    ends += 1;
                    if !seen.contains(&k) {
                        seen.push(k);
                    }
                }
            }
        }
        if ends > best.0 {
            seen.sort_unstable();
            best = (ends, seen.len(), slot, seen);
        }
    }
    best
}

#[derive(Clone, Copy)]
struct Align {
    reach: f64,
    turn: f64,
    ring: f64,
    mark_c: i64,
    mark_e: i64,
    ends: usize,
    kinds: usize,
    slot: i64,
    flat: bool,
}

fn turn_for(w: &Wheel, arm: f64, mark: f64) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, PI);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        let got = 4.0 * w.a as f64 * sweep(w, arm, mid) / PI;
        if (got - mark) * w.side < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn tame(w: &Wheel) -> f64 {
    1.5 * 1.0_f64.min(w.amp)
}

fn arms(t: f64) -> (f64, f64) {
    (2.0 * t / 3.0, 2.0_f64.sqrt() * t / 3.0)
}

fn edge_mark(w: &Wheel, t: f64, x: f64) -> Option<f64> {
    let (arm_c, arm_e) = arms(t);
    let r = ring(w, arm_c, x);
    let (lo, hi) = band(w, arm_e);
    if r <= lo || r >= hi {
        return None;
    }
    Some(mark_at(w, arm_e, r))
}

fn alignments(w: &Wheel, tlo: f64, thi: f64, grid: usize) -> Vec<Align> {
    let sample = curves(w, &carpet(), 0.5 * (tlo + thi)).0;
    let reach = |i: usize| tlo + (thi - tlo) * i as f64 / grid as f64;
    let full = 4.0 * w.b as f64 * w.side;
    let marks: Vec<i64> = if w.side < 0.0 {
        ((full.ceil() as i64 + 1)..0).collect()
    } else {
        (1..(full.floor() as i64)).collect()
    };
    let hold = |t: f64, m: i64| -> f64 {
        let (arm_c, arm_e) = arms(t);
        mark_at(w, arm_e, ring(w, arm_c, turn_for(w, arm_c, m as f64)))
    };
    let mut out: Vec<Align> = Vec::new();
    for &m in &marks {
        let mut back = hold(reach(0), m);
        for cell in 1..=grid {
            let front = hold(reach(cell), m);
            {
                let (d0, d1) = (back, front);
                let (dlo, dhi) = (d0.min(d1), d0.max(d1));
                for n in (dlo.floor() as i64 + 1)..=(dhi.ceil() as i64) {
                    if (n as f64) <= dlo || (n as f64) >= dhi {
                        continue;
                    }
                    let (ends, kinds, slot, _) = meeting(&sample, m, n, false);
                    if ends < 3 {
                        continue;
                    }
                    let (mut u, mut v) = (reach(cell - 1), reach(cell));
                    for _ in 0..90 {
                        let mid = 0.5 * (u + v);
                        if (hold(mid, m) - n as f64) * (d0 - n as f64) > 0.0 {
                            u = mid;
                        } else {
                            v = mid;
                        }
                    }
                    let t = 0.5 * (u + v);
                    let (arm_c, _) = arms(t);
                    let x = turn_for(w, arm_c, m as f64);
                    if edge_mark(w, t, x).is_none() {
                        continue;
                    }
                    out.push(Align {
                        reach: t,
                        turn: x,
                        ring: ring(w, arm_c, x),
                        mark_c: m,
                        mark_e: n,
                        ends,
                        kinds,
                        slot,
                        flat: false,
                    });
                }
            }
            back = front;
        }
    }
    for &m in &marks {
        for (n, outer) in [(0i64, true), (full.round() as i64, false)] {
            let (ends, kinds, slot, _) = meeting(&sample, m, n, true);
            if ends < 3 {
                continue;
            }
            let gap = |t: f64| -> f64 {
                let (arm_c, arm_e) = arms(t);
                let r = ring(w, arm_c, turn_for(w, arm_c, m as f64));
                let (lo, hi) = band(w, arm_e);
                r - if outer { hi } else { lo }
            };
            let mut back = gap(reach(0));
            for cell in 1..=grid {
                let front = gap(reach(cell));
                if back * front < 0.0 {
                    let (mut u, mut v) = (reach(cell - 1), reach(cell));
                    for _ in 0..90 {
                        let mid = 0.5 * (u + v);
                        if gap(mid) * back > 0.0 {
                            u = mid;
                        } else {
                            v = mid;
                        }
                    }
                    let t = 0.5 * (u + v);
                    let (arm_c, _) = arms(t);
                    let x = turn_for(w, arm_c, m as f64);
                    out.push(Align {
                        reach: t,
                        turn: x,
                        ring: ring(w, arm_c, x),
                        mark_c: m,
                        mark_e: n,
                        ends,
                        kinds,
                        slot,
                        flat: true,
                    });
                }
                back = front;
            }
        }
    }
    out.sort_by(|p, q| p.reach.partial_cmp(&q.reach).unwrap());
    out
}

struct Node {
    ring: f64,
    angle: f64,
    curves: Vec<usize>,
    branches: usize,
}

fn arc(w: &Wheel, cs: &[Curve], k: usize, sign: f64, r: f64) -> Option<f64> {
    let (lo, hi) = band(w, cs[k].arm);
    if r < lo || r > hi {
        return None;
    }
    Some(cs[k].oct as f64 / 8.0 + sign * mark_at(w, cs[k].arm, r) / 8.0)
}

fn census(w: &Wheel, cs: &[Curve], grid: usize) -> Vec<Node> {
    let ports: Vec<(usize, f64)> = (0..cs.len()).flat_map(|k| [(k, 1.0), (k, -1.0)]).collect();
    let mut raw: Vec<(f64, f64, usize, usize)> = Vec::new();
    for i in 0..ports.len() {
        for j in (i + 1)..ports.len() {
            if ports[i].0 == ports[j].0 && ports[i].1 == ports[j].1 {
                continue;
            }
            let (bi, bj) = (band(w, cs[ports[i].0].arm), band(w, cs[ports[j].0].arm));
            let (lo, hi) = (bi.0.max(bj.0), bi.1.min(bj.1));
            if hi <= lo {
                continue;
            }
            let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
            let ray = |s: f64| mid - half * (PI * s).cos();
            let gap = |s: f64| {
                let r = ray(s).clamp(lo, hi);
                let p = arc(w, cs, ports[i].0, ports[i].1, r);
                let q = arc(w, cs, ports[j].0, ports[j].1, r);
                match (p, q) {
                    (Some(p), Some(q)) => Some(p - q),
                    _ => None,
                }
            };
            let mut back = match gap(0.0) {
                Some(d) => d,
                None => continue,
            };
            for cell in 1..=grid {
                let s1 = cell as f64 / grid as f64;
                let s0 = (cell - 1) as f64 / grid as f64;
                let front = match gap(s1) {
                    Some(d) => d,
                    None => continue,
                };
                let (dlo, dhi) = (back.min(front), back.max(front));
                for goal in (dlo.floor() as i64 + 1)..=(dhi.ceil() as i64) {
                    let goal = goal as f64;
                    if goal <= dlo || goal >= dhi {
                        continue;
                    }
                    let (mut u, mut v) = (s0, s1);
                    for _ in 0..90 {
                        let mid = 0.5 * (u + v);
                        match gap(mid) {
                            Some(d) => {
                                if (d - goal) * (back - goal) > 0.0 {
                                    u = mid;
                                } else {
                                    v = mid;
                                }
                            }
                            None => break,
                        }
                    }
                    let seat = 0.5 * (u + v);
                    if !(1e-11..=1.0 - 1e-11).contains(&seat) {
                        continue;
                    }
                    let r = ray(seat).clamp(lo, hi);
                    if let Some(p) = arc(w, cs, ports[i].0, ports[i].1, r) {
                        raw.push((r, p - p.floor(), i, j));
                    }
                }
                back = front;
            }
        }
    }
    let reach = cs.iter().map(|c| band(w, c.arm).1).fold(0.0_f64, f64::max);
    raw.sort_by(|p, q| p.0.partial_cmp(&q.0).unwrap());
    let mut out: Vec<Node> = Vec::new();
    let mut pool: Vec<Vec<(f64, f64, usize, usize)>> = Vec::new();
    for row in raw {
        let hit = pool.iter().position(|g: &Vec<(f64, f64, usize, usize)>| {
            (g[0].0 - row.0).abs() < 1e-9 * reach && {
                let d = (g[0].1 - row.1).abs();
                d.min(1.0 - d) < 1e-9
            }
        });
        match hit {
            Some(k) => pool[k].push(row),
            None => pool.push(vec![row]),
        }
    }
    for group in pool {
        let mut seen: Vec<usize> = Vec::new();
        let mut ends: Vec<usize> = Vec::new();
        for row in &group {
            for port in [row.2, row.3] {
                if !ends.contains(&port) {
                    ends.push(port);
                }
                if !seen.contains(&ports[port].0) {
                    seen.push(ports[port].0);
                }
            }
        }
        seen.sort_unstable();
        out.push(Node {
            ring: group[0].0,
            angle: group[0].1,
            branches: ends.len(),
            curves: seen,
        });
    }
    out
}

const CARPET: [u8; 9] = [1, 1, 1, 1, 0, 1, 1, 1, 1];

fn near(px: f64, py: f64, x0: f64, y0: f64, x1: f64, y1: f64) -> f64 {
    let (dx, dy) = (x1 - x0, y1 - y0);
    let len = dx * dx + dy * dy;
    let s = if len > 0.0 {
        (((px - x0) * dx + (py - y0) * dy) / len).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (px - x0 - s * dx).hypot(py - y0 - s * dy)
}

fn probe(w: &Wheel, t: f64, r: f64, angle: f64, samples: usize) -> Vec<f64> {
    let kind = if w.side < 0.0 { "in" } else { "out" };
    let path = track(kind, w.a, w.b, 4, 1).unwrap();
    let pens = pencils(&CARPET, 3, 3, "fill", t, 0.0, 1).unwrap();
    let pts = trace(&path, &pens, samples).unwrap();
    let hit = w.b as f64 * r;
    let th = angle * TAU / w.a as f64;
    let (px, py) = (hit * th.cos(), hit * th.sin());
    (0..pens.len())
        .map(|k| {
            let head = k * samples * 2;
            let mut best = f64::MAX;
            for i in 0..samples - 1 {
                best = best.min(near(
                    px,
                    py,
                    pts[head + 2 * i] as f64,
                    pts[head + 2 * i + 1] as f64,
                    pts[head + 2 * i + 2] as f64,
                    pts[head + 2 * i + 3] as f64,
                ));
            }
            best
        })
        .collect()
}

fn probe_wide(w: &Wheel, t: f64, r: f64, angle: f64, samples: usize) -> Vec<f64> {
    let kind = if w.side < 0.0 { "in" } else { "out" };
    let path = track(kind, w.a, w.b, 4, 1).unwrap();
    let pens = pencils(&CARPET, 3, 3, "fill", t, 0.0, 1).unwrap();
    let hit = w.b as f64 * r;
    let th = angle * TAU / w.a as f64;
    let (px, py) = (hit * th.cos(), hit * th.sin());
    pens.iter()
        .map(|pen| {
            let mut best = f64::MAX;
            let mut back = point(&path, pen, 0.0);
            for i in 1..samples {
                let front = point(&path, pen, path.total * i as f64 / (samples - 1) as f64);
                best = best.min(near(px, py, back.0, back.1, front.0, front.1));
                back = front;
            }
            best
        })
        .collect()
}

fn f32_floor(v: f64) -> f64 {
    let a = v as f32;
    0.5 * (f32::from_bits(a.to_bits() + 1) - a) as f64
}

fn residue(w: &Wheel, arm: f64, x: f64, mark: i64) -> f64 {
    let om = (x.cos(), x.sin());
    let mut left = (1.0, 0.0);
    let mut right = (1.0, 0.0);
    let pow = w.a as i64 + 2 * w.b as i64 * w.side as i64;
    let base = (w.amp + arm * om.0, arm * om.1);
    let flip = (w.amp * om.0 + arm, w.amp * om.1);
    for _ in 0..w.a {
        left = (
            left.0 * base.0 - left.1 * base.1,
            left.0 * base.1 + left.1 * base.0,
        );
        right = (
            right.0 * flip.0 - right.1 * flip.1,
            right.0 * flip.1 + right.1 * flip.0,
        );
    }
    let spin = pow as f64 * x;
    left = (
        left.0 * spin.cos() - left.1 * spin.sin(),
        left.0 * spin.sin() + left.1 * spin.cos(),
    );
    let quarter = mark.rem_euclid(4);
    let unit = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)][quarter as usize];
    let want = (
        right.0 * unit.0 - right.1 * unit.1,
        right.0 * unit.1 + right.1 * unit.0,
    );
    ((left.0 - want.0).hypot(left.1 - want.1)) / (left.0.hypot(left.1)).max(1e-12)
}

fn defect(w: &Wheel, t: f64, m: f64, n: f64) -> f64 {
    let (arm_c, arm_e) = (2.0 * t / 3.0, 2.0_f64.sqrt() * t / 3.0);
    let r = ring(w, arm_c, turn_for(w, arm_c, m));
    mark_at(w, arm_e, r) - n
}

fn deck() -> Vec<(usize, usize, bool)> {
    let mut out = Vec::new();
    for b in 1..=8usize {
        for a in (b + 1)..=13usize {
            let (mut p, mut q) = (a, b);
            while q > 0 {
                let r = p % q;
                p = q;
                q = r;
            }
            if p != 1 {
                continue;
            }
            for inside in [true, false] {
                out.push((a, b, inside));
            }
        }
    }
    out
}

fn reduction() {
    println!("THE REDUCTION");
    let seats = carpet();
    let (mut worst, mut counts) = (0.0_f64, 0usize);
    for &(a, b, inside) in &deck() {
        let w = wheel(a, b, inside);
        let kind = if inside { "in" } else { "out" };
        let path = track(kind, a, b, 4, 1).unwrap();
        let pens = pencils(&CARPET, 3, 3, "fill", 0.83, 0.0, 1).unwrap();
        for pen in &pens {
            for k in 0..29 {
                let psi = TAU * k as f64 / 29.0;
                let (x, y) = point(&path, pen, path.total * psi / TAU);
                let spin = w.b as f64 * psi;
                let seat = (w.b as f64 + w.side * w.a as f64) * psi;
                let zx = b as f64 * (w.amp * spin.cos() + pen.x * seat.cos() - pen.y * seat.sin());
                let zy = b as f64 * (w.amp * spin.sin() + pen.x * seat.sin() + pen.y * seat.cos());
                worst = worst.max((x - zx).hypot(y - zy));
            }
        }
        let (cs, _) = curves(&w, &seats, 0.83);
        if cs.len() == distinct(&path, &pens, true) {
            counts += 1;
        }
    }
    println!("  z(psi) = r e^(i b psi) (A + p e^(i eps a psi)), A = abs(a/b + eps), at reach 0.83 and 29 phases a seat, worst gap to spirograph::point {worst:.3e}");
    println!(
        "  offset classes against spirograph::distinct: {counts} of {} tracks agree",
        deck().len()
    );
    let w = wheel(7, 3, true);
    let path = track("in", 7, 3, 4, 1).unwrap();
    let (alpha, arm) = (0.3_f64, 0.5_f64);
    let seat = |turn: f64| Pencil {
        x: arm * turn.cos(),
        y: arm * turn.sin(),
        seat: (0, 0),
        kind: Kind::Fill,
    };
    let (plain, turned) = (seat(0.0), seat(alpha));
    let spin = w.b as f64 * alpha / w.a as f64;
    let mut gaps = [0.0_f64; 2];
    for (k, lean) in [-1.0_f64, 1.0].iter().enumerate() {
        let shift = lean * alpha / (w.side * w.a as f64);
        let angle = lean * w.side * spin;
        for j in 0..3 {
            let psi = TAU * (0.17 + 0.31 * j as f64);
            let (x0, y0) = point(&path, &plain, path.total * psi / TAU);
            let (x1, y1) = point(
                &path,
                &turned,
                path.total * (psi + shift).rem_euclid(TAU) / TAU,
            );
            let (c, sn) = (angle.cos(), angle.sin());
            gaps[k] = gaps[k].max((x1 - (x0 * c - y0 * sn)).hypot(y1 - (x0 * sn + y0 * c)));
        }
    }
    println!(
        "  turning the seat by {alpha} on 7/3 inside turns the curve by minus b eps alpha over a: worst gap {:.3e}, against {:.3e} for the plus sign",
        gaps[0], gaps[1]
    );
}

fn tally(nodes: &[Node], count: usize) -> (usize, usize, usize, usize, usize) {
    let mut own = vec![0usize; count];
    let mut duo = vec![0usize; count * count];
    let mut branch = 0;
    for node in nodes {
        branch = branch.max(node.branches);
        if node.curves.len() == 1 {
            own[node.curves[0]] += 1;
        } else if node.curves.len() == 2 {
            duo[node.curves[0] * count + node.curves[1]] += 1;
        }
    }
    let pairs: Vec<usize> = (0..count)
        .flat_map(|i| ((i + 1)..count).map(move |j| (i, j)))
        .map(|(i, j)| duo[i * count + j])
        .collect();
    (
        own.iter().cloned().min().unwrap_or(0),
        own.iter().cloned().max().unwrap_or(0),
        pairs.iter().cloned().min().unwrap_or(0),
        pairs.iter().cloned().max().unwrap_or(0),
        branch,
    )
}

fn mark_law() {
    println!("THE MARK");
    println!("  m(x) = 4a(b eps x / a + arg(A + p e^(ix)))/pi runs 0 to 4 b eps while abs(p) < A, and is one to one for abs(p) < min(1, A)");
    println!(
        "  a/b side  reach  abs(p)      A    tame  self/curve  pair/pair  branch   a(b-1)  2ab"
    );
    for &(a, b, inside, t) in &[
        (7usize, 3usize, true, 0.83),
        (7, 3, false, 0.83),
        (5, 2, true, 0.83),
        (11, 4, true, 0.83),
        (5, 1, false, 0.83),
        (13, 7, true, 0.83),
        (7, 5, true, 0.585),
        (7, 5, true, 0.615),
        (7, 5, true, 0.83),
        (4, 3, true, 0.49),
        (4, 3, true, 0.51),
        (9, 8, true, 0.83),
    ] {
        let w = wheel(a, b, inside);
        let (cs, _) = curves(&w, &carpet(), t);
        let nodes = census(&w, &cs, 1500);
        let (s0, s1, p0, p1, branch) = tally(&nodes, cs.len());
        let arm = cs.iter().map(|c| c.arm).fold(0.0_f64, f64::max);
        let tame = if arm < 1.0_f64.min(w.amp) {
            "yes"
        } else {
            "no "
        };
        let side = if inside { "in " } else { "out" };
        println!(
            "  {a}/{b} {side} {t:6.3} {arm:6.4} {:6.4}  {tame}   {s0:2}..{s1:2}      {p0:2}..{p1:2}      {branch}    {:5}  {:4}",
            w.amp,
            a * (b - 1),
            2 * a * b
        );
    }
}

fn carpet_table() {
    println!("THE ALIGNMENT");
    let w = wheel(7, 3, true);
    let (tlo, thi) = (0.05_f64, tame(&w).min(1.4999));
    let found = alignments(&w, tlo, thi, 6000);
    println!("  7/3 inside, carpet fills, the window abs(p) < min(1, A) being the whole of abs(p) < 1 here: {} alignment reaches", found.len());
    let flats = found.iter().filter(|h| h.flat).count();
    let quads = found.iter().filter(|h| !h.flat && h.kinds == 4).count();
    println!("  the scan runs from reach {tlo} to {thi:.4}, which on the corner seats is abs(p) from {:.6} to {:.6}", 2.0 * tlo / 3.0, 2.0 * thi / 3.0);
    println!("  {} of them transversal and {flats} tangential, the tangential ones meeting an edge curve at its own apex", found.len() - flats);
    println!("  every transversal meeting carries four branches; {quads} of the {} carry four distinct curves and {} carry three, one curve bringing two branches as its own self crossing lands on the meeting", found.len() - flats, found.len() - flats - quads);
    println!(
        "  reach          abs(p) corner    ring     m_c  m_e  branches curves  kind   law residual"
    );
    for hit in &found {
        let arm_c = 2.0 * hit.reach / 3.0;
        let arm_e = 2.0_f64.sqrt() * hit.reach / 3.0;
        let turn_e = if hit.flat {
            if hit.mark_e == 0 {
                0.0
            } else {
                PI
            }
        } else {
            ((hit.ring * hit.ring - w.amp * w.amp - arm_e * arm_e) / (2.0 * w.amp * arm_e))
                .clamp(-1.0, 1.0)
                .acos()
        };
        let worst =
            residue(&w, arm_c, hit.turn, hit.mark_c).max(residue(&w, arm_e, turn_e, hit.mark_e));
        println!(
            "  {:.12}   {:.12}  {:.6}  {:3}  {:3}  {:6}  {:5}  {:5}  {:.2e}",
            hit.reach,
            arm_c,
            hit.ring,
            hit.mark_c,
            hit.mark_e,
            hit.ends,
            hit.kinds,
            if hit.flat { "apex" } else { "cross" },
            worst
        );
    }
}

fn node_of(cs: &[Curve], hit: &Align) -> Node {
    let (ends, _, _, seen) = meeting(cs, hit.mark_c, hit.mark_e, hit.flat);
    Node {
        ring: hit.ring,
        angle: hit.slot as f64 / 8.0,
        branches: ends,
        curves: seen,
    }
}

fn pinch(w: &Wheel, t: f64) -> (Vec<Node>, usize, usize) {
    let (cs, _) = curves(w, &carpet(), t);
    let nodes = census(w, &cs, 8000);
    let deep = nodes.iter().filter(|n| n.branches > 2).count();
    let branch = nodes.iter().map(|n| n.branches).max().unwrap_or(0);
    (nodes, deep, branch)
}

fn seat_hits(w: &Wheel, t: f64, node: &Node, samples: usize) -> (f64, f64, usize) {
    seat_gaps(w, t, node, probe(w, t, node.ring, node.angle, samples))
}

fn seat_wide(w: &Wheel, t: f64, node: &Node, samples: usize) -> (f64, f64, usize) {
    seat_gaps(w, t, node, probe_wide(w, t, node.ring, node.angle, samples))
}

fn seat_gaps(w: &Wheel, t: f64, node: &Node, gaps: Vec<f64>) -> (f64, f64, usize) {
    let (cs, owner) = curves(w, &carpet(), t);
    let _ = cs;
    let want: Vec<usize> = (0..owner.len())
        .filter(|&k| node.curves.contains(&owner[k]))
        .collect();
    let near = want.iter().map(|&k| gaps[k]).fold(0.0_f64, f64::max);
    let far = (0..gaps.len())
        .filter(|k| !want.contains(k))
        .map(|k| gaps[k])
        .fold(f64::INFINITY, f64::min);
    (near, far, want.len())
}

fn witness() {
    println!("THE WITNESS");
    let w = wheel(7, 3, true);
    let (mut lo, mut hi) = (0.78_f64, 0.80_f64);
    let sign = defect(&w, lo, -6.0, -9.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if defect(&w, mid, -6.0, -9.0) * sign > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = 0.5 * (lo + hi);
    println!(
        "  7/3 inside, marks (-6, -9): reach {t:.12}, bracket width {:.1e}",
        hi - lo
    );
    println!(
        "  abs(p) corner {:.12}, ring {:.12} wheel radii",
        2.0 * t / 3.0,
        ring(&w, 2.0 * t / 3.0, turn_for(&w, 2.0 * t / 3.0, -6.0))
    );
    let (nodes, deep, branch) = pinch(&w, t);
    println!(
        "  census at the reach: {} classes, {deep} of them with {branch} branches, the rest 2",
        nodes.len()
    );
    let node = nodes.iter().find(|n| n.branches > 2).unwrap();
    let (cs, _) = curves(&w, &carpet(), t);
    let corners = node.curves.iter().filter(|&&k| cs[k].corner).count();
    println!(
        "  the class at ring {:.9}, angle {:.9} of a turn over a: {corners} corner curves and {} edge curves",
        node.ring,
        node.angle,
        node.curves.len() - corners
    );
    let floor = f32_floor(w.b as f64 * node.ring);
    println!("  half an f32 step at that radius is {floor:.3e}, the floor any read of the f32 trace can reach");
    println!("  samples   worst gap in f64 over the meeting seats   the same read off the f32 trace   least gap over the others");
    for samples in [2000usize, 8000, 32000, 128000] {
        let (wide, far, count) = seat_wide(&w, t, node, samples);
        let thin = seat_hits(&w, t, node, samples).0;
        println!("  {samples:7}   {wide:.3e} over {count} seats                {thin:.3e}                    {far:.3e}");
    }
    for shift in [-0.01_f64, 0.01] {
        let control = t + shift;
        let (plain, deep, branch) = pinch(&w, control);
        let (near, far, count) = seat_hits(&w, control, node, 8000);
        println!(
            "  control reach {control:.6}: {} classes, {deep} past 2 branches, deepest {branch}; the same point now {near:.3e} from {count} seats, {far:.3e} from the rest",
            plain.len()
        );
    }
    let all = alignments(&w, 0.05, tame(&w).min(1.4999), 6000);
    let sample = curves(&w, &carpet(), 0.8).0;
    println!("  the tangential family, a corner pair meeting an edge curve at its own apex:");
    for hit in all.iter().filter(|h| h.flat) {
        let node = node_of(&sample, hit);
        let (wide, far, count) = seat_wide(&w, hit.reach, &node, 8000);
        println!(
            "  reach {:.12}, marks ({}, {}), ring {:.9}: {count} seats within {wide:.3e} of the point in f64 and {:.3e} off the f32 trace, floor {:.3e}, the rest no nearer than {far:.3e}",
            hit.reach,
            hit.mark_c,
            hit.mark_e,
            hit.ring,
            seat_hits(&w, hit.reach, &node, 8000).0,
            f32_floor(w.b as f64 * hit.ring)
        );
    }
    let plain = pinch(&w, t - 0.01).0.len();
    println!(
        "  the drop: {plain} classes and so {} nodes at a generic reach, {} classes and {} nodes at the alignment, five double points swallowed by each meeting",
        plain * w.a,
        nodes.len(),
        nodes.len() * w.a
    );
}

fn gcd4(b: usize) -> usize {
    let (mut p, mut q) = (b, 4usize);
    while q > 0 {
        let r = p % q;
        p = q;
        q = r;
    }
    p
}

#[derive(Clone)]
struct Band {
    tracks: usize,
    aligns: usize,
    bent: usize,
    held: usize,
    deep: usize,
    wide: f64,
    thin: f64,
    floor: f64,
    far: f64,
    breach: usize,
}

fn survey() {
    println!("THE SWEEP");
    let mut rows: Vec<Band> = (0..8)
        .map(|_| Band {
            tracks: 0,
            aligns: 0,
            bent: 0,
            held: 0,
            deep: 0,
            wide: 0.0,
            thin: 0.0,
            floor: 0.0,
            far: f64::INFINITY,
            breach: 0,
        })
        .collect();
    let mut book: Vec<(usize, usize, bool, usize, usize)> = Vec::new();
    for &(a, b, inside) in &deck() {
        let w = wheel(a, b, inside);
        let cap = tame(&w).min(1.4999);
        let found = alignments(&w, 0.05, cap, 6000);
        book.push((b, a, inside, found.len(), 4 / gcd4(b)));
        let row = &mut rows[b - 1];
        row.tracks += 1;
        row.aligns += found.len();
        let shape = curves(&w, &carpet(), 0.5 * (0.05 + cap)).0;
        for hit in found.iter() {
            let node = if hit.flat {
                row.bent += 1;
                node_of(&shape, hit)
            } else {
                let (nodes, deep, branch) = pinch(&w, hit.reach);
                if branch < 3 {
                    continue;
                }
                row.deep += deep;
                row.held += 1;
                let seen = nodes.iter().find(|n| n.branches > 2).unwrap();
                Node {
                    ring: seen.ring,
                    angle: seen.angle,
                    branches: seen.branches,
                    curves: seen.curves.clone(),
                }
            };
            let (wide, far, _) = seat_wide(&w, hit.reach, &node, 8000);
            row.wide = row.wide.max(wide);
            row.thin = row.thin.max(seat_hits(&w, hit.reach, &node, 8000).0);
            row.floor = row.floor.max(f32_floor(w.b as f64 * node.ring));
            if far.is_finite() {
                row.far = row.far.min(far);
            }
        }
        let control = match found.len() {
            0 => 0.7,
            1 => found[0].reach + 0.013,
            _ => 0.5 * (found[0].reach + found[1].reach),
        };
        if control > 0.05 && control < 1.4999 {
            let (_, _, branch) = pinch(&w, control);
            if branch > 2 {
                row.breach += 1;
            }
        }
    }
    println!("  b   tracks  alignments  tangential  held in the census  classes past two branches  worst gap f64  worst gap f32  f32 floor  least other gap  control breaches");
    for (k, row) in rows.iter().enumerate() {
        let other = if row.far.is_finite() {
            format!("{:.3e}", row.far)
        } else {
            "no other seat".to_string()
        };
        println!(
            "  {}   {:5}   {:8}   {:8}   {:14}   {:14}   {:.3e}      {:.3e}      {:.3e}  {other:13}     {}",
            k + 1,
            row.tracks,
            row.aligns,
            row.bent,
            row.held,
            row.deep,
            row.wide,
            row.thin,
            row.floor,
            row.breach
        );
    }
    println!("  alignment reaches per track over the window abs(p) < min(1, A), a:count");
    for b in 1..=8usize {
        for inside in [true, false] {
            let line: Vec<String> = book
                .iter()
                .filter(|r| r.0 == b && r.2 == inside)
                .map(|r| format!("{}:{}", r.1, r.3))
                .collect();
            if line.is_empty() {
                continue;
            }
            println!(
                "  b {b} {}  classes {}   {}",
                if inside { "in " } else { "out" },
                4 / gcd4(b),
                line.join(" ")
            );
        }
    }
    let total: usize = rows.iter().map(|r| r.aligns).sum();
    let breach: usize = rows.iter().map(|r| r.breach).sum();
    let held: usize = rows.iter().map(|r| r.held).sum();
    let bent: usize = rows.iter().map(|r| r.bent).sum();
    println!("  {total} alignment reaches over {} circle tracks, {bent} of them tangential; every one is read against mrly_math::spirograph::point and mrly_math::spirograph::trace and holds, {held} of the {} transversal ones also show a class past two branches in the census, and {breach} control reaches show one", deck().len(), total - bent);
    for (a, b, inside) in [
        (7usize, 3usize, true),
        (12, 7, true),
        (13, 7, true),
        (9, 5, true),
        (13, 8, true),
        (13, 7, false),
    ] {
        let probe = wheel(a, b, inside);
        let cap = tame(&probe).min(1.4999);
        let side = if inside { "inside" } else { "outside" };
        println!(
            "  grid check on {a}/{b} {side}, reach up to {cap:.4}: {} at 6000, {} at 24000",
            alignments(&probe, 0.05, cap, 6000).len(),
            alignments(&probe, 0.05, cap, 24000).len()
        );
    }
}

fn main() {
    reduction();
    mark_law();
    carpet_table();
    witness();
    survey();
}
