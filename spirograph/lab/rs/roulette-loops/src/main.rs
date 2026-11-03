use mrly_math::spirograph::{frame, trace, track, Kind, Pencil};
use std::f64::consts::{PI, TAU};

fn rim(a: i64, b: i64, out: bool) -> i64 {
    if out {
        a + b
    } else {
        a - b
    }
}

fn skew(a: i64, b: i64, out: bool) -> i64 {
    if out {
        a + 2 * b
    } else {
        a - 2 * b
    }
}

fn gcd(mut x: i64, mut y: i64) -> i64 {
    while y != 0 {
        let t = x % y;
        x = y;
        y = t;
    }
    x
}

fn bisect(f: &dyn Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let sign = f(lo).signum();
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid).signum() == sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn sweep(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64, steps: usize) -> Vec<f64> {
    let mut out = Vec::new();
    let mut prev = f(lo);
    let mut left = lo;
    for k in 1..=steps {
        let x = lo + (hi - lo) * k as f64 / steps as f64;
        let cur = f(x);
        if prev * cur < 0.0 {
            out.push(bisect(f, left, x));
        }
        prev = cur;
        left = x;
    }
    out
}

fn arm_roots(a: i64, b: i64, out: bool, lam: f64, eps: f64, steps: usize) -> Vec<f64> {
    let (r, bf) = (rim(a, b, out) as f64, b as f64);
    let f = |d: f64| r * (bf * d).sin() - eps * lam * bf * (r * d).sin();
    sweep(&f, 1e-7, PI - 1e-7, steps)
}

fn arms(a: i64, b: i64, out: bool, lam: f64, steps: usize) -> usize {
    arm_roots(a, b, out, lam, 1.0, steps).len() + arm_roots(a, b, out, lam, -1.0, steps).len()
}

fn crossings(a: i64, b: i64, out: bool, lam: f64, steps: usize) -> usize {
    a as usize * arms(a, b, out, lam, steps) / 2
}

fn tangent_angles(a: i64, b: i64, out: bool, steps: usize) -> Vec<f64> {
    let m = skew(a, b, out);
    let (af, mf) = (a as f64, m as f64);
    let g = |d: f64| af * (mf * d).sin() - mf * (af * d).sin();
    let lo = 0.4 / a.max(m.abs()) as f64;
    let mut out_angles = vec![0.0];
    for d in sweep(&g, lo, PI / 2.0, steps) {
        if PI / 2.0 - d > 1e-6 {
            out_angles.push(d);
            out_angles.push(PI - d);
        }
    }
    if g(PI / 2.0).abs() < 1e-9 {
        out_angles.push(PI / 2.0);
    }
    out_angles
}

fn reach_at(a: i64, b: i64, out: bool, d: f64) -> f64 {
    let (r, bf) = (rim(a, b, out) as f64, b as f64);
    let (s, c) = ((r * d).sin(), (r * d).cos());
    if s.abs() > c.abs() {
        (r * (bf * d).sin() / (bf * s)).abs()
    } else {
        ((bf * d).cos() / c).abs()
    }
}

fn steps_of(a: i64, b: i64, out: bool, steps: usize) -> Vec<f64> {
    let mut out_reach: Vec<f64> = tangent_angles(a, b, out, steps)
        .into_iter()
        .map(|d| reach_at(a, b, out, d))
        .collect();
    out_reach.sort_by(|x, y| x.partial_cmp(y).unwrap());
    out_reach
}

fn law(a: i64, b: i64, out: bool, lam: f64, steps: &[f64]) -> i64 {
    let sign = skew(a, b, out).signum();
    let past = steps.iter().filter(|&&t| t < lam).count() as i64;
    a * (b - 1) + sign * a * past
}

fn parity(e: i64) -> i64 {
    if e % 2 == 0 {
        1
    } else {
        -1
    }
}

fn turns(b: i64, r: i64) -> i64 {
    let flip = (r - b).signum();
    let mut marks: Vec<(i64, i64, i64)> = Vec::new();
    for j in 1..r {
        marks.push((j, r, flip * parity(j + 1 + b * j / r)));
    }
    for k in 1..b {
        marks.push((k, b, flip * parity(k + r * k / b)));
    }
    marks.sort_by_key(|&(n, d, _)| (n * (b * r / d), d));
    let mut walk = vec![1_i64];
    walk.extend(marks.iter().map(|&(_, _, s)| s));
    1 + walk.windows(2).filter(|w| w[0] != w[1]).count() as i64
}

fn seat(p: (f64, f64)) -> Pencil {
    Pencil {
        x: p.0,
        y: p.1,
        seat: (0, 0),
        kind: Kind::Fill,
    }
}

fn draw(
    a: i64,
    b: i64,
    out: bool,
    pens: &[(f64, f64)],
    samples: usize,
) -> (Vec<Vec<f64>>, [f64; 4]) {
    let kind = if out { "out" } else { "in" };
    let path = track(kind, a as usize, b as usize, 4, 1).unwrap();
    let pencils: Vec<Pencil> = pens.iter().map(|&p| seat(p)).collect();
    let flat = trace(&path, &pencils, samples).unwrap();
    let curves = pencils
        .iter()
        .enumerate()
        .map(|(k, _)| {
            flat[k * samples * 2..(k + 1) * samples * 2]
                .iter()
                .map(|&v| v as f64)
                .collect()
        })
        .collect();
    (curves, frame(&path, &pencils))
}

fn hits(s: [f64; 4], t: [f64; 4]) -> Option<(f64, f64)> {
    let [ax, ay, bx, by] = s;
    let [cx, cy, dx, dy] = t;
    let (ux, uy) = (bx - ax, by - ay);
    let (vx, vy) = (dx - cx, dy - cy);
    let det = ux * vy - uy * vx;
    if det == 0.0 {
        return None;
    }
    let (wx, wy) = (cx - ax, cy - ay);
    let t = (wx * vy - wy * vx) / det;
    let s = (wx * uy - wy * ux) / det;
    if (0.0..1.0).contains(&t) && (0.0..1.0).contains(&s) {
        Some((ax + t * ux, ay + t * uy))
    } else {
        None
    }
}

fn nodes(
    curves: &[Vec<f64>],
    box_: [f64; 4],
    grid: usize,
    axes: i64,
    tilt: f64,
) -> (usize, usize, f64, usize, usize) {
    let mut segs: Vec<[f64; 4]> = Vec::new();
    let mut owner: Vec<u32> = Vec::new();
    let mut index: Vec<u32> = Vec::new();
    for (k, curve) in curves.iter().enumerate() {
        let n = curve.len() / 2 - 1;
        for i in 0..n {
            segs.push([
                curve[2 * i],
                curve[2 * i + 1],
                curve[2 * i + 2],
                curve[2 * i + 3],
            ]);
            owner.push(k as u32);
            index.push(i as u32);
        }
    }
    let last: Vec<u32> = curves.iter().map(|c| (c.len() / 2 - 2) as u32).collect();
    let (wx, wy) = (box_[2] - box_[0], box_[3] - box_[1]);
    let mut bins: Vec<Vec<u32>> = vec![Vec::new(); grid * grid];
    let put = |v: f64, lo: f64, w: f64| {
        (((v - lo) / w * grid as f64) as isize).clamp(0, grid as isize - 1) as usize
    };
    for (id, s) in segs.iter().enumerate() {
        let (x0, x1) = (
            put(s[0].min(s[2]), box_[0], wx),
            put(s[0].max(s[2]), box_[0], wx),
        );
        let (y0, y1) = (
            put(s[1].min(s[3]), box_[1], wy),
            put(s[1].max(s[3]), box_[1], wy),
        );
        for cx in x0..=x1 {
            for cy in y0..=y1 {
                bins[cx * grid + cy].push(id as u32);
            }
        }
    }
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let (mut own, mut between, mut off) = (0, 0, 0.0_f64);
    let mut mine: Vec<(f64, f64)> = Vec::new();
    let mut yours: Vec<(f64, f64)> = Vec::new();
    for bin in &bins {
        for u in 0..bin.len() {
            for v in u + 1..bin.len() {
                let (i, j) = (bin[u].min(bin[v]) as usize, bin[u].max(bin[v]) as usize);
                if owner[i] == owner[j] {
                    let (p, q) = (index[i], index[j]);
                    if q - p == 1 || (p == 0 && q == last[owner[i] as usize]) {
                        continue;
                    }
                }
                if !seen.insert((i as u64) << 32 | j as u64) {
                    continue;
                }
                if let Some(spot) = hits(segs[i], segs[j]) {
                    if owner[i] == owner[j] {
                        own += 1;
                        off = off.max(stray(spot.0, spot.1, axes, tilt));
                        mine.push(spot);
                    } else {
                        between += 1;
                        yours.push(spot);
                    }
                }
            }
        }
    }
    let tol = (wx.max(wy)) * 2e-3;
    (own, between, off, thin(&mine, tol), thin(&yours, tol))
}

fn thin(spots: &[(f64, f64)], tol: f64) -> usize {
    let mut kept: Vec<(f64, f64)> = Vec::new();
    for &spot in spots {
        if !kept
            .iter()
            .any(|k: &(f64, f64)| (k.0 - spot.0).hypot(k.1 - spot.1) < tol)
        {
            kept.push(spot);
        }
    }
    kept.len()
}

fn stray(x: f64, y: f64, axes: i64, tilt: f64) -> f64 {
    if axes <= 0 {
        return 0.0;
    }
    let (r, theta) = (x.hypot(y), y.atan2(x));
    (0..axes)
        .map(|k| (r * (theta - tilt - PI * k as f64 / axes as f64).sin()).abs())
        .fold(f64::MAX, f64::min)
}

#[derive(Clone, Copy)]
struct Surd {
    a: i128,
    b: i128,
    den: i128,
    root: i128,
}

fn square_part(mut d: i128) -> (i128, i128) {
    let (mut f, mut k) = (1_i128, 2_i128);
    while k * k <= d {
        while d % (k * k) == 0 {
            d /= k * k;
            f *= k;
        }
        k += 1;
    }
    (f, d)
}

impl Surd {
    fn tidy(self) -> Surd {
        let g = gcd(
            gcd(
                self.a.unsigned_abs().min(i64::MAX as u128) as i64,
                self.b.unsigned_abs().min(i64::MAX as u128) as i64,
            ),
            self.den.unsigned_abs().min(i64::MAX as u128) as i64,
        )
        .max(1) as i128;
        let sign = if self.den < 0 { -1 } else { 1 };
        Surd {
            a: sign * self.a / g,
            b: sign * self.b / g,
            den: sign * self.den / g,
            root: self.root,
        }
    }

    fn plus(self, other: Surd) -> Option<Surd> {
        Some(
            Surd {
                a: self
                    .a
                    .checked_mul(other.den)?
                    .checked_add(other.a.checked_mul(self.den)?)?,
                b: self
                    .b
                    .checked_mul(other.den)?
                    .checked_add(other.b.checked_mul(self.den)?)?,
                den: self.den.checked_mul(other.den)?,
                root: self.root,
            }
            .tidy(),
        )
    }

    fn times(self, other: Surd) -> Option<Surd> {
        Some(
            Surd {
                a: self
                    .a
                    .checked_mul(other.a)?
                    .checked_add(self.b.checked_mul(other.b)?.checked_mul(self.root)?)?,
                b: self
                    .a
                    .checked_mul(other.b)?
                    .checked_add(self.b.checked_mul(other.a)?)?,
                den: self.den.checked_mul(other.den)?,
                root: self.root,
            }
            .tidy(),
        )
    }

    fn over(self, other: Surd) -> Option<Surd> {
        let base = other
            .a
            .checked_mul(other.a)?
            .checked_sub(other.b.checked_mul(other.b)?.checked_mul(other.root)?)?;
        if base == 0 {
            return None;
        }
        self.times(Surd {
            a: other.den.checked_mul(other.a)?,
            b: -other.den.checked_mul(other.b)?,
            den: base,
            root: self.root,
        })
    }

    fn value(self) -> f64 {
        (self.a as f64 + self.b as f64 * (self.root as f64).sqrt()) / self.den as f64
    }

    fn text(self) -> String {
        if self.b == 0 {
            return format!("{}/{}", self.a, self.den);
        }
        format!(
            "({} {} {} sqrt {})/{}",
            self.a,
            if self.b < 0 { "-" } else { "+" },
            self.b.abs(),
            self.root,
            self.den
        )
    }
}

fn poly_surd(p: &[i128], u: Surd) -> Option<Surd> {
    let mut acc = Surd {
        a: 0,
        b: 0,
        den: 1,
        root: u.root,
    };
    for &c in p.iter().rev() {
        acc = acc.times(u)?.plus(Surd {
            a: c,
            b: 0,
            den: 1,
            root: u.root,
        })?;
    }
    Some(acc)
}

fn quadratic_roots(q: &[i128]) -> Vec<Surd> {
    if q.len() != 3 {
        return Vec::new();
    }
    let (c, b, a) = (q[0], q[1], q[2]);
    let disc = b * b - 4 * a * c;
    if disc < 0 {
        return Vec::new();
    }
    let (f, d) = square_part(disc);
    [1_i128, -1]
        .iter()
        .map(|&s| {
            Surd {
                a: -b,
                b: s * f,
                den: 2 * a,
                root: d,
            }
            .tidy()
        })
        .collect()
}

fn pair_arms(a: i64, b: i64, out: bool, p: (f64, f64), q: (f64, f64), steps: usize) -> usize {
    let (r, bf) = (rim(a, b, out) as f64, b as f64);
    let sw = if out { 1.0 } else { -1.0 };
    let mid = (0.5 * (p.0 + q.0), 0.5 * (p.1 + q.1));
    let gap = (0.5 * (p.0 - q.0), 0.5 * (p.1 - q.1));
    let big = mid.0 * mid.0 + mid.1 * mid.1;
    let small = gap.0 * gap.0 + gap.1 * gap.1;
    let tilt = mid.1 * gap.0 - mid.0 * gap.1;
    let f = |d: f64| {
        let (s, c) = ((r * d).sin(), (r * d).cos());
        let load = big * s * s + small * c * c - sw * tilt * (2.0 * r * d).sin();
        let w = (bf * d).sin();
        r * r * w * w - bf * bf * load
    };
    sweep(&f, 0.0, TAU, steps).len()
}

fn pair_count(a: i64, b: i64, out: bool, p: (f64, f64), q: (f64, f64), steps: usize) -> usize {
    a as usize * pair_arms(a, b, out, p, q, steps) / 2
}

fn poly_add(p: &[i128], q: &[i128], sp: i128, sq: i128) -> Vec<i128> {
    let n = p.len().max(q.len());
    (0..n)
        .map(|k| sp * p.get(k).copied().unwrap_or(0) + sq * q.get(k).copied().unwrap_or(0))
        .collect()
}

fn poly_shift(p: &[i128]) -> Vec<i128> {
    let mut out = vec![0];
    out.extend_from_slice(p);
    out
}

fn poly_trim(mut p: Vec<i128>) -> Vec<i128> {
    while p.len() > 1 && *p.last().unwrap() == 0 {
        p.pop();
    }
    p
}

fn angle_polys(top: i64) -> Vec<(Vec<i128>, Vec<i128>)> {
    let mut out = vec![(vec![0_i128], vec![1_i128])];
    for k in 0..top {
        let (s, c) = (out[k as usize].0.clone(), out[k as usize].1.clone());
        let next = if k % 2 == 0 {
            let sin = poly_add(&c, &poly_add(&s, &poly_shift(&s), 1, -1), 1, 1);
            let cos = poly_add(&c, &poly_shift(&s), 1, -1);
            (poly_trim(sin), poly_trim(cos))
        } else {
            let sin = poly_add(&s, &c, 1, 1);
            let cos = poly_add(
                &poly_add(&c, &poly_shift(&c), 1, -1),
                &poly_shift(&s),
                1,
                -1,
            );
            (poly_trim(sin), poly_trim(cos))
        };
        out.push(next);
    }
    out
}

fn poly_at(p: &[i128], u: f64) -> f64 {
    p.iter().rev().fold(0.0, |acc, &c| acc * u + c as f64)
}

fn poly_text(p: &[i128]) -> String {
    let mut parts = Vec::new();
    for (k, &c) in p.iter().enumerate().rev() {
        if c == 0 {
            continue;
        }
        parts.push(match k {
            0 => format!("{c}"),
            1 => format!("{c} u"),
            _ => format!("{c} u^{k}"),
        });
    }
    if parts.is_empty() {
        "0".to_string()
    } else {
        parts.join(" + ").replace("+ -", "- ")
    }
}

fn poly_mul(p: &[i128], q: &[i128]) -> Vec<i128> {
    let mut out = vec![0_i128; p.len() + q.len() - 1];
    for (i, &x) in p.iter().enumerate() {
        for (j, &y) in q.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    poly_trim(out)
}

fn cos_square(tab: &[(Vec<i128>, Vec<i128>)], k: i64) -> Vec<i128> {
    let part = &tab[k as usize].1;
    let sq = poly_mul(part, part);
    if k % 2 == 0 {
        sq
    } else {
        poly_mul(&[1, -1], &sq)
    }
}

fn reach_square(a: i64, b: i64, out: bool) -> (Vec<i128>, Vec<i128>) {
    let r = rim(a, b, out);
    let tab = angle_polys(b.max(r));
    (cos_square(&tab, b), cos_square(&tab, r))
}

fn rational_root(p: &[i128], u: f64) -> Option<(i128, i128)> {
    for den in 1_i128..=120 {
        let num = (u * den as f64).round() as i128;
        if num < 0 || gcd(num as i64, den as i64) != 1 {
            continue;
        }
        if (num as f64 / den as f64 - u).abs() > 1e-9 {
            continue;
        }
        let deg = p.len() - 1;
        let mut acc = 0_i128;
        let mut ok = true;
        for (k, &c) in p.iter().enumerate() {
            match num
                .checked_pow(k as u32)
                .and_then(|x| {
                    den.checked_pow((deg - k) as u32)
                        .and_then(|y| x.checked_mul(y))
                })
                .and_then(|z| z.checked_mul(c))
                .and_then(|z| acc.checked_add(z))
            {
                Some(v) => acc = v,
                None => ok = false,
            }
        }
        if ok && acc == 0 {
            return Some((num, den));
        }
    }
    None
}

fn poly_ratio(p: &[i128], num: i128, den: i128) -> Option<(i128, i128)> {
    let deg = p.len() - 1;
    let mut acc = 0_i128;
    for (k, &c) in p.iter().enumerate() {
        let term = num
            .checked_pow(k as u32)?
            .checked_mul(den.checked_pow((deg - k) as u32)?)?
            .checked_mul(c)?;
        acc = acc.checked_add(term)?;
    }
    Some((acc, den.checked_pow(deg as u32)?))
}

fn tangency_poly(a: i64, m: i64) -> Vec<i128> {
    let m = m.abs();
    let tab = angle_polys(a.max(m));
    let (pa, pm) = (&tab[a as usize].0, &tab[m as usize].0);
    poly_trim(poly_add(pm, pa, a as i128, -(m as i128)))
}

fn cases() -> Vec<(i64, i64, bool)> {
    let mut out = Vec::new();
    for b in 1..=6 {
        for a in b + 1..=11 {
            if gcd(a, b) != 1 {
                continue;
            }
            if skew(a, b, false) != 0 {
                out.push((a, b, false));
            }
            out.push((a, b, true));
        }
    }
    out
}

fn side(out: bool) -> f64 {
    if out {
        1.0
    } else {
        -1.0
    }
}

fn tag(out: bool) -> &'static str {
    if out {
        "out"
    } else {
        "in "
    }
}

fn walk(a: i64, b: i64, out: bool, steps: &[f64]) -> String {
    let mut text = format!("{}", a * (b - 1));
    let mut k = 0;
    while k < steps.len() {
        let mut n = 1;
        while k + n < steps.len() && (steps[k + n] - steps[k]).abs() < 1e-9 {
            n += 1;
        }
        k += n;
        text += &format!(
            " |{:.6}| {}",
            steps[k - 1],
            law(a, b, out, steps[k - 1] + 1e-9, steps)
        );
    }
    text
}

fn exact_reach(a: i64, b: i64, out: bool, num: i128, den: i128) -> Option<(i128, i128)> {
    let (top, bot) = reach_square(a, b, out);
    let (tn, td) = poly_ratio(&top, num, den)?;
    let (bn, bd) = poly_ratio(&bot, num, den)?;
    if bn == 0 {
        return None;
    }
    let (mut p, mut q) = (tn.checked_mul(bd)?, bn.checked_mul(td)?);
    if q < 0 {
        p = -p;
        q = -q;
    }
    let g = gcd(p.unsigned_abs().min(i64::MAX as u128) as i64, q as i64).max(1) as i128;
    Some((p / g, q / g))
}

fn axes() {
    println!(
        "THE AXES  every self crossing of one trochoid sits on one of the a mirror lines k pi / a"
    );
    for (a, b, out) in [
        (5, 1, false),
        (7, 2, false),
        (8, 3, false),
        (5, 2, true),
        (7, 3, true),
    ] {
        for lam in [1.4, 2.6, 3.9] {
            for alpha in [0.0_f64, 0.3] {
                let seat = (lam * alpha.cos(), lam * alpha.sin());
                let tilt = -side(out) * b as f64 * alpha / a as f64;
                let (curves, box_) = draw(a, b, out, &[seat], 24001);
                let (own, _, off, _, _) = nodes(&curves, box_, 260, a, tilt);
                let size = (box_[2] - box_[0]).max(box_[3] - box_[1]);
                println!(
                    "  {} {a}/{b} reach {lam} seat angle {alpha}  nodes {own}  worst stray {:.3e} of frame {size:.3}  relative {:.2e}",
                    tag(out),
                    off,
                    off / size
                );
            }
        }
    }
}

fn staircase() {
    println!("THE STAIRCASE  count on each step, tangency reaches between, angles counted with multiplicity");
    for &(a, b, out) in cases().iter().filter(|c| c.0 <= 8) {
        let steps = steps_of(a, b, out, 400_000);
        let m = skew(a, b, out);
        println!(
            "  {} {a}/{b} m {m} angles {} of {}  {}",
            tag(out),
            steps.len(),
            turns(b, rim(a, b, out)),
            walk(a, b, out, &steps)
        );
    }
    println!("THE ALGEBRA  u = sin^2 of the tangency angle, Q(u) = 0 its integer equation, lam^2 rational in u");
    for (a, b, out) in [
        (5, 1, false),
        (7, 2, false),
        (8, 3, false),
        (11, 4, false),
        (3, 1, true),
        (5, 1, true),
        (3, 2, true),
        (5, 2, true),
    ] {
        let m = skew(a, b, out);
        let p = tangency_poly(a, m);
        let q: Vec<i128> = p[1..].to_vec();
        let stem = p[0];
        let (top, bot) = reach_square(a, b, out);
        println!(
            "  {} {a}/{b}  P(0) = {stem}  Q(u) = {}  lam^2 = ({}) / ({})",
            tag(out),
            poly_text(&q),
            poly_text(&top),
            poly_text(&bot)
        );
        for d in tangent_angles(a, b, out, 400_000) {
            if d == 0.0 || d > PI / 2.0 + 1e-9 {
                continue;
            }
            let u = d.sin() * d.sin();
            let lam = reach_at(a, b, out, d);
            let ring = a % 2 == 0 && (u - 1.0).abs() < 1e-12;
            let root = if ring {
                Some((1, 1))
            } else {
                rational_root(&q, u)
            };
            let exact = if ring {
                let r = rim(a, b, out) as i128;
                Some((r * r, (b * b) as i128))
            } else {
                root.and_then(|(n, dn)| exact_reach(a, b, out, n, dn))
            };
            let mut surd = String::new();
            if exact.is_none() {
                for pick in quadratic_roots(&q) {
                    if (pick.value() - u).abs() > 1e-9 {
                        continue;
                    }
                    if let Some(v) = poly_surd(&top, pick)
                        .and_then(|t| poly_surd(&bot, pick).and_then(|w| t.over(w)))
                    {
                        surd = format!("= {} check {:.12}", v.text(), v.value());
                    }
                }
            }
            println!(
                "    u {u:.12}  {}  Q(u) {:.2e}  lam {lam:.12}  lam^2 {:.12} {}",
                match root {
                    Some((n, dn)) => format!("= {n}/{dn}"),
                    None => "irrational".to_string(),
                },
                if ring {
                    0.0
                } else {
                    poly_at(&q, u) / q.iter().map(|c| c.unsigned_abs() as f64).sum::<f64>()
                },
                lam * lam,
                match exact {
                    Some((n, dn)) => format!("= {n}/{dn}"),
                    None => surd,
                }
            );
        }
    }
}

fn signs(top: i64) {
    println!(
        "THE SIGNS  the threshold count from the exact critical values of the crossing integral"
    );
    let (mut pairs, mut wrong) = (0_usize, 0_usize);
    for r in 1..=top {
        for b in 1..=top {
            if b == r || gcd(b, r) != 1 {
                continue;
            }
            pairs += 1;
            if turns(b, r) != (r - b).abs() {
                wrong += 1;
                println!(
                    "    SIGNS b {b} rho {r} count {} want {}",
                    turns(b, r),
                    (r - b).abs()
                );
            }
        }
    }
    println!("  frequency pairs {pairs} to {top} failures {wrong}");
}

fn cutting(b: i64, r: i64) -> Vec<bool> {
    let mut marks: Vec<(i64, bool)> = Vec::new();
    for j in 1..r {
        marks.push((j * b, true));
    }
    for k in 1..b {
        marks.push((k * r, false));
    }
    marks.sort();
    marks.into_iter().map(|(_, v)| v).collect()
}

fn words(top: i64) {
    println!("THE WORDS  equal neighbours in the cutting word of the segment (0,0) to (rho,b)");
    let (mut pairs, mut wrong, mut loose) = (0_usize, 0_usize, 0_usize);
    for r in 1..=top {
        for b in 1..=top {
            if b == r || gcd(b, r) != 1 {
                continue;
            }
            pairs += 1;
            let word = cutting(b, r);
            let equal = word.windows(2).filter(|w| w[0] == w[1]).count() as i64;
            if equal != (r - b).abs() - 1 || equal != turns(b, r) - 1 {
                wrong += 1;
                println!(
                    "    WORDS b {b} rho {r} equal {equal} want {}",
                    (r - b).abs() - 1
                );
            }
            let major = r > b;
            let isolated = word.windows(2).all(|w| w[0] == major || w[1] == major);
            let ends = word.first().map_or(true, |&x| x == major)
                && word.last().map_or(true, |&x| x == major);
            if !isolated || !ends {
                loose += 1;
            }
        }
    }
    println!("  frequency pairs {pairs} to {top} failures {wrong} minority letter loose {loose}");
}

fn lift(b: i64, r: i64, d: f64) -> f64 {
    (r as f64 * (b as f64 * d).sin() / (b as f64 * (r as f64 * d).sin())).abs()
}

fn folds_of(b: i64, r: i64) -> Vec<f64> {
    let mut stops: Vec<f64> = Vec::new();
    for j in 1..r {
        stops.push(PI * j as f64 / r as f64);
    }
    for k in 1..b {
        stops.push(PI * k as f64 / b as f64);
    }
    stops.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let (p, q) = ((r + b) as f64, (r - b) as f64);
    let g = |d: f64| p * (q * d).sin() - q * (p * d).sin();
    let mut found = Vec::new();
    for w in stops.windows(2) {
        if g(w[0]) * g(w[1]) < 0.0 {
            found.push(bisect(&g, w[0], w[1]));
        }
    }
    found
}

fn level_count(pieces: &[(f64, f64)], t: f64) -> i64 {
    pieces.iter().filter(|&&(lo, hi)| lo < t && t < hi).count() as i64
}

fn folds(top: i64) {
    println!("THE FOLDS  roots of |H| = t counted on its monotone pieces, against the staircase");
    let (mut pairs, mut wrong, mut levels, mut closed, mut mins) =
        (0_usize, 0_usize, 0_usize, 0_usize, 0_usize);
    for r in 1..=top {
        for b in 1..=top {
            if b == r || gcd(b, r) != 1 {
                continue;
            }
            pairs += 1;
            let sign = (r - b).signum();
            let cuts = folds_of(b, r);
            let mut stops: Vec<(f64, f64)> = vec![(0.0, 1.0), (PI, 1.0)];
            for j in 1..r {
                stops.push((PI * j as f64 / r as f64, f64::INFINITY));
            }
            for k in 1..b {
                stops.push((PI * k as f64 / b as f64, 0.0));
            }
            let mut reaches = vec![1.0];
            for &d in &cuts {
                let v = lift(b, r, d);
                stops.push((d, v));
                reaches.push(v);
                let e = 1e-3 / r.max(b) as f64;
                let side = |x: f64| lift(b, r, x) - v;
                if sign * ((side(d - e) > 0.0) as i64 - (side(d - e) < 0.0) as i64) < 0
                    || sign * ((side(d + e) > 0.0) as i64 - (side(d + e) < 0.0) as i64) < 0
                {
                    mins += 1;
                }
            }
            stops.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            let pieces: Vec<(f64, f64)> = stops
                .windows(2)
                .map(|w| (w[0].1.min(w[1].1), w[0].1.max(w[1].1)))
                .collect();
            reaches.sort_by(|x, y| x.partial_cmp(y).unwrap());
            if cuts.len() as i64 + 1 != (r - b).abs() {
                wrong += 1;
                println!(
                    "    FOLDS b {b} rho {r} angles {} want {}",
                    cuts.len() + 1,
                    (r - b).abs()
                );
            }
            let mut probes = vec![(0.5 * reaches[0], 0_i64)];
            for k in 1..reaches.len() {
                if reaches[k] - reaches[k - 1] > 1e-7 * reaches[k] {
                    probes.push((0.5 * (reaches[k] + reaches[k - 1]), k as i64));
                }
            }
            probes.push((2.0 * reaches[reaches.len() - 1] + 1.0, reaches.len() as i64));
            for (t, below) in probes {
                levels += 1;
                if level_count(&pieces, t) != 2 * (b - 1) + 2 * sign * below {
                    wrong += 1;
                    println!(
                        "    FOLDS b {b} rho {r} level {t} count {} want {}",
                        level_count(&pieces, t),
                        2 * (b - 1) + 2 * sign * below
                    );
                }
            }
            if 2 * (b - 1) + 2 * sign * reaches.len() as i64 == 2 * (r - 1) {
                closed += 1;
            }
        }
    }
    println!("  frequency pairs {pairs} to {top} levels {levels} failures {wrong} ends closed {closed} folds on the wrong side {mins}");
}

fn census(top: i64) {
    println!("THE CENSUS  thresholds counted against min(|m|, a), the closed form against the crossing equation");
    let (mut cases, mut wrong, mut ends) = (0_usize, 0_usize, 0_usize);
    let mut fractions = 0_usize;
    for a in 2..=top {
        for b in 1..a {
            if gcd(a, b) != 1 {
                continue;
            }
            fractions += 1;
            for out in [false, true] {
                let m = skew(a, b, out);
                if m == 0 {
                    continue;
                }
                cases += 1;
                let steps = steps_of(a, b, out, 600_000);
                if steps.len() as i64 != turns(b, rim(a, b, out)) {
                    wrong += 1;
                    println!(
                        "    STEPS {} {a}/{b} m {m} found {} want {}",
                        tag(out),
                        steps.len(),
                        turns(b, rim(a, b, out))
                    );
                }
                if law(a, b, out, steps.last().unwrap() * 4.0 + 4.0, &steps)
                    != a * (rim(a, b, out) - 1)
                {
                    ends += 1;
                    println!("    LIMIT {} {a}/{b} m {m}", tag(out));
                }
                let mut marks = vec![steps[0] * 0.5];
                for k in 1..steps.len() {
                    if steps[k] - steps[k - 1] > 1e-6 {
                        marks.push((steps[k] * steps[k - 1]).sqrt());
                    }
                }
                marks.push(steps.last().unwrap() * 1.7 + 0.3);
                for lam in marks {
                    let root = crossings(a, b, out, lam, 60_000) as i64;
                    let want = law(a, b, out, lam, &steps);
                    if root != want {
                        wrong += 1;
                        println!(
                            "    LAW {} {a}/{b} reach {lam:.6} equation {root} law {want}",
                            tag(out)
                        );
                    }
                }
            }
        }
    }
    println!("  fractions {fractions} cases {cases} step-count or law failures {wrong} limit failures {ends}");
}

fn spot(a: i64, b: i64, out: bool, lam: f64, phi: f64) -> (f64, f64) {
    let (r, bf, sw) = (rim(a, b, out) as f64, b as f64, side(out));
    let (u, v) = (bf * phi, sw * r * phi);
    (
        r * u.cos() + lam * bf * v.cos(),
        r * u.sin() + lam * bf * v.sin(),
    )
}

fn mirror(top: i64) {
    println!("THE RECIPROCITY  swapping the wheel and the rim frequency inverts every threshold");
    let (mut pairs, mut worst) = (0_usize, 0.0_f64);
    for b in 1..top {
        for r in 1..top {
            if b == r || b + r > top || gcd(b, r) != 1 {
                continue;
            }
            pairs += 1;
            let here = steps_of(b + r, b, false, 400_000);
            let there = steps_of(b + r, r, false, 400_000);
            for (x, y) in here.iter().zip(there.iter().rev()) {
                worst = worst.max((x * y - 1.0).abs());
            }
        }
    }
    println!("  wheel and rim swaps {pairs} to {top}  worst deviation of the product from one {worst:.2e}");
}

fn contact() {
    println!("THE CONTACT  on a tangency reach the meetings are the transversal ones plus a/2 per tangency angle, each a tacnode");
    for (a, b, out) in [
        (5, 1, false),
        (7, 2, false),
        (11, 4, false),
        (3, 1, true),
        (5, 2, true),
    ] {
        let steps = steps_of(a, b, out, 400_000);
        let angles = tangent_angles(a, b, out, 400_000);
        let mut hits: Vec<f64> = Vec::new();
        for &t in steps.iter().filter(|&&t| t > 1.0 + 1e-9) {
            if !hits.iter().any(|h| (h - t).abs() < 1e-7) {
                hits.push(t);
            }
        }
        for hit in hits {
            let here: Vec<f64> = angles
                .iter()
                .copied()
                .filter(|&d| d > 1e-9 && (reach_at(a, b, out, d) - hit).abs() < 1e-7)
                .collect();
            let simple = arms(a, b, out, hit, 400_000);
            let under = law(a, b, out, hit * (1.0 - 1e-7), &steps);
            let (r, bf) = (rim(a, b, out) as f64, b as f64);
            let d0 = here[0];
            let eps = if (r * (bf * d0).sin() - hit * bf * (r * d0).sin()).abs() < 1e-6 {
                1.0
            } else {
                -1.0
            };
            let sigma = if eps * side(out) > 0.0 {
                PI / a as f64
            } else {
                0.0
            };
            let (p0, p1) = (
                spot(a, b, out, hit, sigma + d0),
                spot(a, b, out, hit, sigma - d0),
            );
            println!(
                "  {} {a}/{b} reach {hit:.9}  angles {}  simple roots {simple}  meetings {}  under {under}  over {}  contact gap {:.2e} at radius {:.6}",
                tag(out),
                here.len(),
                a as usize * (simple + here.len()) / 2,
                law(a, b, out, hit * (1.0 + 1e-7), &steps),
                (p0.0 - p1.0).hypot(p0.1 - p1.1),
                p0.0.hypot(p0.1)
            );
        }
    }
    println!("THE BIRTHS  a root is born at each end of (0, pi) as the reach passes 1, so the jump of a is two births of a/2");
    for (a, b, out) in [(3, 1, false), (5, 1, false), (7, 2, false), (3, 1, true)] {
        for lam in [0.98, 1.02] {
            let mut ends = Vec::new();
            for eps in [1.0_f64, -1.0] {
                let roots = arm_roots(a, b, out, lam, eps, 400_000);
                ends.push(format!(
                    "e{}:{}",
                    if eps > 0.0 { "+" } else { "-" },
                    roots
                        .iter()
                        .map(|d| format!("{d:.4}"))
                        .collect::<Vec<String>>()
                        .join(",")
                ));
            }
            println!(
                "  {} {a}/{b} reach {lam}  roots {}  count {}",
                tag(out),
                ends.join(" "),
                crossings(a, b, out, lam, 400_000)
            );
        }
    }
}

fn centre() {
    println!("THE CENTRE  at reach rho/b the curve runs through the origin with a branches and the point count drops by C(a,2) - 1");
    for (a, b, out) in [
        (3, 1, false),
        (5, 1, false),
        (5, 2, false),
        (5, 3, false),
        (5, 4, false),
        (7, 2, false),
        (7, 4, false),
        (3, 1, true),
        (3, 2, true),
        (5, 2, true),
        (4, 1, false),
    ] {
        let steps = steps_of(a, b, out, 400_000);
        let hub = rim(a, b, out) as f64 / b as f64;
        let mut row = Vec::new();
        for (name, lam) in [("under", hub * 0.97), ("hub", hub), ("over", hub * 1.03)] {
            let (curves, box_) = draw(a, b, out, &[(lam, 0.0)], 24001);
            let (own, _, _, spots, _) = nodes(&curves, box_, 260, 0, 0.0);
            row.push(format!(
                "{name} pairs {own} points {spots} law {}",
                law(a, b, out, lam, &steps)
            ));
        }
        println!(
            "  {} {a}/{b} rho/b {hub:.6}  {}  hub law less C(a,2)-1 {}{}",
            tag(out),
            row.join("  "),
            law(a, b, out, hub, &steps) - a * (a - 1) / 2 + 1,
            if a % 2 == 0 {
                "  even a, the hub is a threshold too"
            } else {
                ""
            }
        );
    }
}

fn hunt() {
    println!("THE SWEEP  reach 0.5 to 4, trace against the closed form, two sample counts");
    let (mut tests, mut wrong, mut worst_gap) = (0_usize, 0_usize, 0.0_f64);
    let mut blur = 0.0_f64;
    for (a, b, out) in cases() {
        let steps = steps_of(a, b, out, 400_000);
        let inband: Vec<f64> = steps
            .iter()
            .copied()
            .filter(|&t| (0.4987..=4.0014).contains(&t))
            .collect();
        let mut seen: Vec<(f64, f64, i64, i64)> = Vec::new();
        let mut prev: Option<(f64, i64)> = None;
        for k in 0..=202 {
            let lam = 0.4987 + 0.01734 * k as f64;
            let want = law(a, b, out, lam, &steps);
            let root = crossings(a, b, out, lam, 40_000) as i64;
            if root != want {
                wrong += 1;
                println!(
                    "    ROOTS {} {a}/{b} reach {lam:.4} equation {root} law {want}",
                    tag(out)
                );
            }
            let mut got = 0;
            for samples in [4001, 12001] {
                let (curves, box_) = draw(a, b, out, &[(lam, 0.0)], samples);
                got = nodes(&curves, box_, 220, 0, 0.0).0 as i64;
                tests += 1;
                if got != want {
                    wrong += 1;
                    let gap = inband
                        .iter()
                        .map(|t| (t - lam).abs())
                        .fold(f64::MAX, f64::min);
                    blur = blur.max(gap);
                    println!("    MISS {} {a}/{b} reach {lam:.4} samples {samples} trace {got} law {want} gap {gap:.4}", tag(out));
                }
            }
            if let Some((last, count)) = prev {
                if got != count {
                    seen.push((last, lam, count, got));
                    worst_gap = worst_gap.max(lam - last);
                }
            }
            prev = Some((lam, got));
        }
        let jumps: Vec<String> = seen
            .iter()
            .map(|(x, y, c0, c1)| {
                let inside = inband.iter().filter(|t| **t > *x && **t <= *y).count();
                format!("({x:.3},{y:.3}) {c0}->{c1} holds {inside}")
            })
            .collect();
        let caught: usize = seen
            .iter()
            .map(|(x, y, _, _)| inband.iter().filter(|t| **t > *x && **t <= *y).count())
            .sum();
        println!(
            "  {} {a}/{b}  predicted {:?}  jumps {}  caught {caught} of {}",
            tag(out),
            inband
                .iter()
                .map(|t| (t * 1e4).round() / 1e4)
                .collect::<Vec<f64>>(),
            jumps.join("  "),
            inband.len()
        );
    }
    println!("  tests {tests} disagreements {wrong} widest jump bracket {worst_gap:.4} widest blur {blur:.4}");
}

fn pairs() {
    println!("THE PAIRS  two pencils on one wheel, crossings between the two curves");
    for (a, b, out) in [(5, 1, false), (7, 2, false), (8, 3, false), (5, 2, true)] {
        for lam in [0.6, 1.3, 2.4, 3.7] {
            let p = (lam, 0.0);
            for (name, q) in [
                ("even", (lam * 0.6, lam * 0.8)),
                ("short", (lam * 0.45, lam * 0.35)),
                ("centre", (0.0, 0.0)),
            ] {
                let (curves, box_) = draw(a, b, out, &[p, q], 12001);
                let (_, between, _, _, spots) = nodes(&curves, box_, 220, 0, 0.0);
                println!(
                    "  {} {a}/{b} reach {lam} seat {name}  trace pairs {between} points {spots}  law {}  2ab {}",
                    tag(out),
                    pair_count(a, b, out, p, q, 400_000),
                    2 * a * b
                );
            }
        }
    }
    println!("THE PAIR POINTS  the pair law counts parameter pairs, and a point count needs no crossing of one curve to sit on the other");
    for (a, b, out, lam, name) in [
        (3, 1, false, 5.0_f64.sqrt() - 1.0, "witness"),
        (3, 1, false, (5.0_f64.sqrt() - 1.0) * 0.98, "under"),
        (3, 1, false, (5.0_f64.sqrt() - 1.0) * 1.02, "over"),
        (5, 1, false, 2.6, "generic"),
    ] {
        let (p, q) = ((lam, 0.0), (0.0, 0.0));
        let (curves, box_) = draw(a, b, out, &[p, q], 24001);
        let (_, between, _, _, spots) = nodes(&curves, box_, 260, 0, 0.0);
        println!(
            "  {} {a}/{b} reach {lam:.9} {name}  law pairs {}  trace pairs {between}  trace points {spots}  2ab {}",
            tag(out),
            pair_count(a, b, out, p, q, 400_000),
            2 * a * b
        );
    }
    println!(
        "THE PAIR THRESHOLD  two seats at one reach, half angle nu apart, the reach where 2ab dies"
    );
    for (a, b, out) in [(5, 1, false), (7, 2, false), (8, 3, false), (5, 2, true)] {
        let mut row = Vec::new();
        for nu in [1.0_f64, 0.5, 0.2, 0.05, 0.01, 0.001, 0.0001] {
            let alive = |lam: f64| {
                let p = (lam * nu.cos(), lam * nu.sin());
                let q = (lam * nu.cos(), -lam * nu.sin());
                pair_count(a, b, out, p, q, 20_000) as i64 == 2 * a * b
            };
            let (mut lo, mut hi) = (1.0, f64::NAN);
            for k in 0..=600 {
                let lam = 1.0 + 1e-6 * 1e7_f64.powf(k as f64 / 600.0);
                if alive(lam) {
                    lo = lam;
                } else {
                    hi = lam;
                    break;
                }
            }
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                if alive(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            row.push(format!("{nu}:{hi:.6}"));
        }
        println!(
            "  {} {a}/{b} 2ab {}  first departure by nu  {}",
            tag(out),
            2 * a * b,
            row.join("  ")
        );
    }
}

fn main() {
    if std::env::args().nth(1).as_deref() == Some("words") {
        words(220);
        folds(100);
        return;
    }
    axes();
    staircase();
    signs(220);
    words(220);
    folds(100);
    census(24);
    mirror(26);
    contact();
    centre();
    pairs();
    hunt();
}
