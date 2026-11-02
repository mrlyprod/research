use std::env;
use std::time::Instant;

// WIDE

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct W([u64; 3]);

impl W {
    const ZERO: W = W([0, 0, 0]);
    const ONE: W = W([1, 0, 0]);

    fn int(v: i64) -> W {
        let s = if v < 0 { u64::MAX } else { 0 };
        W([v as u64, s, s])
    }

    fn add(self, o: W) -> W {
        let (a0, c0) = self.0[0].overflowing_add(o.0[0]);
        let (a1, c1) = self.0[1].overflowing_add(o.0[1]);
        let (a1, c2) = a1.overflowing_add(c0 as u64);
        let a2 = self.0[2]
            .wrapping_add(o.0[2])
            .wrapping_add((c1 || c2) as u64);
        W([a0, a1, a2])
    }

    fn neg(self) -> W {
        W([!self.0[0], !self.0[1], !self.0[2]]).add(W::ONE)
    }

    fn sub(self, o: W) -> W {
        self.add(o.neg())
    }

    fn dbl(self) -> W {
        self.add(self)
    }

    fn negative(self) -> bool {
        self.0[2] >> 63 == 1
    }

    fn abs(self) -> W {
        if self.negative() {
            self.neg()
        } else {
            self
        }
    }

    fn is_zero(self) -> bool {
        self == W::ZERO
    }

    fn bits(self) -> u32 {
        for i in (0..3).rev() {
            if self.0[i] != 0 {
                return 64 * i as u32 + 64 - self.0[i].leading_zeros();
            }
        }
        0
    }

    fn bit(self, i: u32) -> bool {
        (self.0[(i / 64) as usize] >> (i % 64)) & 1 == 1
    }

    fn below(self, o: W) -> bool {
        for i in (0..3).rev() {
            if self.0[i] != o.0[i] {
                return self.0[i] < o.0[i];
            }
        }
        false
    }

    fn divmod(self, d: u64) -> (W, u64) {
        let mut r: u128 = 0;
        let mut q = [0u64; 3];
        for i in (0..3).rev() {
            let cur = (r << 64) | self.0[i] as u128;
            q[i] = (cur / d as u128) as u64;
            r = cur % d as u128;
        }
        (W(q), r as u64)
    }

    fn mul(self, k: u64) -> W {
        let mut out = [0u64; 3];
        let mut carry: u128 = 0;
        for (i, limb) in self.0.iter().enumerate() {
            let cur = *limb as u128 * k as u128 + carry;
            out[i] = cur as u64;
            carry = cur >> 64;
        }
        W(out)
    }

    fn shr(self, k: u32) -> W {
        let mut x = self;
        for _ in 0..k {
            x = W([
                (x.0[0] >> 1) | (x.0[1] << 63),
                (x.0[1] >> 1) | (x.0[2] << 63),
                x.0[2] >> 1,
            ]);
        }
        x
    }

    fn trailing(self) -> u32 {
        for i in 0..3 {
            if self.0[i] != 0 {
                return 64 * i as u32 + self.0[i].trailing_zeros();
            }
        }
        192
    }

    fn float(self) -> f64 {
        self.0[0] as f64 + self.0[1] as f64 * 2f64.powi(64) + self.0[2] as f64 * 2f64.powi(128)
    }

    fn dec(self) -> String {
        if self.is_zero() {
            return "0".into();
        }
        let mut parts = vec![];
        let mut x = self;
        while !x.is_zero() {
            let (q, r) = x.divmod(10_000_000_000_000_000_000);
            parts.push(r);
            x = q;
        }
        let mut s = parts.pop().unwrap().to_string();
        for p in parts.iter().rev() {
            s += &format!("{p:019}");
        }
        s
    }

    fn pow(b: u64, e: u32) -> W {
        let mut x = W::ONE;
        for _ in 0..e {
            x = x.mul(b);
        }
        x
    }
}

// LAYERS

fn tm(h: W) -> (W, u32) {
    let l = h.bits();
    let (mut x, mut y) = (W::int(3), W::int(-1));
    for i in (0..l).rev() {
        let s = x.add(y).neg();
        if h.bit(i) {
            x = s;
            y = y.dbl();
        } else {
            x = x.dbl();
            y = s;
        }
    }
    (x, l)
}

fn t3(h: W) -> Option<(bool, u32)> {
    let (mut h, mut m, mut neg) = (h, 0, false);
    loop {
        if h.is_zero() {
            return Some((neg, m));
        }
        if h == W::ONE {
            return None;
        }
        let (q, r) = h.divmod(3);
        h = if r == 1 { q.add(W::ONE) } else { q };
        if r != 0 {
            m += 1;
            neg = !neg;
        }
    }
}

fn runs(h: W) -> u32 {
    let l = h.bits();
    1 + (0..l.saturating_sub(1))
        .filter(|&i| h.bit(i) != h.bit(i + 1))
        .count() as u32
}

fn weight3(h: W) -> u32 {
    let (mut h, mut w) = (h, 0);
    while !h.is_zero() {
        let (q, r) = h.divmod(3);
        h = if r == 2 { q.add(W::ONE) } else { q };
        w += (r != 0) as u32;
    }
    w
}

struct Exact {
    num: W,
    two: u32,
    three: u32,
}

impl Exact {
    fn of(h: W) -> Exact {
        let (x, l) = tm(h);
        let Some((neg, m)) = t3(h) else {
            return Exact {
                num: W::ZERO,
                two: 0,
                three: 0,
            };
        };
        let mut num = if neg { x.neg() } else { x };
        let mut two = l;
        let mut three = m + 1;
        if num.is_zero() {
            return Exact {
                num,
                two: 0,
                three: 0,
            };
        }
        let z = num.abs().trailing().min(two);
        let sign = num.negative();
        let mut mag = num.abs().shr(z);
        two -= z;
        while three > 0 {
            let (q, r) = mag.divmod(3);
            if r != 0 {
                break;
            }
            mag = q;
            three -= 1;
        }
        num = if sign { mag.neg() } else { mag };
        Exact { num, two, three }
    }

    fn value(&self) -> f64 {
        let v = self.num.abs().float() / 2f64.powi(self.two as i32) / 3f64.powi(self.three as i32);
        if self.num.negative() {
            -v
        } else {
            v
        }
    }

    fn show(&self) -> String {
        let s = if self.num.negative() { "-" } else { "" };
        format!(
            "{s}{}/(2^{} 3^{})",
            self.num.abs().dec(),
            self.two,
            self.three
        )
    }
}

fn tm64(h: u64) -> f64 {
    let l = 64 - h.leading_zeros();
    let (mut x, mut y) = (3i64, -1i64);
    for i in (0..l).rev() {
        let s = -(x + y);
        if (h >> i) & 1 == 1 {
            x = s;
            y *= 2;
        } else {
            x *= 2;
            y = s;
        }
    }
    x as f64 / 3.0 / 2f64.powi(l as i32)
}

fn t364(h: u64) -> f64 {
    let (mut h, mut v) = (h, 1.0f64);
    loop {
        if h == 0 {
            return v;
        }
        if h == 1 {
            return 0.0;
        }
        let r = h % 3;
        h = if r == 1 { h / 3 + 1 } else { h / 3 };
        if r != 0 {
            v /= -3.0;
        }
    }
}

fn runs64(h: u64) -> u32 {
    let l = 64 - h.leading_zeros();
    if l <= 1 {
        return 1;
    }
    1 + ((h ^ (h >> 1)) & ((1u64 << (l - 1)) - 1)).count_ones()
}

fn weight64(h: u64) -> u32 {
    let (mut h, mut w) = (h, 0);
    while h != 0 {
        let r = h % 3;
        h = if r == 2 { h / 3 + 1 } else { h / 3 };
        w += (r != 0) as u32;
    }
    w
}

// RATIONALS

#[derive(Clone, Copy, Debug, PartialEq)]
struct Q(i128, i128);

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

fn q(n: i128, d: i128) -> Q {
    let g = gcd(n, d).max(1) * d.signum();
    Q(n / g, d / g)
}

impl Q {
    fn add(self, o: Q) -> Q {
        q(self.0 * o.1 + o.0 * self.1, self.1 * o.1)
    }
    fn mul(self, o: Q) -> Q {
        q(self.0 * o.0, self.1 * o.1)
    }
    fn neg(self) -> Q {
        Q(-self.0, self.1)
    }
    fn inv(self) -> Q {
        q(self.1, self.0)
    }
    fn show(self) -> String {
        if self.1 == 1 {
            format!("{}", self.0)
        } else {
            format!("{}/{}", self.0, self.1)
        }
    }
}

type Poly = Vec<Q>;

fn padd(a: &Poly, b: &Poly) -> Poly {
    (0..a.len().max(b.len()))
        .map(|i| {
            a.get(i)
                .copied()
                .unwrap_or(Q(0, 1))
                .add(b.get(i).copied().unwrap_or(Q(0, 1)))
        })
        .collect()
}

fn pmul(a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![Q(0, 1); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = out[i + j].add(x.mul(*y));
        }
    }
    out
}

fn pscale(a: &Poly, k: Q) -> Poly {
    a.iter().map(|x| x.mul(k)).collect()
}

fn trim(mut a: Poly) -> Poly {
    while a.len() > 1 && a.last().unwrap().0 == 0 {
        a.pop();
    }
    a
}

fn at(a: &Poly, t: Q) -> Q {
    a.iter().rev().fold(Q(0, 1), |acc, c| acc.mul(t).add(*c))
}

fn is_zero(a: &Poly) -> bool {
    a.iter().all(|c| c.0 == 0)
}

fn pder(a: &Poly) -> Poly {
    let d: Poly = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.mul(Q(i as i128, 1)))
        .collect();
    if d.is_empty() {
        vec![Q(0, 1)]
    } else {
        d
    }
}

fn prem(a: &Poly, b: &Poly) -> Poly {
    let b = trim(b.clone());
    let lead = b.last().unwrap().inv();
    let mut r = trim(a.clone());
    while !is_zero(&r) && r.len() >= b.len() {
        let k = r.last().unwrap().mul(lead);
        let s = r.len() - b.len();
        for (i, c) in b.iter().enumerate() {
            r[s + i] = r[s + i].add(c.mul(k).neg());
        }
        r.pop();
        if r.is_empty() {
            return vec![Q(0, 1)];
        }
        r = trim(r);
    }
    r
}

fn sturm(p: &Poly, lo: Q, hi: Q) -> usize {
    let mut seq = vec![trim(p.clone()), trim(pder(p))];
    while !is_zero(seq.last().unwrap()) {
        let n = seq.len();
        let r = prem(&seq[n - 2], &seq[n - 1]);
        seq.push(pscale(&r, Q(-1, 1)));
    }
    seq.pop();
    let changes = |x: Q| {
        let signs: Vec<i128> = seq
            .iter()
            .map(|s| at(s, x).0.signum())
            .filter(|&s| s != 0)
            .collect();
        signs.windows(2).filter(|w| w[0] != w[1]).count()
    };
    assert!(at(p, lo).0 != 0);
    changes(lo) - changes(hi)
}

fn pshow(a: &Poly) -> String {
    a.iter()
        .enumerate()
        .rev()
        .map(|(i, c)| format!("{} t^{i}", c.show()))
        .collect::<Vec<_>>()
        .join(" + ")
}

// LEMMA

const FORM: [i128; 3] = [20, 10, 41];
const RHO2: Q = Q(103, 200);
const KBOUND: f64 = 1.0093;
const RHO: f64 = 0.717636;

fn lemma() {
    let [al, be, ga] = FORM;
    let t: Poly = vec![Q(0, 1), Q(1, 1)];
    let c: Poly = vec![Q(-1, 3), Q(1, 3)];
    let t2 = pmul(&t, &t);
    let tc = pmul(&t, &c);
    let c2 = pmul(&c, &c);
    let a = padd(&vec![RHO2.mul(Q(al, 1))], &pscale(&t2, Q(-ga, 1)));
    let b = padd(
        &padd(&vec![RHO2.mul(Q(2 * be, 1))], &pscale(&t, Q(-2 * be, 1))),
        &pscale(&tc, Q(-2 * ga, 1)),
    );
    let d = padd(
        &padd(
            &vec![RHO2.mul(Q(ga, 1)).add(Q(-al, 1))],
            &pscale(&c, Q(-2 * be, 1)),
        ),
        &pscale(&c2, Q(-ga, 1)),
    );
    let delta = trim(padd(
        &pscale(&pmul(&a, &d), Q(4, 1)),
        &pscale(&pmul(&b, &b), Q(-1, 1)),
    ));
    let ends = [Q(-1, 2), Q(1, 4)];
    println!(
        "form Q(x, y) = {al} x^2 + {} xy + {ga} y^2, rho^2 = {}",
        2 * be,
        RHO2.show()
    );
    let mut ok = true;
    for (name, p) in [("a", trim(a)), ("d", trim(d)), ("Delta", delta)] {
        let lead = *p.last().unwrap();
        let concave = p.len() == 3 && lead.0 < 0;
        let vals: Vec<Q> = ends.iter().map(|&e| at(&p, e)).collect();
        let pos = vals.iter().all(|v| v.0 > 0);
        let roots = sturm(&p, ends[0], ends[1]);
        ok &= concave && pos && roots == 0;
        println!(
            "  {name}(t) = {}; concave {concave}; at -1/2 and 1/4: {} and {}; real roots on [-1/2, 1/4] by an exact Sturm sequence: {roots}",
            pshow(&p),
            vals[0].show(),
            vals[1].show()
        );
    }
    println!(
        "  certificate on [-1/2, 1/4]: {}",
        if ok { "holds" } else { "FAILS" }
    );
    assert!(ok);
    let mut worst = 0f64;
    for i in 0..=3000 {
        let s = -0.5 + 0.75 * i as f64 / 3000.0;
        worst = worst.max(opnorm(s));
    }
    println!(
        "  sup of the form norm of A(t) on a 3001-point grid: {worst:.6} against rho = {:.6}",
        (103f64 / 200.0).sqrt()
    );
    let k2: f64 = 41.0 * 161.0 / 720.0 / 9.0;
    println!(
        "  K = sqrt(41/720) Q(1, -1/3)^(1/2) = sqrt(41 . 161/720)/3 = {:.6}, odd r: sqrt(20 . 161/720)/3 = {:.6}",
        k2.sqrt(),
        (20.0f64 * 161.0 / 720.0 / 9.0).sqrt()
    );
    let lam = (3.0 + 73f64.sqrt()) / 32.0;
    println!(
        "  lower: A(1/4) A(-1/2) A(-1/2) has char poly x^2 - 3/16 x - 1/16, radius (3 + sqrt 73)/32 = {lam:.6}, per run {:.6}",
        lam.cbrt()
    );
}

fn opnorm(t: f64) -> f64 {
    let p = [[20.0, 10.0], [10.0, 41.0]];
    let a = [[0.0, 1.0], [t, (t - 1.0) / 3.0]];
    let mut best = 0f64;
    for k in 0..4000 {
        let th = std::f64::consts::PI * k as f64 / 4000.0;
        let u = [th.cos(), th.sin()];
        let v = [
            a[0][0] * u[0] + a[0][1] * u[1],
            a[1][0] * u[0] + a[1][1] * u[1],
        ];
        let qf = |w: [f64; 2]| {
            p[0][0] * w[0] * w[0] + 2.0 * p[0][1] * w[0] * w[1] + p[1][1] * w[1] * w[1]
        };
        best = best.max((qf(v) / qf(u)).sqrt());
    }
    best
}

// RUNS

fn decay() {
    const BITS: u32 = 25;
    let mut top = vec![0f64; 30];
    let mut arg = vec![0u64; 30];
    let mut over = 0u64;
    let mut worst = 0f64;
    let mut sums = vec![0f64; BITS as usize + 1];
    let mut sumz = vec![0f64; BITS as usize + 1];
    let mut sum3 = vec![0f64; BITS as usize + 1];
    let mut tern = 0u64;
    for h in 1u64..(1 << BITS) {
        let c = tm64(h).abs();
        let r = runs64(h) as usize;
        if c > top[r] {
            top[r] = c;
            arg[r] = h;
        }
        let ratio = c / (KBOUND * RHO.powi(r as i32));
        worst = worst.max(ratio);
        if ratio > 1.0 {
            over += 1;
        }
        let c3 = t364(h).abs();
        if c3 > 3f64.powf(-(weight64(h) as f64) / 2.0) * (1.0 + 1e-12) {
            tern += 1;
        }
        let l = (64 - h.leading_zeros()) as usize;
        sums[l] += c * c;
        sumz[l] += (c * c3) * (c * c3);
        sum3[l] += c3 * c3;
    }
    for l in 1..=BITS as usize {
        sums[l] += sums[l - 1];
        sumz[l] += sumz[l - 1];
        sum3[l] += sum3[l - 1];
    }
    println!("all h < 2^{BITS}: |C_t(h)| <= {KBOUND} . {RHO}^r(h) fails at {over} h, largest ratio {worst:.6}; |C3| <= 3^(-w_3/2) fails at {tern} h");
    println!("  r: max |C_t| over h < 2^{BITS} with r runs, at h, and max/rho^r");
    for r in 1..=BITS as usize {
        println!(
            "  {r:2} {:.6e} at {:>9} ratio {:.4}",
            top[r],
            arg[r],
            top[r] / RHO.powi(r as i32)
        );
    }
    println!("  Cesaro over 1 <= h < 2^L: L, mean C_t^2, its bound K^2 ((1 + rho^2)/2)^L, mean C3^2, mean C_z^2, the partial sum of C_z^2");
    for l in (5..=BITS as usize).step_by(2) {
        let n = 2f64.powi(l as i32);
        println!(
            "  {l:2} {:.6e} {:.6e} {:.6e} {:.6e} {:.9}",
            sums[l] / n,
            KBOUND * KBOUND * ((1.0 + 103.0 / 200.0) / 2f64).powi(l as i32),
            sum3[l] / n,
            sumz[l] / n,
            sumz[l]
        );
    }
    println!(
        "  slope of log2 of the C3^2 mean per bit, L = 9..25: {:.4}",
        (sum3[25] / 2f64.powi(25) / (sum3[9] / 2f64.powi(9))).log2() / 16.0
    );
    println!(
        "  slopes of log2 of the C_t^2 mean per bit, L = 9..25: {:.4}",
        (sums[25] / 2f64.powi(25) / (sums[9] / 2f64.powi(9))).log2() / 16.0
    );
    let theta = (400f64 / 303.0).log2();
    let ternary = (9f64 / 5.0).ln() / 3f64.ln();
    println!(
        "  proved exponents: theta = log2(400/303) = {theta:.6}, log3(9/5) = {ternary:.6}, 1 - theta - log3(9/5) = {:.6}",
        1.0 - theta - ternary
    );
    println!("family: runs (1, 1, 2) repeated from the top, |C_t| exact, per-run root and three-run ratio");
    let mut prev: Vec<f64> = vec![];
    for j in 1..=30u32 {
        let pattern = [1u32, 1, 2];
        let mut h = W::ZERO;
        let mut digit = 1u64;
        for i in 0..j {
            for _ in 0..pattern[(i % 3) as usize] {
                h = h.dbl().add(W::int(digit as i64));
            }
            digit ^= 1;
        }
        let (x, l) = tm(h);
        let v = x.abs().float() / 3.0 / 2f64.powi(l as i32);
        assert_eq!(runs(h), j);
        prev.push(v);
        if j % 3 == 0 {
            let three = if j > 3 {
                (v / prev[j as usize - 4]).cbrt()
            } else {
                f64::NAN
            };
            println!(
                "  r {j:2} h {:>20} |C_t| {v:.6e} root {:.6} ratio {three:.6}",
                h.dec(),
                v.powf(1.0 / j as f64)
            );
        }
    }
    println!("obstruction: C_t(2^k - 1) = (-1/2)^k + (1 - (-1/2)^k)/9, one run and k ones");
    for k in [1u32, 2, 3, 4, 6, 8, 12, 16, 24, 40, 60] {
        let h = W::pow(2, k).sub(W::ONE);
        let (x, l) = tm(h);
        let v = (if x.negative() {
            -x.abs().float()
        } else {
            x.float()
        }) / 3.0
            / 2f64.powi(l as i32);
        let s = (-0.5f64).powi(k as i32);
        let f = s + (1.0 - s) / 9.0;
        assert!((v - f).abs() < 1e-15);
        print!(" k {k}: {v:.9}");
    }
    println!();
}

// WINDOWS

fn windows() {
    const TOP: u32 = 24;
    const EPS: f64 = 1e-3;
    let mut tail = vec![(0f64, 0u64); TOP as usize + 2];
    let (mut many, mut last, mut near) = (0u64, 0u64, 0u64);
    println!("k, max |C_z| on [2^k, 2^(k+1)) at h, r_2(h), w_3(h), exact");
    for k in 0..=TOP {
        let (mut best, mut arg) = (0f64, 0u64);
        for h in (1u64 << k)..(1u64 << (k + 1)) {
            let c3 = t364(h);
            if c3 == 0.0 || (c3.abs() <= best && c3.abs() < 0.98 * EPS) {
                continue;
            }
            let v = (tm64(h) * c3).abs();
            if v >= EPS {
                many += 1;
                last = h;
            }
            if (v - EPS).abs() <= 0.02 * EPS {
                near += 1;
            }
            if v > best {
                best = v;
                arg = h;
            }
        }
        tail[k as usize] = (best, arg);
        if arg == 0 {
            println!("  {k:2} 0, C_z vanishes on the window");
            continue;
        }
        let e = Exact::of(W([arg, 0, 0]));
        assert!((e.value().abs() - best).abs() <= 1e-12 * best);
        println!(
            "  {k:2} {best:.6e} at {arg:>9} r {} w {} {}",
            runs64(arg),
            weight64(arg),
            e.show()
        );
    }
    let mut sup = 0f64;
    for k in (0..=TOP as usize).rev() {
        sup = sup.max(tail[k].0);
        if k % 4 == 0 {
            print!(" sup on [2^{k}, 2^25): {sup:.4e};");
        }
    }
    println!();
    let beyond = tail[22..=TOP as usize]
        .iter()
        .map(|t| t.0)
        .fold(0f64, f64::max);
    let e = Exact::of(W([last, 0, 0]));
    assert!(e.value().abs() >= EPS);
    println!(
        "  1 <= h < 2^25 with |C_z| >= {EPS:e}: {many} lags, the last at {last} r {} w {} C_z {:.4e} = {}; {near} lags within 0.02 of the threshold; sup on [2^22, 2^25): {beyond:.4e}",
        runs64(last),
        weight64(last),
        e.value(),
        e.show()
    );
}

// CHANGES

fn changes_of(h: u64, b: u64) -> u32 {
    let (mut h, mut prev, mut n) = (h, None, 0);
    while h != 0 {
        let d = h % b;
        if prev.is_some_and(|p| p != d) {
            n += 1;
        }
        prev = Some(d);
        h /= b;
    }
    n
}

fn signed_lengths(n: usize, b: u64) -> Vec<u32> {
    let mut f = vec![0u32; n + 1];
    for v in 1..=n {
        let (q, r) = ((v as u64 / b) as usize, (v as u64 % b) as u32);
        f[v] = if r == 0 {
            f[q]
        } else if v == 1 {
            1
        } else {
            (r + f[q]).min(b as u32 - r + f[q + 1])
        };
    }
    f
}

fn searched_lengths(m: i64, b: i64) -> Vec<u8> {
    let mut steps = Vec::new();
    let mut p = 1i64;
    while p <= m {
        steps.push(p);
        steps.push(-p);
        p *= b;
    }
    let mut dist = vec![u8::MAX; (2 * m + 1) as usize];
    let mut queue = vec![0i64];
    dist[m as usize] = 0;
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head];
        head += 1;
        let dv = dist[(v + m) as usize];
        for &s in &steps {
            let w = v + s;
            if w.abs() <= m && dist[(w + m) as usize] == u8::MAX {
                dist[(w + m) as usize] = dv + 1;
                queue.push(w);
            }
        }
    }
    dist
}

fn changes() {
    const N: u64 = 531441;
    let m = 3 * N as i64;
    println!("BC(h, b) <= 2 l_b(h) - 1 at every 1 <= h < 3^12 = {N}");
    for b in [2u64, 3] {
        let f = signed_lengths(N as usize, b);
        let g = searched_lengths(m, b as i64);
        let (mut fail, mut tight, mut worst, mut below, mut above) =
            (0u64, 0u64, i64::MIN, 0u64, 0u64);
        let mut first = Vec::new();
        for h in 1..N {
            let l = f[h as usize];
            let s = g[(h as i64 + m) as usize] as u32;
            below += (s < l) as u64;
            above += (s > l) as u64;
            let bc = changes_of(h, b);
            let gap = bc as i64 - (2 * l as i64 - 1);
            worst = worst.max(gap);
            fail += (gap > 0) as u64;
            if gap == 0 {
                tight += 1;
                if first.len() < 6 {
                    first.push(h);
                }
            }
            if b == 2 {
                assert_eq!(bc + 1, runs64(h));
                assert!(l <= runs64(h) + 1);
            } else {
                assert_eq!(l, weight64(h));
            }
        }
        println!("  b = {b}: failures {fail}, max of BC - (2 l - 1) {worst}, equality at {tight} lags, the first {first:?}");
        println!("    l_b by the digit recursion against breadth-first search on [-{m}, {m}]: {below} below, {above} above");
    }
    println!("  b = 2: BC(h, 2) = r_2(h) - 1 and l_2(h) <= r_2(h) + 1 at every lag; b = 3: l_3(h) = w_3(h) at every lag");
    let rho = (103.0f64 / 200.0).sqrt();
    let k = (41.0f64 * 161.0 / 720.0).sqrt() / 3.0;
    println!(
        "  rho^4 = 10609/40000 and 3 . 10609 = {} < 40000, so rho < 3^(-1/4): rho = {rho:.9}, 3^(-1/4) = {:.9}",
        3 * 10609,
        3f64.powf(-0.25)
    );
    println!(
        "  K = {k:.9}, K 3^(1/4) = {:.9}, log(3)/4 = {:.9}, log(3)/8 = {:.9}",
        k * 3f64.powf(0.25),
        3f64.ln() / 4.0,
        3f64.ln() / 8.0
    );
}

// SPARSE

struct Census {
    rmax: u32,
    cut: W,
    p64: Vec<u64>,
    p: Vec<W>,
    found: Vec<W>,
    survivors: u64,
    nodes: u64,
}

impl Census {
    fn pass(&self, x: u64) -> bool {
        ((x ^ (x >> 1)) & (u64::MAX >> 1)).count_ones() <= self.rmax
    }

    fn exact(&mut self, terms: &[(u32, bool)]) {
        self.survivors += 1;
        let mut h = W::ZERO;
        for &(b, neg) in terms {
            h = if neg {
                h.sub(self.p[b as usize])
            } else {
                h.add(self.p[b as usize])
            };
        }
        if h.negative() || h.is_zero() || !h.below(self.cut) || runs(h) > self.rmax {
            return;
        }
        assert_eq!(weight3(h), terms.len() as u32);
        self.found.push(h);
    }

    fn grow(&mut self, terms: &mut Vec<(u32, bool)>, low: u64, wmax: usize) {
        let last = terms.last().unwrap().0;
        if terms.len() + 1 == wmax {
            for b in 0..last {
                let pb = self.p64[b as usize];
                self.nodes += 2;
                for (neg, x) in [(false, low.wrapping_add(pb)), (true, low.wrapping_sub(pb))] {
                    if self.pass(x) {
                        terms.push((b, neg));
                        self.exact(terms);
                        terms.pop();
                    }
                }
            }
            return;
        }
        for b in 0..last {
            let pb = self.p64[b as usize];
            for (neg, x) in [(false, low.wrapping_add(pb)), (true, low.wrapping_sub(pb))] {
                self.nodes += 1;
                terms.push((b, neg));
                if self.pass(x) {
                    self.exact(terms);
                }
                self.grow(terms, x, wmax);
                terms.pop();
            }
        }
    }
}

fn sparse(rmax: u32, wmax: usize, digits: u32) {
    let clock = Instant::now();
    let cut = W::pow(10, digits);
    let mut top = 0u32;
    while W::pow(3, top + 1).below(cut.dbl()) {
        top += 1;
    }
    let p: Vec<W> = (0..=top).map(|b| W::pow(3, b)).collect();
    let p64 = p.iter().map(|x| x.0[0]).collect();
    let mut c = Census {
        rmax,
        cut,
        p64,
        p,
        found: vec![],
        survivors: 0,
        nodes: 0,
    };
    for b in 0..=top {
        let mut terms = vec![(b, false)];
        let low = c.p64[b as usize];
        c.nodes += 1;
        if c.pass(low) {
            c.exact(&terms);
        }
        if wmax > 1 {
            c.grow(&mut terms, low, wmax);
        }
    }
    let mut found = c.found.clone();
    found.sort_by(|a, b| {
        if a.below(*b) {
            std::cmp::Ordering::Less
        } else if b.below(*a) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    let before = found.len();
    found.dedup();
    assert_eq!(before, found.len());
    println!(
        "all h < 10^{digits} with r_2(h) <= {rmax} and w_3(h) <= {wmax}: ternary top exponent {top}, {} balanced words, {} passed the low-64-bit run filter, {} found, {:.1} s",
        c.nodes,
        c.survivors,
        found.len(),
        clock.elapsed().as_secs_f64()
    );
    let mut line = String::new();
    let shown = if found.len() <= 20000 {
        0
    } else {
        found.len() - 48
    };
    if shown > 0 {
        println!("  the list has more than 20000 elements; its last 48:");
    }
    for (i, h) in found.iter().enumerate().skip(shown) {
        line += &format!(" {}", h.dec());
        if i % 12 == 11 {
            println!(" {line}");
            line.clear();
        }
    }
    if !line.is_empty() {
        println!(" {line}");
    }
    let mut blocks = vec![0usize; 1 + (digits as usize * 10) / 24];
    for h in &found {
        blocks[(h.bits() as usize - 1) / 8] += 1;
    }
    let shape: Vec<String> = blocks.iter().map(|n| n.to_string()).collect();
    println!(
        "  elements per bit-length block [2^(8j), 2^(8j+8)), j = 0, 1, ...: {}",
        shape.join(" ")
    );
    let largest = *found.last().unwrap();
    println!(
        "  largest {} = 2-runs {}, w_3 {}",
        largest.dec(),
        runs(largest),
        weight3(largest)
    );
    let mut rows: Vec<(f64, W)> = found
        .iter()
        .map(|&h| (Exact::of(h).value().abs(), h))
        .collect();
    rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let zero = rows.iter().filter(|r| r.0 == 0.0).count();
    println!("  C_z = 0 at {zero} of them; the 8 largest elements and the max |C_z| over elements above each threshold:");
    for h in found.iter().rev().take(8) {
        let e = Exact::of(*h);
        println!(
            "    {} r {} w {} C_z {:.6e} = {}",
            h.dec(),
            runs(*h),
            weight3(*h),
            e.value(),
            e.show()
        );
    }
    for k in [8u32, 12, 16, 20, 24, 28, 32] {
        let floor = W::pow(2, k);
        if let Some(r) = rows.iter().find(|r| !r.1.below(floor)) {
            println!("    h >= 2^{k}: max |C_z| {:.6e} at {}", r.0, r.1.dec());
        } else {
            println!("    h >= 2^{k}: none");
        }
    }
    let near = W::pow(10, digits - 2).mul(98);
    println!(
        "  within 0.02 of the cut, h >= 0.98 . 10^{digits}: {}",
        found.iter().filter(|h| !h.below(near)).count()
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let clock = Instant::now();
    match args.get(1).map(String::as_str) {
        Some("lemma") => lemma(),
        Some("decay") => decay(),
        Some("windows") => windows(),
        Some("changes") => changes(),
        Some("sparse") => {
            let num = |i: usize, d: u32| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
            sparse(num(2, 6), num(3, 6) as usize, num(4, 40));
        }
        _ => {
            eprintln!("verbs: lemma, decay, windows, changes, sparse [rmax wmax digits]");
            std::process::exit(2);
        }
    }
    eprintln!("{:.2} s", clock.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layers_match_the_block_recursions() {
        for h in 1u64..3000 {
            let w = W([h, 0, 0]);
            let (x, l) = tm(w);
            let exact = (if x.negative() {
                -x.abs().float()
            } else {
                x.float()
            }) / 3.0
                / 2f64.powi(l as i32);
            let rec = if h % 2 == 0 {
                tm64(h / 2)
            } else {
                -(tm64(h / 2) + tm64(h / 2 + 1)) / 2.0
            };
            assert!((exact - rec).abs() < 1e-12 && (exact - tm64(h)).abs() < 1e-12);
            let c3 = match t3(w) {
                Some((neg, m)) => (if neg { -1.0 } else { 1.0 }) / 3f64.powi(m as i32),
                None => 0.0,
            };
            let rec3 = match h % 3 {
                0 => t364(h / 3),
                1 => -t364(h / 3 + 1) / 3.0,
                _ => -t364(h / 3) / 3.0,
            };
            assert!((c3 - rec3).abs() < 1e-12 && (c3 - t364(h)).abs() < 1e-12);
            assert_eq!(runs(w), runs64(h));
            assert_eq!(weight3(w), weight64(h));
        }
        assert_eq!(tm64(1), -1.0 / 3.0);
        assert_eq!(t364(2), -1.0 / 3.0);
    }

    #[test]
    fn wide_arithmetic_round_trips() {
        let big = W::pow(10, 40);
        assert_eq!(big.dec(), format!("1{}", "0".repeat(40)));
        assert_eq!(big.sub(W::ONE).add(W::ONE), big);
        assert_eq!(W::pow(3, 84).divmod(3).0, W::pow(3, 83));
        assert_eq!(weight3(W::pow(3, 84).sub(W::pow(3, 10))), 2);
        assert_eq!(runs(W::pow(2, 130).sub(W::ONE)), 1);
    }
}
