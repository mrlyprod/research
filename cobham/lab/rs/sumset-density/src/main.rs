use std::collections::HashSet;
use std::env;
use std::time::Instant;

// BITSET

struct Bits {
    top: u64,
    w: Vec<u64>,
}

impl Bits {
    fn new(top: u64) -> Bits {
        let words = (top / 64 + 1) as usize;
        Bits {
            top,
            w: vec![0; words],
        }
    }

    fn set(&mut self, x: u64) {
        self.w[(x / 64) as usize] |= 1u64 << (x % 64);
    }

    fn get(&self, x: u64) -> bool {
        (self.w[(x / 64) as usize] >> (x % 64)) & 1 == 1
    }

    fn mask(&mut self) {
        let r = self.top % 64;
        let last = self.w.len() - 1;
        if r != 63 {
            self.w[last] &= (1u64 << (r + 1)) - 1;
        }
    }

    fn or_shift(&mut self, sh: u64) {
        if sh > self.top {
            return;
        }
        let wsh = (sh / 64) as usize;
        let bsh = (sh % 64) as u32;
        let len = self.w.len();
        for i in (wsh..len).rev() {
            let src = i - wsh;
            let mut v = self.w[src] << bsh;
            if bsh > 0 && src > 0 {
                v |= self.w[src - 1] >> (64 - bsh);
            }
            self.w[i] |= v;
        }
        self.mask();
    }

    fn count(&self) -> u64 {
        self.w.iter().map(|x| x.count_ones() as u64).sum()
    }

    fn count_to(&self, d: u64) -> u64 {
        let full = (d / 64) as usize;
        let head: u64 = self.w[..full].iter().map(|x| x.count_ones() as u64).sum();
        let r = d % 64;
        let tail = if r == 63 {
            self.w[full]
        } else {
            self.w[full] & ((1u64 << (r + 1)) - 1)
        };
        head + tail.count_ones() as u64
    }
}

// DESIGNS

fn powers(base: u64, top: u64) -> Vec<u64> {
    let mut v = vec![1u64];
    while v
        .last()
        .unwrap()
        .checked_mul(base)
        .is_some_and(|p| p <= top)
    {
        v.push(v.last().unwrap() * base);
    }
    v
}

fn members(base: u64, top: u64) -> Vec<u64> {
    let mut list = vec![0u64];
    for p in powers(base, top) {
        let n = list.len();
        for i in 0..n {
            let x = list[i] + p;
            if x <= top {
                list.push(x);
            }
        }
    }
    list.sort_unstable();
    list
}

fn sumset(top: u64, direct: u64, shifted: u64) -> Bits {
    let mut s = Bits::new(top);
    for a in members(direct, top) {
        s.set(a);
    }
    for p in powers(shifted, top) {
        s.or_shift(p);
    }
    s
}

fn sumset_direct(top: u64) -> Bits {
    let mut s = Bits::new(top);
    let a = members(3, top);
    let b = members(4, top);
    for &x in &a {
        for &y in &b {
            if x + y <= top {
                s.set(x + y);
            }
        }
    }
    s
}

// SCAN

#[derive(Clone, Copy)]
struct Ext {
    c: u64,
    x: u64,
}

impl Ext {
    fn none() -> Ext {
        Ext { c: 0, x: 0 }
    }

    fn above(&self, c: u64, x: u64) -> bool {
        self.x == 0 || (c as u128) * (self.x as u128) > (self.c as u128) * (x as u128)
    }

    fn below(&self, c: u64, x: u64) -> bool {
        self.x == 0 || (c as u128) * (self.x as u128) < (self.c as u128) * (x as u128)
    }
}

#[derive(Clone, Copy)]
struct Window {
    lo: u64,
    hi: u64,
    max: Ext,
    min: Ext,
}

impl Window {
    fn new(lo: u64, hi: u64) -> Window {
        Window {
            lo,
            hi,
            max: Ext::none(),
            min: Ext::none(),
        }
    }

    fn feed(&mut self, c: u64, x: u64, up: bool, down: bool) {
        if up && self.max.above(c, x) {
            self.max = Ext { c, x };
        }
        if down && self.min.below(c, x) {
            self.min = Ext { c, x };
        }
    }
}

struct Scan {
    marks: Vec<(u64, u64)>,
    windows3: Vec<Window>,
    windows2: Vec<Window>,
    global: Window,
}

fn window_list(base: u64, top: u64) -> Vec<Window> {
    let p = powers(base, top);
    let mut v = Vec::new();
    for i in 0..p.len() {
        let lo = p[i];
        let hi = if i + 1 < p.len() { p[i + 1] - 1 } else { top };
        v.push(Window::new(lo, hi));
    }
    v
}

fn scan(s: &Bits, marks: &[u64]) -> Scan {
    let top = s.top;
    let windows3 = window_list(3, top);
    let windows2 = window_list(2, top);
    let mut special: HashSet<u64> = HashSet::new();
    for &m in marks {
        special.insert(m / 64);
    }
    for w in windows3.iter().chain(windows2.iter()) {
        special.insert(w.lo / 64);
        special.insert(w.hi / 64);
        if w.lo > 0 {
            special.insert((w.lo - 1) / 64);
        }
    }
    special.insert(0);
    special.insert(top / 64);
    let mut markset: Vec<u64> = marks.to_vec();
    markset.sort_unstable();
    markset.dedup();
    let mut sc = Scan {
        marks: Vec::new(),
        windows3,
        windows2,
        global: Window::new(1, top),
    };
    let mut i3 = 0usize;
    let mut i2 = 0usize;
    let mut mi = 0usize;
    let mut c: u64 = 0;
    for (wi, &word) in s.w.iter().enumerate() {
        let base = (wi as u64) * 64;
        if !special.contains(&(wi as u64)) {
            if word == 0 {
                let x = base + 63;
                sc.windows3[i3].feed(c, x, false, true);
                sc.windows2[i2].feed(c, x, false, true);
                sc.global.feed(c, x, false, true);
            } else if word == u64::MAX {
                c += 64;
                let x = base + 63;
                sc.windows3[i3].feed(c, x, true, false);
                sc.windows2[i2].feed(c, x, true, false);
                sc.global.feed(c, x, true, false);
            } else {
                for j in 0..64u64 {
                    let x = base + j;
                    let on = (word >> j) & 1 == 1;
                    if on {
                        c += 1;
                    }
                    sc.windows3[i3].feed(c, x, on, !on);
                    sc.windows2[i2].feed(c, x, on, !on);
                    sc.global.feed(c, x, on, !on);
                }
            }
            continue;
        }
        for j in 0..64u64 {
            let x = base + j;
            if x == 0 {
                continue;
            }
            if x > top {
                break;
            }
            if (word >> j) & 1 == 1 {
                c += 1;
            }
            while i3 + 1 < sc.windows3.len() && x >= sc.windows3[i3 + 1].lo {
                i3 += 1;
            }
            while i2 + 1 < sc.windows2.len() && x >= sc.windows2[i2 + 1].lo {
                i2 += 1;
            }
            sc.windows3[i3].feed(c, x, true, true);
            sc.windows2[i2].feed(c, x, true, true);
            sc.global.feed(c, x, true, true);
            while mi < markset.len() && markset[mi] == x {
                sc.marks.push((x, c));
                mi += 1;
            }
        }
    }
    sc
}

// PRINT

fn dens(c: u64, x: u64) -> String {
    let q = (c as u128) * 1_000_000u128 / (x as u128);
    format!("{}.{:06}", q / 1_000_000, q % 1_000_000)
}

fn expo(c: u64, x: u64) -> String {
    if x < 2 || c == 0 {
        return "-".to_string();
    }
    let e = (c as f64).ln() / (x as f64).ln();
    let t = (e * 1_000_000.0).floor() / 1_000_000.0;
    format!("{:.6}", t)
}

fn print_window(tag: &str, i: usize, w: &Window) {
    println!(
        "{} {:>2} [{}, {}] max {} at x = {} c = {} min {} at x = {} c = {}",
        tag,
        i,
        w.lo,
        w.hi,
        dens(w.max.c, w.max.x),
        w.max.x,
        w.max.c,
        dens(w.min.c, w.min.x),
        w.min.x,
        w.min.c
    );
}

fn centres(top: u64) -> Vec<(u32, u32, u64, bool)> {
    let p3 = powers(3, top * 3);
    let p4 = powers(4, top * 4);
    let mut v = Vec::new();
    for (s, &four) in p4.iter().enumerate().skip(1) {
        for (r, &three) in p3.iter().enumerate().skip(1) {
            let d = (three - 1) / 2 + (four - 1) / 3;
            if d > top {
                continue;
            }
            let clean = three > d && four > d;
            if three < 3 * four && four < 3 * three {
                v.push((r as u32, s as u32, d, clean));
            }
        }
    }
    v.sort_by_key(|t| t.2);
    v
}

fn density(k: u32) {
    let t0 = Instant::now();
    let top = 3u64.pow(k);
    let s = sumset(top, 3, 4);
    let built = t0.elapsed().as_secs_f64();
    let p3 = powers(3, top);
    let p4 = powers(4, top);
    let mut marks: Vec<u64> = Vec::new();
    marks.extend(p3.iter().copied());
    marks.extend(p4.iter().copied());
    marks.extend(p3.iter().map(|p| p / 2));
    marks.extend(p4.iter().map(|p| p / 3));
    let cs = centres(top);
    marks.extend(cs.iter().map(|t| t.2));
    marks.retain(|&x| x >= 1);
    let sc = scan(&s, &marks);
    let at = |x: u64| -> u64 { sc.marks.iter().find(|m| m.0 == x).map(|m| m.1).unwrap() };
    println!("top 3^{} = {}", k, top);
    println!("members of A up to top {}", members(3, top).len());
    println!("members of B up to top {}", members(4, top).len());
    println!("card(S meet [1, top]) {}", s.count() - 1);
    for (i, &x) in p3.iter().enumerate() {
        let c = at(x);
        println!(
            "3^{:<2} x = {} c = {} D = {} exponent {}",
            i,
            x,
            c,
            dens(c, x),
            expo(c, x)
        );
    }
    for (i, &x) in p4.iter().enumerate() {
        let c = at(x);
        println!(
            "4^{:<2} x = {} c = {} D = {} exponent {}",
            i,
            x,
            c,
            dens(c, x),
            expo(c, x)
        );
    }
    for (i, &p) in p3.iter().enumerate().skip(1) {
        let x = p / 2;
        let c = at(x);
        println!("3^{:<2}/2 x = {} c = {} D = {}", i, x, c, dens(c, x));
    }
    for (i, &p) in p4.iter().enumerate().skip(1) {
        let x = p / 3;
        let c = at(x);
        println!("4^{:<2}/3 x = {} c = {} D = {}", i, x, c, dens(c, x));
    }
    for &(r, sx, d, clean) in &cs {
        let c = at(d);
        println!(
            "centre r = {:<2} s = {:<2} 4^s/3^r = {} d = {} c = {} D = {} {}",
            r,
            sx,
            ratio(r, sx),
            d,
            c,
            dens(c, d),
            if clean { "clean" } else { "mixed" }
        );
    }
    for (i, w) in sc.windows3.iter().enumerate() {
        print_window("window3", i, w);
    }
    for (i, w) in sc.windows2.iter().enumerate() {
        print_window("window2", i, w);
    }
    print_window("global", 0, &sc.global);
    println!(
        "built in {:.2} s, scanned in {:.2} s",
        built,
        t0.elapsed().as_secs_f64() - built
    );
}

// ENERGY

fn zeros3(mut t: i64, k: u32) -> Option<u32> {
    let mut z = 0;
    for _ in 0..k {
        let r = t.rem_euclid(3);
        t = (t - r) / 3;
        if r == 2 {
            t += 1;
        } else if r == 0 {
            z += 1;
        }
    }
    if t == 0 {
        Some(z)
    } else {
        None
    }
}

fn energy(k: u32, m: u32) -> u128 {
    let p4: Vec<i64> = (0..m).map(|i| 4i64.pow(i)).collect();
    let mut e: u128 = 0;
    let mut digits = vec![-1i64; m as usize];
    loop {
        let mut t = 0i64;
        let mut z4 = 0u32;
        for i in 0..m as usize {
            t += digits[i] * p4[i];
            if digits[i] == 0 {
                z4 += 1;
            }
        }
        if let Some(z3) = zeros3(t, k) {
            e += 1u128 << (z3 + z4);
        }
        let mut i = 0usize;
        loop {
            if i == m as usize {
                return e;
            }
            if digits[i] < 1 {
                digits[i] += 1;
                break;
            }
            digits[i] = -1;
            i += 1;
        }
    }
}

#[cfg(test)]
fn energy_direct(k: u32, m: u32) -> u128 {
    let a = members(3, 3u64.pow(k) - 1);
    let b = members(4, 4u64.pow(m) - 1);
    let d = (3u64.pow(k) - 1) / 2 + (4u64.pow(m) - 1) / 3;
    let mut r = vec![0u64; d as usize + 1];
    for &x in &a {
        for &y in &b {
            r[(x + y) as usize] += 1;
        }
    }
    r.iter().map(|&v| (v as u128) * (v as u128)).sum()
}

fn six(num: u128, den: u128) -> String {
    let q = num * 1_000_000 / den;
    format!("{}.{:06}", q / 1_000_000, q % 1_000_000)
}

fn six_up(num: u128, den: u128) -> String {
    let q = (num * 1_000_000).div_ceil(den);
    format!("{}.{:06}", q / 1_000_000, q % 1_000_000)
}

fn ratio(r: u32, m: u32) -> String {
    six_up(4u128.pow(m) * 1_000_000, 3u128.pow(r) * 1_000_000)
}

fn energies(kmax: u32) {
    let t0 = Instant::now();
    let top = 3u64.pow(kmax);
    let s = sumset(top, 3, 4);
    let mut first: Vec<Option<(u32, f64)>> = vec![None, None];
    let mut last: Option<(u32, f64)> = None;
    for (r, m, d, clean) in centres(top) {
        if r < 4 {
            continue;
        }
        let e = energy(r, m);
        let pairs = 1u128 << (r + m);
        let card = s.count_to(d) as u128;
        let range = d as u128 + 1;
        println!(
            "k = {:<2} m = {:<2} {} 4^m/3^k = {} d = {} card = {} fill {} energy {} random {} ratio {} bound {}",
            r,
            m,
            if clean { "clean" } else { "mixed" },
            ratio(r, m),
            d,
            card,
            six(card, range),
            e,
            six(pairs * pairs, range),
            six_up(e * range, pairs * pairs),
            six(pairs * pairs, e * range)
        );
        let q = (e as f64) * (range as f64) / (pairs as f64) / (pairs as f64);
        for (i, start) in [6u32, 11].iter().enumerate() {
            if r >= *start && first[i].is_none() {
                first[i] = Some((r, q));
            }
        }
        if r >= 6 {
            last = Some((r, q));
        }
    }
    for f in first {
        if let (Some((k0, q0)), Some((k1, q1))) = (f, last) {
            let eta = (q1 / q0).ln() / 3f64.ln() / ((k1 - k0) as f64);
            println!(
                "growth of Q from k = {} to k = {} as 3^(eta k): eta = {:.6}",
                k0,
                k1,
                (eta * 1_000_000.0).ceil() / 1_000_000.0
            );
        }
    }
    println!("{:.2} s", t0.elapsed().as_secs_f64());
}

// LADDER

const CHUNK: u32 = 8;
const WIDTH: u64 = 6561;

fn chunk_tables() -> (Vec<u8>, Vec<u64>) {
    let mut zt = vec![0u8; WIDTH as usize];
    let mut ct = vec![0u64; WIDTH as usize];
    for r in 0..WIDTH {
        let mut x = r;
        let mut z = 0u8;
        for _ in 0..CHUNK {
            let d = x % 3;
            x /= 3;
            if d == 2 {
                x += 1;
            } else if d == 0 {
                z += 1;
            }
        }
        zt[r as usize] = z;
        ct[r as usize] = x;
    }
    (zt, ct)
}

fn low_strings(h: u32) -> Vec<(i64, u32)> {
    let mut v = vec![(0i64, 0u32)];
    let mut p = 1i64;
    for _ in 0..h {
        let mut next = Vec::with_capacity(v.len() * 3);
        for &(x, z) in &v {
            next.push((x - p, z));
            next.push((x, z + 1));
            next.push((x + p, z));
        }
        v = next;
        p *= 4;
    }
    v
}

fn energy_fast(k: u32, m: u32) -> u128 {
    let (zt, ct) = chunk_tables();
    let chunks = k.div_ceil(CHUNK);
    let pad = chunks * CHUNK - k;
    let half = (3u64.pow(k) - 1) / 2;
    let h = m.min(12);
    let low = low_strings(h);
    let reach = (4i64.pow(h) - 1) / 3;
    let top = 3u64.pow(m - h);
    let step = 4i64.pow(h);
    let next = std::sync::atomic::AtomicU64::new(0);
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    std::thread::scope(|sc| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                sc.spawn(|| {
                    let mut acc: u128 = 0;
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= top {
                            return acc;
                        }
                        let mut j = i;
                        let mut u = 0i64;
                        let mut zu = 0u32;
                        let mut p = 1i64;
                        for _ in h..m {
                            let d = (j % 3) as i64 - 1;
                            j /= 3;
                            u += d * p;
                            if d == 0 {
                                zu += 1;
                            }
                            p *= 4;
                        }
                        let base = u * step;
                        if base.unsigned_abs() > half + reach as u64 {
                            continue;
                        }
                        let mut part: u64 = 0;
                        for &(v, zv) in &low {
                            let mut t = (base + v).unsigned_abs();
                            if t > half {
                                continue;
                            }
                            let mut z = 0u32;
                            for _ in 0..chunks {
                                let r = (t % WIDTH) as usize;
                                z += zt[r] as u32;
                                t = t / WIDTH + ct[r];
                            }
                            part += 1u64 << (z - pad + zu + zv);
                        }
                        acc += part as u128;
                    }
                })
            })
            .collect();
        handles.into_iter().map(|t| t.join().unwrap()).sum()
    })
}

fn ladder_pairs(kmax: u32) -> Vec<(u32, u32)> {
    let mut v = Vec::new();
    for k in 6..=kmax {
        let p3 = 3u128.pow(k);
        for m in 1..=kmax {
            let p4 = 4u128.pow(m);
            if p3 < 3 * p4 && p4 < 4 * p3 {
                v.push((k, m));
            }
        }
    }
    v
}

fn least_squares(pts: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let slope = sxy / sxx;
    let icpt = my - slope * mx;
    let rss: f64 = pts.iter().map(|p| (p.1 - icpt - slope * p.0).powi(2)).sum();
    let se = (rss / (n - 2.0) / sxx).sqrt();
    (slope, icpt, se, (rss / n).sqrt())
}

fn fits(tag: &str, rows: &[(u32, f64)], k0: u32) {
    let pts: Vec<(u32, f64)> = rows.iter().filter(|r| r.0 >= k0).cloned().collect();
    if pts.len() < 3 {
        return;
    }
    let s = 2f64.ln() / 3f64.ln() - 0.5;
    let power: Vec<(f64, f64)> = pts
        .iter()
        .map(|&(k, q)| (k as f64, q.ln() / 3f64.ln()))
        .collect();
    let (eta, c, se, _) = least_squares(&power);
    let prss: f64 = pts
        .iter()
        .map(|&(k, q)| (q - 3f64.powf(c + eta * k as f64)).powi(2))
        .sum::<f64>();
    let sat: Vec<(f64, f64)> = pts
        .iter()
        .map(|&(k, q)| (3f64.powf(-s * k as f64), q))
        .collect();
    let (b, a, _, srms) = least_squares(&sat);
    let rss_q = |e: f64| {
        let g: Vec<f64> = pts.iter().map(|&(k, _)| 3f64.powf(e * k as f64)).collect();
        let c = pts.iter().zip(&g).map(|(p, g)| p.1 * g).sum::<f64>()
            / g.iter().map(|g| g * g).sum::<f64>();
        (
            pts.iter()
                .zip(&g)
                .map(|(p, g)| (p.1 - c * g).powi(2))
                .sum::<f64>(),
            c,
        )
    };
    let (mut lo, mut hi) = (-0.1f64, 0.1f64);
    for _ in 0..200 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if rss_q(m1).0 < rss_q(m2).0 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let eq = (lo + hi) / 2.0;
    let qrms = (rss_q(eq).0 / pts.len() as f64).sqrt();
    println!(
        "fit {} k >= {} over {} pairs: power in log_3 Q: eta = {:.6} se {:.6} rms in Q {:.6}; power fitted in Q: eta = {:.6} rms in Q {:.6}; saturation a + b 3^(-s k), s = {:.6}: a = {:.6} b = {:.6} rms in Q {:.6}",
        tag,
        k0,
        pts.len(),
        eta,
        se,
        (prss / pts.len() as f64).sqrt(),
        eq,
        qrms,
        s,
        a,
        b,
        srms
    );
}

fn ladder(kmax: u32, cache: Option<&str>) {
    let t0 = Instant::now();
    let mut known: std::collections::HashMap<(u32, u32), u128> = std::collections::HashMap::new();
    if let Some(path) = cache {
        if let Ok(text) = std::fs::read_to_string(path) {
            for line in text.lines() {
                let f: Vec<&str> = line.split_whitespace().collect();
                if f.len() == 3 {
                    known.insert(
                        (f[0].parse().unwrap(), f[1].parse().unwrap()),
                        f[2].parse().unwrap(),
                    );
                }
            }
        }
    }
    let mut all = Vec::new();
    let mut chain = Vec::new();
    let mut clean = Vec::new();
    let mut worst = (0u32, 0u32, 0u128, 1u128, 0f64);
    let mut least = (0u32, 0u32, 0u128, 1u128, f64::MAX);
    let mut whole = Vec::new();
    let mut chain_lo = (0u128, 1u128, f64::MAX);
    let mut chain_hi = (0u128, 1u128, 0f64);
    for (k, m) in ladder_pairs(kmax) {
        let t1 = Instant::now();
        let e = match known.get(&(k, m)) {
            Some(&e) => e,
            None => {
                let e = energy_fast(k, m);
                if let Some(path) = cache {
                    use std::io::Write;
                    let mut f = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(path)
                        .unwrap();
                    writeln!(f, "{} {} {}", k, m, e).unwrap();
                }
                e
            }
        };
        let pairs = 1u128 << (k + m);
        let d = (3u128.pow(k) - 1) / 2 + (4u128.pow(m) - 1) / 3;
        let range = d + 1;
        let is_clean = 3u128.pow(k) > d && 4u128.pow(m) > d;
        let copy = 2 * 4u128.pow(m) < 3u128.pow(k) + 5 || 3u128.pow(k + 1) < 4u128.pow(m) + 5;
        let on_chain = 4u128.pow(m) >= 3u128.pow(k) && 4u128.pow(m) < 4 * 3u128.pow(k);
        let q = (e as f64) * (range as f64) / (pairs as f64) / (pairs as f64);
        println!(
            "k = {:<2} m = {:<2} {}{} 4^m/3^k = {} energy {} Q {} bound {} {:.2} s",
            k,
            m,
            if on_chain { "chain " } else { "      " },
            if copy {
                "copy "
            } else if is_clean {
                "clean"
            } else {
                "mixed"
            },
            ratio(k, m),
            e,
            six_up(e * range, pairs * pairs),
            six(pairs * pairs, e * range),
            t1.elapsed().as_secs_f64()
        );
        let (num, den) = (e * range, pairs * pairs);
        if on_chain {
            whole.push(q);
            if q < chain_lo.2 {
                chain_lo = (num, den, q);
            }
            if q > chain_hi.2 {
                chain_hi = (num, den, q);
            }
        }
        if copy {
            continue;
        }
        all.push((k, q));
        if on_chain {
            chain.push((k, q));
        }
        if is_clean {
            clean.push((k, q));
        }
        if q > worst.4 {
            worst = (k, m, num, den, q);
        }
        if q < least.4 {
            least = (k, m, num, den, q);
        }
    }
    println!(
        "{} pairs without copies, {} of them on the chain; Q in [{}, {}], the lower end truncated and the upper rounded up, least at ({}, {}), largest at ({}, {}); the whole chain, copies included, {} levels, Q in [{}, {}], mean {:.6}",
        all.len(),
        chain.len(),
        six(least.2, least.3),
        six_up(worst.2, worst.3),
        least.0,
        least.1,
        worst.0,
        worst.1,
        whole.len(),
        six(chain_lo.0, chain_lo.1),
        six_up(chain_hi.0, chain_hi.1),
        whole.iter().sum::<f64>() / whole.len() as f64
    );
    for k0 in [6u32, 11, 16] {
        fits("all", &all, k0);
        fits("chain", &chain, k0);
        fits("clean", &clean, k0);
    }
    println!("{:.2} s", t0.elapsed().as_secs_f64());
}

fn constants() {
    let alpha = 2f64.ln() / 3f64.ln();
    let gamma = 1.0 / (2.0 - 2.0 * alpha);
    let eb = 6f64.powf(1.0 / (2.0 * gamma)) * gamma / (gamma - 1.0);
    let t = 30f64.powf(alpha);
    let c0 = 2.0 * t * eb + 3.0 * t * 6f64.sqrt();
    let up = |x: f64| (x * 1e6).ceil() / 1e6;
    println!(
        "alpha = {:.6} gamma = {:.6} E|b - b'|^(alpha - 1) <= {:.6} 30^alpha = {:.6} C_0 = 2 30^alpha E + 3 30^alpha 6^(1/2) = {:.6}, each rounded up",
        alpha,
        gamma,
        up(eb),
        up(t),
        up(c0)
    );
}

// PHASES

fn window_energy(k: u32, strings: &[(i64, u32)], sigma: f64, zt: &[u8], ct: &[u64]) -> u128 {
    let chunks = k.div_ceil(CHUNK);
    let pad = chunks * CHUNK - k;
    let half = (3u64.pow(k) - 1) / 2;
    let mut acc: u128 = 0;
    for &(v, zv) in strings {
        let mut t = (-(sigma * v as f64)).round().abs() as u64;
        if t > half {
            continue;
        }
        let mut z = 0u32;
        for _ in 0..chunks {
            let r = (t % WIDTH) as usize;
            z += zt[r] as u32;
            t = t / WIDTH + ct[r];
        }
        acc += 1u128 << (z - pad + zv);
    }
    acc
}

fn tent_energy(k: u32, strings: &[(i64, u32)], sigma: f64, zt: &[u8], ct: &[u64]) -> f64 {
    let chunks = k.div_ceil(CHUNK);
    let pad = chunks * CHUNK - k;
    let half = (3u64.pow(k) - 1) / 2;
    let mut acc = 0f64;
    for &(v, zv) in strings {
        let y = -(sigma * v as f64);
        let lo = y.floor();
        for c in [lo, lo + 1.0] {
            let w = 1.0 - (c - y).abs();
            let mut t = c.abs() as u64;
            if w <= 0.0 || t > half {
                continue;
            }
            let mut z = 0u32;
            for _ in 0..chunks {
                let r = (t % WIDTH) as usize;
                z += zt[r] as u32;
                t = t / WIDTH + ct[r];
            }
            acc += w * (1u64 << (z - pad + zv)) as f64;
        }
    }
    acc
}

fn phases(kmax: u32, grid: usize) {
    let t0 = Instant::now();
    let (zt, ct) = chunk_tables();
    let mut ranks = Vec::new();
    let mut qranks = Vec::new();
    let mut locals = Vec::new();
    let mut tents = Vec::new();
    let mut lratios = Vec::new();
    let mut tratios = Vec::new();
    for k in 8..=kmax {
        let mut m = 1u32;
        while 4u128.pow(m) < 3u128.pow(k) {
            m += 1;
        }
        let tau_k = 4f64.powi(m as i32) / 3f64.powi(k as i32);
        let strings = low_strings(m);
        let scale = 3f64.powi(k as i32) / 4f64.powi((k + m) as i32);
        let exact = energy_fast(k, m);
        let h0 = exact as f64 * scale;
        let next = std::sync::atomic::AtomicUsize::new(0);
        let threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4);
        let mut values: Vec<(usize, f64, f64)> = std::thread::scope(|sc| {
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    sc.spawn(|| {
                        let mut out = Vec::new();
                        loop {
                            let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if i >= grid {
                                return out;
                            }
                            let tau = 4f64.powf(i as f64 / grid as f64);
                            let w = window_energy(k, &strings, tau / tau_k, &zt, &ct);
                            let e = tent_energy(k, &strings, tau / tau_k, &zt, &ct);
                            out.push((i, w as f64 * scale, e * scale));
                        }
                    })
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect()
        });
        values.sort_by_key(|v| v.0);
        let below = values.iter().filter(|v| v.1 < h0).count();
        let flat = |i: usize, h: f64| h * (0.5 + 4f64.powf(i as f64 / grid as f64) / 3.0);
        let q0 = h0 * (0.5 + tau_k / 3.0);
        let qbelow = values.iter().filter(|v| flat(v.0, v.1) < q0).count();
        let ik = (tau_k.ln() / 4f64.ln() * grid as f64) as i64;
        let near: Vec<&(usize, f64, f64)> = values
            .iter()
            .filter(|v| (v.0 as i64 - ik).abs() <= (grid / 50) as i64)
            .collect();
        let local = near.iter().filter(|v| v.1 < h0).count() as f64 / near.len() as f64;
        let tent = near.iter().filter(|v| v.2 < h0).count() as f64 / near.len() as f64;
        let lmean = near.iter().map(|v| v.1).sum::<f64>() / near.len() as f64;
        let tmean = near.iter().map(|v| v.2).sum::<f64>() / near.len() as f64;
        tents.push(tent);
        lratios.push(h0 / lmean);
        tratios.push(h0 / tmean);
        qranks.push(qbelow as f64 / grid as f64);
        locals.push(local);
        let mean = values.iter().map(|v| v.1).sum::<f64>() / grid as f64;
        let lebmean = values
            .iter()
            .map(|v| v.1 * 4f64.powf(v.0 as f64 / grid as f64))
            .sum::<f64>()
            * 4f64.ln()
            / grid as f64;
        let (imax, hmax) =
            values.iter().fold(
                (0usize, 0f64),
                |a, v| if v.1 > a.1 { (v.0, v.1) } else { a },
            );
        let rank = below as f64 / grid as f64;
        ranks.push(rank);
        println!(
            "k = {:<2} m = {:<2} tau_k = {:.6} h_k(tau_k) = {:.6} rank {:.4} flat rank {:.4} local rank {:.4} tent local rank {:.4} local ratio {:.4} tent local ratio {:.4} log-mean {:.6} integral over [1, 4] {:.6} max {:.6} at tau = {:.6} {:.2} s",
            k,
            m,
            tau_k,
            h0,
            rank,
            qranks[qranks.len() - 1],
            local,
            tent,
            lratios[lratios.len() - 1],
            tratios[tratios.len() - 1],
            mean,
            lebmean,
            hmax,
            4f64.powf(imax as f64 / grid as f64),
            t0.elapsed().as_secs_f64()
        );
    }
    let n = ranks.len() as f64;
    let mr = ranks.iter().sum::<f64>() / n;
    let mq = qranks.iter().sum::<f64>() / n;
    let ml = locals.iter().sum::<f64>() / n;
    let mt = tents.iter().sum::<f64>() / n;
    println!(
        "mean rank {:.4}, mean flat rank {:.4}, mean local rank {:.4}, mean tent local rank {:.4} over {} levels, grid {}",
        mr,
        mq,
        ml,
        mt,
        ranks.len(),
        grid
    );
}

// CONTROL

const OEIS_A367090: [u64; 58] = [
    62, 63, 143, 144, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221,
    222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
    241, 242, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479,
    480,
];

fn same(a: &Bits, b: &Bits) -> bool {
    a.top == b.top && a.w == b.w
}

fn symmetric(s: &Bits, d: u64) -> bool {
    (0..=d).all(|x| s.get(x) == s.get(d - x))
}

fn complement(s: &Bits, n: usize) -> Vec<u64> {
    (1..=s.top).filter(|&x| !s.get(x)).take(n).collect()
}

fn control() {
    let t0 = Instant::now();
    let top = 3u64.pow(13);
    let direct = sumset_direct(top);
    let shifted = sumset(top, 3, 4);
    println!(
        "double loop against shift-or at 3^13: {}",
        if same(&direct, &shifted) {
            "agree"
        } else {
            "DIFFER"
        }
    );
    let top = 3u64.pow(17);
    let ab = sumset(top, 3, 4);
    let ba = sumset(top, 4, 3);
    println!(
        "A direct with B shifted against B direct with A shifted at 3^17: {}",
        if same(&ab, &ba) { "agree" } else { "DIFFER" }
    );
    let gaps = complement(&ab, OEIS_A367090.len());
    println!(
        "first {} non-members against A367090: {}",
        OEIS_A367090.len(),
        if gaps == OEIS_A367090 {
            "agree"
        } else {
            "DIFFER"
        }
    );
    for (r, s, d, clean) in centres(top) {
        println!(
            "centre r = {:<2} s = {:<2} d = {} {} symmetric {}",
            r,
            s,
            d,
            if clean { "clean" } else { "mixed" },
            symmetric(&ab, d)
        );
    }
    println!("{:.2} s", t0.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("density") => density(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(17)),
        Some("control") => control(),
        Some("energy") => energies(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(17)),
        Some("ladder") => ladder(
            args.get(2).and_then(|s| s.parse().ok()).unwrap_or(17),
            args.get(3).map(|s| s.as_str()),
        ),
        Some("constants") => constants(),
        Some("phases") => phases(
            args.get(2).and_then(|s| s.parse().ok()).unwrap_or(12),
            args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024),
        ),
        _ => {
            println!("verbs: density K, energy K, ladder K [cache], phases K G, constants, control")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shift_or_matches_double_loop() {
        let top = 3u64.pow(9);
        assert!(same(&sumset_direct(top), &sumset(top, 3, 4)));
        assert!(same(&sumset_direct(top), &sumset(top, 4, 3)));
    }

    #[test]
    fn first_gaps_are_a367090() {
        let s = sumset(3u64.pow(7), 3, 4);
        assert_eq!(complement(&s, OEIS_A367090.len()), OEIS_A367090.to_vec());
    }

    #[test]
    fn member_counts_are_powers_of_two() {
        assert_eq!(members(3, 3u64.pow(10) - 1).len(), 1 << 10);
        assert_eq!(members(4, 4u64.pow(6) - 1).len(), 1 << 6);
    }

    #[test]
    fn scan_counts_match_popcount() {
        let top = 3u64.pow(9);
        let s = sumset(top, 3, 4);
        let sc = scan(&s, &[top]);
        assert_eq!(sc.marks, vec![(top, s.count() - 1)]);
    }

    #[test]
    fn scan_extremes_match_direct_search() {
        let top = 3u64.pow(9);
        let s = sumset(top, 3, 4);
        let sc = scan(&s, &[]);
        for w in sc.windows3.iter().chain(sc.windows2.iter()) {
            let mut c = (1..w.lo).filter(|&x| s.get(x)).count() as u64;
            let mut best_max = Ext::none();
            let mut best_min = Ext::none();
            for x in w.lo..=w.hi {
                if s.get(x) {
                    c += 1;
                }
                if best_max.above(c, x) {
                    best_max = Ext { c, x };
                }
                if best_min.below(c, x) {
                    best_min = Ext { c, x };
                }
            }
            assert_eq!((w.max.c, w.max.x), (best_max.c, best_max.x));
            assert_eq!((w.min.c, w.min.x), (best_min.c, best_min.x));
        }
    }

    #[test]
    fn energy_matches_representation_histogram() {
        for (k, m) in [(4, 3), (6, 5), (9, 7)] {
            assert_eq!(energy(k, m), energy_direct(k, m));
        }
    }

    #[test]
    fn chunked_energy_matches_digit_energy() {
        for k in 1..=11 {
            for m in 1..=14 {
                assert_eq!(energy_fast(k, m), energy(k, m), "k = {} m = {}", k, m);
            }
        }
        assert_eq!(energy_fast(9, 7), energy_direct(9, 7));
    }

    #[test]
    fn gap_copies_double_the_energy() {
        for (k, m) in [(6, 4), (7, 5), (11, 8), (12, 9)] {
            assert!(2 * 4u128.pow(m) < 3u128.pow(k) + 5);
            assert_eq!(energy_fast(k, m), 2 * energy_fast(k - 1, m));
        }
        assert_ne!(energy_fast(13, 10), 2 * energy_fast(12, 10));
        for (k, m) in [(9, 8), (14, 12)] {
            assert!(3u128.pow(k + 1) < 4u128.pow(m) + 5);
            assert_eq!(energy_fast(k, m), 2 * energy_fast(k, m - 1));
        }
    }

    #[test]
    fn window_energy_at_the_lattice_is_the_energy() {
        let (zt, ct) = chunk_tables();
        for (k, m) in [(6, 5), (9, 7), (10, 8), (12, 10)] {
            assert_eq!(
                window_energy(k, &low_strings(m), 1.0, &zt, &ct),
                energy(k, m)
            );
            assert_eq!(
                tent_energy(k, &low_strings(m), 1.0, &zt, &ct),
                energy(k, m) as f64
            );
        }
    }

    #[test]
    fn prefix_count_matches_bit_scan() {
        let s = sumset(3u64.pow(9), 3, 4);
        for d in [0u64, 63, 64, 705, 15302, 19683] {
            assert_eq!(s.count_to(d), (0..=d).filter(|&x| s.get(x)).count() as u64);
        }
    }

    #[test]
    fn clean_centres_are_symmetric() {
        let top = 3u64.pow(11);
        let s = sumset(top, 3, 4);
        for (_, _, d, clean) in centres(top) {
            if clean {
                assert!(symmetric(&s, d));
            }
        }
    }
}
