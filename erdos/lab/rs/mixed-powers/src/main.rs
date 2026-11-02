use std::collections::HashMap;
use std::env;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

// SETS

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn sigma(bases: &[u64]) -> (u128, u128) {
    let den = bases.iter().fold(1u128, |l, &d| {
        let e = (d - 1) as u128;
        l / gcd(l, e) * e
    });
    let num = bases.iter().map(|&d| den / (d - 1) as u128).sum();
    (num, den)
}

fn coprime(bases: &[u64]) -> bool {
    bases.iter().fold(0u128, |g, &d| gcd(g, d as u128)) == 1
}

fn family(bases: &[u64]) -> bool {
    let (num, den) = sigma(bases);
    !bases.is_empty() && num > den && coprime(bases)
}

fn subsets(r: u64) -> Vec<Vec<u64>> {
    let pool: Vec<u64> = (3..=r).collect();
    (1u64..(1 << pool.len()))
        .map(|mask| {
            pool.iter()
                .enumerate()
                .filter(|(i, _)| mask >> i & 1 == 1)
                .map(|(_, &d)| d)
                .collect()
        })
        .collect()
}

fn order(sets: &mut [Vec<u64>]) {
    sets.sort_by(|a, b| {
        a.iter()
            .max()
            .cmp(&b.iter().max())
            .then(a.len().cmp(&b.len()))
            .then(a.cmp(b))
    });
}

fn minimal(r: u64) -> Vec<Vec<u64>> {
    let mut out: Vec<Vec<u64>> = subsets(r)
        .into_iter()
        .filter(|s| family(s))
        .filter(|s| {
            (0..s.len()).all(|i| {
                let mut t = s.clone();
                t.remove(i);
                !family(&t)
            })
        })
        .collect();
    order(&mut out);
    out
}

fn unit(r: u64) -> Vec<Vec<u64>> {
    let mut out: Vec<Vec<u64>> = subsets(r)
        .into_iter()
        .filter(|s| {
            let (num, den) = sigma(s);
            num == den && coprime(s)
        })
        .collect();
    order(&mut out);
    out
}

// ELEMENTS

fn elements(bases: &[u64], k: u32, cap: u128) -> Vec<u128> {
    let mut v = vec![];
    for &d in bases {
        let mut x = (d as u128).pow(k);
        while x <= cap {
            v.push(x);
            x *= d as u128;
        }
    }
    v.sort();
    v
}

// BITS

struct Bits {
    w: Vec<u64>,
}

impl Bits {
    fn new(bits: usize) -> Bits {
        let mut w = vec![0u64; bits / 64 + 2];
        w[0] = 1;
        Bits { w }
    }

    fn fold(&mut self, a: usize, top: usize) {
        let q = a / 64;
        let r = (a % 64) as u32;
        for i in (q..=top / 64).rev() {
            let v = if r == 0 {
                self.w[i - q]
            } else {
                let lo = if i > q {
                    self.w[i - q - 1] >> (64 - r)
                } else {
                    0
                };
                (self.w[i - q] << r) | lo
            };
            self.w[i] |= v;
        }
    }

    fn hole(&self, h: usize) -> Option<usize> {
        let mut i = h / 64;
        let keep = h % 64;
        let mut word = !self.w[i];
        if keep < 63 {
            word &= (1u64 << (keep + 1)) - 1;
        }
        loop {
            if word != 0 {
                return Some(i * 64 + 63 - word.leading_zeros() as usize);
            }
            if i == 0 {
                return None;
            }
            i -= 1;
            word = !self.w[i];
        }
    }

    fn gaps(&self, below: usize) -> usize {
        (0..below)
            .filter(|&x| self.w[x / 64] >> (x % 64) & 1 == 0)
            .count()
    }

    fn get(&self, x: usize) -> bool {
        self.w[x / 64] >> (x % 64) & 1 == 1
    }

    fn each_zero(&self, lo: usize, hi: usize, mut f: impl FnMut(usize) -> bool) -> bool {
        if lo > hi {
            return false;
        }
        for i in lo / 64..=hi / 64 {
            let mut word = !self.w[i];
            while word != 0 {
                let x = i * 64 + word.trailing_zeros() as usize;
                word &= word - 1;
                if x >= lo && x <= hi && f(x) {
                    return true;
                }
            }
        }
        false
    }

    fn hole_in(&self, lo: usize, h: usize) -> Option<usize> {
        if lo > h {
            return None;
        }
        let mut i = h / 64;
        let keep = h % 64;
        let mut word = !self.w[i];
        if keep < 63 {
            word &= (1u64 << (keep + 1)) - 1;
        }
        loop {
            if word != 0 {
                let x = i * 64 + 63 - word.leading_zeros() as usize;
                return if x >= lo { Some(x) } else { None };
            }
            if i == lo / 64 {
                return None;
            }
            i -= 1;
            word = !self.w[i];
        }
    }
}

// CERTIFY

#[derive(Debug, Clone, PartialEq)]
enum Verdict {
    Found {
        f: u128,
        gaps: usize,
        n: usize,
        last: u128,
        sum: u128,
        reach: u128,
    },
    Cap,
}

fn far() -> u128 {
    10u128.pow(30)
}

fn grows(a: &[u128], n: usize, s: u128, t: u128, surplus: (u128, u128, u128)) -> Option<u128> {
    let (gain, den, c) = surplus;
    let mut sm = s;
    for &x in &a[n + 1..] {
        if gain * x >= c + den * (2 * t).saturating_sub(1) {
            return Some(x);
        }
        if x + 2 * t > sm + 1 {
            return None;
        }
        sm += x;
    }
    None
}

fn certify(bases: &[u64], k: u32, cap: usize) -> Verdict {
    let (num, den) = sigma(bases);
    if num <= den {
        return Verdict::Cap;
    }
    let c: u128 = bases
        .iter()
        .map(|&d| (d as u128).pow(k) * (den / (d - 1) as u128))
        .sum();
    let a = elements(bases, k, far());
    let mut bits = Bits::new(cap + 64);
    let mut s: u128 = 0;
    for n in 0..a.len() - 1 {
        if s + a[n] > cap as u128 {
            return Verdict::Cap;
        }
        bits.fold(a[n] as usize, (s + a[n]) as usize);
        s += a[n];
        let t = bits.hole((s / 2) as usize).map_or(0, |x| x as u128 + 1);
        if t > a[n + 1] || t == 0 {
            continue;
        }
        if let Some(reach) = grows(&a, n, s, t, (num - den, den, c)) {
            return Verdict::Found {
                f: t - 1,
                gaps: bits.gaps(t as usize),
                n: n + 1,
                last: a[n],
                sum: s,
                reach,
            };
        }
    }
    Verdict::Cap
}

// WIDE

#[derive(Debug, Clone, PartialEq, Default)]
struct Lanes {
    seed: usize,
    window: usize,
    pair: usize,
    last: u8,
}

fn lane(a: &[u128], pre: &[u128], z: u128) -> u8 {
    let c = a.partition_point(|&t| t <= z);
    if c == 0 {
        0
    } else if z > pre[c - 1] {
        1
    } else {
        2
    }
}

fn wide(bases: &[u64], k: u32, bits: u32) -> (Verdict, Lanes) {
    wide_with(bases, k, bits, |_, _, _| {})
}

fn wide_with(
    bases: &[u64],
    k: u32,
    bits: u32,
    mut log: impl FnMut(usize, u128, u128),
) -> (Verdict, Lanes) {
    let (num, den) = sigma(bases);
    if num <= den {
        return (Verdict::Cap, Lanes::default());
    }
    let c: u128 = bases
        .iter()
        .map(|&d| (d as u128).pow(k) * (den / (d - 1) as u128))
        .sum();
    let a = elements(bases, k, far());
    let mut pre = vec![0u128];
    for &x in &a {
        pre.push(pre[pre.len() - 1] + x);
    }
    let w = 1usize << bits;
    let wu = w as u128;
    let mut b = Bits::new(w + 64);
    let limit = 1usize << 25;
    let mut sparse: Vec<u64> = vec![];
    let mut next: Vec<u64> = vec![];
    let mut s: u128 = 0;
    let mut seen = 0usize;
    let mut top: Option<usize> = None;
    for n in 0..a.len() - 1 {
        let x = a[n];
        let snew = s + x;
        if snew >> 63 != 0 {
            return (Verdict::Cap, Lanes::default());
        }
        if snew / 2 > wu {
            let half = snew / 2;
            if x > s && half > s.max(wu) && half - s.max(wu) > limit as u128 {
                return (Verdict::Cap, Lanes::default());
            }
            let mem = |z: i128| -> bool {
                if z < 0 || z as u128 > s {
                    return false;
                }
                let z = z as u128;
                if z <= wu {
                    b.get(z as usize)
                } else if z >= s - wu {
                    b.get((s - z) as usize)
                } else if 2 * z <= s {
                    sparse.binary_search(&(z as u64)).is_err()
                } else {
                    sparse.binary_search(&((s - z) as u64)).is_err()
                }
            };
            next.clear();
            next.reserve(limit + 1);
            let mut keep = |z: u128| -> bool {
                if z > wu && z <= half && !mem(z as i128 - x as i128) {
                    next.push(z as u64);
                }
                next.len() > limit
            };
            let mut spill = false;
            for &g in &sparse {
                spill |= keep(g as u128) || keep(s - g as u128);
            }
            if s > wu {
                let d = s as i128 - x as i128;
                let glo = if d > 0 { ((d + 1) / 2) as usize } else { 0 };
                let ghi = w.min((s - wu - 1) as usize);
                spill |= b.each_zero(glo, ghi, |g| keep(s - g as u128));
            }
            if x > s {
                for z in s.max(wu) + 1..=half {
                    spill |= keep(z);
                }
            }
            if spill {
                return (Verdict::Cap, Lanes::default());
            }
            next.sort_unstable();
            next.dedup();
            std::mem::swap(&mut sparse, &mut next);
        }
        if x <= wu {
            b.fold(x as usize, w.min(snew as usize));
        }
        s = snew;
        let h = w.min((s / 2) as usize);
        let fresh = if h >= seen { b.hole_in(seen, h) } else { None };
        top = fresh.or_else(|| top.and_then(|t| b.hole_in(0, t)));
        seen = seen.max(h + 1);
        let t = match sparse.last() {
            Some(&z) => z as u128 + 1,
            None => top.map_or(0, |z| z as u128 + 1),
        };
        log(n + 1, s, t);
        if t > a[n + 1] || t == 0 {
            continue;
        }
        if let Some(reach) = grows(&a, n, s, t, (num - den, den, c)) {
            let mut lanes = Lanes::default();
            let mut count = |z: u128| match lane(&a, &pre, z) {
                0 => lanes.seed += 1,
                1 => lanes.window += 1,
                _ => lanes.pair += 1,
            };
            b.each_zero(0, w.min((t - 1) as usize), |z| {
                count(z as u128);
                false
            });
            for &z in &sparse {
                count(z as u128);
            }
            lanes.last = lane(&a, &pre, t - 1);
            let v = Verdict::Found {
                f: t - 1,
                gaps: lanes.seed + lanes.window + lanes.pair,
                n: n + 1,
                last: a[n],
                sum: s,
                reach,
            };
            return (v, lanes);
        }
    }
    (Verdict::Cap, Lanes::default())
}

// CONTROL

fn knapsack(bases: &[u64], k: u32, b: usize) -> Vec<bool> {
    let mut r = vec![false; b + 1];
    r[0] = true;
    for a in elements(bases, k, b as u128) {
        let a = a as usize;
        for x in (a..=b).rev() {
            if r[x - a] {
                r[x] = true;
            }
        }
    }
    r
}

fn control() {
    let cells: [(&[u64], u32); 9] = [
        (&[3, 4, 5], 1),
        (&[3, 4, 5], 2),
        (&[3, 4, 5], 3),
        (&[3, 4, 6], 2),
        (&[3, 5, 6, 7], 3),
        (&[4, 5, 6, 7, 8], 3),
        (&[3, 4, 7, 8], 3),
        (&[3, 5, 6, 9], 2),
        (&[3, 4, 9, 10], 2),
    ];
    for (bases, k) in cells {
        let clock = Instant::now();
        let v = certify(bases, k, 1 << 30);
        let Verdict::Found { f, gaps, .. } = v else {
            println!("{:?} k {} cap", bases, k);
            continue;
        };
        let b = (4 * f as usize).max(1000);
        let r = knapsack(bases, k, b);
        let last = (0..=b).rev().find(|&x| !r[x]).unwrap();
        let count = (0..=b).filter(|&x| !r[x]).count();
        println!(
            "{:?} k {} certificate F {} gaps {} | knapsack to {} F {} gaps {} | {} | {:.2} s",
            bases,
            k,
            f,
            gaps,
            b,
            last,
            count,
            if last as u128 == f && count == gaps {
                "agree"
            } else {
                "DISAGREE"
            },
            clock.elapsed().as_secs_f64()
        );
    }
    for bases in [&[3u64, 4][..], &[3, 6, 9, 12, 15, 21][..]] {
        println!("{:?} k 1 {:?} to 2^24", bases, certify(bases, 1, 1 << 24));
    }
}

fn wcontrol(bits: u32) {
    let clock = Instant::now();
    let (mut agree, mut differ, mut capped, mut past) = (0, 0, 0, 0);
    for set in minimal(10) {
        for k in 1..=4 {
            let full = certify(&set, k, 1 << 30);
            if full == Verdict::Cap {
                continue;
            }
            let (v, _) = wide(&set, k, bits);
            match &v {
                Verdict::Cap => capped += 1,
                Verdict::Found { sum, .. } => {
                    past += (*sum > 2u128 << bits) as usize;
                    if v == full {
                        agree += 1;
                    } else {
                        differ += 1;
                        println!("{} k {} full {:?} wide {:?}", show(&set), k, full, v);
                    }
                }
            }
        }
    }
    println!(
        "wide at 2^{} against the full certificate at 2^30, minimal sets in [3, 10], k = 1..4: {} agree, {} differ, {} capped, {} certified with S_n past twice the window; {:.1} s",
        bits, agree, differ, capped, past, clock.elapsed().as_secs_f64()
    );
}

fn split(bases: &[u64], k: u32, x0: u128, count: u128) {
    let clock = Instant::now();
    let terms = elements(bases, k, x0 + count);
    let (l, r) = terms.split_at(terms.len() / 2);
    let sums = |t: &[u128]| -> Vec<u128> {
        let mut v = vec![0u128];
        for &x in t {
            for i in 0..v.len() {
                v.push(v[i] + x);
            }
        }
        v.sort();
        v
    };
    let (ls, rs) = (sums(l), sums(r));
    let hit = |x: u128| -> bool {
        let (mut i, mut j) = (0usize, rs.len());
        while i < ls.len() && j > 0 {
            let v = ls[i] + rs[j - 1];
            if v == x {
                return true;
            }
            if v < x {
                i += 1;
            } else {
                j -= 1;
            }
        }
        false
    };
    let out: Vec<u128> = (x0..x0 + count).filter(|&x| !hit(x)).collect();
    println!(
        "{} k {}: meet in the middle over {} terms, the non-sums in [{}, {}] are {:?}; {:.1} s",
        show(bases),
        k,
        terms.len(),
        x0,
        x0 + count - 1,
        out,
        clock.elapsed().as_secs_f64()
    );
}

// CENSUS

fn show(bases: &[u64]) -> String {
    format!(
        "{{{}}}",
        bases
            .iter()
            .map(|d| d.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

type Row = Vec<(u32, Verdict, Lanes, f64)>;

fn census(r: u64, kmax: u32, bits: u32, threads: usize, deep: bool) {
    let clock = Instant::now();
    let sets = minimal(r);
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<(usize, Row)>> = Mutex::new(vec![]);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= sets.len() {
                    break;
                }
                let mut row = vec![];
                for k in 1..=kmax {
                    let t = Instant::now();
                    let (v, lanes) = if deep {
                        wide(&sets[i], k, bits)
                    } else {
                        (certify(&sets[i], k, 1usize << bits), Lanes::default())
                    };
                    let stop = v == Verdict::Cap;
                    row.push((k, v, lanes, t.elapsed().as_secs_f64()));
                    if stop {
                        break;
                    }
                }
                rows.lock().unwrap().push((i, row));
            });
        }
    });
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by_key(|x| x.0);
    let mut depth = vec![0usize; kmax as usize + 1];
    for (i, row) in &rows {
        let (num, den) = sigma(&sets[*i]);
        let mut line = format!(
            "{} sigma {}/{}",
            show(&sets[*i]),
            num / gcd(num, den),
            den / gcd(num, den)
        );
        let mut d = 0;
        for (k, v, lanes, t) in row {
            match v {
                Verdict::Found {
                    f,
                    gaps,
                    n,
                    last,
                    sum,
                    reach,
                } => {
                    d = *k;
                    line += &format!(
                        " | k {} F {} gaps {} n {} last {} sum {} reach {} {:.2}s",
                        k, f, gaps, n, last, sum, reach, t
                    );
                    if deep {
                        line += &format!(
                            " lanes {}/{}/{} F {}",
                            lanes.seed,
                            lanes.window,
                            lanes.pair,
                            ["seed", "window", "pair"][lanes.last as usize]
                        );
                    }
                }
                Verdict::Cap => line += &format!(" | k {} cap 2^{}", k, bits),
            }
        }
        depth[d as usize] += 1;
        println!("{}", line);
    }
    println!("minimal sets {} below {}", sets.len(), r + 1);
    for (d, count) in depth.iter().enumerate() {
        if *count > 0 {
            println!("certified to k = {}: {} sets", d, count);
        }
    }
    if deep {
        for k in 1..=kmax {
            let mut last = [0usize; 3];
            let mut total = [0usize; 3];
            for (_, row) in &rows {
                for (kk, v, lanes, _) in row {
                    if *kk == k && *v != Verdict::Cap {
                        last[lanes.last as usize] += 1;
                        total[0] += lanes.seed;
                        total[1] += lanes.window;
                        total[2] += lanes.pair;
                    }
                }
            }
            println!(
                "k = {}: F below the least term in {} cells, by the window route in {}, by the pair route in {}; non-sums seed/window/pair {}/{}/{}",
                k, last[0], last[1], last[2], total[0], total[1], total[2]
            );
        }
    }
    let floor = rows
        .iter()
        .map(|(_, row)| row.iter().filter(|x| x.1 != Verdict::Cap).count())
        .min()
        .unwrap_or(0);
    println!(
        "every D in [3, {}] with sigma > 1 and gcd 1 is complete at k = 1..{}",
        r, floor
    );
    for s in unit(r) {
        println!("sigma = 1, gcd 1, outside the certificate: {}", show(&s));
    }
    println!("{:.1} s", clock.elapsed().as_secs_f64());
}

fn cell(bases: &[u64], k: u32, bits: u32) {
    let clock = Instant::now();
    let v = certify(bases, k, 1usize << bits);
    println!(
        "{} k {} {:?} {:.1} s",
        show(bases),
        k,
        v,
        clock.elapsed().as_secs_f64()
    );
}

fn wcell(bases: &[u64], k: u32, bits: u32) {
    let clock = Instant::now();
    let (v, lanes) = wide(bases, k, bits);
    println!(
        "{} k {} wide 2^{} {:?} {:?} {:.1} s",
        show(bases),
        k,
        bits,
        v,
        lanes,
        clock.elapsed().as_secs_f64()
    );
}

fn sets(bases: &[u64], k: u32, b: usize) {
    let mut terms = elements(bases, k, b as u128);
    terms.dedup();
    let mut r = vec![false; b + 1];
    r[0] = true;
    for a in terms {
        let a = a as usize;
        for x in (a..=b).rev() {
            if r[x - a] {
                r[x] = true;
            }
        }
    }
    let last = (0..=b).rev().find(|&x| !r[x]).unwrap();
    println!(
        "{} k {} as a set: largest non-sum below {} is {}",
        show(bases),
        k,
        b,
        last
    );
}

// WINDOWS

fn windows(bases: &[u64], k: u32, lo: u128, hi: u128) {
    let a = elements(bases, k, far());
    let mut sums = vec![];
    let mut s: u128 = 0;
    for &x in &a {
        s += x;
        sums.push(s);
    }
    let open: Vec<usize> = (0..a.len() - 2)
        .filter(|&n| a[n + 1] >= lo && a[n + 1] <= hi && a[n + 2] > sums[n])
        .collect();
    let seen = (0..a.len() - 2)
        .filter(|&n| a[n + 1] >= lo && a[n + 1] <= hi)
        .count();
    let mut chains = 0;
    for &m in &open {
        for &n in &open {
            if n > m && sums[n].saturating_sub(a[n + 1]) < a[m + 2] && sums[m] < a[n + 2] - a[n + 1]
            {
                chains += 1;
            }
        }
    }
    println!(
        "{} k {}: {} of {} terms a_(n+1) in [{}, {}] open a window a_(n+2) > S_n; {} pairs m < n of open windows where (S_m, a_(m+2)) meets (S_n - a_(n+1), a_(n+2) - a_(n+1))",
        show(bases), k, open.len(), seen, lo, hi, chains
    );
}

fn graham(t: u128, p: u128, q: u128, bits: u32) {
    let cap = 1usize << bits;
    let mut terms = vec![];
    let (mut num, mut den) = (t, 1u128);
    loop {
        num *= p;
        den *= q;
        let x = num / den;
        if x as usize > cap {
            break;
        }
        if x > 0 {
            terms.push(x as usize);
        }
    }
    terms.dedup();
    let mut bits_ = Bits::new(terms.iter().sum::<usize>() + 64);
    let mut s = 0usize;
    let mut sums = vec![];
    for &x in &terms {
        bits_.fold(x, s + x);
        s += x;
        sums.push(s);
    }
    let mut chained = 0;
    let mut windows = 0;
    for n in 0..terms.len().saturating_sub(2) {
        if terms[n + 2] > sums[n] && terms[n + 2] <= cap {
            windows += 1;
            if (sums[n] + 1..terms[n + 2]).any(|x| bits_.w[x / 64] >> (x % 64) & 1 == 0) {
                chained += 1;
            }
        }
    }
    let last = bits_
        .hole(cap / 2)
        .map_or("none".to_string(), |x| x.to_string());
    println!(
        "floor({} ({}/{})^n), n >= 1, {} terms to 2^{}: {} windows a_(n+2) > S_n, {} of them holding a non-sum; largest non-sum up to 2^{} is {}",
        t, p, q, terms.len(), bits, windows, chained, bits - 1, last
    );
}

fn fibonacci(cap: usize) -> Vec<usize> {
    let mut v = vec![1usize, 3];
    while v[v.len() - 1] + v[v.len() - 2] < cap {
        let n = v.len();
        v.push(v[n - 1] + v[n - 2] + 1);
    }
    v
}

fn missing(terms: &[usize]) -> Bits {
    let mut bits = Bits::new(terms.iter().sum::<usize>() + 64);
    let mut s = 0;
    for &x in terms {
        bits.fold(x, s + x);
        s += x;
    }
    bits
}

fn refute(bits_: u32) {
    let terms = fibonacci(1 << bits_);
    let bits = missing(&terms);
    let mut s = 0i128;
    let mut least = i128::MAX;
    let (mut windows, mut held) = (0, 0);
    let mut sums = vec![];
    for &x in &terms {
        s += x as i128;
        sums.push(s as usize);
    }
    for n in 0..terms.len() - 1 {
        least = least.min(2 * (sums[n] as i128 - terms[n + 1] as i128) - terms[n + 1] as i128);
    }
    for n in 0..terms.len() - 2 {
        if terms[n + 2] > sums[n] {
            windows += 1;
            if (sums[n] + 1..terms[n + 2]).any(|x| bits.w[x / 64] >> (x % 64) & 1 == 0) {
                held += 1;
            }
        }
    }
    let low: Vec<usize> = (1..1_000_000)
        .filter(|&x| bits.w[x / 64] >> (x % 64) & 1 == 0)
        .collect();
    println!(
        "2 F_m - 1, m >= 2, {} terms to 2^{}: least 2 (S_n - a_(n+1)) - a_(n+1) = {}; {} windows, {} holding a non-sum; {} non-sums below 10^6: {:?}",
        terms.len(), bits_, least, windows, held, low.len(), low
    );
}

fn kick(n: usize) -> usize {
    2 - (n / 2) % 2
}

fn chainless_terms(cap: usize) -> Vec<usize> {
    let mut a = vec![4usize, 5, 6];
    loop {
        let n = a.len() - 1;
        let v = a[..n].iter().sum::<usize>() - kick(n) - kick(n + 1);
        if v > cap {
            return a;
        }
        a.push(v);
    }
}

fn forced(bits_: u32) {
    let cap = 1usize << bits_;
    for p in [6usize, 8, 10] {
        for c in [2usize, 3] {
            let mut a = vec![4usize, 5, 6];
            loop {
                let n = a.len() - 1;
                let s: usize = a[..n].iter().sum();
                let v = if n >= 6 && n % p == 0 {
                    s + c
                } else {
                    s - kick(n) - kick(n + 1)
                };
                if v > cap {
                    break;
                }
                a.push(v);
            }
            let bits = missing(&a);
            let mut sums = vec![0usize];
            for &x in &a {
                sums.push(sums[sums.len() - 1] + x);
            }
            let least = (1..a.len())
                .map(|n| 2 * sums[n] as i64 - 3 * a[n] as i64)
                .min()
                .unwrap();
            let open: Vec<usize> = (2..a.len() - 1).filter(|&n| a[n + 1] > sums[n]).collect();
            let mut chains = 0;
            for &m in &open {
                for &n in &open {
                    if n > m && sums[n] - a[n] < a[m + 1] && sums[m] < a[n + 1] - a[n] {
                        chains += 1;
                    }
                }
            }
            let top = a[a.len() - 1];
            let mut total = 0;
            bits.each_zero(1, top - 1, |_| {
                total += 1;
                false
            });
            let mut late = 0;
            bits.each_zero(a[a.len() - 13], top - 1, |_| {
                late += 1;
                false
            });
            println!(
                "E with a_(n+2) = S_n + {} at n >= 6, n = 0 mod {}: {} terms to 2^{}, {} open windows at n >= 2, {} chained pairs, least 2 (S_n - a_(n+1)) - a_(n+1) = {}; {} non-sums below the last term, {} in the last 12 term intervals",
                c, p, a.len(), bits_, open.len(), chains, least, total, late
            );
        }
    }
}

fn chainless(bits_: u32) {
    let a = chainless_terms(1 << bits_);
    let bits = missing(&a);
    let mut sums = vec![0usize];
    for &x in &a {
        sums.push(sums[sums.len() - 1] + x);
    }
    let top = a[a.len() - 1];
    let early = (0..2).filter(|&n| a[n + 1] > sums[n]).count();
    let late = (2..a.len() - 1).filter(|&n| a[n + 1] > sums[n]).count();
    let least = (1..a.len())
        .map(|n| 2 * sums[n] as i64 - 3 * a[n] as i64)
        .min()
        .unwrap();
    let mut kicks: Vec<i64> = (2..a.len() - 2)
        .map(|n| a[n + 2] as i64 - a[n + 1] as i64 - a[n] as i64)
        .collect();
    kicks.sort();
    kicks.dedup();
    let chain: Vec<usize> = (3..a.len() - 1)
        .filter(|&n| a[n] + kick(n) < a[n + 1])
        .map(|n| a[n] + kick(n))
        .collect();
    let held = chain.iter().filter(|&&x| !bits.get(x)).count();
    let wide_a: Vec<u128> = a.iter().map(|&x| x as u128).collect();
    let pre: Vec<u128> = sums.iter().map(|&x| x as u128).collect();
    let mut lanes = [0usize; 3];
    bits.each_zero(1, top - 1, |z| {
        lanes[lane(&wide_a, &pre, z as u128) as usize] += 1;
        false
    });
    let low: Vec<usize> = (1..200).filter(|&x| !bits.get(x)).collect();
    println!(
        "a_1..a_3 = 4, 5, 6 and a_(n+2) = S_n - y_n - y_(n+1), y_n = 2 - (floor(n/2) mod 2), {} terms to 2^{}: {:?} ...; windows a_(n+2) > S_n at n < 2: {}, at n >= 2: {}; least 2 (S_n - a_(n+1)) - a_(n+1) = {}; a_(n+3) - a_(n+2) - a_(n+1) over n >= 2 takes {:?}; chain a_(n+1) + y_n, n >= 3: {} of {} non-sums, first {:?}; non-sums below {} by lane seed/window/pair: {}/{}/{}; non-sums below 200: {:?}",
        a.len(), bits_, &a[..12], early, late, least, kicks, held, chain.len(), &chain[..8], top, lanes[0], lanes[1], lanes[2], low
    );
}

fn windowed_terms(cap: usize) -> (Vec<usize>, Vec<i64>) {
    let mut a = vec![4usize, 5, 6];
    let mut u = vec![0i64, 0];
    let rule = |n: usize, a: &[usize], u: &[i64]| -> i64 {
        let x = |m: usize| a[m] as i64 + u[m];
        match n % 5 {
            0 if n >= 5 => 2,
            1 if n >= 6 => -x(n - 2),
            2 if n >= 7 => x(n - 3),
            _ => 1,
        }
    };
    loop {
        let n = a.len() - 1;
        while u.len() <= n + 1 {
            let v = rule(u.len(), &a, &u);
            u.push(v);
        }
        let v = a[..n].iter().sum::<usize>() as i64 - u[n] - u[n + 1];
        if v as usize > cap {
            return (a, u);
        }
        a.push(v as usize);
    }
}

fn windowed(bits_: u32) {
    let (a, u) = windowed_terms(1 << bits_);
    let bits = missing(&a);
    let mut sums = vec![0usize];
    for &x in &a {
        sums.push(sums[sums.len() - 1] + x);
    }
    let len = a.len();
    let top = a[len - 1];
    let sorted = a.windows(2).all(|w| w[0] <= w[1]);
    let open: Vec<usize> = (2..len - 1).filter(|&n| a[n + 1] > sums[n]).collect();
    let fifth = open.iter().all(|&n| n % 5 == 0)
        && (2..len - 1).filter(|&n| n % 5 == 0).count() == open.len();
    let size = open
        .iter()
        .map(|&n| 1000 * (a[n + 1] - sums[n]) / a[n + 1])
        .min()
        .unwrap();
    let mut chains = 0;
    for &m in &open {
        for &n in &open {
            if n > m && sums[n] - a[n] < a[m + 1] && sums[m] < a[n + 1] - a[n] {
                chains += 1;
            }
        }
    }
    let held = open
        .iter()
        .filter(|&&n| {
            let mut hit = false;
            bits.each_zero(sums[n] + 1, a[n + 1] - 1, |_| {
                hit = true;
                true
            });
            hit
        })
        .count();
    let least = (3..len)
        .map(|n| 12 * (sums[n] as i64 - a[n] as i64) - a[n] as i64)
        .min()
        .unwrap();
    let steep = (1..len)
        .map(|n| a[n] as i64 - 3 * a[n - 1] as i64)
        .max()
        .unwrap();
    let mut kicks: Vec<i64> = (3..len - 2)
        .filter(|&n| n % 5 == 3)
        .map(|n| a[n + 2] as i64 - a[n + 1] as i64 - a[n] as i64)
        .collect();
    kicks.sort();
    kicks.dedup();
    let chain: Vec<usize> = (3..len - 1)
        .filter(|&n| ((a[n] as i64 + u[n]) as usize) < a[n + 1])
        .map(|n| (a[n] as i64 + u[n]) as usize)
        .collect();
    let gone = chain.iter().filter(|&&x| !bits.get(x)).count();
    let wide_a: Vec<u128> = a.iter().map(|&x| x as u128).collect();
    let pre: Vec<u128> = sums.iter().map(|&x| x as u128).collect();
    let mut lanes = [0usize; 3];
    bits.each_zero(1, top - 1, |z| {
        lanes[lane(&wide_a, &pre, z as u128) as usize] += 1;
        false
    });
    println!(
        "V: a_1..a_3 = 4, 5, 6, a_(n+2) = S_n - u_n - u_(n+1), u_n = 2, -x_(n-2), x_(n-3) at n = 0, 1, 2 mod 5 from n = 5, 6, 7, else 1, x_n = a_(n+1) + u_n; {} terms to 2^{}: {:?} ...; sorted {}; windows a_(n+2) > S_n at n >= 2: {} at {:?} ..., exactly n = 0 mod 5: {}; least (a_(n+2) - S_n)/a_(n+2) over them {}/1000; chained pairs {}; windows holding a non-sum {} of {}; least 12 (S_n - a_(n+1)) - a_(n+1) over n >= 3: {}; largest a_(n+1) - 3 a_n: {}; a_(n+3) - a_(n+2) - a_(n+1) at n = 3 mod 5, n >= 3, takes {:?}; chain x_n, 3 <= n, x_n below a_(n+2): {} of {} non-sums, first {:?}; non-sums below {}: {}, by lane seed/window/pair: {}/{}/{}",
        len, bits_, &a[..14], sorted, open.len(), &open[..open.len().min(6)], fifth, size, chains, held, open.len(), least, steep, kicks, gone, chain.len(), &chain[..8], top, lanes.iter().sum::<usize>(), lanes[0], lanes[1], lanes[2]
    );
}

// DEPTH

fn base_of(bases: &[u64], x: u128) -> u64 {
    for &d in bases {
        let mut y = x;
        while y.is_multiple_of(d as u128) {
            y /= d as u128;
        }
        if y == 1 {
            return d;
        }
    }
    0
}

fn trace(bases: &[u64], k: u32, bits: u32) {
    let a = elements(bases, k, far());
    let mut rows = vec![];
    let (v, _) = wide_with(bases, k, bits, |n, s, t| rows.push((n, s, t)));
    println!("{} k {} wide 2^{} {:?}", show(bases), k, bits, v);
    for (i, &(n, s, t)) in rows.iter().enumerate() {
        let next = a[n];
        let rise = rows.get(i + 1).is_some_and(|r| r.2 > t);
        println!(
            "n {} a_n {} = {}^. S_n {} a_(n+1) {} R_n {} T_n {} T_n/a_(n+1) {}.{:04}{}",
            n,
            a[n - 1],
            base_of(bases, a[n - 1]),
            s,
            next,
            s as i128 - next as i128,
            t,
            t / next,
            t % next * 10000 / next,
            if rise { " rises" } else { "" }
        );
    }
}

#[derive(Default, Clone)]
struct Jumps {
    rises: usize,
    pauses: usize,
    pauses_at_last: usize,
    last: usize,
    run: usize,
    open: usize,
    bad: usize,
    birth: usize,
    y: i128,
    deep: bool,
    fs: (u128, u128),
    peak: (u128, u128, usize),
}

fn frac(p: (u128, u128), places: u32) -> String {
    let m = 10u128.pow(places);
    format!(
        "{}.{:0w$}",
        p.0 / p.1,
        p.0 % p.1 * m / p.1,
        w = places as usize
    )
}

fn jumps_of(
    bases: &[u64],
    k: u32,
    a: &[u128],
    rows: &[(usize, u128, u128)],
    f: Option<u128>,
) -> Jumps {
    let mut j = Jumps {
        peak: (0, 1, 0),
        fs: (1, 1),
        ..Default::default()
    };
    let first = bases.iter().map(|&d| (d as u128).pow(k)).max().unwrap();
    let all = a.partition_point(|&x| x < first) + 1;
    let mut run = 0;
    for i in 0..rows.len() {
        let (n, s, t) = rows[i];
        let next = a[n];
        let d = s as i128 - next as i128;
        let ti = t as i128;
        if n >= all && t * j.peak.1 > j.peak.0 * next {
            j.peak = (t, next, n);
        }
        if i > 0 && t > rows[i - 1].2 && d > 2 * ti - 2 {
            j.bad += 1;
        }
        if let Some(f) = f {
            let s1 = if i + 1 < rows.len() {
                rows[i + 1].1
            } else {
                s + next
            };
            if s < 2 * f && 2 * f <= s1 {
                j.birth = n;
                j.y = f as i128 - next as i128;
                j.deep = ti - 1 > d;
                j.fs = (f, s);
            }
        }
        if i + 1 == rows.len() {
            break;
        }
        let (_, s1, t1) = rows[i + 1];
        let rise = t1 > t;
        if n >= all && !rise {
            j.pauses += 1;
        }
        if rise {
            j.last = n;
            j.pauses_at_last = j.pauses;
            j.rises += 1;
            run += 1;
            j.run = j.run.max(run);
            if t1 + t < s + 2 || t1 > s1 / 2 + 1 || d > 2 * ti - 2 {
                j.bad += 1;
            }
        } else {
            run = 0;
        }
        if d < ti - 1 && 2 * (ti - 1) < s as i128 && !rise {
            j.bad += 1;
        }
    }
    j.open = run;
    j
}

type Cell = (usize, u32, Verdict, Lanes, Jumps);

fn jumps(r: u64, kmax: u32, bits: u32, threads: usize) {
    let clock = Instant::now();
    let sets = minimal(r);
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Cell>> = Mutex::new(vec![]);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= sets.len() {
                    break;
                }
                for k in 1..=kmax {
                    let a = elements(&sets[i], k, far());
                    let mut rows = vec![];
                    let (v, lanes) = wide_with(&sets[i], k, bits, |n, s, t| rows.push((n, s, t)));
                    let f = match v {
                        Verdict::Found { f, .. } => Some(f),
                        Verdict::Cap => None,
                    };
                    let j = jumps_of(&sets[i], k, &a, &rows, f);
                    let stop = v == Verdict::Cap;
                    out.lock().unwrap().push((i, k, v, lanes, j));
                    if stop {
                        break;
                    }
                }
            });
        }
    });
    let mut out = out.into_inner().unwrap();
    out.sort_by_key(|x| (x.0, x.1));
    let mut bad = 0;
    let mut worst: Option<(u128, u128, String)> = None;
    let mut low: Option<(u128, u128, String)> = None;
    let mut gaps = vec![0usize; 64];
    let mut born = 0;
    let mut caps = (0usize, usize::MAX, 0usize, 0usize);
    let mut tally = vec![[0usize; 7]; kmax as usize + 1];
    for (i, k, v, lanes, j) in &out {
        let name = format!("{} k {}", show(&sets[*i]), k);
        bad += j.bad;
        match v {
            Verdict::Found { f, n, .. } => {
                if worst
                    .as_ref()
                    .is_none_or(|w| j.peak.0 * w.1 > w.0 * j.peak.1)
                {
                    worst = Some((j.peak.0, j.peak.1, name.clone()));
                }
                if *k >= 2 && low.as_ref().is_none_or(|w| j.peak.0 * w.1 < w.0 * j.peak.1) {
                    low = Some((j.peak.0, j.peak.1, name.clone()));
                }
                gaps[(n - j.last).min(63)] += 1;
                born += (j.birth == j.last) as usize;
                let row = &mut tally[*k as usize];
                row[0] += 1;
                row[1] += (j.y >= 0) as usize;
                row[2] += (lanes.last == 2) as usize;
                row[3] += (lanes.last == 1 && j.deep) as usize;
                row[5] += j.deep as usize;
                row[6] += (j.pauses_at_last == 0) as usize;
                row[4] +=
                    (j.y >= 0 && lanes.last != 2) as usize + (lanes.last == 1 && !j.deep) as usize;
                println!(
                    "{} F {} N {} rises {} longest run {}, steps without a rise from the entry of the last base to the last rise at m {}: {}; F born at b {}, F - a_(b+1) {}, T_b - 1 {} R_b, F/S_b {}, F by {}; peak T_n/a_(n+1) {} at n {}; bad {}",
                    name, f, n, j.rises, j.run, j.last, j.pauses_at_last, j.birth, j.y, if j.deep { ">" } else { "<=" }, frac(j.fs, 4), ["seed", "window", "pair"][lanes.last as usize], frac((j.peak.0, j.peak.1), 4), j.peak.2, j.bad
                );
            }
            Verdict::Cap => {
                caps.0 += 1;
                caps.1 = caps.1.min(j.open);
                caps.2 = caps.2.max(j.open);
                caps.3 += (j.open == j.rises) as usize;
                println!(
                "{} cap 2^{}: rises {} longest run {}, run open at the cap {}; peak T_n/a_(n+1) {} at n {}; bad {}",
                name, bits, j.rises, j.run, j.open, frac((j.peak.0, j.peak.1), 4), j.peak.2, j.bad
                )
            }
        }
        bad += tally[*k as usize][4];
        tally[*k as usize][4] = 0;
    }
    for (k, row) in tally.iter().enumerate().skip(1) {
        println!(
            "k = {}: {} certified cells; T rises at every step from the entry of the last base to the last rise in {}; T_b - 1 > R_b at the birth b of F in {}; F >= a_(b+1) in {}; F by the pair route in {}; F by the window route, each with T_b - 1 > R_b, in {}",
            k, row[0], row[6], row[5], row[1], row[2], row[3]
        );
    }
    if let Some(w) = worst {
        println!(
            "largest peak T_n/a_(n+1) after every base enters: {} at {}",
            frac((w.0, w.1), 4),
            w.2
        );
    }
    if let Some(w) = low {
        println!(
            "least peak T_n/a_(n+1) after every base enters over k >= 2: {} at {}",
            frac((w.0, w.1), 4),
            w.2
        );
    }
    let close: Vec<String> = gaps
        .iter()
        .enumerate()
        .filter(|x| *x.1 > 0)
        .map(|(g, c)| format!("{} steps in {}", g, c))
        .collect();
    println!("the certificate closes after the last rise of T: {}; F is born at the last rise in {} cells", close.join(", "), born);
    println!(
        "capped cells {}: every rise in one run open at the cap in {}, open runs of {} to {} rises",
        caps.0, caps.3, caps.1, caps.2
    );
    println!(
        "cells run {}; violations of the jump law and the route lemma {}; {:.1} s",
        out.len(),
        bad,
        clock.elapsed().as_secs_f64()
    );
}

// TERNARY

const WORDS: usize = 128;

#[derive(Clone, Copy)]
struct Sums {
    w: [u64; WORDS],
    top: usize,
}

impl Sums {
    fn zero() -> Sums {
        let mut w = [0u64; WORDS];
        w[0] = 1;
        Sums { w, top: 0 }
    }

    fn meets(&self, a: usize) -> bool {
        let q = a / 64;
        let r = (a % 64) as u32;
        for i in q..=(self.top + a) / 64 {
            let lo = if r > 0 && i > q {
                self.w[i - q - 1] >> (64 - r)
            } else {
                0
            };
            let v = if r == 0 {
                self.w[i - q]
            } else {
                (self.w[i - q] << r) | lo
            };
            if v & self.w[i] != 0 {
                return true;
            }
        }
        false
    }

    fn with(&self, a: usize) -> Sums {
        let mut out = *self;
        for shift in [a, 2 * a] {
            let q = shift / 64;
            let r = (shift % 64) as u32;
            for i in q..=(self.top + shift) / 64 {
                let lo = if r > 0 && i > q {
                    self.w[i - q - 1] >> (64 - r)
                } else {
                    0
                };
                let v = if r == 0 {
                    self.w[i - q]
                } else {
                    (self.w[i - q] << r) | lo
                };
                out.w[i] |= v;
            }
        }
        out.top = self.top + 2 * a;
        out
    }
}

struct Hunt<'a> {
    floor: &'a [usize],
    near: usize,
    need: u128,
    nodes: u64,
    found: Option<Vec<usize>>,
}

impl Hunt<'_> {
    fn dfs(&mut self, sums: &Sums, chosen: &mut Vec<usize>, sq: u128, left: usize) {
        self.nodes += 1;
        if self.found.is_some() {
            return;
        }
        if left == 0 {
            self.found = Some(chosen.clone());
            return;
        }
        let top = *chosen.last().unwrap();
        let mut lo = self.floor[left - 1] + 1;
        if left >= 2 {
            lo = lo.max(self.near);
        }
        for a in (lo..top).rev() {
            let rest = (left - 1) as u128 * ((a - 1) as u128).pow(2);
            if sq + (a as u128).pow(2) + rest < self.need {
                break;
            }
            if sums.meets(a) || sums.meets(2 * a) {
                continue;
            }
            chosen.push(a);
            let next = sums.with(a);
            self.dfs(&next, chosen, sq + (a as u128).pow(2), left - 1);
            chosen.pop();
            if self.found.is_some() {
                return;
            }
        }
    }
}

fn hunt(n: usize, top: usize, floor: &[usize], threads: usize) -> (Option<Vec<usize>>, u64) {
    near(n, top, floor, threads, 0)
}

fn near(
    n: usize,
    top: usize,
    floor: &[usize],
    threads: usize,
    band: usize,
) -> (Option<Vec<usize>>, u64) {
    let need = (9u128.pow(n as u32) - 1).div_ceil(8);
    let root = Sums::zero().with(top);
    let lo = floor[n - 2] + 1;
    let seconds: Vec<usize> = (lo..top).rev().collect();
    let next = AtomicUsize::new(0);
    let out: Mutex<(Option<Vec<usize>>, u64)> = Mutex::new((None, 0));
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= seconds.len() || out.lock().unwrap().0.is_some() {
                    break;
                }
                let a = seconds[i];
                let sq = (top as u128).pow(2) + (a as u128).pow(2);
                if sq + (n as u128 - 2) * ((a - 1) as u128).pow(2) < need
                    || root.meets(a)
                    || root.meets(2 * a)
                {
                    continue;
                }
                if band > 0 && a + band < top {
                    continue;
                }
                let mut h = Hunt {
                    floor,
                    near: if band > 0 { top - band } else { 0 },
                    need,
                    nodes: 0,
                    found: None,
                };
                let mut chosen = vec![top, a];
                h.dfs(&root.with(a), &mut chosen, sq, n - 2);
                let mut o = out.lock().unwrap();
                o.1 += h.nodes;
                if h.found.is_some() && o.0.is_none() {
                    o.0 = h.found;
                }
            });
        }
    });
    out.into_inner().unwrap()
}

fn admissible(a: &[usize]) -> bool {
    let n = a.len();
    let mut seen = std::collections::HashSet::new();
    for c in 0..3usize.pow(n as u32) {
        let (mut x, mut s) = (c, 0);
        for &v in a {
            s += (x % 3) * v;
            x /= 3;
        }
        if !seen.insert(s) {
            return false;
        }
    }
    true
}

fn probe(n: usize, top: usize, threads: usize) {
    let floor = [0usize, 1, 3, 8, 22, 60, 168];
    let clock = Instant::now();
    let (found, nodes) = hunt(n, top, &floor[..n], threads);
    println!(
        "n {} top {}: {:?}, {} nodes, {:.1} s",
        n,
        top,
        found,
        nodes,
        clock.elapsed().as_secs_f64()
    );
}

fn band(lo: usize, hi: usize, width: usize, threads: usize) {
    let floor = [0usize, 1, 3, 8, 22, 60, 168];
    let clock = Instant::now();
    for top in lo..=hi {
        let (found, nodes) = near(7, top, &floor, threads, width);
        if let Some(mut set) = found {
            set.reverse();
            assert!(admissible(&set));
            println!(
                "band {}: first top {} with six terms within {} of it: {:?}, {} nodes, {:.1} s",
                width,
                top,
                width,
                set,
                nodes,
                clock.elapsed().as_secs_f64()
            );
            return;
        }
    }
    println!(
        "band {}: none with top in [{}, {}], {:.1} s",
        width,
        lo,
        hi,
        clock.elapsed().as_secs_f64()
    );
}

fn offsets(lo: usize, hi: usize) {
    let families: [[usize; 6]; 3] = [
        [0, 1, 3, 8, 22, 60],
        [0, 2, 6, 9, 23, 61],
        [0, 2, 5, 7, 21, 60],
    ];
    for fam in families {
        let mut best: Option<Vec<usize>> = None;
        'outer: for top in lo..=hi {
            for x in fam[5] + 1..top {
                let mut set: Vec<usize> = fam.iter().map(|o| top - o).collect();
                set.push(top - x);
                set.sort();
                if admissible(&set) {
                    best = Some(set);
                    break 'outer;
                }
            }
        }
        println!("offsets {:?} plus one: {:?}", fam, best);
    }
}

fn ternary(nmax: usize, hi: usize, threads: usize) {
    let mut floor = vec![0usize, 1];
    println!("g_3(1) = 1");
    for n in 2..=nmax {
        let clock = Instant::now();
        let mut total = 0;
        let mut hit = None;
        for top in floor[n - 1] + 1..=hi {
            let (found, nodes) = hunt(n, top, &floor, threads);
            total += nodes;
            if let Some(set) = found {
                hit = Some((top, set));
                break;
            }
        }
        match hit {
            Some((top, mut set)) => {
                set.reverse();
                assert!(admissible(&set));
                println!(
                    "g_3({}) = {}, witness {:?}, {} nodes, {:.1} s",
                    n,
                    top,
                    set,
                    total,
                    clock.elapsed().as_secs_f64()
                );
                floor.push(top);
            }
            None => {
                println!(
                    "g_3({}) > {}, {} nodes, {:.1} s",
                    n,
                    hi,
                    total,
                    clock.elapsed().as_secs_f64()
                );
                return;
            }
        }
    }
}

// TIGHT

struct Wide {
    w: usize,
    limit: usize,
    b: Bits,
    sparse: Vec<u128>,
    s: u128,
    seen: usize,
    top: Option<usize>,
}

impl Wide {
    fn new(bits: u32, holes: u32) -> Wide {
        let w = 1usize << bits;
        Wide {
            w,
            limit: 1usize << holes,
            b: Bits::new(w + 64),
            sparse: vec![],
            s: 0,
            seen: 0,
            top: None,
        }
    }

    fn push(&mut self, x: u128) -> Option<u128> {
        let w = self.w;
        let wu = w as u128;
        let limit = self.limit;
        let s = self.s;
        let snew = s + x;
        if snew >> 63 != 0 {
            return None;
        }
        if snew / 2 > wu {
            let half = snew / 2;
            if x > s && half > s.max(wu) && half - s.max(wu) > limit as u128 {
                return None;
            }
            let b = &self.b;
            let sparse = &self.sparse;
            let mem = |z: i128| -> bool {
                if z < 0 || z as u128 > s {
                    return false;
                }
                let z = z as u128;
                if z <= wu {
                    b.get(z as usize)
                } else if z >= s - wu {
                    b.get((s - z) as usize)
                } else if 2 * z <= s {
                    sparse.binary_search(&z).is_err()
                } else {
                    sparse.binary_search(&(s - z)).is_err()
                }
            };
            let mut next: Vec<u128> = vec![];
            let mut keep = |z: u128| -> bool {
                if z > wu && z <= half && !mem(z as i128 - x as i128) {
                    next.push(z);
                }
                next.len() > limit
            };
            let mut spill = false;
            for &g in sparse {
                spill |= keep(g) || keep(s - g);
            }
            if s > wu {
                let d = s as i128 - x as i128;
                let glo = if d > 0 { ((d + 1) / 2) as usize } else { 0 };
                let ghi = w.min((s - wu - 1) as usize);
                spill |= b.each_zero(glo, ghi, |g| keep(s - g as u128));
            }
            if x > s {
                for z in s.max(wu) + 1..=half {
                    spill |= keep(z);
                }
            }
            if spill {
                return None;
            }
            next.sort();
            next.dedup();
            self.sparse = next;
        }
        if x <= wu {
            self.b.fold(x as usize, w.min(snew as usize));
        }
        self.s = snew;
        let h = w.min((snew / 2) as usize);
        let fresh = if h >= self.seen {
            self.b.hole_in(self.seen, h)
        } else {
            None
        };
        let b = &self.b;
        self.top = fresh.or_else(|| self.top.and_then(|t| b.hole_in(0, t)));
        self.seen = self.seen.max(h + 1);
        Some(match self.sparse.last() {
            Some(&z) => z + 1,
            None => self.top.map_or(0, |z| z as u128 + 1),
        })
    }

    fn gaps(&self, t: u128) -> usize {
        let mut count = 0;
        self.b.each_zero(0, self.w.min((t - 1) as usize), |_| {
            count += 1;
            false
        });
        count + self.sparse.len()
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Tight {
    f: u128,
    gaps: usize,
    n: usize,
    last: u128,
    sum: u128,
    l: u128,
    big: u32,
    checked: usize,
    slack: i128,
    margin: f64,
    tries: usize,
}

fn ellison(l: u128) -> (u32, f64) {
    let rate = 4f64.ln() - 0.2;
    let need = ((3 * l) as f64).ln();
    let big = 8u32.max((need / rate).ceil() as u32 + 1);
    (big, (big + 1) as f64 * rate - need)
}

fn tight_cell(k: u32, bits: u32, holes: u32) -> Result<Tight, (usize, u128, usize)> {
    let c = (3 * 3u128.pow(k) + 2 * 4u128.pow(k) + 7u128.pow(k)) / 6;
    let a = elements(&[3, 4, 7], k, far());
    let mut wide = Wide::new(bits, holes);
    let mut tries = 0;
    for n in 0..a.len() - 1 {
        let Some(t) = wide.push(a[n]) else {
            return Err((n, wide.s, wide.sparse.len()));
        };
        if t == 0 || t > a[n + 1] {
            continue;
        }
        tries += 1;
        let l = c + 2 * t - 1;
        let (big, margin) = ellison(l);
        let h = 4u128.pow(big);
        assert!(h < far());
        let mut sm = wide.s;
        let mut slack = i128::MAX;
        let mut checked = 0;
        let mut holds = true;
        for &x in &a[n + 1..] {
            if x > h {
                break;
            }
            let room = sm as i128 + 1 - 2 * t as i128 - x as i128;
            if room < 0 {
                holds = false;
                break;
            }
            slack = slack.min(room);
            checked += 1;
            sm += x;
        }
        if holds {
            return Ok(Tight {
                f: t - 1,
                gaps: wide.gaps(t),
                n: n + 1,
                last: a[n],
                sum: wide.s,
                l,
                big,
                checked,
                slack,
                margin,
                tries,
            });
        }
    }
    Err((a.len(), wide.s, wide.sparse.len()))
}

fn tight(kmax: u32, bits: u32, holes: u32) {
    let clock = Instant::now();
    for k in 1..=kmax {
        let t = Instant::now();
        let c = (3 * 3u128.pow(k) + 2 * 4u128.pow(k) + 7u128.pow(k)) / 6;
        match tight_cell(k, bits, holes) {
            Ok(r) => {
                println!(
                    "{{3,4,7}} k {} C_k {} F {} gaps {} n {} last {} S_n {} L {} B {} 4^B {} terms checked {} least slack {} Ellison margin {:.3} tries {} {:.2}s",
                    k,
                    c,
                    r.f,
                    r.gaps,
                    r.n,
                    r.last,
                    r.sum,
                    r.l,
                    r.big,
                    4u128.pow(r.big),
                    r.checked,
                    r.slack,
                    r.margin,
                    r.tries,
                    t.elapsed().as_secs_f64()
                );
            }
            Err((n, s, len)) => {
                println!(
                    "{{3,4,7}} k {} cap W = 2^{}, 2^{} holes: folding term {} with S_n {} and {} middle holes {:.2}s",
                    k,
                    bits,
                    holes,
                    n + 1,
                    s,
                    len,
                    t.elapsed().as_secs_f64()
                );
                break;
            }
        }
    }
    println!("{:.1} s", clock.elapsed().as_secs_f64());
}

// PAIRS

const EXCEPT: [u32; 16] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 19, 27];

fn root(d: u64) -> u64 {
    (2..=d)
        .find(|&r| {
            let mut x = r;
            while x < d {
                x *= r;
            }
            x == d
        })
        .unwrap()
}

fn smooth(d: u64) -> bool {
    let mut x = d;
    while x.is_multiple_of(2) {
        x /= 2;
    }
    while x.is_multiple_of(3) {
        x /= 3;
    }
    x == 1
}

fn ones(r: u64, size: usize) -> Vec<Vec<u64>> {
    let den = (3..=r).fold(1u128, |l, d| {
        let e = (d - 1) as u128;
        l / gcd(l, e) * e
    });
    fn walk(
        d0: u64,
        r: u64,
        left: usize,
        den: u128,
        sum: u128,
        cur: &mut Vec<u64>,
        out: &mut Vec<Vec<u64>>,
    ) {
        for d in d0..=r {
            let w = den / (d - 1) as u128;
            if sum + w > den {
                continue;
            }
            if sum + left as u128 * w < den {
                break;
            }
            cur.push(d);
            if sum + w == den {
                out.push(cur.clone());
            } else if left > 1 {
                walk(d + 1, r, left - 1, den, sum + w, cur, out);
            }
            cur.pop();
        }
    }
    let mut out = vec![];
    walk(3, r, size, den, 0, &mut vec![], &mut out);
    let mut out: Vec<Vec<u64>> = out
        .into_iter()
        .filter(|s| coprime(s))
        .filter(|s| {
            let roots: Vec<u64> = s.iter().map(|&d| root(d)).collect();
            (0..s.len()).all(|i| (0..i).all(|j| roots[i] != roots[j]))
        })
        .collect();
    order(&mut out);
    out
}

fn pair_of(bases: &[u64]) -> Option<(u64, u64)> {
    let s: Vec<u64> = bases.iter().copied().filter(|&d| smooth(d)).collect();
    (s.len() >= 2).then(|| (s[0], s[1]))
}

fn rho() -> (u128, u32, u32) {
    let mut best: (u128, u32, u32) = (1, 0, 0);
    for &x in &EXCEPT {
        let p = 1u128 << x;
        for y in 1..=20u32 {
            let c = p.abs_diff(3u128.pow(y));
            if best.1 == 0 || c << best.1 < best.0 << x {
                best = (c, x, y);
            }
        }
    }
    best
}

fn horizon(m: u64, l: u128) -> u128 {
    let (c, x, _) = rho();
    let need = (m as u128 - 1) * l;
    let h1 = (need << x).div_ceil(c);
    let theta = 1.0 - 1.0 / (10.0 * 2f64.ln());
    let h2 = ((need as f64).ln() / theta).exp() * (1.0 + 1e-9);
    h1.max(h2.ceil() as u128 + 1)
}

#[derive(Debug, Clone, PartialEq)]
struct Paired {
    e: u64,
    pf: u64,
    c: u128,
    f: u128,
    gaps: usize,
    n: usize,
    last: u128,
    sum: u128,
    l: u128,
    h: u128,
    checked: usize,
    slack: i128,
    tries: usize,
}

enum Stop {
    Cap(usize, u128, usize),
    Far(usize, u128),
}

#[derive(Clone, Copy, PartialEq)]
enum Rule {
    Smooth,
    Logs,
}

fn paired_cell(bases: &[u64], k: u32, bits: u32, holes: u32, rule: Rule) -> Result<Paired, Stop> {
    let near = match rule {
        Rule::Smooth => vec![],
        Rule::Logs => nears(bases),
    };
    let reach = |l: u128| -> (u128, u64, u64) {
        match rule {
            Rule::Smooth => {
                let (e, f) = pair_of(bases).unwrap();
                (horizon(f, l), e, f)
            }
            Rule::Logs => log_horizon(&near, l),
        }
    };
    let c = 1 + bases
        .iter()
        .map(|&d| ((d as u128).pow(k) - 1) / (d as u128 - 1))
        .sum::<u128>();
    let a = elements(bases, k, far());
    let mut wide = Wide::new(bits, holes);
    let mut tries = 0;
    for n in 0..a.len() - 1 {
        let Some(t) = wide.push(a[n]) else {
            return Err(Stop::Cap(n, wide.s, wide.sparse.len()));
        };
        if t == 0 || t > a[n + 1] {
            continue;
        }
        tries += 1;
        let l = c + 2 * t - 1;
        let (h, e, pf) = reach(l);
        if h >= far() {
            return Err(Stop::Far(n, h));
        }
        let mut sm = wide.s;
        let mut slack = i128::MAX;
        let mut checked = 0;
        let mut holds = true;
        for &x in &a[n + 1..] {
            if x > h {
                break;
            }
            let room = sm as i128 + 1 - 2 * t as i128 - x as i128;
            if room < 0 {
                holds = false;
                break;
            }
            slack = slack.min(room);
            checked += 1;
            sm += x;
        }
        if holds {
            return Ok(Paired {
                e,
                pf,
                c,
                f: t - 1,
                gaps: wide.gaps(t),
                n: n + 1,
                last: a[n],
                sum: wide.s,
                l,
                h,
                checked,
                slack,
                tries,
            });
        }
    }
    Err(Stop::Far(a.len(), far()))
}

fn paired_row(
    bases: &[u64],
    k: u32,
    bits: u32,
    holes: u32,
    rule: Rule,
) -> (String, Option<Paired>) {
    let t = Instant::now();
    match paired_cell(bases, k, bits, holes, rule) {
        Ok(r) => (
            format!(
                "{} k {} pair {{{},{}}} C_k {} F {} gaps {} n {} last {} S_n {} L {} H {} terms checked {} least slack {} tries {} {:.2}s",
                show(bases), k, r.e, r.pf, r.c, r.f, r.gaps, r.n, r.last, r.sum, r.l, r.h, r.checked, r.slack, r.tries,
                t.elapsed().as_secs_f64()
            ),
            Some(r),
        ),
        Err(Stop::Cap(n, s, len)) => (
            format!(
                "{} k {} cap W = 2^{}, 2^{} holes: folding term {} with S_n {} and {} middle holes {:.2}s",
                show(bases), k, bits, holes, n + 1, s, len, t.elapsed().as_secs_f64()
            ),
            None,
        ),
        Err(Stop::Far(n, h)) => (
            format!("{} k {} horizon {} at term {} passes 10^30 {:.2}s", show(bases), k, h, n + 1, t.elapsed().as_secs_f64()),
            None,
        ),
    }
}

fn unit_sets(r: u64, size: usize) {
    let all = ones(r, size);
    let (with, without): (Vec<_>, Vec<_>) = all.iter().partition(|s| pair_of(s).is_some());
    let (c, x, y) = rho();
    println!(
        "sigma = 1, gcd = 1, at most {} bases in [3, {}], no two powers of one integer: {} sets; {} hold two {{2,3}}-smooth bases, {} do not",
        size, r, all.len(), with.len(), without.len()
    );
    println!(
        "least |2^x - 3^y|/2^x over x in Ellison's S: {}/2^{} at (x, y) = ({}, {})",
        c, x, x, y
    );
    for s in &with {
        let (e, f) = pair_of(s).unwrap();
        println!("pair {{{},{}}} {}", e, f, show(s));
    }
    for s in &without {
        println!("no pair {}", show(s));
    }
}

fn chosen(r: u64, size: usize, rule: Rule, all: bool) -> Vec<Vec<u64>> {
    ones(r, size)
        .into_iter()
        .filter(|s| (all && rule == Rule::Logs) || pair_of(s).is_some() == (rule == Rule::Smooth))
        .collect()
}

fn pcontrol(r: u64, size: usize, kmax: u32, rule: Rule, all: bool) {
    let clock = Instant::now();
    let sets = chosen(r, size, rule, all);
    let (mut agree, mut differ) = (0, vec![]);
    for s in &sets {
        for k in 1..=kmax {
            let Ok(p) = paired_cell(s, k, 31, 25, rule) else {
                differ.push((s.clone(), k));
                continue;
            };
            let f = p.f as usize;
            let sums = knapsack(s, k, 4 * f + 4);
            let ok = !sums[f]
                && sums[f + 1..].iter().all(|&x| x)
                && (1..=f).filter(|&x| !sums[x]).count() == p.gaps;
            if ok {
                agree += 1;
            } else {
                differ.push((s.clone(), k));
            }
        }
    }
    println!(
        "{} sets at k = 1..{}: a plain knapsack to 4 F + 4 agrees on F and the count of non-sums in {} cells, differs in {:?}; {:.1} s",
        sets.len(),
        kmax,
        agree,
        differ,
        clock.elapsed().as_secs_f64()
    );
}

#[allow(clippy::too_many_arguments)]
fn pairs(
    r: u64,
    size: usize,
    kmax: u32,
    bits: u32,
    holes: u32,
    threads: usize,
    rule: Rule,
    all: bool,
) {
    let clock = Instant::now();
    let sets = chosen(r, size, rule, all);
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<(usize, Vec<(String, Option<Paired>)>)>> = Mutex::new(vec![]);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= sets.len() {
                    break;
                }
                let mut row = vec![];
                for k in 1..=kmax {
                    let (line, cell) = paired_row(&sets[i], k, bits, holes, rule);
                    let ok = cell.is_some();
                    row.push((line, cell));
                    if !ok {
                        break;
                    }
                }
                rows.lock().unwrap().push((i, row));
            });
        }
    });
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by_key(|r| r.0);
    let mut closed = vec![0usize; kmax as usize + 1];
    let mut top: Vec<Option<(u128, usize, u128, usize)>> = vec![None; kmax as usize + 1];
    for (i, row) in &rows {
        for (j, (line, cell)) in row.iter().enumerate() {
            println!("{}", line);
            if let Some(p) = cell {
                closed[j + 1] += 1;
                let t = top[j + 1].get_or_insert((p.f, *i, p.h, *i));
                if p.f > t.0 {
                    (t.0, t.1) = (p.f, *i);
                }
                if p.h > t.2 {
                    (t.2, t.3) = (p.h, *i);
                }
            }
        }
    }
    for k in 1..=kmax as usize {
        if let Some((f, i, h, j)) = top[k] {
            println!(
                "k {}: {} closed, largest F {} at {}, largest H {} at {}",
                k,
                closed[k],
                f,
                show(&sets[i]),
                h,
                show(&sets[j])
            );
        }
    }
    println!(
        "{} sets; closed at k = 1..{}: {:?}; {:.1} s",
        sets.len(),
        kmax,
        &closed[1..],
        clock.elapsed().as_secs_f64()
    );
}

// LOGS

const GAP: u32 = 80;

#[derive(Clone)]
struct Big(Vec<u64>);

impl Big {
    fn one() -> Big {
        Big(vec![1])
    }

    fn times(&self, d: u64) -> Big {
        let mut out = Vec::with_capacity(self.0.len() + 1);
        let mut carry = 0u128;
        for &x in &self.0 {
            let p = x as u128 * d as u128 + carry;
            out.push(p as u64);
            carry = p >> 64;
        }
        if carry > 0 {
            out.push(carry as u64);
        }
        Big(out)
    }

    fn rank(&self, o: &Big) -> std::cmp::Ordering {
        self.0
            .len()
            .cmp(&o.0.len())
            .then_with(|| self.0.iter().rev().cmp(o.0.iter().rev()))
    }

    fn small(&self) -> Option<u128> {
        match self.0.len() {
            1 => Some(self.0[0] as u128),
            2 => Some(self.0[0] as u128 | (self.0[1] as u128) << 64),
            _ => None,
        }
    }
}

fn gap(u: &Big, v: &Big) -> Option<u128> {
    let (hi, lo) = if u.rank(v).is_ge() { (u, v) } else { (v, u) };
    let mut out = Vec::with_capacity(hi.0.len());
    let mut borrow = false;
    for (i, &x) in hi.0.iter().enumerate() {
        let (d1, b1) = x.overflowing_sub(lo.0.get(i).copied().unwrap_or(0));
        let (d2, b2) = d1.overflowing_sub(borrow as u64);
        out.push(d2);
        borrow = b1 || b2;
    }
    while out.len() > 1 && out.last() == Some(&0) {
        out.pop();
    }
    Big(out).small()
}

struct Near {
    e: u64,
    f: u64,
    top: u32,
    rows: Vec<(u128, Option<u128>, u32, u32)>,
}

fn ceiling(e: u64, f: u64) -> u32 {
    let (le, lf) = ((e as f64).ln(), (f as f64).ln());
    (((GAP + 1) as f64 * 2f64.ln() + 2520.0 * le * lf) / le).ceil() as u32 + 1
}

fn closest(e: u64, f: u64) -> Near {
    let top = ceiling(e, f);
    let (mut u, mut v, mut b) = (Big::one(), Big::one(), 0u32);
    let mut rows = vec![];
    for a in 1..=top {
        u = u.times(e);
        let mut w = v.times(f);
        while w.rank(&u).is_le() {
            v = w;
            w = v.times(f);
            b += 1;
        }
        for (j, p) in [(b, &v), (b + 1, &w)] {
            if j == 0 {
                continue;
            }
            let g = gap(&u, p);
            assert_ne!(g, Some(0));
            if let Some(g) = g.filter(|&g| g < 1 << GAP) {
                let low = if u.rank(p).is_le() { &u } else { p };
                rows.push((g, low.small(), a, j));
            }
        }
    }
    Near { e, f, top, rows }
}

fn nears(bases: &[u64]) -> Vec<Near> {
    let mut out = vec![];
    for (i, &e) in bases.iter().enumerate() {
        for &f in &bases[i + 1..] {
            out.push(closest(e.min(f), e.max(f)));
        }
    }
    out
}

fn reach(p: &Near, g: u128) -> u128 {
    assert!(g < 1 << GAP);
    p.rows
        .iter()
        .filter(|r| r.0 < g)
        .map(|r| r.1.unwrap_or(u128::MAX))
        .max()
        .unwrap_or(0)
}

fn log_horizon(near: &[Near], l: u128) -> (u128, u64, u64) {
    near.iter()
        .map(|p| (reach(p, (p.f as u128 - 1) * l), p.e, p.f))
        .min()
        .unwrap()
}

fn near_pair(e: u64, f: u64, bits: u32) {
    let clock = Instant::now();
    let p = closest(e.min(f), e.max(f));
    let (e, f) = (p.e, p.f);
    let huge = p.rows.iter().filter(|r| r.1.is_none()).count();
    println!(
        "pair {{{},{}}}: Laurent ceiling 2520 ln {} ln {} = {:.3}, exponents a <= {} of {}, {} neighbour pairs with 0 < |{}^a - {}^b| < 2^{}, {} with min(...) >= 2^128",
        e, f, e, f, 2520.0 * (e as f64).ln() * (f as f64).ln(), p.top, e, p.rows.len(), e, f, GAP, huge
    );
    let mut rows: Vec<_> = p.rows.iter().filter(|r| r.0 < 1 << bits).collect();
    rows.sort();
    for r in &rows {
        println!(
            "a {} b {} gap {} min {}",
            r.2,
            r.3,
            r.0,
            r.1.map_or("huge".to_string(), |m| m.to_string())
        );
    }
    let top = p
        .rows
        .iter()
        .max_by_key(|r| r.1.unwrap_or(u128::MAX))
        .unwrap();
    let ratio = p
        .rows
        .iter()
        .max_by(|x, y| {
            (x.1.unwrap_or(u128::MAX) as f64 / x.0 as f64)
                .total_cmp(&(y.1.unwrap_or(u128::MAX) as f64 / y.0 as f64))
        })
        .unwrap();
    let m = ratio.1.unwrap_or(u128::MAX);
    let exact = p
        .rows
        .iter()
        .all(|r| r.1.unwrap_or(u128::MAX) <= m.div_ceil(ratio.0) * r.0);
    println!(
        "largest min over the {} pairs: {} at (a, b) = ({}, {}), gap {}; largest min/gap: {}/{} at (a, b) = ({}, {}), ceil {}, bounding every pair: {}; {:.2} s",
        p.rows.len(), top.1.unwrap_or(u128::MAX), top.2, top.3, top.0, m, ratio.0, ratio.2, ratio.3, m.div_ceil(ratio.0), exact, clock.elapsed().as_secs_f64()
    );
}

// SWEEP

fn word(w: &[u64], bit: usize) -> u64 {
    let q = bit >> 6;
    let r = bit & 63;
    if r == 0 {
        w[q]
    } else {
        (w[q] >> r) | (w[q + 1] << (64 - r))
    }
}

fn bit(w: &[u64], i: usize) -> bool {
    (w[i >> 6] >> (i & 63)) & 1 == 1
}

fn owed_span(r: usize, m: usize) -> usize {
    if r == 0 {
        0
    } else {
        2 * (r * (m - 1) - r * (r - 1) / 2)
    }
}

fn owed_sq(j: usize, a: usize) -> u128 {
    (1..=j).map(|i| ((a - i) as u128).pow(2)).sum()
}

fn root_ceil(v: u128) -> usize {
    let mut s = (v as f64).sqrt() as u128;
    while s * s < v {
        s += 1;
    }
    while s > 0 && (s - 1) * (s - 1) >= v {
        s -= 1;
    }
    s as usize
}

fn evens(x: u64) -> u64 {
    let mut x = x & 0x5555_5555_5555_5555;
    x = (x | (x >> 1)) & 0x3333_3333_3333_3333;
    x = (x | (x >> 2)) & 0x0f0f_0f0f_0f0f_0f0f;
    x = (x | (x >> 4)) & 0x00ff_00ff_00ff_00ff;
    x = (x | (x >> 8)) & 0x0000_ffff_0000_ffff;
    (x | (x >> 16)) & 0x0000_0000_ffff_ffff
}

fn widen(src: &[u64], dst: &mut [u64], off: usize, a: usize, b: usize) {
    for i in (off - b) >> 6..=(off + b) >> 6 {
        let g = i << 6;
        dst[i] = src[i]
            | word(src, g - a)
            | word(src, g + a)
            | word(src, g - 2 * a)
            | word(src, g + 2 * a);
    }
}

fn kernel(t: f64, h: f64) -> f64 {
    let u = t / h;
    h * (((-u * u / 30.0 + 2.0 / 3.0) * u - 4.0 / 3.0) * u * u + 16.0 / 15.0)
}

fn fejer(known: &[usize], r: usize, top: usize, atoms: f64) -> bool {
    let pi = std::f64::consts::PI;
    let two_h = 2.0 * pi / (3.0 * top as f64);
    let h = two_h / 2.0;
    let k = known.len();
    let all = k + r;
    let mut b = [0f64; 16];
    for (i, &a) in known.iter().enumerate() {
        b[i] = a as f64;
    }
    for j in 0..r {
        b[k + j] = (top - j) as f64;
    }
    let total: f64 = b[..all].iter().sum();
    let grid = ((((2.0 * two_h * total).ceil() as usize).max(64) + 1) / 2) * 2;
    let step = two_h / grid as f64;
    let (mut rc, mut rs, mut c, mut s) = ([0f64; 16], [0f64; 16], [1f64; 16], [0f64; 16]);
    let (mut pc, mut ps) = ([1f64; 16], [0f64; 16]);
    for i in 0..all {
        rc[i] = (step * b[i]).cos();
        rs[i] = (step * b[i]).sin();
    }
    let (mut simpson, mut tail) = (0.0, 0.0);
    for q in 0..=grid {
        let t = step * q as f64;
        if q > 0 {
            if q % 256 == 0 {
                for i in 0..all {
                    c[i] = (t * b[i]).cos();
                    s[i] = (t * b[i]).sin();
                }
            } else {
                for i in 0..all {
                    let (x, y) = (c[i], s[i]);
                    c[i] = x * rc[i] - y * rs[i];
                    s[i] = x * rs[i] + y * rc[i];
                }
            }
        }
        let mut kp = 1.0;
        let mut neg = false;
        let mut sup = 1.0;
        for i in 0..k {
            let f = (1.0 + 2.0 * c[i]) / 3.0;
            kp *= f;
            if q > 0 {
                let fp = (1.0 + 2.0 * pc[i]) / 3.0;
                neg |= f < 1e-12 || fp < 1e-12;
                sup *= if s[i] * ps[i] > 0.0 {
                    f.abs().max(fp.abs())
                } else {
                    1.0
                };
            }
        }
        let mut um = 1.0;
        for i in k..all {
            um *= (1.0 + 2.0 * c[i]) / 3.0;
        }
        let w = if q == 0 || q == grid {
            1.0
        } else if q % 2 == 1 {
            4.0
        } else {
            2.0
        };
        simpson += w * kernel(t, h) * kp * um;
        if neg {
            tail += step * kernel(t - step, h) * sup * (1.0 - um);
        }
        pc[..k].copy_from_slice(&c[..k]);
        ps[..k].copy_from_slice(&s[..k]);
    }
    simpson *= step / 3.0;
    let d4 = h * 1.07 * total.powi(4)
        + 4.0 * total.powi(3)
        + 6.0 * 2.67 / h * total.powi(2)
        + 16.0 / (h * h) * total
        + 8.0 / (h * h * h);
    let err = two_h * step.powi(4) * d4 / 180.0;
    2.0 * (simpson - err) - 2.0 * tail > 2.0 * pi * kernel(0.0, h) / atoms * (1.0 + 1e-9) + 1e-18
}

const MIRROR: usize = 4096;

struct Sweep<'a> {
    n: usize,
    floor: &'a [usize],
    need: u128,
    atoms: f64,
    off: usize,
    bufs: Vec<Vec<u64>>,
    even: Vec<u64>,
    odd: Vec<u64>,
    single: Vec<u64>,
    chosen: Vec<usize>,
    nodes: [u64; 8],
    tests: u64,
    cut: bool,
    counting: bool,
    sets: Vec<Vec<usize>>,
    found: Option<Vec<usize>>,
    collect: usize,
    out: Vec<Vec<usize>>,
    stop: &'a AtomicBool,
}

impl<'a> Sweep<'a> {
    fn new(
        n: usize,
        floor: &'a [usize],
        top: usize,
        cut: bool,
        counting: bool,
        stop: &'a AtomicBool,
    ) -> Sweep<'a> {
        let off = (owed_span(n - 1, top).max(2 * top) + 256) & !63;
        let words = (2 * off + 6 * top + 512) / 64;
        Sweep {
            n,
            floor,
            need: (9u128.pow(n as u32) - 1) / 8,
            atoms: 3f64.powi(n as i32),
            off,
            bufs: vec![vec![0u64; words]; n + 1],
            even: vec![0u64; words],
            odd: vec![0u64; words],
            single: vec![0u64; (MIRROR + 3 * top) / 64 + 4],
            chosen: vec![],
            nodes: [0; 8],
            tests: 0,
            cut,
            counting,
            sets: vec![],
            found: None,
            collect: 0,
            out: vec![],
            stop,
        }
    }

    fn root(&mut self, prefix: &[usize]) -> u128 {
        let off = self.off;
        let top = prefix[0];
        let d = &mut self.bufs[0];
        d.iter_mut().for_each(|w| *w = 0);
        for x in [off - 2 * top, off - top, off, off + top, off + 2 * top] {
            d[x >> 6] |= 1 << (x & 63);
        }
        self.chosen = vec![top];
        let mut sq = (top as u128).pow(2);
        for (i, &a) in prefix.iter().enumerate().skip(1) {
            let (lo, hi) = self.bufs.split_at_mut(i);
            widen(&lo[i - 1], &mut hi[0], off, a, owed_span(self.n - i - 1, a));
            self.chosen.push(a);
            sq += (a as u128).pow(2);
        }
        sq
    }

    fn done(&self) -> bool {
        self.found.is_some() || self.stop.load(Ordering::Relaxed)
    }

    fn lift(&mut self, r: usize, lo: usize, m: usize) -> usize {
        if !self.cut || r < 3 || lo >= m {
            return lo;
        }
        let (mut l, mut h) = (lo - 1, m);
        while h - l > 1 {
            let mid = (l + h + 1) / 2;
            self.tests += 1;
            if fejer(&self.chosen, r, mid, self.atoms) {
                l = mid;
            } else {
                h = mid;
            }
        }
        l + 1
    }

    fn dfs(&mut self, depth: usize, r: usize, sq: u128) {
        self.nodes[depth] += 1;
        if self.done() {
            return;
        }
        if self.collect > 0 && self.chosen.len() == self.collect {
            self.out.push(self.chosen.clone());
            return;
        }
        if r == 2 && self.collect == 0 {
            self.pairs(depth, sq);
            return;
        }
        let m = *self.chosen.last().unwrap();
        let off = self.off;
        let mut lo = self.floor[r].max(r);
        while lo < m && sq + (lo as u128).pow(2) + owed_sq(r - 1, lo) < self.need {
            lo += 1;
        }
        lo = self.lift(r, lo, m);
        for a in (lo..m).rev() {
            if bit(&self.bufs[depth], off + a) || bit(&self.bufs[depth], off + 2 * a) {
                continue;
            }
            if r == 1 {
                let mut set = self.chosen.clone();
                set.push(a);
                if self.counting {
                    self.sets.push(set);
                    continue;
                }
                self.found = Some(set);
                return;
            }
            {
                let (lo_, hi_) = self.bufs.split_at_mut(depth + 1);
                widen(&lo_[depth], &mut hi_[0], off, a, owed_span(r - 1, a));
            }
            self.chosen.push(a);
            self.dfs(depth + 1, r - 1, sq + (a as u128).pow(2));
            self.chosen.pop();
            if self.done() {
                return;
            }
        }
    }

    fn pairs(&mut self, depth: usize, sq: u128) {
        let m = *self.chosen.last().unwrap();
        let off = self.off;
        let h = off / 2;
        let b = owed_span(2, m);
        {
            let d = &self.bufs[depth];
            for j in (off - b) >> 7..=((off + b) >> 7) + 1 {
                self.even[j] = evens(d[2 * j]) | (evens(d[2 * j + 1]) << 32);
                self.odd[j] = evens(d[2 * j] >> 1) | (evens(d[2 * j + 1] >> 1) << 32);
            }
            for j in (MIRROR - m) >> 6..=(MIRROR + 2 * m) / 64 + 1 {
                self.single[j] =
                    !(word(d, off + 64 * j - MIRROR) | word(&self.even, h + 64 * j - MIRROR));
            }
        }
        let mut lo = self.floor[2].max(2);
        while lo < m && sq + (lo as u128).pow(2) + ((lo - 1) as u128).pow(2) < self.need {
            lo += 1;
        }
        let d = &self.bufs[depth];
        let (e, o, u) = (&self.even, &self.odd, &self.single);
        for x in (lo..m).rev() {
            if !bit(u, MIRROR + x) {
                continue;
            }
            self.nodes[depth + 1] += 1;
            let base = sq + (x as u128).pow(2);
            let ylo = if base >= self.need {
                1
            } else {
                root_ceil(self.need - base).max(1)
            };
            let half = if x % 2 == 0 { e } else { o };
            let (sp, sm) = if x % 2 == 0 {
                (x / 2, x / 2)
            } else {
                ((x - 1) / 2, (x + 1) / 2)
            };
            for chunk in ylo >> 6..=(x - 1) >> 6 {
                let y0 = chunk << 6;
                let bad = !(u[(MIRROR >> 6) + chunk]
                    & word(u, MIRROR + y0 - x)
                    & word(u, MIRROR + y0 + x))
                    | word(d, off + y0 + 2 * x)
                    | word(d, off + y0 - 2 * x)
                    | word(half, h + y0 + sp)
                    | word(half, h + y0 - sm);
                let from = ylo.max(y0) - y0;
                let to = (x - 1).min(y0 + 63) - y0;
                let mut good = !bad & (!0u64 >> (63 - to)) & (!0u64 << from);
                while good != 0 {
                    let y = y0 + good.trailing_zeros() as usize;
                    good &= good - 1;
                    let mut set = self.chosen.clone();
                    set.push(x);
                    set.push(y);
                    if self.counting {
                        self.sets.push(set);
                        continue;
                    }
                    self.found = Some(set);
                    return;
                }
            }
        }
    }
}

struct Tally {
    found: Option<Vec<usize>>,
    sets: Vec<Vec<usize>>,
    nodes: [u64; 8],
    tests: u64,
    tasks: usize,
}

fn sweep_top(
    n: usize,
    top: usize,
    floor: &[usize],
    threads: usize,
    cut: bool,
    counting: bool,
) -> Tally {
    let stop = AtomicBool::new(false);
    let t = if n >= 5 {
        3
    } else if n == 4 {
        2
    } else {
        1
    };
    let tasks: Vec<Vec<usize>> = if t == 1 {
        vec![vec![top]]
    } else {
        let mut s = Sweep::new(n, floor, top, cut, counting, &stop);
        let sq = s.root(&[top]);
        s.collect = t;
        s.dfs(0, n - 1, sq);
        s.out
    };
    let next = AtomicUsize::new(0);
    let out = Mutex::new(Tally {
        found: None,
        sets: vec![],
        nodes: [0; 8],
        tests: 0,
        tasks: tasks.len(),
    });
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut s = Sweep::new(n, floor, top, cut, counting, &stop);
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= tasks.len() || stop.load(Ordering::Relaxed) {
                        break;
                    }
                    let sq = s.root(&tasks[i]);
                    s.dfs(tasks[i].len() - 1, n - tasks[i].len(), sq);
                    if s.found.is_some() {
                        stop.store(true, Ordering::Relaxed);
                        let mut o = out.lock().unwrap();
                        if o.found.is_none() {
                            o.found = s.found.take();
                        }
                    }
                }
                let mut o = out.lock().unwrap();
                for k in 0..8 {
                    o.nodes[k] += s.nodes[k];
                }
                o.tests += s.tests;
                o.sets.append(&mut s.sets);
            });
        }
    });
    let mut t = out.into_inner().unwrap();
    for set in t.sets.iter_mut() {
        set.reverse();
    }
    t.sets.sort();
    if let Some(set) = t.found.as_mut() {
        set.reverse();
    }
    t
}

fn sweep(n: usize, lo: usize, hi: usize, threads: usize) {
    assert!((2..=7).contains(&n) && lo > [0, 1, 3, 8, 22, 60, 168][n - 1]);
    let floor = [0usize, 1, 3, 8, 22, 60, 168];
    let clock = Instant::now();
    for top in lo..=hi {
        let t0 = Instant::now();
        let t = sweep_top(n, top, &floor, threads, true, false);
        println!(
            "top {}: {}, {} tasks, nodes {:?}, {} cut tests, {:.2} s, {:.1} s in all",
            top,
            t.found
                .as_ref()
                .map_or("none".to_string(), |s| format!("{:?}", s)),
            t.tasks,
            &t.nodes[..n],
            t.tests,
            t0.elapsed().as_secs_f64(),
            clock.elapsed().as_secs_f64()
        );
        if let Some(set) = t.found {
            assert!(admissible(&set));
            println!(
                "sweep: least top in [{}, {}] of an admissible {}-set is {}, witness {:?}",
                lo, hi, n, top, set
            );
            return;
        }
    }
    println!(
        "sweep: no admissible {}-set with top in [{}, {}], {:.1} s",
        n,
        lo,
        hi,
        clock.elapsed().as_secs_f64()
    );
}

fn spared(set: &[usize], atoms: f64) -> bool {
    let mut down = set.to_vec();
    down.reverse();
    let n = down.len();
    (1..n - 1).all(|k| (down[k]..down[k - 1]).all(|a| !fejer(&down[..k], n - k, a, atoms)))
}

fn tally(n: usize, lo: usize, hi: usize, threads: usize, plain: bool) {
    assert!((2..=7).contains(&n));
    let floor = [0usize, 1, 3, 8, 22, 60, 168];
    let atoms = 3f64.powi(n as i32);
    let (mut total, mut fired, mut differ) = (0usize, 0usize, 0usize);
    let clock = Instant::now();
    for top in lo..=hi {
        let t = sweep_top(n, top, &floor, threads, true, true);
        for set in &t.sets {
            assert!(admissible(set));
            if !spared(set, atoms) {
                fired += 1;
            }
        }
        total += t.sets.len();
        let line = if plain {
            let p = sweep_top(n, top, &floor, threads, false, true);
            if p.sets != t.sets {
                differ += 1;
            }
            format!("{} with the cut, {} without", t.sets.len(), p.sets.len())
        } else {
            format!("{} with the cut", t.sets.len())
        };
        let shown: Vec<&Vec<usize>> = t.sets.iter().take(4).collect();
        println!("n {} top {}: {}; {:?}", n, top, line, shown);
    }
    println!(
        "tally: {} admissible {}-sets with top in [{}, {}], cut fired on a prefix of {}, lists differ at {} tops, {:.1} s",
        total, n, lo, hi, fired, differ, clock.elapsed().as_secs_f64()
    );
}

// ROUTE

fn route(r: u64, lo: u128) {
    for bases in minimal(r) {
        let a = elements(&bases, 1, far());
        let mut below: u128 = 0;
        let mut best = (u128::MAX, 0u128);
        for &x in &a {
            let key = below * 1_000_000 / x;
            if x >= lo && key < best.0 {
                best = (key, x);
            }
            below += x;
        }
        println!(
            "{} least sum(M below N)/N over N in M, {} <= N <= 10^30: {}.{:06} at N = {}",
            show(&bases),
            lo,
            best.0 / 1_000_000,
            best.0 % 1_000_000,
            best.1
        );
    }
}

// CANTOR

const DEAD: u32 = u32::MAX;

struct Carries {
    next: Vec<[u32; 2]>,
    capped: bool,
}

fn carries(mults: &[u128], cap: usize) -> Carries {
    let r = mults.len();
    assert!(r <= 4 && mults.iter().all(|&m| m < 1 << 126));
    let mut index: HashMap<[u128; 4], u32> = HashMap::new();
    let mut states: Vec<[u128; 4]> = vec![[0; 4]];
    index.insert([0; 4], 0);
    let mut next = vec![];
    while next.len() < states.len() {
        if states.len() > cap {
            return Carries { next, capped: true };
        }
        let at = states[next.len()];
        let mut out = [DEAD; 2];
        for (d, slot) in out.iter_mut().enumerate() {
            let mut to = [0u128; 4];
            let mut ok = true;
            for i in 0..r {
                let t = at[i] + mults[i] * d as u128;
                if t % 3 == 2 {
                    ok = false;
                    break;
                }
                to[i] = t / 3;
            }
            if ok {
                let fresh = states.len() as u32;
                let id = *index.entry(to).or_insert(fresh);
                if id == fresh {
                    states.push(to);
                }
                *slot = id;
            }
        }
        next.push(out);
    }
    Carries {
        next,
        capped: false,
    }
}

fn components(next: &[[u32; 2]]) -> (Vec<u32>, usize) {
    let n = next.len();
    let mut index = vec![u32::MAX; n];
    let mut low = vec![0u32; n];
    let mut on = vec![false; n];
    let mut comp = vec![u32::MAX; n];
    let mut stack: Vec<u32> = vec![];
    let mut call: Vec<(u32, u8)> = vec![];
    let mut counter = 0u32;
    let mut comps = 0usize;
    for root in 0..n {
        if index[root] != u32::MAX {
            continue;
        }
        index[root] = counter;
        low[root] = counter;
        counter += 1;
        stack.push(root as u32);
        on[root] = true;
        call.push((root as u32, 0));
        while let Some(top) = call.last_mut() {
            let v = top.0 as usize;
            if top.1 < 2 {
                let w = next[v][top.1 as usize];
                top.1 += 1;
                if w == DEAD {
                    continue;
                }
                let w = w as usize;
                if index[w] == u32::MAX {
                    index[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push(w as u32);
                    on[w] = true;
                    call.push((w as u32, 0));
                } else if on[w] {
                    low[v] = low[v].min(index[w]);
                }
            } else {
                call.pop();
                if let Some(parent) = call.last() {
                    let p = parent.0 as usize;
                    low[p] = low[p].min(low[v]);
                }
                if low[v] == index[v] {
                    loop {
                        let w = stack.pop().unwrap() as usize;
                        on[w] = false;
                        comp[w] = comps as u32;
                        if w == v {
                            break;
                        }
                    }
                    comps += 1;
                }
            }
        }
    }
    (comp, comps)
}

fn radius(next: &[[u32; 2]], comp: &[u32], members: &[u32]) -> (f64, f64) {
    let c = comp[members[0] as usize];
    let pos: HashMap<u32, usize> = members.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let succ: Vec<Vec<usize>> = members
        .iter()
        .map(|&v| {
            next[v as usize]
                .iter()
                .filter(|&&w| w != DEAD && comp[w as usize] == c)
                .map(|w| pos[w])
                .collect()
        })
        .collect();
    let mut v = vec![1.0f64; members.len()];
    let (mut lo, mut hi) = (0.0f64, f64::MAX);
    for _ in 0..100000 {
        let w: Vec<f64> = (0..v.len())
            .map(|i| v[i] + succ[i].iter().map(|&j| v[j]).sum::<f64>())
            .collect();
        lo = f64::MAX;
        hi = 0.0;
        for i in 0..v.len() {
            lo = lo.min(w[i] / v[i]);
            hi = hi.max(w[i] / v[i]);
        }
        let top = w.iter().cloned().fold(0.0, f64::max);
        v = w.iter().map(|x| x / top).collect();
        if hi - lo < 1e-13 {
            break;
        }
    }
    (lo - 1.0, hi - 1.0)
}

struct Shape {
    states: usize,
    capped: bool,
    loops: usize,
    rich: usize,
    lo: f64,
    hi: f64,
    size: usize,
    hub: Option<u32>,
    next: Vec<[u32; 2]>,
    comp: Vec<u32>,
}

fn shape(mults: &[u128], cap: usize) -> Shape {
    let Carries { next, capped } = carries(mults, cap);
    let mut out = Shape {
        states: next.len(),
        capped,
        loops: 0,
        rich: 0,
        lo: 1.0,
        hi: 1.0,
        size: 0,
        hub: None,
        next: vec![],
        comp: vec![],
    };
    if capped {
        return out;
    }
    let (comp, comps) = components(&next);
    let mut nodes = vec![0usize; comps];
    let mut edges = vec![0usize; comps];
    let mut members: Vec<Vec<u32>> = vec![vec![]; comps];
    for (v, out) in next.iter().enumerate() {
        let c = comp[v] as usize;
        nodes[c] += 1;
        members[c].push(v as u32);
        edges[c] += out
            .iter()
            .filter(|&&w| w != DEAD && comp[w as usize] as usize == c)
            .count();
    }
    for c in 0..comps {
        if edges[c] == 0 {
            continue;
        }
        out.loops += 1;
        if edges[c] > nodes[c] {
            out.rich += 1;
            let (lo, hi) = radius(&next, &comp, &members[c]);
            if lo > out.lo || out.hub.is_none() {
                out.lo = lo;
                out.hi = hi;
                out.size = nodes[c];
                out.hub = members[c].iter().copied().find(|&v| {
                    next[v as usize]
                        .iter()
                        .all(|&w| w != DEAD && comp[w as usize] as usize == c)
                });
            }
        }
    }
    out.next = next;
    out.comp = comp;
    out
}

fn walk(next: &[[u32; 2]], comp: Option<&[u32]>, from: u32, to: u32) -> Vec<u8> {
    let mut parent: HashMap<u32, (u32, u8)> = HashMap::new();
    let mut queue = std::collections::VecDeque::from([from]);
    parent.insert(from, (DEAD, 0));
    while let Some(v) = queue.pop_front() {
        if v == to {
            break;
        }
        for d in 0..2u8 {
            let w = next[v as usize][d as usize];
            if w == DEAD || parent.contains_key(&w) {
                continue;
            }
            if let Some(c) = comp {
                if c[w as usize] != c[from as usize] {
                    continue;
                }
            }
            parent.insert(w, (v, d));
            queue.push_back(w);
        }
    }
    let mut word = vec![];
    let mut at = to;
    while at != from {
        let (p, d) = parent[&at];
        word.push(d);
        at = p;
    }
    word.reverse();
    word
}

fn hub_cycles(s: &Shape) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let hub = s.hub?;
    let prefix = walk(&s.next, None, 0, hub);
    let cycle = |d: u8| {
        let mut w = vec![d];
        w.extend(walk(
            &s.next,
            Some(&s.comp),
            s.next[hub as usize][d as usize],
            hub,
        ));
        w
    };
    Some((prefix, cycle(0), cycle(1)))
}

fn run(m: u128, start: u128, word: &[u8]) -> Option<u128> {
    let mut n = start;
    for &d in word {
        let t = n + m * d as u128;
        if t % 3 == 2 {
            return None;
        }
        n = t / 3;
    }
    Some(n)
}

fn closes(mults: &[u128], prefix: &[u8], cycles: &[&[u8]]) -> bool {
    mults.iter().all(|&m| match run(m, 0, prefix) {
        Some(hub) => cycles.iter().all(|c| run(m, hub, c) == Some(hub)),
        None => false,
    })
}

fn moran(p: usize, q: usize) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..200 {
        let s = (lo + hi) / 2.0;
        if 3f64.powf(-(p as f64) * s) + 3f64.powf(-(q as f64) * s) > 1.0 {
            lo = s;
        } else {
            hi = s;
        }
    }
    lo
}

fn digits(w: &[u8]) -> String {
    w.iter().map(|d| char::from(b'0' + d)).collect()
}

fn ternary_of(mut m: u128) -> String {
    let mut s = vec![];
    while m > 0 {
        s.push(char::from(b'0' + (m % 3) as u8));
        m /= 3;
    }
    s.iter().rev().collect()
}

fn paths(next: &[[u32; 2]], len: usize) -> Vec<u128> {
    let mut at = vec![0u128; next.len()];
    at[0] = 1;
    let mut out = vec![];
    for _ in 0..len {
        let mut to = vec![0u128; next.len()];
        for (v, &c) in at.iter().enumerate() {
            if c > 0 {
                for &w in &next[v] {
                    if w != DEAD {
                        to[w as usize] += c;
                    }
                }
            }
        }
        out.push(to.iter().sum());
        at = to;
    }
    out
}

fn cantor(exps: &[u32], cap: usize) {
    let mults: Vec<u128> = exps.iter().map(|&e| 1u128 << e).collect();
    let name = exps
        .iter()
        .map(|e| format!("2^{}", e))
        .collect::<Vec<_>>()
        .join(", ");
    carry_set(&mults, &name, cap)
}

fn carry_set(mults: &[u128], name: &str, cap: usize) {
    let t = Instant::now();
    let s = shape(mults, cap);
    for &m in mults {
        println!("{} = ({})_3", m, ternary_of(m));
    }
    if s.capped {
        println!(
            "C(1, {}): capped at {} states, {:.2}s",
            name,
            cap,
            t.elapsed().as_secs_f64()
        );
        return;
    }
    println!(
        "C(1, {}): {} states, {} cyclic components, {} with more edges than nodes, radius in [{:.9}, {:.9}], dim {:.6}, component size {}, {:.2}s",
        name,
        s.states,
        s.loops,
        s.rich,
        s.lo,
        s.hi,
        s.lo.max(1.0).ln() / 3f64.ln(),
        s.size,
        t.elapsed().as_secs_f64()
    );
    if let Some((prefix, u, v)) = hub_cycles(&s) {
        println!(
            "witness: prefix {} then any word in cycles {} and {}, low digit first; direct carries close: {}; dim >= {:.6}",
            digits(&prefix),
            digits(&u),
            digits(&v),
            closes(mults, &prefix, &[&u, &v]),
            moran(u.len(), v.len())
        );
    }
    let counts = paths(&s.next, 28);
    println!("allowed prefixes of length 1..28: {:?}", counts);
}

fn subsets_of(pool: &[u32], r: usize) -> Vec<Vec<u32>> {
    if r == 0 {
        return vec![vec![]];
    }
    let mut out = vec![];
    for (i, &x) in pool.iter().enumerate() {
        for mut rest in subsets_of(&pool[i + 1..], r - 1) {
            rest.insert(0, x);
            out.push(rest);
        }
    }
    out
}

type Reading = (Vec<u32>, usize, bool, usize, f64, f64, usize);

fn tuples(r: usize, top: u32, threads: usize, cap: usize) {
    let clock = Instant::now();
    let pool: Vec<u32> = (1..=top / 2).map(|j| 2 * j).collect();
    let jobs = subsets_of(&pool, r);
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<Reading>> = Mutex::new(vec![]);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= jobs.len() {
                    break;
                }
                let mults: Vec<u128> = jobs[i].iter().map(|&e| 1u128 << e).collect();
                let s = shape(&mults, cap);
                rows.lock().unwrap().push((
                    jobs[i].clone(),
                    s.states,
                    s.capped,
                    s.rich,
                    s.lo,
                    s.hi,
                    s.loops,
                ));
            });
        }
    });
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let mut most = (0usize, vec![]);
    let (mut zero, mut capped, mut origin) = (0usize, 0usize, 0usize);
    for (e, states, cap_hit, rich, lo, hi, loops) in &rows {
        if *states > most.0 {
            most = (*states, e.clone());
        }
        if *cap_hit {
            capped += 1;
            println!("capped {:?}", e);
        } else if *rich > 0 {
            println!(
                "positive {:?}: {} states, radius in [{:.9}, {:.9}], dim {:.6}",
                e,
                states,
                lo,
                hi,
                lo.ln() / 3f64.ln()
            );
        } else {
            zero += 1;
            if *loops == 1 {
                origin += 1;
            } else {
                println!("countable {:?}: {} states, {} cycles", e, states, loops);
            }
        }
    }
    println!(
        "{} sets C(1, 2^m_1, ..., 2^m_{}) with even 2 <= m_1 < ... <= {}: {} of dimension zero, {} of them equal to {{0}}, {} positive, {} capped at {} states; most states {} at {:?}; {:.2}s",
        rows.len(),
        r,
        top,
        zero,
        origin,
        rows.len() - zero - capped,
        capped,
        cap,
        most.0,
        most.1,
        clock.elapsed().as_secs_f64()
    );
}

fn death(units: &[u128], dmax: u32) -> Option<u32> {
    let mut first = vec![];
    for &u in units {
        if u % 3 == 2 {
            return Some(1);
        }
        first.push(u / 3);
    }
    let mut stack = vec![(first, 1u32)];
    let mut deepest = 1;
    while let Some((carry, len)) = stack.pop() {
        deepest = deepest.max(len);
        if len == dmax {
            return None;
        }
        for d in 0..2u128 {
            let mut to = Vec::with_capacity(units.len());
            let mut ok = true;
            for (n, u) in carry.iter().zip(units) {
                let t = n + u * d;
                if t % 3 == 2 {
                    ok = false;
                    break;
                }
                to.push(t / 3);
            }
            if ok {
                stack.push((to, len + 1));
            }
        }
    }
    Some(deepest + 1)
}

fn residue_power(e: u32, modulus: u128) -> u128 {
    let mut x = 1u128 % modulus;
    for _ in 0..e {
        x = 2 * x % modulus;
    }
    x
}

struct Deaths {
    sets: usize,
    hist: Vec<usize>,
    open: Vec<Vec<u32>>,
    worst: (u32, Vec<u32>),
}

fn deaths(r: usize, top: u32, dmax: u32, threads: usize) {
    let clock = Instant::now();
    assert!(dmax <= 79 && r >= 1);
    let modulus = 3u128.pow(dmax);
    let pool: Vec<u32> = (1..=top / 2).map(|j| 2 * j).collect();
    let units: Vec<u128> = (0..=top).map(|e| residue_power(e, modulus)).collect();
    let next = AtomicUsize::new(0);
    let total: Mutex<Deaths> = Mutex::new(Deaths {
        sets: 0,
        hist: vec![0; dmax as usize + 2],
        open: vec![],
        worst: (0, vec![]),
    });
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= pool.len() {
                    break;
                }
                let mut mine = Deaths {
                    sets: 0,
                    hist: vec![0; dmax as usize + 2],
                    open: vec![],
                    worst: (0, vec![]),
                };
                for mut rest in subsets_of(&pool[i + 1..], r - 1) {
                    rest.insert(0, pool[i]);
                    let us: Vec<u128> = rest.iter().map(|&e| units[e as usize]).collect();
                    mine.sets += 1;
                    match death(&us, dmax) {
                        Some(d) => {
                            mine.hist[d as usize] += 1;
                            if d > mine.worst.0 {
                                mine.worst = (d, rest);
                            }
                        }
                        None => mine.open.push(rest),
                    }
                }
                let mut all = total.lock().unwrap();
                all.sets += mine.sets;
                for (x, y) in all.hist.iter_mut().zip(&mine.hist) {
                    *x += y;
                }
                all.open.extend(mine.open);
                if mine.worst.0 > all.worst.0
                    || (mine.worst.0 == all.worst.0 && mine.worst.1 < all.worst.1)
                {
                    all.worst = mine.worst;
                }
            });
        }
    });
    let mut all = total.into_inner().unwrap();
    all.open.sort();
    for e in all.open.iter().take(40) {
        println!("survives to length {}: {:?}", dmax, e);
    }
    let shown: Vec<String> = all
        .hist
        .iter()
        .enumerate()
        .filter(|x| *x.1 > 0)
        .map(|(d, c)| format!("{}:{}", d, c))
        .collect();
    println!(
        "{} sets of {} even exponents in [2, {}]: {} with every admissible prefix of lowest digit 1 dead by length {}, so C = {{0}}, {} surviving; deepest death {} at {:?}; death lengths {}; {:.2}s",
        all.sets,
        r,
        top,
        all.sets - all.open.len(),
        dmax,
        all.open.len(),
        all.worst.0,
        all.worst.1,
        shown.join(" "),
        clock.elapsed().as_secs_f64()
    );
}

fn dead_classes(r: usize, d: u32, threads: usize) -> (usize, usize) {
    let modulus = 3u128.pow(d);
    let period = 3usize.pow(d - 1);
    let mut units = vec![1u128; period];
    for i in 1..period {
        units[i] = units[i - 1] * 4 % modulus;
    }
    let next = AtomicUsize::new(0);
    let dead = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut local = 0usize;
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= period {
                        break;
                    }
                    let mut idx = vec![0usize; r];
                    idx[0] = i;
                    for k in 0..period.pow(r as u32 - 1) {
                        let mut rest = k;
                        for slot in idx.iter_mut().skip(1) {
                            *slot = rest % period;
                            rest /= period;
                        }
                        let us: Vec<u128> = idx.iter().map(|&j| units[j]).collect();
                        if death(&us, d).is_some() {
                            local += 1;
                        }
                    }
                }
                dead.fetch_add(local, Ordering::SeqCst);
            });
        }
    });
    (dead.into_inner(), period.pow(r as u32))
}

fn classes(r: usize, d: u32, threads: usize) {
    let clock = Instant::now();
    let (dead, all) = dead_classes(r, d, threads);
    println!(
        "{} of the {} ordered {}-tuples of exponent classes mod 2 3^{} die by length {}, {} survive against the first-moment bound 2^{}: lower density {:.6} of {}-tuples of even exponents with C = {{0}}, truncated; {:.2}s",
        dead,
        all,
        r,
        d - 1,
        d,
        all - dead,
        (r as u32 + 1) * (d - 1),
        (dead as f64 / all as f64 * 1e6).floor() / 1e6,
        r,
        clock.elapsed().as_secs_f64()
    );
}

// FOUR

fn on(bits: &[u64], x: usize) -> bool {
    bits[x >> 6] >> (x & 63) & 1 == 1
}

fn extend4(sums: &[usize], bits: &mut [u64], a: usize) -> Option<(Vec<usize>, Vec<usize>)> {
    let fresh: Vec<usize> = sums
        .iter()
        .map(|&x| x + a)
        .filter(|&t| !on(bits, t))
        .collect();
    for &t in &fresh {
        bits[t >> 6] |= 1 << (t & 63);
    }
    let mut all: Vec<usize> = sums.iter().chain(fresh.iter()).copied().collect();
    all.sort_unstable();
    let mut free = true;
    'out: for &t in &fresh {
        for &u in &all {
            if u > t {
                let diff = u - t;
                for k in 1..=3 {
                    if diff % k == 0 {
                        let d = diff / k;
                        let i = 3 - k;
                        if t >= i * d && (0..4).all(|l| on(bits, t - i * d + l * d)) {
                            free = false;
                            break 'out;
                        }
                    }
                }
            } else if u < t && (t - u) % 3 == 0 {
                let d = (t - u) / 3;
                if on(bits, u + d) && on(bits, u + 2 * d) {
                    free = false;
                    break 'out;
                }
            }
        }
    }
    if free {
        Some((all, fresh))
    } else {
        for &t in &fresh {
            bits[t >> 6] &= !(1 << (t & 63));
        }
        None
    }
}

fn plain_free(set: &[usize]) -> bool {
    let mut sums = vec![0usize];
    for &a in set {
        let more: Vec<usize> = sums.iter().map(|&x| x + a).collect();
        sums.extend(more);
    }
    sums.sort_unstable();
    sums.dedup();
    let seen: std::collections::HashSet<usize> = sums.iter().copied().collect();
    for (i, &x) in sums.iter().enumerate() {
        for &y in &sums[i + 1..] {
            let d = y - x;
            if seen.contains(&(y + d)) && seen.contains(&(y + 2 * d)) {
                return false;
            }
        }
    }
    true
}

fn descend(
    sums: &[usize],
    bits: &mut [u64],
    cands: &[usize],
    r: usize,
    floor: &[usize],
    chosen: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
    all: bool,
    nodes: &mut u64,
) {
    for (at, &b) in cands.iter().enumerate().rev() {
        if b < floor[r].max(1) {
            break;
        }
        *nodes += 1;
        if r == 1 {
            chosen.push(b);
            out.push(chosen.clone());
            chosen.pop();
            if !all {
                return;
            }
            continue;
        }
        let (next, fresh) = extend4(sums, bits, b).unwrap();
        let mut kept = vec![];
        for &y in &cands[..at] {
            if let Some((_, back)) = extend4(&next, bits, y) {
                for x in back {
                    bits[x >> 6] &= !(1 << (x & 63));
                }
                kept.push(y);
            }
        }
        let room = (1..r).all(|j| kept.iter().filter(|&&y| y >= floor[j].max(1)).count() >= r - j);
        if room {
            chosen.push(b);
            descend(&next, bits, &kept, r - 1, floor, chosen, out, all, nodes);
            chosen.pop();
        }
        for x in fresh {
            bits[x >> 6] &= !(1 << (x & 63));
        }
        if !all && !out.is_empty() {
            return;
        }
    }
}

fn four_top(n: usize, top: usize, floor: &[usize], all: bool) -> (Vec<Vec<usize>>, u64) {
    let mut bits = vec![0u64; (n * top) / 64 + 2];
    bits[0] |= 1;
    let sums = extend4(&[0], &mut bits, top).unwrap().0;
    let mut chosen = vec![top];
    let mut out = vec![];
    let mut nodes = 0;
    if n == 1 {
        out.push(chosen.clone());
    } else {
        let mut cands = vec![];
        for y in 1..top {
            if let Some((_, back)) = extend4(&sums, &mut bits, y) {
                for x in back {
                    bits[x >> 6] &= !(1 << (x & 63));
                }
                cands.push(y);
            }
        }
        descend(
            &sums,
            &mut bits,
            &cands,
            n - 1,
            floor,
            &mut chosen,
            &mut out,
            all,
            &mut nodes,
        );
    }
    for set in out.iter_mut() {
        set.reverse();
    }
    (out, nodes)
}

fn four(nmax: usize, hi: usize, threads: usize) {
    let mut floor = vec![0usize];
    for n in 1..=nmax {
        let clock = Instant::now();
        let lo = floor[n - 1] + 1;
        let nodes = AtomicUsize::new(0);
        let next = AtomicUsize::new(lo);
        let best = AtomicUsize::new(usize::MAX);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| loop {
                    let top = next.fetch_add(1, Ordering::SeqCst);
                    if top > hi || top >= best.load(Ordering::SeqCst) {
                        break;
                    }
                    let (found, k) = four_top(n, top, &floor, false);
                    nodes.fetch_add(k as usize, Ordering::SeqCst);
                    if !found.is_empty() {
                        best.fetch_min(top, Ordering::SeqCst);
                    }
                });
            }
        });
        let g = best.into_inner();
        if g == usize::MAX {
            println!(
                "g_4({}) > {}, {} nodes, {:.1} s",
                n,
                hi,
                nodes.into_inner(),
                clock.elapsed().as_secs_f64()
            );
            return;
        }
        let (sets, _) = four_top(n, g, &floor, true);
        assert!(sets.iter().all(|s| plain_free(s)));
        println!(
            "g_4({}) = {}, {} extremal sets {:?}, {} nodes, {:.1} s",
            n,
            g,
            sets.len(),
            sets,
            nodes.into_inner(),
            clock.elapsed().as_secs_f64()
        );
        floor.push(g);
    }
}

fn four_plain(n: usize, hi: usize) {
    fn walk(set: &mut Vec<usize>, n: usize, best: &mut usize, wit: &mut Vec<usize>) {
        if set.len() == n {
            *best = *set.last().unwrap();
            *wit = set.clone();
            return;
        }
        let from = set.last().map_or(1, |&x| x + 1);
        let mut a = from;
        while a < *best {
            set.push(a);
            if plain_free(set) {
                walk(set, n, best, wit);
            }
            set.pop();
            a += 1;
        }
    }
    let clock = Instant::now();
    let mut best = hi + 1;
    let mut wit = vec![];
    walk(&mut vec![], n, &mut best, &mut wit);
    if wit.is_empty() {
        println!(
            "no {}-set in [1, {}], {:.1} s",
            n,
            hi,
            clock.elapsed().as_secs_f64()
        );
    } else {
        println!("least largest element of a {}-set with no 4-term progression in its subset sums: {}, witness {:?}, {:.1} s", n, best, wit, clock.elapsed().as_secs_f64());
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let num = |i: usize, d: u64| args.get(i).map(|s| s.parse().unwrap()).unwrap_or(d);
    match args.get(1).map(|s| s.as_str()).unwrap_or("") {
        "census" => census(num(2, 10), num(3, 4) as u32, num(4, 31) as u32, num(5, 4) as usize, false),
        "deep" => census(num(2, 10), num(3, 5) as u32, num(4, 31) as u32, num(5, 4) as usize, true),
        "wcell" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            wcell(&bases, num(3, 1) as u32, num(4, 31) as u32)
        }
        "chainless" => chainless(num(2, 27) as u32),
        "forced" => forced(num(2, 26) as u32),
        "windowed" => windowed(num(2, 27) as u32),
        "jumps" => jumps(num(2, 10), num(3, 5) as u32, num(4, 31) as u32, num(5, 4) as usize),
        "trace" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            trace(&bases, num(3, 1) as u32, num(4, 31) as u32)
        }
        "wcontrol" => wcontrol(num(2, 16) as u32),
        "split" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            let x0: u128 = args[4].parse().unwrap();
            split(&bases, num(3, 1) as u32, x0, num(5, 64) as u128)
        }
        "control" => control(),
        "cell" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            cell(&bases, num(3, 1) as u32, num(4, 31) as u32)
        }
        "set" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            sets(&bases, num(3, 1) as u32, num(4, 20000) as usize)
        }
        "windows" => {
            let bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            windows(&bases, num(3, 1) as u32, 10u128.pow(num(4, 6) as u32), 10u128.pow(num(5, 30) as u32))
        }
        "graham" => graham(num(2, 2) as u128, num(3, 5) as u128, num(4, 3) as u128, num(5, 22) as u32),
        "refute" => refute(num(2, 27) as u32),
        "ternary" => ternary(num(2, 6) as usize, num(3, 200) as usize, num(4, 8) as usize),
        "probe" => probe(num(2, 7) as usize, num(3, 420) as usize, num(4, 8) as usize),
        "offsets" => offsets(num(2, 420) as usize, num(3, 504) as usize),
        "band" => band(num(2, 419) as usize, num(3, 504) as usize, num(4, 70) as usize, num(5, 8) as usize),
        "sweep" => sweep(num(2, 7) as usize, num(3, 169) as usize, num(4, 473) as usize, num(5, 8) as usize),
        "tally" => tally(num(2, 6) as usize, num(3, 168) as usize, num(4, 176) as usize, num(5, 8) as usize, num(6, 1) == 1),
        "route" => route(num(2, 10), 10u128.pow(num(3, 12) as u32)),
        "tight" => tight(num(2, 4) as u32, num(3, 31) as u32, num(4, 25) as u32),
        "ones" => unit_sets(num(2, 47), num(3, 6) as usize),
        "pair" | "lpair" => {
            let mut bases: Vec<u64> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            bases.sort();
            let rule = if args[1] == "pair" { Rule::Smooth } else { Rule::Logs };
            for k in 1..=num(3, 1) as u32 {
                let (line, cell) = paired_row(&bases, k, num(4, 31) as u32, num(5, 25) as u32, rule);
                println!("{}", line);
                if cell.is_none() {
                    break;
                }
            }
        }
        "pcontrol" => pcontrol(num(2, 47), num(3, 6) as usize, num(4, 2) as u32, Rule::Smooth, false),
        "lcontrol" => pcontrol(num(2, 47), num(3, 6) as usize, num(4, 2) as u32, Rule::Logs, num(5, 0) == 1),
        "pairs" => pairs(num(2, 47), num(3, 6) as usize, num(4, 4) as u32, num(5, 31) as u32, num(6, 25) as u32, num(7, 4) as usize, Rule::Smooth, false),
        "logs" => pairs(num(2, 47), num(3, 6) as usize, num(4, 4) as u32, num(5, 31) as u32, num(6, 25) as u32, num(7, 4) as usize, Rule::Logs, num(8, 0) == 1),
        "near" => near_pair(num(2, 3), num(3, 5), num(4, 24) as u32),
        "cantor" => {
            let exps: Vec<u32> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            cantor(&exps, num(3, 1 << 24) as usize)
        }
        "carry" => {
            let mults: Vec<u128> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            carry_set(&mults, &args[2], num(3, 1 << 24) as usize)
        }
        "deaths" => deaths(num(2, 3) as usize, num(3, 400) as u32, num(4, 79) as u32, num(5, 8) as usize),
        "classes" => classes(num(2, 3) as usize, num(3, 5) as u32, num(4, 8) as usize),
        "tuples" => tuples(num(2, 2) as usize, num(3, 40) as u32, num(4, 8) as usize, num(5, 1 << 24) as usize),
        "four" => four(num(2, 7) as usize, num(3, 5000) as usize, num(4, 8) as usize),
        "fourtop" => {
            let floor: Vec<usize> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
            let clock = Instant::now();
            let (sets, nodes) = four_top(floor.len(), num(3, 100) as usize, &floor, num(4, 1) == 1);
            println!("{} sets {:?}, {} nodes, {:.1} s", sets.len(), sets, nodes, clock.elapsed().as_secs_f64());
        }
        "fourplain" => four_plain(num(2, 5) as usize, num(3, 100) as usize),
        _ => println!("verbs census R K BITS THREADS, deep R K BITS THREADS, jumps R K BITS THREADS, trace D K BITS, windowed BITS, cell D K BITS, wcell D K BITS, chainless BITS, forced BITS, wcontrol BITS, split D K X COUNT, set D K B, windows D K E F, graham T P Q BITS, refute BITS, ternary N HI THREADS, probe N TOP THREADS, band LO HI W THREADS, offsets LO HI, sweep N LO HI THREADS, tally N LO HI THREADS PLAIN, control, route R E, tight K BITS HOLES, ones R SIZE, pair D K BITS HOLES, pairs R SIZE K BITS HOLES THREADS, pcontrol R SIZE K, near E F BITS, lpair D K BITS HOLES, logs R SIZE K BITS HOLES THREADS ALL, lcontrol R SIZE K ALL, cantor E1,E2,... CAP, carry M1,M2,... CAP, tuples R TOP THREADS CAP, deaths R TOP D THREADS, classes R D THREADS, four N HI THREADS, fourplain N HI"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_four_term_values_to_five() {
        let mut floor = vec![0usize];
        for (n, g) in [1usize, 3, 5, 14, 40].into_iter().enumerate() {
            assert!(four_top(n + 1, g - 1, &floor, false).0.is_empty() || g == 1);
            assert!(!four_top(n + 1, g, &floor, false).0.is_empty());
            floor.push(g);
        }
        assert!(!plain_free(&[1, 2]));
        assert!(plain_free(&[2, 29, 45, 74, 77, 79]));
        assert!(plain_free(&[1, 6, 87, 135, 222, 231, 237]));
        assert!(!plain_free(&[1, 2, 29, 45, 74, 77, 79]));
    }

    #[test]
    fn the_unit_sets_below_eight_are_one() {
        assert_eq!(unit(8), vec![vec![3, 4, 7]]);
        assert_eq!(sigma(&[3, 4, 7]), (6, 6));
    }

    #[test]
    fn the_minimal_sets_below_nine() {
        let want: Vec<Vec<u64>> = vec![
            vec![3, 4, 5],
            vec![3, 4, 6],
            vec![3, 4, 7, 8],
            vec![3, 5, 6, 7],
            vec![3, 5, 6, 8],
            vec![3, 5, 7, 8],
            vec![3, 6, 7, 8],
            vec![4, 5, 6, 7, 8],
        ];
        let mut got = minimal(8);
        got.sort();
        assert_eq!(got, want);
    }

    #[test]
    fn the_bit_array_matches_a_plain_array() {
        let mut bits = Bits::new(4096);
        let mut plain = vec![false; 4096];
        plain[0] = true;
        let mut s = 0;
        for a in [3usize, 64, 5, 130, 77, 1000, 128] {
            bits.fold(a, s + a);
            for x in (a..=s + a).rev() {
                if plain[x - a] {
                    plain[x] = true;
                }
            }
            s += a;
            for h in [1usize, 63, 64, 65, s / 2, s] {
                let want = (0..=h).rev().find(|&x| !plain[x]);
                assert_eq!(bits.hole(h), want);
            }
        }
    }

    #[test]
    fn three_four_five_at_k_one_misses_seventy_nine_last() {
        let Verdict::Found { f, .. } = certify(&[3, 4, 5], 1, 1 << 20) else {
            panic!()
        };
        assert_eq!(f, 79);
        let r = knapsack(&[3, 4, 5], 1, 5000);
        assert!(!r[79]);
        assert!((80..=5000).all(|x| r[x]));
    }

    #[test]
    fn surplus_and_a_full_spectrum_leave_a_hole_chain() {
        let terms = fibonacci(1 << 20);
        assert_eq!(&terms[..8], &[1, 3, 5, 9, 15, 25, 41, 67]);
        let bits = missing(&terms);
        for x in [
            2usize, 7, 22, 63, 172, 459, 1212, 3185, 11, 36, 103, 280, 745, 1964, 5157,
        ] {
            assert_eq!(bits.w[x / 64] >> (x % 64) & 1, 0);
        }
        let mut s = 0i64;
        for n in 0..terms.len() - 1 {
            s += terms[n] as i64;
            assert!(2 * (s - terms[n + 1] as i64) - terms[n + 1] as i64 >= -9);
        }
    }

    #[test]
    fn the_seven_set_under_four_hundred_seventy_five_is_admissible() {
        assert!(admissible(&[302, 409, 447, 459, 465, 466, 474]));
        assert!(!admissible(&[302, 409, 447, 459, 465, 467, 474]));
        let floor = [0usize, 1, 3, 8, 22];
        assert_eq!(hunt(5, 59, &floor, 2).0, None);
        assert!(hunt(5, 60, &floor, 2).0.is_some());
    }

    #[test]
    fn the_cut_spares_every_prefix_of_the_extremal_sets() {
        for set in [
            vec![19, 52, 57, 59, 60],
            vec![38, 52, 57, 59, 60],
            vec![107, 145, 159, 162, 164, 168],
            vec![107, 145, 159, 162, 166, 168],
            vec![302, 409, 447, 459, 465, 466, 474],
        ] {
            assert!(admissible(&set));
            assert!(spared(&set, 3f64.powi(set.len() as i32)));
        }
    }

    #[test]
    fn the_cut_keeps_every_set_the_plain_search_finds() {
        let floor = [0usize, 1, 3, 8, 22, 60, 168];
        for top in 168..=172 {
            let cut = sweep_top(6, top, &floor, 2, true, true).sets;
            assert_eq!(cut, sweep_top(6, top, &floor, 2, false, true).sets);
        }
        assert_eq!(sweep_top(6, 168, &floor, 2, true, true).sets.len(), 2);
        assert!(sweep_top(6, 167, &floor, 2, true, false).found.is_none());
    }

    #[test]
    fn the_wide_certificate_matches_the_full_one_past_its_window() {
        let mut truncated = 0;
        for (bases, k) in [
            (&[3u64, 4, 5][..], 2),
            (&[3, 4, 6], 2),
            (&[3, 5, 6, 7], 2),
            (&[3, 4, 7, 8], 2),
        ] {
            let full = certify(bases, k, 1 << 26);
            let (v, lanes) = wide(bases, k, 20);
            assert_eq!(v, full);
            let Verdict::Found { gaps, sum, .. } = v else {
                panic!()
            };
            assert_eq!(gaps, lanes.seed + lanes.window + lanes.pair);
            truncated += (sum > 1 << 21) as usize;
        }
        assert!(truncated > 0);
    }

    #[test]
    fn a_multiset_with_surplus_full_spectrum_and_no_window_is_incomplete() {
        let a = chainless_terms(1 << 22);
        assert_eq!(&a[..10], &[4, 5, 6, 7, 12, 18, 31, 50, 80, 129]);
        let bits = missing(&a);
        let mut sums = vec![0usize];
        for &x in &a {
            sums.push(sums[sums.len() - 1] + x);
        }
        for n in 2..a.len() - 1 {
            assert!(a[n + 1] < sums[n]);
        }
        for n in 3..a.len() - 1 {
            assert!(a[n] + kick(n) < a[n + 1] && !bits.get(a[n] + kick(n)));
        }
    }

    #[test]
    fn windows_at_every_fifth_step_leave_the_chain_open() {
        let (a, u) = windowed_terms(1 << 22);
        assert_eq!(&a[..10], &[4, 5, 6, 7, 13, 19, 47, 54, 86, 153]);
        let bits = missing(&a);
        let mut sums = vec![0usize];
        for &x in &a {
            sums.push(sums[sums.len() - 1] + x);
        }
        for n in 2..a.len() - 1 {
            assert_eq!(a[n + 1] > sums[n], n % 5 == 0);
            assert!(12 * sums[n + 1] + 2 >= 13 * a[n + 1]);
        }
        for n in 3..a.len() - 1 {
            let x = (a[n] as i64 + u[n]) as usize;
            assert!(x < a[n + 1] && !bits.get(x));
        }
    }

    #[test]
    fn the_deepest_hole_rises_only_past_the_middle() {
        for (bases, k) in [(&[3u64, 4, 5][..], 3), (&[3, 4, 6], 2), (&[3, 5, 6, 7], 2)] {
            let a = elements(bases, k, far());
            let mut rows = vec![];
            let (v, lanes) = wide_with(bases, k, 20, |n, s, t| rows.push((n, s, t)));
            let Verdict::Found { f, .. } = v else {
                panic!()
            };
            let j = jumps_of(bases, k, &a, &rows, Some(f));
            assert_eq!(j.bad, 0);
            assert!(j.rises > 10);
            assert!(j.y < 0 || lanes.last == 2);
            assert!(lanes.last != 1 || j.deep);
        }
    }

    #[test]
    fn the_unit_set_certifies_through_the_power_gap() {
        for k in [1, 2] {
            let r = tight_cell(k, 20, 25).unwrap();
            let f = r.f as usize;
            let sums = knapsack(&[3, 4, 7], k, 4 * f);
            assert!(!sums[f] && sums[f + 1..].iter().all(|&x| x));
            assert_eq!((1..=f).filter(|&x| !sums[x]).count(), r.gaps);
            assert!(r.margin > 0.0 && r.big >= 8);
        }
        assert_eq!(tight_cell(1, 16, 25).unwrap().f, 581);
    }

    #[test]
    fn the_smooth_pair_reproduces_the_tight_certificate() {
        for k in [1, 2, 3] {
            let t = tight_cell(k, 20, 25).unwrap();
            let p = paired_cell(&[3, 4, 7], k, 20, 25, Rule::Smooth)
                .ok()
                .unwrap();
            assert_eq!((p.f, p.gaps, p.n, p.slack), (t.f, t.gaps, t.n, t.slack));
        }
    }

    #[test]
    fn the_close_powers_match_machine_integers() {
        let p = closest(3, 5);
        let mut direct = vec![];
        for a in 1..=80u32 {
            let u = 3u128.pow(a);
            let b = (1..=55u32)
                .filter(|&b| 5u128.pow(b) <= u)
                .max()
                .unwrap_or(0);
            for j in [b, b + 1] {
                let v = 5u128.pow(j);
                if j >= 1 && u.abs_diff(v) < 1 << GAP {
                    direct.push((u.abs_diff(v), Some(u.min(v)), a, j));
                }
            }
        }
        let small: Vec<_> = p.rows.iter().filter(|r| r.2 <= 80).cloned().collect();
        assert_eq!(small, direct);
        assert!(p.rows.iter().all(|r| r.2 <= 80));
    }

    #[test]
    fn the_two_log_rule_matches_a_knapsack() {
        for k in [1, 2] {
            let p = paired_cell(&[3, 5, 7, 13], k, 20, 25, Rule::Logs)
                .ok()
                .unwrap();
            let f = p.f as usize;
            let sums = knapsack(&[3, 5, 7, 13], k, 4 * f + 4);
            assert!(!sums[f] && sums[f + 1..].iter().all(|&x| x));
            assert_eq!((1..=f).filter(|&x| !sums[x]).count(), p.gaps);
        }
        assert_eq!(
            paired_cell(&[3, 5, 7, 13], 1, 20, 25, Rule::Logs)
                .ok()
                .unwrap()
                .f,
            112
        );
    }

    #[test]
    fn a_set_below_the_threshold_never_certifies() {
        assert_eq!(certify(&[3, 4], 1, 1 << 22), Verdict::Cap);
        assert_eq!(certify(&[3, 6, 9, 12, 15, 21], 1, 1 << 22), Verdict::Cap);
    }

    #[test]
    fn the_carry_automaton_counts_prefixes_by_multiplication() {
        for mults in [vec![256u128], vec![4, 256], vec![1 << 12, 1 << 36, 1 << 38]] {
            let s = shape(&mults, 1 << 20);
            let mut level = vec![0u128];
            let mut direct = vec![];
            for n in 0..28u32 {
                let p = 3u128.pow(n);
                let mut next = vec![];
                for &x in &level {
                    for y in [x, x + p] {
                        if mults.iter().all(|&m| (m * y / p) % 3 != 2) {
                            next.push(y);
                        }
                    }
                }
                direct.push(next.len() as u128);
                level = next;
            }
            assert_eq!(paths(&s.next, 28), direct);
        }
        let s = shape(&[4, 256], 1 << 20);
        let root = s.lo;
        assert!((root.powi(6) - root.powi(5) - 1.0).abs() < 1e-9 && s.hi - s.lo < 1e-9);
    }

    #[test]
    fn the_prefix_test_agrees_with_the_automaton() {
        for r in [2, 3] {
            for e in subsets_of(&[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24], r) {
                let mults: Vec<u128> = e.iter().map(|&x| 1u128 << x).collect();
                let s = shape(&mults, 1 << 22);
                assert_eq!(
                    death(&mults, 79).is_some(),
                    s.loops == 1 && s.rich == 0,
                    "{:?}",
                    e
                );
            }
        }
        assert_eq!(death(&[4, 256, 1 << 72], 79), Some(23));
    }

    #[test]
    fn the_surviving_classes_obey_the_first_moment() {
        for (r, d) in [(2, 4), (2, 5), (3, 3), (3, 4)] {
            let (dead, all) = dead_classes(r, d, 4);
            assert!(all - dead <= 1 << ((r as u32 + 1) * (d - 1)));
        }
    }
}
