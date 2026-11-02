use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

// THE DESIGN

const CHUNK_BITS: usize = 16;
const THREADS: usize = 4;
const SIEVE_LEVEL: usize = 15;

fn read3(mut k: u64) -> u64 {
    let mut n = 0u64;
    let mut p = 1u64;
    while k > 0 {
        if k & 1 == 1 {
            n += p;
        }
        p *= 3;
        k >>= 1;
    }
    n
}

// THE FACTORING

fn small_primes(limit: u64) -> Vec<u64> {
    let mut sieve = vec![true; limit as usize + 1];
    let mut out = Vec::new();
    for p in 2..=limit as usize {
        if sieve[p] {
            out.push(p as u64);
            let mut q = p * p;
            while q <= limit as usize {
                sieve[q] = false;
                q += p;
            }
        }
    }
    out
}

fn mu_sieve(limit: usize) -> Vec<i8> {
    let mut mu = vec![1i8; limit + 1];
    let mut composite = vec![false; limit + 1];
    mu[0] = 0;
    for p in 2..=limit {
        if composite[p] {
            continue;
        }
        let mut q = p;
        while q <= limit {
            composite[q] = true;
            mu[q] = -mu[q];
            q += p;
        }
        let mut q = p * p;
        while q <= limit {
            mu[q] = 0;
            q += p * p;
        }
    }
    mu
}

fn mulmod(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1u64;
    a %= n;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, n);
        }
        a = mulmod(a, a, n);
        e >>= 1;
    }
    r
}

const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &WITNESSES {
        if n % p == 0 {
            return n == p;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'outer: for &a in &WITNESSES {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    r
}

fn mobius(mut n: u64, primes: &[u64]) -> i8 {
    let mut mu = 1i8;
    for &p in primes {
        if p * p > n {
            break;
        }
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return 0;
            }
            mu = -mu;
        }
    }
    if n == 1 {
        return mu;
    }
    let bound = *primes.last().unwrap();
    if n < bound * bound {
        return -mu;
    }
    let r = isqrt(n);
    if r * r == n {
        return 0;
    }
    if is_prime(n) {
        -mu
    } else {
        mu
    }
}

// THE WALK

#[derive(Clone, Copy, Debug, PartialEq)]
struct Rung {
    level: u32,
    m: i64,
    a: i64,
    max: u64,
}

#[derive(Clone, Copy, Default)]
struct Chunk {
    sum: i64,
    odd: i64,
    lo: i64,
    hi: i64,
}

fn chunk(j: u64, primes: &[u64]) -> Chunk {
    let mut c = Chunk::default();
    for k in (j << CHUNK_BITS) + 1..=(j + 1) << CHUNK_BITS {
        let v = mobius(read3(k), primes) as i64;
        c.sum += v;
        if k & 1 == 1 {
            c.odd += v;
        }
        c.lo = c.lo.min(c.sum);
        c.hi = c.hi.max(c.sum);
    }
    c
}

fn walk(level: usize, primes: &[u64]) -> Vec<Rung> {
    let total = 1u64 << level;
    let first = total.min(1 << CHUNK_BITS);
    let mut rungs = Vec::new();
    let (mut sum, mut odd, mut max) = (0i64, 0i64, 0u64);
    for k in 1..=first {
        let v = mobius(read3(k), primes) as i64;
        sum += v;
        if k & 1 == 1 {
            odd += v;
        }
        max = max.max(sum.unsigned_abs());
        if k.is_power_of_two() {
            rungs.push(Rung {
                level: k.trailing_zeros(),
                m: sum,
                a: odd,
                max,
            });
        }
    }
    if total > first {
        let chunks = (total >> CHUNK_BITS) as usize;
        let next = AtomicUsize::new(1);
        let mut results = vec![Chunk::default(); chunks];
        let parts: Vec<Vec<(usize, Chunk)>> = std::thread::scope(|s| {
            let workers: Vec<_> = (0..THREADS)
                .map(|_| {
                    s.spawn(|| {
                        let mut out = Vec::new();
                        loop {
                            let j = next.fetch_add(1, Ordering::Relaxed);
                            if j >= chunks {
                                return out;
                            }
                            out.push((j, chunk(j as u64, primes)));
                        }
                    })
                })
                .collect();
            workers.into_iter().map(|w| w.join().unwrap()).collect()
        });
        for (j, c) in parts.into_iter().flatten() {
            results[j] = c;
        }
        for (j, c) in results.iter().enumerate().skip(1) {
            max = max
                .max((sum + c.lo).unsigned_abs())
                .max((sum + c.hi).unsigned_abs());
            sum += c.sum;
            odd += c.odd;
            let end = ((j + 1) as u64) << CHUNK_BITS;
            if end.is_power_of_two() {
                rungs.push(Rung {
                    level: end.trailing_zeros(),
                    m: sum,
                    a: odd,
                    max,
                });
            }
        }
    }
    rungs
}

// THE READINGS

fn theta(max: u64, level: u32) -> f64 {
    (max as f64).ln() / (level as f64 * 3f64.ln())
}

fn fit(rungs: &[Rung], i: usize) -> Option<(f64, f64, f64)> {
    if i < 5 {
        return None;
    }
    let window = &rungs[i - 4..=i];
    let xs: Vec<f64> = window.iter().map(|r| r.level as f64 * 3f64.ln()).collect();
    let ys: Vec<f64> = window.iter().map(|r| (r.max as f64).ln()).collect();
    let n = xs.len() as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let thetas: Vec<f64> = window.iter().map(|r| theta(r.max, r.level)).collect();
    let lo = thetas.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = thetas.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    Some((sxy / sxx, lo, hi))
}

fn recursion_holds(rungs: &[Rung]) -> bool {
    let mut running = 0i64;
    rungs.iter().enumerate().all(|(i, r)| {
        running += r.m;
        let prev = if i == 0 { 0 } else { rungs[i - 1].a };
        r.m == r.a - prev && r.a == running
    })
}

fn trial_bound(level: usize) -> u64 {
    let top = 3u128.pow(level as u32);
    let mut b = 3u64.pow((level / 3) as u32);
    while (b as u128).pow(3) <= top {
        b += 1;
    }
    b
}

// THE LEVEL

fn mu_level(level: usize, primes: &[u64]) -> Vec<i8> {
    let total = 1usize << level;
    let block = total.min(1 << CHUNK_BITS);
    let threads = std::thread::available_parallelism().map_or(4, |t| t.get());
    let mut a = vec![0i8; total];
    let mut lanes: Vec<Vec<(usize, &mut [i8])>> = (0..threads).map(|_| Vec::new()).collect();
    for (j, c) in a.chunks_mut(block).enumerate() {
        lanes[j % threads].push((j, c));
    }
    std::thread::scope(|s| {
        for lane in lanes {
            s.spawn(move || {
                for (j, c) in lane {
                    for (i, v) in c.iter_mut().enumerate() {
                        let k = (j * block + i) as u64;
                        *v = if k == 0 { 0 } else { mobius(read3(k), primes) };
                    }
                }
            });
        }
    });
    a
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct LevelRow {
    level: u32,
    q: i64,
    q_cop: i64,
    m: i64,
    max: u64,
    sp: i128,
    so: i128,
    su: i128,
    half: i128,
    cross: i128,
}

fn level_rows(a: &[i8]) -> Vec<LevelRow> {
    let mut rows = Vec::new();
    let (mut p, mut o, mut o_half) = (0i64, 0i64, 0i64);
    let (mut q, mut q_cop, mut max) = (0i64, 0i64, 0u64);
    let (mut sp, mut so, mut su, mut half, mut cross) = (0i128, 0i128, 0i128, 0i128, 0i128);
    for (k, &v) in a.iter().enumerate() {
        let v = v as i64;
        p += v;
        q += v * v;
        if k & 1 == 1 {
            o += v;
            q_cop += v * v;
        }
        if k >= 1 && k & 1 == 0 && (k >> 1) & 1 == 1 {
            o_half += a[k >> 1] as i64;
        }
        assert_eq!(p, o - o_half, "recursion at k = {k}");
        max = max.max(p.unsigned_abs());
        sp += (p * p) as i128;
        so += (o * o) as i128;
        if k & 1 == 1 {
            su += (o * o) as i128;
        }
        half += (o_half * o_half) as i128;
        cross += (o * o_half) as i128;
        if (k + 1).is_power_of_two() && k >= 1 {
            rows.push(LevelRow {
                level: (k + 1).trailing_zeros(),
                q,
                q_cop,
                m: p,
                max,
                sp,
                so,
                su,
                half,
                cross,
            });
        }
    }
    rows
}

fn fwht(h: &mut [i32]) {
    let mut len = 1;
    while len < h.len() {
        for block in h.chunks_mut(2 * len) {
            let (lo, hi) = block.split_at_mut(len);
            for (x, y) in lo.iter_mut().zip(hi.iter_mut()) {
                let (u, w) = (*x, *y);
                *x = u + w;
                *y = u - w;
            }
        }
        len *= 2;
    }
}

fn ladder(a: &[i8], low: bool) -> Vec<i128> {
    let n = a.len().trailing_zeros() as usize;
    let mut v: Vec<i32> = a.iter().map(|&x| x as i32).collect();
    let mut out = vec![0i128; n + 1];
    for j in (0..=n).rev() {
        let width = 1usize << j;
        out[j] = v[..width].iter().map(|&x| (x as i128) * (x as i128)).sum();
        if j > 0 {
            let w = width / 2;
            for r in 0..w {
                v[r] = if low {
                    v[r] + v[r + w]
                } else {
                    v[2 * r] + v[2 * r + 1]
                };
            }
        }
    }
    out
}

struct Spectrum {
    low: Vec<i128>,
    top: Vec<i128>,
    below: u64,
    ties: u64,
    max: i64,
    max_index: usize,
}

fn spectrum(a: &[i8], q: i64, m: i64) -> Spectrum {
    let n = a.len().trailing_zeros() as usize;
    let low = ladder(a, true);
    let top = ladder(a, false);
    let mut h: Vec<i32> = a.iter().map(|&x| x as i32).collect();
    fwht(&mut h);
    assert_eq!(h[0] as i64, m, "trivial twist is the meter");
    let sq = |x: i32| (x as i128) * (x as i128);
    let total: i128 = h.iter().map(|&x| sq(x)).sum();
    assert_eq!(total, (q as i128) << n, "Plancherel");
    assert_eq!(low[n], q as i128);
    assert_eq!(top[n], q as i128);
    assert_eq!(low[0], (m as i128) * (m as i128));
    for j in 0..=n {
        let lo: i128 = h[..1 << j].iter().map(|&x| sq(x)).sum();
        let hi: i128 = h.iter().step_by(1 << (n - j)).map(|&x| sq(x)).sum();
        assert_eq!(lo, low[j] << j, "low ladder at j = {j}");
        assert_eq!(hi, top[j] << j, "top ladder at j = {j}");
    }
    let h0 = sq(h[0]);
    let below = h.iter().filter(|&&x| sq(x) < h0).count() as u64;
    let ties = h.iter().filter(|&&x| sq(x) == h0).count() as u64;
    let (max_index, max) = h
        .iter()
        .enumerate()
        .map(|(i, &x)| (i, (x as i64).abs()))
        .max_by_key(|&(_, x)| x)
        .unwrap();
    Spectrum {
        low,
        top,
        below,
        ties,
        max,
        max_index,
    }
}

fn level_verb(level: usize) {
    let start = Instant::now();
    let primes = small_primes(trial_bound(level));
    let alpha = 2f64.ln() / 3f64.ln();
    let a = mu_level(level, &primes);
    println!("level verb: mu on the 2^{level} elements of B_{level}, the base 3 {{0,1}} strings below 3^{level}, at {:.1}s", start.elapsed().as_secs_f64());
    let rows = level_rows(&a);
    println!("check recursion: M_F(y) = A(y) - A(y/3) at every one of the 2^{level} points y of B_{level}: exact");
    println!("columns: n, Q = #squarefree in B_n, Q_cop = those coprime to 3, M = M_F(3^n - 1), Mmax, W/Q = mean of M_F(y)^2 over y in B_n over Q, V/Q_cop = mean of A(3w + 1)^2 over w in B_(n-1) over Q_cop, rho = scale correlation of the law of cosines, read = log W / log 3^n");
    let mut window: Vec<(f64, f64)> = Vec::new();
    for r in &rows {
        let n = r.level as usize;
        let size = (1u64 << n) as f64;
        let w = r.sp as f64 / size;
        let v = r.su as f64 / (size / 2.0);
        assert_eq!(r.sp, r.so + r.half - 2 * r.cross, "law of cosines");
        let rho = if r.so > 0 && r.half > 0 {
            r.cross as f64 / ((r.so as f64) * (r.half as f64)).sqrt()
        } else {
            0.0
        };
        let read = if n >= 2 && w > 0.0 {
            w.ln() / (n as f64 * 3f64.ln())
        } else {
            0.0
        };
        println!(
            "row n={n} Q={} Q_cop={} M={} Mmax={} W={:.3} W/Q={:.6} V/Q_cop={:.6} rho={:.6} read={:.6}",
            r.q, r.q_cop, r.m, r.max, w, w / r.q as f64, v / r.q_cop as f64, rho, read
        );
        if n >= 13 {
            window.push((w / r.q as f64, rho));
        }
    }
    let bracket = |xs: Vec<f64>| {
        let lo = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    };
    let (wl, wh) = bracket(window.iter().map(|x| x.0).collect());
    let (rl, rh) = bracket(window.iter().map(|x| x.1).collect());
    println!("bracket n = 13..{level}: W/Q in [{wl:.6}, {wh:.6}], rho in [{rl:.6}, {rh:.6}]");
    println!("check law of cosines: 2^n W_n = sum O(k)^2 + sum O(k/2)^2 - 2 sum O(k) O(k/2) exact at every n = 1..{level}; alpha = {alpha:.6}");
    println!("columns: n, rank = #{{I : |M^I| < |M_F|}}/2^n, ties, max_I |M^I| / Q^(1/2), its I as a bit string over the digits n-1..0, the top digit first, low = E_j/Q at j = 1, n/2, top = E_j/Q at j = 1, n/2");
    let mut ranks: Vec<f64> = Vec::new();
    for r in rows.iter().filter(|r| r.level >= 2) {
        let n = r.level as usize;
        let s = spectrum(&a[..1 << n], r.q, r.m);
        let size = (1u64 << n) as f64;
        if n >= 12 {
            ranks.push(s.below as f64 / size);
        }
        println!(
            "walsh n={n} rank={:.6} ties={} max/sqrtQ={:.4} I={:0width$b} low1={:.6} lowhalf={:.6} top1={:.6} tophalf={:.6}",
            s.below as f64 / size,
            s.ties,
            s.max as f64 / (r.q as f64).sqrt(),
            s.max_index,
            s.low[1] as f64 / r.q as f64,
            s.low[n / 2] as f64 / r.q as f64,
            s.top[1] as f64 / r.q as f64,
            s.top[n / 2] as f64 / r.q as f64,
            width = n
        );
        if n == level {
            let fmt = |e: &[i128]| {
                e.iter()
                    .map(|&x| format!("{:.4}", x as f64 / r.q as f64))
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            println!("ladder low n={n} E_j/Q at j = 0..{n}: {}", fmt(&s.low));
            println!("ladder top n={n} E_j/Q at j = 0..{n}: {}", fmt(&s.top));
            let near = |e: &[i128]| {
                (0..=n)
                    .rev()
                    .take_while(|&j| (e[j] as f64 / r.q as f64 - 1.0).abs() <= 0.03)
                    .last()
                    .unwrap_or(n)
            };
            println!("ladder n={n}: E_j/Q within 0.03 of 1 at every j from {} to {n} on the low ladder and from {} to {n} on the top ladder", near(&s.low), near(&s.top));
        }
    }
    let mean = ranks.iter().sum::<f64>() / ranks.len() as f64;
    let (kl, kh) = bracket(ranks.clone());
    println!(
        "ranks n = 12..{level}: {} rungs, mean {mean:.6}, range [{kl:.6}, {kh:.6}], {} below 0.1",
        ranks.len(),
        ranks.iter().filter(|&&x| x < 0.1).count()
    );
    println!("check Plancherel: 2^-n sum_I (M^I)^2 = Q_n and both ladders 2^-j sum = E_j at every j <= n and every n = 2..{level}: exact");
    println!("level verb done at {:.1}s", start.elapsed().as_secs_f64());
}

// THE MULTIPLICATIVE FAKES

fn factor(mut n: u64, primes: &[u64]) -> Option<Vec<u64>> {
    let mut out = Vec::new();
    for &p in primes {
        if p * p > n {
            break;
        }
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return None;
            }
            out.push(p);
        }
    }
    if n > 1 {
        out.push(n);
    }
    Some(out)
}

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn coprime_factors(level: usize, primes: &[u64]) -> Vec<Option<Vec<u64>>> {
    let half = 1usize << (level - 1);
    let threads = std::thread::available_parallelism().map_or(4, |t| t.get());
    let mut out: Vec<Option<Vec<u64>>> = vec![None; half];
    let block = half.div_ceil(threads * 16).max(1);
    let mut lanes: Vec<Vec<(usize, &mut [Option<Vec<u64>>])>> =
        (0..threads).map(|_| Vec::new()).collect();
    for (j, c) in out.chunks_mut(block).enumerate() {
        lanes[j % threads].push((j, c));
    }
    std::thread::scope(|s| {
        for lane in lanes {
            s.spawn(move || {
                for (j, c) in lane {
                    for (i, v) in c.iter_mut().enumerate() {
                        let k = 2 * (j * block + i) as u64 + 1;
                        *v = factor(read3(k), primes);
                    }
                }
            });
        }
    });
    out
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct FakeRow {
    level: u32,
    squarefree: u64,
    forced: u64,
    max_mu: u64,
    max_greedy: u64,
    max_random: u64,
}

fn lift(odd: &[i8], k: usize) -> i64 {
    if k & 1 == 1 {
        odd[k >> 1] as i64
    } else if k > 0 && (k >> 1) & 1 == 1 {
        -(odd[k >> 2] as i64)
    } else {
        0
    }
}

fn walk_rungs(odd: &[i8], level: usize) -> Vec<(i64, u64)> {
    let total = 1usize << level;
    let (mut p, mut max) = (0i64, 0u64);
    let mut out = Vec::new();
    for k in 0..total {
        p += lift(odd, k);
        max = max.max(p.unsigned_abs());
        if (k + 1).is_power_of_two() && k >= 1 {
            out.push((p, max));
        }
    }
    out
}

fn walk_max(odd: &[i8], level: usize) -> Vec<u64> {
    walk_rungs(odd, level).into_iter().map(|(_, m)| m).collect()
}

const PRETENDERS: [(&str, u64); 10] = [
    ("principal", 1),
    ("chi4+", 4),
    ("chi4-", 4),
    ("chi8", 8),
    ("chi-8", 8),
    ("chi12", 12),
    ("chi5", 5),
    ("chi7", 7),
    ("chi11", 11),
    ("chi13", 13),
];

fn legendre(a: u64, q: u64) -> i8 {
    match powmod(a % q, (q - 1) / 2, q) {
        0 => 0,
        1 => 1,
        _ => -1,
    }
}

fn character(name: &str, m: u64) -> i8 {
    match name {
        "principal" => 1,
        "chi4+" | "chi4-" => match m % 4 {
            1 => 1,
            3 => -1,
            _ => 0,
        },
        "chi8" => match m % 8 {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        },
        "chi-8" => match m % 8 {
            1 | 3 => 1,
            5 | 7 => -1,
            _ => 0,
        },
        "chi12" => match m % 12 {
            1 | 11 => 1,
            5 | 7 => -1,
            _ => 0,
        },
        "chi5" => legendre(m, 5),
        "chi7" => legendre(m, 7),
        "chi11" => legendre(m, 11),
        "chi13" => legendre(m, 13),
        _ => unreachable!(),
    }
}

fn pretender(name: &str, q: u64, p: u64) -> i8 {
    if q % p == 0 {
        if name == "chi4-" {
            -1
        } else {
            1
        }
    } else {
        character(name, p)
    }
}

fn pretend_verb(level: usize) {
    let start = Instant::now();
    let primes = small_primes(isqrt(3u64.pow(level as u32)) + 1);
    let facs = coprime_factors(level, &primes);
    let mut per_level = vec![0u64; level];
    for (i, f) in facs.iter().enumerate() {
        if f.is_some() {
            let k = 2 * i + 1;
            per_level[(usize::BITS - k.leading_zeros()) as usize - 1] += 1;
        }
    }
    let mut cumulative = vec![0u64; level];
    let mut acc = 0u64;
    for n in 0..level {
        acc += per_level[n];
        cumulative[n] = acc;
    }
    println!("pretend verb: the pretentious members f = mu^2 psi on the base 3 {{0,1}} design to 3^{level}; psi(3) = -1 always; principal has psi(p) = 1 at every other prime; chiq has psi(p) = chi(p) off q for the real primitive character chi mod q and psi(p) = +1 at the primes dividing q, chi4- taking psi(2) = -1");
    println!("columns: member, n, M = M_f(3^n - 1) the meter at the last element below 3^n, max = max |M_f(y)| over y < 3^n, Q_A = coprime squarefree elements below 3^n, max/sqrt(Q_A)");
    for (name, q) in PRETENDERS {
        let vals: Vec<i8> = facs
            .iter()
            .map(|f| match f {
                Some(ps) => ps.iter().fold(1i8, |s, &p| s * pretender(name, q, p)),
                None => 0,
            })
            .collect();
        let rungs = walk_rungs(&vals, level);
        for (n, &(m, max)) in rungs.iter().enumerate() {
            println!(
                "row member={name} n={} M={m} max={max} Q_A={} max/sqrt(Q_A)={:.4}",
                n + 1,
                cumulative[n],
                max as f64 / (cumulative[n] as f64).sqrt()
            );
        }
        if name == "principal" {
            for n in 0..level {
                assert_eq!(rungs[n].0, per_level[n] as i64);
            }
            println!("principal: M_f(3^n - 1) equals the coprime squarefree count of level n, elements in [3^(n-1), 3^n), at every n <= {level}");
        }
    }
    let names: Vec<&str> = PRETENDERS[1..].iter().map(|(s, _)| *s).collect();
    let sf: Vec<i8> = facs.iter().map(|f| f.is_some() as i8).collect();
    let (plain, squarefree) = plain_sums(level, &sf);
    println!("plain: sum of chi(m) over m in B_n, and squarefree: the same over the squarefree m in B_n, listed for n = 1..{level}");
    for (c, name) in names.iter().enumerate() {
        println!("plain member={name} {:?}", plain[c]);
        println!("squarefree member={name} {:?}", squarefree[c]);
    }
    assert!(chi4_closed_form(&plain[0]));
    println!("chi4: the plain sum over B_n equals 2^((n-1)/2) at odd n and 0 at even n, asserted at every n <= {level}");
    for (c, q) in [(5usize, 5u64), (6, 7), (7, 11)] {
        let o = order_of_three(q);
        let bound = ((q - 1) << (o - 1)) as i64;
        assert!(plain[c].iter().all(|v| v.abs() <= bound));
        println!(
            "{}: the plain sums stay within (q - 1) 2^(o - 1) = {bound}, o = ord_q(3) = {o}, at every n <= {level}, largest {}",
            names[c],
            plain[c].iter().map(|v| v.abs()).max().unwrap()
        );
    }
    println!("pretend verb done at {:.1}s", start.elapsed().as_secs_f64());
}

fn order_of_three(q: u64) -> u32 {
    let (mut o, mut t) = (1u32, 3 % q);
    while t != 1 {
        t = t * 3 % q;
        o += 1;
    }
    o
}

fn plain_sums(level: usize, sf: &[i8]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let names: Vec<&str> = PRETENDERS[1..].iter().map(|(s, _)| *s).collect();
    let mut plain = vec![vec![0i64; level]; names.len()];
    let mut squarefree = vec![vec![0i64; level]; names.len()];
    let mut run = vec![0i64; names.len()];
    let mut run_sf = vec![0i64; names.len()];
    for k in 1..(1usize << level) {
        let m = read3(k as u64);
        let w = lift(sf, k).unsigned_abs() as i64;
        for (c, name) in names.iter().enumerate() {
            let v = character(name, m) as i64;
            run[c] += v;
            run_sf[c] += v * w;
        }
        if (k + 1).is_power_of_two() {
            let n = (usize::BITS - k.leading_zeros()) as usize - 1;
            for c in 0..names.len() {
                plain[c][n] = run[c];
                squarefree[c][n] = run_sf[c];
            }
        }
    }
    (plain, squarefree)
}

fn chi4_closed_form(plain: &[i64]) -> bool {
    plain.iter().enumerate().all(|(i, &v)| {
        let n = i + 1;
        v == if n % 2 == 1 { 1i64 << ((n - 1) / 2) } else { 0 }
    })
}

fn fakes(level: usize, primes: &[u64]) -> Vec<FakeRow> {
    let facs = coprime_factors(level, primes);
    let half = facs.len();
    let mut mu = vec![0i8; half];
    let mut greedy = vec![0i8; half];
    let mut random = vec![0i8; half];
    let mut psi: std::collections::HashMap<u64, i8> = std::collections::HashMap::new();
    let mut a = 0i64;
    let mut rows = vec![FakeRow::default(); level];
    for (i, f) in facs.iter().enumerate() {
        let Some(ps) = f else { continue };
        let k = 2 * i + 1;
        let n = (usize::BITS - k.leading_zeros()) as usize;
        let row = &mut rows[n - 1];
        row.squarefree += 1;
        mu[i] = if ps.len() % 2 == 0 { 1 } else { -1 };
        random[i] = ps
            .iter()
            .fold(1i8, |s, &p| if splitmix(p) >> 63 == 0 { s } else { -s });
        let free: Vec<u64> = ps
            .iter()
            .copied()
            .filter(|p| !psi.contains_key(p))
            .collect();
        let known: i8 = ps
            .iter()
            .filter_map(|p| psi.get(p))
            .fold(1i8, |s, &v| s * v);
        let want: i8 = if a > 0 { -1 } else { 1 };
        greedy[i] = if let Some((&last, rest)) = free.split_last() {
            for &p in rest {
                psi.insert(p, 1);
            }
            psi.insert(last, want * known);
            want
        } else {
            row.forced += 1;
            known
        };
        a += greedy[i] as i64;
    }
    let (m0, m1, m2) = (
        walk_max(&mu, level),
        walk_max(&greedy, level),
        walk_max(&random, level),
    );
    for (n, row) in rows.iter_mut().enumerate() {
        row.level = n as u32 + 1;
        row.max_mu = m0[n];
        row.max_greedy = m1[n];
        row.max_random = m2[n];
    }
    rows
}

fn fake_verb(level: usize) {
    let start = Instant::now();
    let primes = small_primes(isqrt(3u64.pow(level as u32)) + 1);
    let rows = fakes(level, &primes);
    println!("fake verb: the class f = mu^2 psi, psi completely multiplicative with values +-1 and psi(3) = -1, on the base 3 {{0,1}} design to 3^{level}; every member has the support of mu and the recursion f(3m) = -f(m); mu is the member psi = -1 at every prime");
    println!("greedy: along the coprime squarefree elements in ascending order, psi is fixed at every new prime, +1 at all but the largest new prime of the element and the largest chosen to step A toward 0; an element with no new prime is forced");
    println!("random: psi(p) = sign bit of splitmix64(p)");
    println!("columns: n, squarefree = coprime squarefree elements in [3^(n-1), 3^n), forced = those with no new prime under the greedy, max |M_f(y)| over y < 3^n for mu, greedy, random");
    for r in &rows {
        let frac = if r.squarefree > 0 {
            r.forced as f64 / r.squarefree as f64
        } else {
            0.0
        };
        println!(
            "row n={} squarefree={} forced={} forced/squarefree={:.6} max_mu={} max_greedy={} max_random={}",
            r.level, r.squarefree, r.forced, frac, r.max_mu, r.max_greedy, r.max_random
        );
    }
    println!("fake verb done at {:.1}s", start.elapsed().as_secs_f64());
}

fn main() {
    let start = Instant::now();
    let first = std::env::args().nth(1);
    if first.as_deref() == Some("fake") {
        let level: usize = std::env::args()
            .nth(2)
            .map(|s| s.parse().expect("a level"))
            .unwrap_or(16);
        assert!((2..=24).contains(&level));
        fake_verb(level);
        return;
    }
    if first.as_deref() == Some("pretend") {
        let level: usize = std::env::args()
            .nth(2)
            .map(|s| s.parse().expect("a level"))
            .unwrap_or(16);
        assert!((2..=24).contains(&level));
        pretend_verb(level);
        return;
    }
    if first.as_deref() == Some("level") {
        let level: usize = std::env::args()
            .nth(2)
            .map(|s| s.parse().expect("a level"))
            .unwrap_or(20);
        assert!((2..=27).contains(&level));
        level_verb(level);
        return;
    }
    let level: usize = first.map(|s| s.parse().expect("a level")).unwrap_or(24);
    assert!((1..=30).contains(&level));
    let bound = trial_bound(level);
    let primes = small_primes(bound);
    let top = 3u64.pow(level as u32);
    let alpha = 2f64.ln() / 3f64.ln();
    println!("design base 3 digits {{0,1}}, level {level}, 3^{level} = {top}, 2^{level} = {} elements, alpha = {alpha:.6}, alpha/2 = {:.6}", 1u64 << level, alpha / 2.0);
    println!("factoring trial division by the {} primes to {bound}, {bound}^3 > 3^{level}, cofactor 1, p, p^2 or pq, Miller-Rabin on the witnesses {:?}", primes.len(), WITNESSES);
    let sieve_top = 3usize.pow(SIEVE_LEVEL as u32);
    let mu = mu_sieve(sieve_top);
    let mismatches = (1..1u64 << SIEVE_LEVEL)
        .filter(|&k| mobius(read3(k), &primes) != mu[read3(k) as usize])
        .count();
    println!("check sieve: {} design elements below 3^{SIEVE_LEVEL} factored against a linear Mobius sieve, {mismatches} mismatches", (1u64 << SIEVE_LEVEL) - 1);
    assert_eq!(mismatches, 0);
    println!("sieve done at {:.1}s", start.elapsed().as_secs_f64());
    let rungs = walk(level, &primes);
    println!(
        "walk done at {:.1}s on {THREADS} threads",
        start.elapsed().as_secs_f64()
    );
    println!("columns: n, M_F(3^n), A(3^n), Mmax = max |M_F(y)| over y <= 3^n, theta = log Mmax / log 3^n, theta/alpha = log Mmax / log 2^n, fit = least-squares slope of log Mmax on log x over n-4..n, bracket = [min, max] of theta over n-4..n, |A|/(n+1) = the floor the recursion puts under max_(j <= n) |M_F(3^j)|");
    for (i, r) in rungs.iter().enumerate() {
        let th = theta(r.max, r.level.max(1));
        let (f, lo, hi) = match fit(&rungs, i) {
            Some((f, lo, hi)) => (format!("{f:.4}"), format!("{lo:.4}"), format!("{hi:.4}")),
            None => ("-".into(), "-".into(), "-".into()),
        };
        println!(
            "row n={} M={} A={} Mmax={} theta={:.6} theta/alpha={:.4} fit={} bracket=[{}, {}] |A|/(n+1)={:.3}",
            r.level, r.m, r.a, r.max, th, th / alpha, f, lo, hi, r.a.unsigned_abs() as f64 / (r.level as f64 + 1.0)
        );
    }
    let ok = recursion_holds(&rungs);
    println!("check recursion: M_F(3^n) = A(3^n) - A(3^(n-1)) and A(3^n) = sum_(j <= n) M_F(3^j) at every n = 0..{level}: {}", if ok { "exact" } else { "BROKEN" });
    assert!(ok);
    let last = rungs.last().unwrap();
    println!("box: 3^{level}, M_F = {}, A = {}, Mmax = {}, theta = {:.6} against alpha/2 = {:.6}, {:.1}s", last.m, last.a, last.max, theta(last.max, last.level), alpha / 2.0, start.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn primes() -> Vec<u64> {
        small_primes(6600)
    }

    #[test]
    fn trial_bound_cubes_past_the_top() {
        assert_eq!(trial_bound(24), 6562);
        assert_eq!(trial_bound(27), 19684);
        assert!((trial_bound(26) as u128).pow(3) > 3u128.pow(26));
        assert!(((trial_bound(26) - 1) as u128).pow(3) <= 3u128.pow(26));
    }

    #[test]
    fn elements_ascend_and_start_right() {
        let first: Vec<u64> = (1..8).map(read3).collect();
        assert_eq!(first, vec![1, 3, 4, 9, 10, 12, 13]);
        assert!((1..1u64 << 12).all(|k| read3(k) < read3(k + 1)));
        assert_eq!(read3(1 << 24), 3u64.pow(24));
    }

    #[test]
    fn mu_matches_sieve() {
        let mu = mu_sieve(200_000);
        let p = primes();
        assert!((1..=200_000u64).all(|n| mobius(n, &p) == mu[n as usize]));
    }

    #[test]
    fn two_large_primes_resolved() {
        let p = primes();
        assert_eq!(mobius(6607 * 6619, &p), 1);
        assert_eq!(mobius(6607 * 6607, &p), 0);
        assert_eq!(mobius(2 * 6607 * 6619, &p), -1);
        assert_eq!(mobius(281_474_976_710_597, &p), -1);
        assert_eq!(mobius(999_983 * 999_979, &p), 1);
    }

    #[test]
    fn strong_pseudoprimes_rejected() {
        for n in [
            3_215_031_751u64,
            2_152_302_898_747,
            3_474_749_660_383,
            341_550_071_728_321,
        ] {
            assert!(!is_prime(n));
        }
        for n in [2u64, 3, 5, 6607, 999_983, 281_474_976_710_597] {
            assert!(is_prime(n));
        }
    }

    #[test]
    fn isqrt_exact() {
        for r in [1u64, 2, 6600, 65_535, 524_287, 524_288] {
            assert_eq!(isqrt(r * r), r);
            assert_eq!(isqrt(r * r + 1), r);
            assert_eq!(isqrt(r * r - 1), r - 1);
        }
    }

    #[test]
    fn walk_matches_mobius_census_to_16() {
        let rungs = walk(16, &primes());
        let census = [
            (4, -2, 3),
            (6, 2, 5),
            (8, 2, 8),
            (10, 5, 13),
            (12, 56, 61),
            (14, 11, 105),
            (16, 149, 173),
        ];
        for (level, m, max) in census {
            let r = rungs.iter().find(|r| r.level == level).unwrap();
            assert_eq!((r.m, r.max), (m, max as u64), "level {level}");
        }
    }

    #[test]
    fn walk_recursion_exact_to_16() {
        assert!(recursion_holds(&walk(16, &primes())));
    }

    #[test]
    fn level_rows_match_walk_to_16() {
        let p = primes();
        let rows = level_rows(&mu_level(16, &p));
        let rungs = walk(16, &p);
        for r in rows.iter().filter(|r| r.level >= 2) {
            let g = rungs.iter().find(|g| g.level == r.level).unwrap();
            assert_eq!((r.m, r.max), (g.m, g.max), "level {}", r.level);
            assert_eq!(r.sp, r.so + r.half - 2 * r.cross);
        }
    }

    #[test]
    fn walsh_spectrum_matches_brute_force_at_10() {
        let a = mu_level(10, &primes());
        let mut h: Vec<i32> = a.iter().map(|&x| x as i32).collect();
        fwht(&mut h);
        for (i, &v) in h.iter().enumerate() {
            let brute: i32 = a
                .iter()
                .enumerate()
                .map(|(k, &x)| {
                    if (i & k).count_ones() % 2 == 0 {
                        x as i32
                    } else {
                        -(x as i32)
                    }
                })
                .sum();
            assert_eq!(v, brute);
        }
        let rows = level_rows(&a);
        for r in rows.iter().filter(|r| r.level >= 2) {
            spectrum(&a[..1 << r.level], r.q, r.m);
        }
    }

    #[test]
    fn fakes_carry_the_true_meter_at_14() {
        let p = small_primes(isqrt(3u64.pow(14)) + 1);
        let rows = fakes(14, &p);
        let level = level_rows(&mu_level(14, &primes()));
        for (f, l) in rows.iter().zip(&level) {
            assert_eq!((f.level, f.max_mu), (l.level, l.max));
        }
        let sq: u64 = rows.iter().map(|r| r.squarefree).sum();
        assert_eq!(sq as i64, level.last().unwrap().q_cop);
    }

    #[test]
    fn parallel_walk_matches_sequential() {
        let p = primes();
        let par = walk(18, &p);
        let (mut sum, mut odd, mut max) = (0i64, 0i64, 0u64);
        let mut seq = Vec::new();
        for k in 1..=1u64 << 18 {
            let v = mobius(read3(k), &p) as i64;
            sum += v;
            if k & 1 == 1 {
                odd += v;
            }
            max = max.max(sum.unsigned_abs());
            if k.is_power_of_two() {
                seq.push(Rung {
                    level: k.trailing_zeros(),
                    m: sum,
                    a: odd,
                    max,
                });
            }
        }
        assert_eq!(par, seq);
    }

    #[test]
    fn plain_character_sums_at_level_12() {
        let level = 12;
        let primes = small_primes(isqrt(3u64.pow(level as u32)) + 1);
        let facs = coprime_factors(level, &primes);
        let sf: Vec<i8> = facs.iter().map(|f| f.is_some() as i8).collect();
        let (plain, squarefree) = plain_sums(level, &sf);
        assert!(chi4_closed_form(&plain[0]));
        for (c, q) in [(5usize, 5u64), (6, 7), (7, 11)] {
            let bound = ((q - 1) << (order_of_three(q) - 1)) as i64;
            assert!(plain[c].iter().all(|v| v.abs() <= bound));
        }
        let literal: i64 = (1..1u64 << level)
            .map(read3)
            .filter(|&m| mobius(m, &primes) != 0)
            .map(|m| character("chi5", m) as i64)
            .sum();
        assert_eq!(squarefree[5][level - 1], literal);
    }
}
