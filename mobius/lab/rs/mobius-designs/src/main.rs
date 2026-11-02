mod census;
mod factor;

use census::{digit_gcd, families, run_family, sweep, Family, Outcome};
use factor::{mobius, mu_sieve, small_primes};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const SIEVE_LIMIT: usize = 129_140_163;
const KEMPNER_LIMIT: u64 = 100_000_000;

// READOUT

fn theta(m: i64, a: u64) -> Option<f64> {
    if m == 0 || a < 2 {
        return None;
    }
    Some((m.unsigned_abs() as f64).ln() / (a as f64).ln())
}

fn theta_str(m: i64, a: u64) -> String {
    match theta(m, a) {
        Some(t) => format!("{t:.4}"),
        None => "-".to_string(),
    }
}

fn row(tag: &str, q: u64, label: &str, lev: usize, a: u64, m: i64, mx: u64) {
    println!(
        "{tag} q={q} F={label} l={lev} A={a} M={m} theta={} Mmax={mx} thetamax={}",
        theta_str(m, a),
        theta_str(mx as i64, a)
    );
}

fn drift(meter: &[u64], counts: &[u64], last: usize) -> String {
    let l = meter.len() - 1;
    let lo = if l > last { l - last + 1 } else { 1 };
    let mut values = Vec::new();
    for lev in lo..=l {
        if let Some(t) = theta(meter[lev] as i64, counts[lev]) {
            values.push(t);
        }
    }
    if values.len() < 2 {
        return "-".to_string();
    }
    let max = values.iter().cloned().fold(f64::MIN, f64::max);
    let min = values.iter().cloned().fold(f64::MAX, f64::min);
    format!("{:.4}", max - min)
}

// CONTROLS

struct ControlColumn {
    q: u64,
    levels: usize,
    meter: Vec<i64>,
    mmax: Vec<u64>,
}

type Controls = (
    Vec<ControlColumn>,
    Vec<Vec<i64>>,
    Vec<Vec<u64>>,
    Vec<Vec<u64>>,
);

fn controls_and_kempner(mu: &[i8]) -> Controls {
    let mut ctrl: Vec<ControlColumn> = [(3u64, 17usize), (4, 13), (5, 11), (10, 8)]
        .iter()
        .map(|&(q, levels)| ControlColumn {
            q,
            levels,
            meter: vec![0i64; levels + 1],
            mmax: vec![0u64; levels + 1],
        })
        .collect();
    let mut marks: Vec<(u64, usize, usize)> = Vec::new();
    for (ci, c) in ctrl.iter().enumerate() {
        let mut p = 1u64;
        for lev in 1..=c.levels {
            p *= c.q;
            marks.push((p, ci, lev));
        }
    }
    marks.sort();
    let mut kem_m = vec![vec![0i64; 9]; 10];
    let mut kem_a = vec![vec![0u64; 9]; 10];
    let mut kem_x = vec![vec![0u64; 9]; 10];
    let mut kem_run = [0i64; 10];
    let mut kem_cnt = [0u64; 10];
    let mut kem_peak = [0u64; 10];
    let mut run = 0i64;
    let mut peak = 0u64;
    let mut next = 0usize;
    for n in 1..=SIEVE_LIMIT as u64 {
        let m = mu[n as usize] as i64;
        run += m;
        peak = peak.max(run.unsigned_abs());
        if n <= KEMPNER_LIMIT {
            let mut mask = 0u16;
            let mut x = n;
            while x > 0 {
                mask |= 1 << (x % 10);
                x /= 10;
            }
            for d in 0..10 {
                if mask >> d & 1 == 0 {
                    kem_run[d] += m;
                    kem_cnt[d] += 1;
                    kem_peak[d] = kem_peak[d].max(kem_run[d].unsigned_abs());
                }
            }
        }
        while next < marks.len() && marks[next].0 == n {
            let (_, ci, lev) = marks[next];
            ctrl[ci].meter[lev] = run;
            ctrl[ci].mmax[lev] = peak;
            if ctrl[ci].q == 10 {
                for d in 0..10 {
                    kem_m[d][lev] = kem_run[d];
                    kem_a[d][lev] = kem_cnt[d];
                    kem_x[d][lev] = kem_peak[d];
                }
            }
            next += 1;
        }
    }
    assert_eq!(next, marks.len());
    (ctrl, kem_m, kem_a, kem_x)
}

// CROSSCHECK

fn crosscheck(mu: &[i8], primes: &[u64]) {
    let lmax = 16usize;
    let digits = [1u64, 2];
    let mut by_factor = vec![0i64; lmax + 1];
    let mut by_sieve = vec![0i64; lmax + 1];
    sweep(3, &digits, lmax, &mut |v, len| {
        by_factor[len] += mobius(v, primes) as i64;
        by_sieve[len] += mu[v as usize] as i64;
    });
    let mut cf = 0i64;
    let mut cs = 0i64;
    for lev in 1..=lmax {
        cf += by_factor[lev];
        cs += by_sieve[lev];
        assert_eq!(cf, cs, "method mismatch at level {lev}");
    }
    println!("crosscheck q=3 F=12 L=16 factorization equals sieve at every level, M(3^16)={cf}");
}

// MAIN

fn main() {
    println!("mobius-designs generator: cargo run --release --manifest-path research/mobius/lab/rs/mobius-designs/Cargo.toml");
    let primes = small_primes(1024);
    let mu = mu_sieve(SIEVE_LIMIT);
    let (ctrl, kem_m, kem_a, kem_x) = controls_and_kempner(&mu);
    for c in &ctrl {
        let mut p = 1u64;
        for lev in 1..=c.levels {
            p *= c.q;
            row("control", c.q, "full", lev, p, c.meter[lev], c.mmax[lev]);
        }
    }
    for d in 0..10 {
        for lev in 1..=8usize {
            row(
                "kempner",
                10,
                &format!("X{d}"),
                lev,
                kem_a[d][lev],
                kem_m[d][lev],
                kem_x[d][lev],
            );
        }
    }
    crosscheck(&mu, &primes);
    drop(mu);
    let jobs = families();
    let results: Vec<Mutex<Option<Outcome>>> = jobs.iter().map(|_| Mutex::new(None)).collect();
    let cursor = AtomicUsize::new(0);
    std::thread::scope(|s| {
        for _ in 0..6 {
            s.spawn(|| loop {
                let i = cursor.fetch_add(1, Ordering::Relaxed);
                if i >= jobs.len() {
                    break;
                }
                let out = run_family(&jobs[i], &primes);
                *results[i].lock().unwrap() = Some(out);
            });
        }
    });
    let outcomes: Vec<Outcome> = results
        .into_iter()
        .map(|r| r.into_inner().unwrap().unwrap())
        .collect();
    for (fam, out) in jobs.iter().zip(outcomes.iter()) {
        for lev in 1..=fam.lmax {
            row(
                "row",
                fam.q,
                &fam.label,
                lev,
                out.counts[lev],
                out.meter[lev],
                out.mmax[lev],
            );
        }
    }
    verify_scalings(&jobs, &outcomes);
    for (fam, out) in jobs.iter().zip(outcomes.iter()) {
        let l = fam.lmax;
        println!(
            "slope q={} F={} L={l} theta={} thetamax={} driftmax5={}",
            fam.q,
            fam.label,
            theta_str(out.meter[l], out.counts[l]),
            theta_str(out.mmax[l] as i64, out.counts[l]),
            drift(&out.mmax, &out.counts, 5)
        );
    }
    for q in [3u64, 4, 5] {
        let mut entries: Vec<(String, Option<f64>)> = jobs
            .iter()
            .zip(outcomes.iter())
            .filter(|(f, _)| f.q == q)
            .map(|(f, o)| {
                (
                    f.label.clone(),
                    theta(o.mmax[f.lmax] as i64, o.counts[f.lmax]),
                )
            })
            .collect();
        entries.sort_by(|a, b| {
            a.1.unwrap_or(f64::MIN)
                .partial_cmp(&b.1.unwrap_or(f64::MIN))
                .unwrap()
        });
        let text: Vec<String> = entries
            .iter()
            .map(|(l, t)| match t {
                Some(t) => format!("F={l} {t:.4}"),
                None => format!("F={l} -"),
            })
            .collect();
        println!(
            "distribution q={q} final-level thetamax ascending: {}",
            text.join(" | ")
        );
    }
    let mut kem: Vec<(usize, f64)> = (0..10)
        .filter_map(|d| theta(kem_x[d][8] as i64, kem_a[d][8]).map(|t| (d, t)))
        .collect();
    kem.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let text: Vec<String> = kem.iter().map(|(d, t)| format!("X{d} {t:.4}")).collect();
    println!(
        "distribution q=10 kempner thetamax at l=8 ascending: {}",
        text.join(" | ")
    );
    let mut tmax_all: Vec<f64> = Vec::new();
    let mut cut_all: Vec<f64> = Vec::new();
    let mut drift_all: Vec<f64> = Vec::new();
    for (fam, out) in jobs.iter().zip(outcomes.iter()) {
        let l = fam.lmax;
        if let Some(t) = theta(out.mmax[l] as i64, out.counts[l]) {
            tmax_all.push(t);
            if let Some(c) = theta(out.meter[l], out.counts[l]) {
                cut_all.push(c);
            }
            let vals: Vec<f64> = (l - 4..=l)
                .filter_map(|v| theta(out.mmax[v] as i64, out.counts[v]))
                .collect();
            let lo = vals.iter().cloned().fold(f64::MAX, f64::min);
            let hi = vals.iter().cloned().fold(f64::MIN, f64::max);
            drift_all.push(hi - lo);
        }
    }
    for d in 0..10 {
        if let Some(t) = theta(kem_x[d][8] as i64, kem_a[d][8]) {
            tmax_all.push(t);
        }
        if let Some(c) = theta(kem_m[d][8], kem_a[d][8]) {
            cut_all.push(c);
        }
    }
    let mut ctl_all: Vec<f64> = Vec::new();
    for c in &ctrl {
        let mut p = 1u64;
        for _ in 0..c.levels {
            p *= c.q;
        }
        if let Some(t) = theta(c.mmax[c.levels] as i64, p) {
            ctl_all.push(t);
        }
    }
    let band = |v: &[f64]| {
        let lo = v.iter().cloned().fold(f64::MAX, f64::min);
        let hi = v.iter().cloned().fold(f64::MIN, f64::max);
        (lo, hi)
    };
    let (tl, th) = band(&tmax_all);
    let (cl, ch) = band(&cut_all);
    let (dl, dh) = band(&drift_all);
    let (gl, gh) = band(&ctl_all);
    let dev = tmax_all
        .iter()
        .map(|t| (t - 0.5).abs())
        .fold(0.0f64, f64::max);
    println!(
        "band families={} thetamax {tl:.4}..{th:.4} maxdev {dev:.4} driftmax5 {dl:.4}..{dh:.4} cut {cl:.4}..{ch:.4} controls {gl:.4}..{gh:.4}",
        tmax_all.len()
    );
    println!("mobius-designs all checks pass");
}

// SCALING LAW

fn verify_scalings(jobs: &[Family], outcomes: &[Outcome]) {
    for (fam, out) in jobs.iter().zip(outcomes.iter()) {
        let g = digit_gcd(&fam.digits);
        if g < 2 {
            continue;
        }
        let base_digits: Vec<u64> = fam.digits.iter().map(|d| d / g).collect();
        let (bi, base) = jobs
            .iter()
            .enumerate()
            .find(|(_, f)| f.q == fam.q && f.digits == base_digits)
            .unwrap();
        let bout = &outcomes[bi];
        let (_, expected) = bout.twisted.iter().find(|(a, _)| *a == g).unwrap();
        let has01 = base.digits.contains(&0) && base.digits.contains(&1);
        for (lev, &want) in expected.iter().enumerate().take(fam.lmax + 1).skip(1) {
            assert_eq!(
                out.meter[lev], want,
                "scaling meter q={} F={} l={lev}",
                fam.q, fam.label
            );
            let base_a = bout.counts[lev] - if has01 { 1 } else { 0 };
            assert_eq!(
                out.counts[lev], base_a,
                "scaling count q={} F={} l={lev}",
                fam.q, fam.label
            );
        }
        println!(
            "identity q={} F={} equals the a={g} twist of F={} at every level 1..{}",
            fam.q, fam.label, base.label, fam.lmax
        );
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;
    use factor::{is_prime, isqrt};

    #[test]
    fn mu_matches_sieve() {
        let primes = small_primes(1024);
        let mu = mu_sieve(50_000);
        for n in 1..=50_000u64 {
            assert_eq!(mobius(n, &primes), mu[n as usize], "mu({n})");
        }
    }

    #[test]
    fn strong_pseudoprimes_rejected() {
        for n in [
            2047u64,
            1373653,
            25326001,
            3215031751,
            3474749660383,
            341550071728321,
        ] {
            assert!(!is_prime(n), "{n} is composite");
        }
        for n in [
            2u64,
            61,
            1_000_000_007,
            999_999_999_989,
            2_305_843_009_213_693_951,
        ] {
            assert!(is_prime(n), "{n} is prime");
        }
    }

    #[test]
    fn isqrt_exact() {
        for n in [0u64, 1, 2, 3, 4, 24, 25, 26, 999999999999999999] {
            let r = isqrt(n);
            assert!(r as u128 * r as u128 <= n as u128);
            assert!((r as u128 + 1) * (r as u128 + 1) > n as u128);
        }
    }

    #[test]
    fn methods_agree_on_family() {
        let primes = small_primes(1024);
        let mu = mu_sieve(6561);
        let mut a = vec![0i64; 9];
        let mut b = vec![0i64; 9];
        sweep(3, &[1, 2], 8, &mut |v, len| {
            a[len] += mobius(v, &primes) as i64;
            b[len] += mu[v as usize] as i64;
        });
        assert_eq!(a, b);
    }

    #[test]
    fn scaled_family_vanishes() {
        let primes = small_primes(1024);
        let fam = census::make_family_depth(5, vec![0, 4], 8);
        let out = run_family(&fam, &primes);
        for lev in 1..=fam.lmax {
            assert_eq!(out.meter[lev], 0);
        }
    }

    #[test]
    fn scaling_identity_small() {
        let primes = small_primes(1024);
        for (q, scaled, base, a) in [
            (3u64, vec![0u64, 2], vec![0u64, 1], 2u64),
            (4, vec![0, 3], vec![0, 1], 3),
        ] {
            let fs = census::make_family_depth(q, scaled, 8);
            let fb = census::make_family_depth(q, base, 8);
            let os = run_family(&fs, &primes);
            let ob = run_family(&fb, &primes);
            let (_, expected) = ob.twisted.iter().find(|(c, _)| *c == a).unwrap();
            for (lev, &want) in expected.iter().enumerate().take(9).skip(1) {
                assert_eq!(os.meter[lev], want);
            }
        }
    }

    #[test]
    fn counting_closed_form() {
        let primes = small_primes(1024);
        let f01 = census::make_family_depth(3, vec![0, 1], 8);
        let o01 = run_family(&f01, &primes);
        let f12 = census::make_family_depth(3, vec![1, 2], 8);
        let o12 = run_family(&f12, &primes);
        for lev in 1..=8usize {
            assert_eq!(o01.counts[lev], 1 << lev);
            assert_eq!(o12.counts[lev], (1 << (lev + 1)) - 2);
        }
    }

    #[test]
    fn ordered_enumeration_ascending() {
        for digits in [vec![0u64, 1], vec![0, 2], vec![1, 2]] {
            let mut prev = 0u64;
            for len in 1..=6 {
                census::sweep_length(3, &digits, len, &mut |v| {
                    assert!(v > prev);
                    prev = v;
                });
            }
        }
    }

    #[test]
    fn ordered_walk_matches_prefix() {
        let primes = small_primes(1024);
        let mu = mu_sieve(19683);
        let fam = census::make_family_depth(3, vec![0, 1], 8);
        let out = run_family(&fam, &primes);
        let mut cum = [0i64; 9];
        let mut cnt = [0u64; 9];
        sweep(3, &[0, 1], 8, &mut |v, len| {
            cum[len] += mu[v as usize] as i64;
            cnt[len] += 1;
        });
        let mut run = 0i64;
        let mut tot = 0u64;
        for lev in 1..=8usize {
            run += cum[lev];
            tot += cnt[lev];
            let boundary = if lev == 1 { mu[3] as i64 } else { 0 };
            assert_eq!(out.meter[lev], run + boundary);
            assert_eq!(out.counts[lev], tot + 1);
            assert!(out.mmax[lev] >= out.meter[lev].unsigned_abs());
            assert!(out.mmax[lev] >= out.mmax[lev - 1]);
        }
    }

    #[test]
    fn mertens_prefix_known() {
        let mu = mu_sieve(100_000);
        let mut run = 0i64;
        let mut got = Vec::new();
        for (n, &sign) in mu.iter().enumerate().take(100_001).skip(1) {
            run += sign as i64;
            if n == 10 || n == 100 || n == 1000 || n == 10_000 || n == 100_000 {
                got.push(run);
            }
        }
        assert_eq!(got, vec![-1, 1, 2, -23, -48]);
    }
}
