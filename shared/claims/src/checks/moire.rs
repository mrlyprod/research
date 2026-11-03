use mrly_math::moire::pairs::{correlation, sampled, witness};
use mrly_num::factor::gcd;
use mrly_num::prime::is_prime;

const SAMPLED: [(usize, usize); 7] = [(3, 5), (3, 9), (5, 15), (9, 15), (5, 7), (7, 21), (9, 21)];

pub fn the_moire_correlation_law_for_odd() -> Result<(), String> {
    let (mut worst, mut at) = (0.0f64, (0, 0));
    for m in (3..=99).step_by(2) {
        for n in (m + 2..=99).step_by(2) {
            let exact = correlation(m, n);
            if exact != correlation(n, m) {
                return Err(format!("the law is not symmetric at {m} and {n}"));
            }
            if (exact == 0.0) != (gcd(m as u128, n as u128) == 1) {
                return Err(format!(
                    "the zero set is not the coprime set at {m} and {n}"
                ));
            }
            let gap = (exact - closed(m, n)).abs();
            if gap > worst {
                (worst, at) = (gap, (m, n));
            }
        }
    }
    if worst > 5.6e-17 {
        return Err(format!(
            "the closed form parts from the integration by {worst} at {} and {}",
            at.0, at.1
        ));
    }
    for (m, n) in SAMPLED {
        let gap = (correlation(m, n) - sampled(m, n).map_err(|e| e.to_string())?).abs();
        if gap > 1e-12 {
            return Err(format!(
                "the sampled layers part from the law by {gap} at {m} and {n}"
            ));
        }
    }
    Ok(())
}

fn closed(m: usize, n: usize) -> f64 {
    let factor = gcd(m as u128, n as u128) as i128;
    let (m, n) = (m as i128, n as i128);
    let apart = (m - 1) * (n - 1);
    let echo = factor * factor - 1;
    let spread = |side: i128| 4 * side * side - (side - 1) * (side - 1);
    (echo * (2 * apart + echo)) as f64 / (apart as f64 * ((spread(m) * spread(n)) as f64).sqrt())
}

pub fn the_stack_is_an_exact_prime() -> Result<(), String> {
    let (mut primes, mut composites) = (0, 0);
    let (mut least, mut at) = (1.0f64, 0);
    for n in (3..=199).step_by(2) {
        let trial = witness(n).map_err(|_| format!("the stack cannot try {n}"))?;
        if trial.prime != is_prime(n) {
            return Err(format!("the detector misreads {n}"));
        }
        if trial.prime {
            primes += 1;
            continue;
        }
        composites += 1;
        if trial.max <= 0.0 {
            return Err(format!("the composite {n} sits at {}", trial.max));
        }
        if trial.max < least {
            (least, at) = (trial.max, n);
        }
    }
    if (primes, composites) != (45, 54) {
        return Err(format!(
            "the stack reads {primes} primes and {composites} composites"
        ));
    }
    if at != 169 {
        return Err(format!("the least composite correlation sits at {at}"));
    }
    if (least - 0.0517383).abs() >= 1e-7 {
        return Err(format!("the least composite correlation reads {least}"));
    }
    Ok(())
}
