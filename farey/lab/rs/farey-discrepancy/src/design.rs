use mrly_num::factor::mobius_sieve;
use mrly_num::factor::totients;

pub struct Design {
    pub name: &'static str,
    pub base: u64,
    pub mask: u32,
}

pub const BASE3: Design = Design {
    name: "base 3, digits {0,1}",
    base: 3,
    mask: 0b11,
};

pub const BASE10: Design = Design {
    name: "base 10, digit 9 missing",
    base: 10,
    mask: 0b01_1111_1111,
};

pub const LADDER3: [usize; 10] = [9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147];

pub const LADDER10: [usize; 5] = [10, 100, 1000, 10000, 100000];

impl Design {
    pub fn alpha(&self) -> f64 {
        (self.mask.count_ones() as f64).ln() / (self.base as f64).ln()
    }

    pub fn holds(&self, mut n: u64) -> bool {
        if n == 0 {
            return true;
        }
        while n > 0 {
            if self.mask & (1u32 << (n % self.base)) == 0 {
                return false;
            }
            n /= self.base;
        }
        true
    }

    pub fn flags(&self, q: usize) -> Vec<bool> {
        (0..=q).map(|n| self.holds(n as u64)).collect()
    }
}

// THE INDEPENDENT COUNT

pub fn card_both(q: usize, keep: &[bool], mu: &[i8]) -> u64 {
    let mut total = 0i64;
    for (d, &sign) in mu.iter().enumerate().take(q + 1).skip(1) {
        if sign == 0 {
            continue;
        }
        let mut run = 0i64;
        let mut acc = 0i64;
        let mut m = d;
        while m <= q {
            if keep[m] {
                run += 1;
                acc += run;
            }
            m += d;
        }
        total += sign as i64 * acc;
    }
    total as u64
}

pub fn card_den(q: usize, keep: &[bool], phi: &[u64]) -> u64 {
    (1..=q).filter(|&b| keep[b]).map(|b| phi[b]).sum()
}

// THE METER

#[derive(Clone, Copy)]
pub struct Lane {
    total: u64,
    seen: u64,
    s1: f64,
    s2: f64,
    last: f64,
    gap: f64,
    gap_lo: f64,
}

impl Lane {
    fn new(total: u64) -> Self {
        Lane {
            total,
            seen: 0,
            s1: 0.0,
            s2: 0.0,
            last: 0.0,
            gap: 0.0,
            gap_lo: 0.0,
        }
    }

    fn push(&mut self, value: f64) {
        self.seen += 1;
        let delta = value - self.seen as f64 / self.total as f64;
        self.s1 += delta.abs();
        self.s2 += delta * delta;
        if value - self.last > self.gap {
            self.gap = value - self.last;
            self.gap_lo = self.last;
        }
        self.last = value;
    }
}

pub fn sweep(q: usize, keep: &[bool]) -> (Lane, Lane, Lane, u64, u64, u64) {
    let mu = mobius_sieve(q);
    let phi = totients(q);
    let all = vec![true; q + 1];
    let n_both = card_both(q, keep, &mu);
    let n_den = card_den(q, keep, &phi);
    let n_full = card_den(q, &all, &phi);
    let mut both = Lane::new(n_both);
    let mut den = Lane::new(n_den);
    let mut full = Lane::new(n_full);
    let order = q as u64;
    let (mut a, mut b, mut c, mut d) = (0u64, 1u64, 1u64, order);
    while c <= order {
        let k = (order + b) / d;
        (a, b, c, d) = (c, d, k * c - a, k * d - b);
        let value = a as f64 / b as f64;
        full.push(value);
        if keep[b as usize] {
            den.push(value);
            if keep[a as usize] {
                both.push(value);
            }
        }
    }
    (both, den, full, n_both, n_den, n_full)
}

// THE TABLE

fn verdict(ok: bool) -> &'static str {
    if ok {
        "PASS"
    } else {
        "FAIL"
    }
}

fn cell(now: f64, was: Option<f64>, span: f64) -> (String, String) {
    match was {
        Some(old) => (
            format!("{:.3}", now / old),
            format!("{:+.3}", (now / old).ln() / span),
        ),
        None => ("-".to_string(), "-".to_string()),
    }
}

struct Row {
    q: usize,
    nodes: u64,
    sieve: u64,
    s1: f64,
    s2: f64,
    gap: f64,
    gap_lo: f64,
}

fn header(title: &str, rule: &str) {
    println!("{title}");
    println!("  {rule}");
}

fn meter_rows(rows: &[Row]) -> (f64, f64, f64) {
    println!("      Q         card        sieve   chk          S2   r S2    e_2          S1   r S1    e_1");
    let mut prev: Option<&Row> = None;
    let (mut em, mut e2, mut e1) = (0.0, 0.0, 0.0);
    for row in rows {
        let span = prev
            .map(|p| (row.q as f64 / p.q as f64).ln())
            .unwrap_or(1.0);
        let (r2, x2) = cell(row.s2, prev.map(|p| p.s2), span);
        let (r1, x1) = cell(row.s1, prev.map(|p| p.s1), span);
        if let Some(p) = prev {
            em = (row.nodes as f64 / p.nodes as f64).ln() / span;
            e2 = (row.s2 / p.s2).ln() / span;
            e1 = (row.s1 / p.s1).ln() / span;
        }
        println!(
            "  {:>5}  {:>11}  {:>11}  {:>4}  {:>10.3e}  {:>5}  {:>6}  {:>10.3e}  {:>5}  {:>6}",
            row.q,
            row.nodes,
            row.sieve,
            verdict(row.nodes == row.sieve),
            row.s2,
            r2,
            x2,
            row.s1,
            r1,
            x1
        );
        prev = Some(row);
    }
    (em, e2, e1)
}

fn scale_rows(rows: &[Row], alpha: f64, e: f64) {
    println!(
        "      Q         card  exp card   S2 Q^{:.3}  S1/Q^{:.3}      S1/card      S2/card   widest gap   from",
        e - alpha,
        alpha / 2.0
    );
    let mut prev: Option<&Row> = None;
    for row in rows {
        let xm = match prev {
            Some(p) => {
                let span = (row.q as f64 / p.q as f64).ln();
                format!("{:+.3}", (row.nodes as f64 / p.nodes as f64).ln() / span)
            }
            None => "-".to_string(),
        };
        let q = row.q as f64;
        println!(
            "  {:>5}  {:>11}  {:>8}  {:>10.4e}  {:>10.4e}  {:>11.4e}  {:>11.4e}  {:>11.5}  {:>7.5}",
            row.q,
            row.nodes,
            xm,
            row.s2 * q.powf(e - alpha),
            row.s1 / q.powf(alpha / 2.0),
            row.s1 / row.nodes as f64,
            row.s2 / row.nodes as f64,
            row.gap,
            row.gap_lo
        );
        prev = Some(row);
    }
}

pub fn ladder(design: &Design, qs: &[usize]) {
    let alpha = design.alpha();
    let top = *qs.last().unwrap();
    let keep = design.flags(top);
    let mut both: Vec<Row> = Vec::new();
    let mut den: Vec<Row> = Vec::new();
    let mut full: Vec<Row> = Vec::new();
    for &q in qs {
        let (lb, ld, lf, nb, nd, nf) = sweep(q, &keep[..=q]);
        for (lane, sieve, sink) in [(lb, nb, &mut both), (ld, nd, &mut den), (lf, nf, &mut full)] {
            sink.push(Row {
                q,
                nodes: lane.seen,
                sieve,
                s1: lane.s1,
                s2: lane.s2,
                gap: lane.gap,
                gap_lo: lane.gap_lo,
            });
        }
    }
    println!();
    println!("{}   alpha = {:.6}", design.name.to_uppercase(), alpha);
    println!();
    header(
        "  CONVENTION strict",
        "F_Q(S_F) = { a/b reduced : 0 < a <= b <= Q, a in S_F, b in S_F }",
    );
    let read_both = meter_rows(&both);
    println!();
    scale_rows(&both, alpha, 2.0 * alpha);
    println!();
    header(
        "  CONVENTION denominator",
        "F_Q(S_F) = { a/b reduced : 0 < a <= b <= Q, b in S_F }",
    );
    let read_den = meter_rows(&den);
    println!();
    scale_rows(&den, alpha, 1.0 + alpha);
    println!();
    header(
        "  CONTROL full set",
        "F_Q = { a/b reduced : 0 < a <= b <= Q }",
    );
    let read_full = meter_rows(&full);
    println!();
    scale_rows(&full, 1.0, 2.0);
    println!();
    println!("  THE TRANSPLANTED SHAPE   D_Q = #{{b in S_F, b <= Q}} ~ Q^a and card ~ Q^e.");
    println!("  Square-root cancellation in the denominators is a count error of order sqrt(D_Q),");
    println!(
        "  which puts e_2 at a - e, and Cauchy-Schwarz on S1 <= sqrt(card S2) caps e_1 at a/2."
    );
    println!("  Franel-Landau are the case a = 1, e = 2. Both are caps on a limsup, never values:");
    println!(
        "  the control's own e_2 and e_1 wander across this ladder, so one rung proves nothing."
    );
    println!("  lane              a      e   a - e     e_2     a/2     e_1  exp card");
    shape("strict", alpha, 2.0 * alpha, read_both);
    shape("denominator", alpha, 1.0 + alpha, read_den);
    shape("full set", 1.0, 2.0, read_full);
}

fn shape(lane: &str, alpha: f64, e: f64, read: (f64, f64, f64)) {
    let (em, e2, e1) = read;
    println!(
        "  {:<12}  {:>5.3}  {:>5.3}  {:>+6.3}  {:>+6.3}  {:>6.3}  {:>+6.3}  {:>+8.3}",
        lane,
        alpha,
        e,
        alpha - e,
        e2,
        alpha / 2.0,
        e1,
        em
    );
}

pub fn run() {
    println!("THE RESTRICTED FAREY METER");
    println!("  S_F is the design: every base digit of the integer drawn from the digit set.");
    println!("  rho_1 < ... < rho_m ascending, delta_j = rho_j - j/m, S2 = sum delta^2, S1 = sum |delta|.");
    println!("  nodes counts the enumeration, sieve the restricted totient sum built without it.");
    ladder(&BASE3, &LADDER3);
    ladder(&BASE10, &LADDER10);
    println!();
    println!("CUT   base 3 stops at 3^11 = 177147 and base 10 at 10^5. Rungs are powers of the");
    println!(
        "      base so the design's set is self-similar at every rung; a rung inside a decade"
    );
    println!("      truncates the top digit and the meter jumps. The full-set control costs one");
    println!(
        "      Stern-Brocot step per node of F_Q, 9.6e9 steps at the base 3 top: the reach is"
    );
    println!("      the control's, the restricted lanes riding the same walk for free.");
}

#[cfg(test)]
mod tests {
    use super::*;
    use mrly_num::factor::gcd;

    const FULL: Design = Design {
        name: "full set, the control",
        base: 10,
        mask: 0b11_1111_1111,
    };

    fn brute_both(q: usize, design: &Design) -> u64 {
        let mut count = 0u64;
        for b in 1..=q as u64 {
            if !design.holds(b) {
                continue;
            }
            for a in 1..=b {
                if design.holds(a) && gcd(a as u128, b as u128) == 1 {
                    count += 1;
                }
            }
        }
        count
    }

    fn brute_den(q: usize, design: &Design) -> u64 {
        let mut count = 0u64;
        for b in 1..=q as u64 {
            if !design.holds(b) {
                continue;
            }
            for a in 1..=b {
                if gcd(a as u128, b as u128) == 1 {
                    count += 1;
                }
            }
        }
        count
    }

    #[test]
    fn restricted_count_matches_brute_force() {
        for (design, q) in [(&BASE3, 243usize), (&BASE10, 200), (&FULL, 120)] {
            let keep = design.flags(q);
            let mu = mobius_sieve(q);
            assert_eq!(card_both(q, &keep, &mu), brute_both(q, design));
        }
    }

    #[test]
    fn restricted_totient_sum_matches_brute_force() {
        for (design, q) in [(&BASE3, 243usize), (&BASE10, 200), (&FULL, 120)] {
            let keep = design.flags(q);
            let phi = totients(q);
            assert_eq!(card_den(q, &keep, &phi), brute_den(q, design));
        }
    }
}
