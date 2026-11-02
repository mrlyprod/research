use mrly_num::series::{zeta, EULER};
use std::f64::consts::{LN_2, PI};

const LAYERS: usize = 500_000;
const BRUTE: usize = 6_001;
const SPAN: usize = 999_999;
const CUT: usize = 10_000_001;
const POINTS: [f64; 3] = [0.6180339887498949, 0.4142135623730951, 0.7320508075688772];
const LINES: [(i64, i64, i64); 20] = [
    (1, -1, 0),
    (1, 1, 1),
    (3, -3, 1),
    (3, -3, 2),
    (5, -5, 2),
    (1, -3, 0),
    (3, 1, 1),
    (3, 5, 2),
    (3, -5, 1),
    (5, 7, 3),
    (3, -9, 2),
    (9, -3, 1),
    (7, -3, 4),
    (1, 2, 0),
    (2, 1, 1),
    (2, -2, 1),
    (2, 3, 1),
    (3, 4, 1),
    (4, -6, 1),
    (2, -4, 1),
];

fn parity(x: f64) -> i64 {
    if (x.floor() as i64).rem_euclid(2) == 0 {
        1
    } else {
        -1
    }
}

fn divisor(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        divisor(b, a % b)
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn predicted(alpha: i64, beta: i64, gamma: i64) -> f64 {
    if alpha % 2 == 0 || beta % 2 == 0 {
        return 0.0;
    }
    let sign = if (alpha < 0) == (beta < 0) { -1.0 } else { 1.0 };
    let phase = if gamma % 2 == 0 { 1.0 } else { -1.0 };
    sign * phase / (alpha * beta).abs() as f64
}

fn read(alpha: i64, beta: i64, gamma: i64, point: f64) -> (f64, f64, f64) {
    let other = (gamma as f64 - alpha as f64 * point) / beta as f64;
    let (mut pair, mut carpet, mut void) = (0i64, 0i64, 0i64);
    for n in (1..2 * LAYERS).step_by(2) {
        let scale = n as f64;
        let (left, right) = (parity(scale * point), parity(scale * other));
        pair += left * right;
        carpet += i64::from(left == 1 || right == 1);
        void += i64::from(left == right);
    }
    let layers = LAYERS as f64;
    (
        pair as f64 / layers,
        carpet as f64 / layers - 0.75,
        void as f64 / layers - 0.5,
    )
}

pub fn rays() {
    println!("pair term, carpet and void ink excess at generic points of a u + b v = c, {LAYERS} odd layers, three points");
    println!("  predicted pair -sgn(ab) (-1)^c/|ab| for odd a and b, else 0; carpet excess -pair/4, void excess +pair/2");
    let (mut odd, mut even, mut carpet, mut void) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for (alpha, beta, gamma) in LINES {
        assert_eq!(divisor(divisor(alpha, beta), gamma), 1);
        let want = predicted(alpha, beta, gamma);
        let reads: Vec<(f64, f64, f64)> = POINTS
            .iter()
            .map(|point| read(alpha, beta, gamma, *point))
            .collect();
        let gap = reads.iter().map(|r| (r.0 - want).abs()).fold(0.0, f64::max);
        carpet = carpet.max(
            reads
                .iter()
                .map(|r| (r.1 + want / 4.0).abs())
                .fold(0.0, f64::max),
        );
        void = void.max(
            reads
                .iter()
                .map(|r| (r.2 - want / 2.0).abs())
                .fold(0.0, f64::max),
        );
        if want == 0.0 {
            even = even.max(gap);
        } else {
            odd = odd.max(gap);
        }
        println!(
            "  {alpha:>2} u {beta:+} v = {gamma:>2}: predicted {want:+.6}  pair {:+.6} {:+.6} {:+.6}  carpet {:+.6}  void {:+.6}",
            reads[0].0, reads[1].0, reads[2].0, reads[0].1, reads[0].2
        );
    }
    println!(
        "  worst gap {odd:.2e} on the odd-odd lines, {even:.2e} on the rest; carpet excess {carpet:.2e} and void excess {void:.2e} from their prediction"
    );
}

pub fn crosshairs() {
    let (mut matched, mut total) = (0, 0);
    for b in 1..=24usize {
        for a in 0..=b {
            if gcd(a, b) != 1 {
                continue;
            }
            let sum: i64 = (1..2 * b)
                .step_by(2)
                .map(|n| if (n * a / b) % 2 == 0 { 1 } else { -1 })
                .sum();
            let want = match (b % 2, a % 2) {
                (0, _) => 0,
                (_, 0) => 1,
                _ => -1,
            };
            total += 1;
            matched += usize::from(sum == want);
            if b <= 7 && a < b {
                println!("  u = {a}/{b}: sum of s(n a/b) over the {b} odd residues {sum:+}, predicted {want:+}");
            }
        }
    }
    println!("crosshair u = a/b, exact period average of s(n a/b) over odd n: {matched} of {total} reduced a/b with b <= 24 equal (-1)^a for odd b and 0 for even b, over b");
}

pub fn solid() {
    let mut flawed = 0;
    for point in POINTS {
        for n in (1..200_000usize).step_by(2) {
            let scale = n as f64;
            if parity(scale * point) * parity(scale * (1.0 - point)) != 1 {
                flawed += 1;
            }
        }
    }
    println!("u + v = 1 is solid void ink: s(nu) s(nv) = 1 at {flawed} of {} layers tried, u = v is s(nu)^2 = 1", 3 * 100_000);
}

struct Tables {
    factor_two: Vec<f64>,
    factor_four: Vec<f64>,
    totient: Vec<u32>,
    odd_one: Vec<f64>,
    odd_two: Vec<f64>,
}

fn tables(limit: usize) -> Tables {
    let mut smallest = vec![0u32; limit + 1];
    for p in 2..=limit {
        if smallest[p] == 0 {
            for multiple in (p..=limit).step_by(p) {
                if smallest[multiple] == 0 {
                    smallest[multiple] = p as u32;
                }
            }
        }
    }
    let mut factor_two = vec![1.0f64; limit + 1];
    let mut factor_four = vec![1.0f64; limit + 1];
    let mut totient = vec![1u32; limit + 1];
    for d in 2..=limit {
        let prime = smallest[d] as usize;
        let (mut rest, mut power) = (d, 1usize);
        while rest % prime == 0 {
            rest /= prime;
            power *= prime;
        }
        let p = prime as f64;
        factor_two[d] = factor_two[rest] * (1.0 - 1.0 / (p * p));
        factor_four[d] = factor_four[rest] * (1.0 - 1.0 / p.powi(4));
        totient[d] = totient[rest] * (power / prime * (prime - 1)) as u32;
    }
    let mut odd_one = vec![0.0f64; limit + 1];
    let mut odd_two = vec![0.0f64; limit + 1];
    for k in 1..=limit {
        let (mut one, mut two) = (odd_one[k - 1], odd_two[k - 1]);
        if k % 2 == 1 {
            one += 1.0 / k as f64;
            two += 1.0 / (k as f64 * k as f64);
        }
        odd_one[k] = one;
        odd_two[k] = two;
    }
    Tables {
        factor_two,
        factor_four,
        totient,
        odd_one,
        odd_two,
    }
}

impl Tables {
    fn odd(&self, order: i32, reach: usize) -> f64 {
        if order == 1 {
            self.odd_one[reach]
        } else {
            self.odd_two[reach]
        }
    }

    fn gcd_sum(&self, power: i32, left: i32, right: i32, limit: usize) -> f64 {
        (1..=limit)
            .step_by(2)
            .map(|d| {
                let reach = limit / d;
                let factor = if power == 2 {
                    self.factor_two[d]
                } else {
                    self.factor_four[d]
                };
                (d as f64).powi(power - left - right)
                    * factor
                    * self.odd(left, reach)
                    * self.odd(right, reach)
            })
            .sum()
    }
}

fn brute(limit: usize) -> (f64, f64, f64) {
    let (mut s2, mut s4, mut s21) = (0.0, 0.0, 0.0);
    for m in (1..=limit).step_by(2) {
        for n in (1..=limit).step_by(2) {
            let g = gcd(m, n) as f64;
            let (mf, nf) = (m as f64, n as f64);
            let ratio = g * g / (mf * nf);
            s2 += ratio;
            s4 += ratio * ratio;
            s21 += g * g / (mf * mf * nf);
        }
    }
    (s2, s4, s21)
}

fn odd_euler() -> (f64, f64, f64) {
    let (mut h, mut h2) = (0.0, 0.0);
    let (mut odd_sum, mut pair, mut cube) = (0.0, 0.0, 0.0);
    for c in (1..=CUT).step_by(2) {
        let c = c as f64;
        h += 1.0 / c;
        h2 += 1.0 / (c * c);
        odd_sum += h / (c * c);
        pair += (2.0 * h - 1.0 / c) / (c * c);
        cube += (2.0 * h2 - 1.0 / (c * c)) / (c * c * c);
    }
    let reach = CUT as f64 + 1.0;
    let tail = (reach.ln() + 1.0 + EULER + LN_2) / (2.0 * reach);
    (odd_sum + tail / 2.0, pair + tail, cube)
}

fn all_euler() -> f64 {
    let (mut h, mut total) = (0.0, 0.0);
    for n in 1..=CUT {
        let n = n as f64;
        if n % 2.0 == 1.0 {
            h += 1.0 / n;
        }
        total += h / (n * n);
    }
    let reach = CUT as f64 + 0.5;
    total + (reach.ln() + 1.0 + EULER + LN_2) / (2.0 * reach)
}

pub fn sums() {
    let zeta3 = zeta(3.0, 2_000_000).expect("three is above one");
    let zeta5 = zeta(5.0, 2_000_000).expect("five is above one");
    let small = tables(BRUTE);
    let (b2, b4, b21) = brute(BRUTE);
    let f2 = small.gcd_sum(2, 1, 1, BRUTE);
    let f4 = small.gcd_sum(4, 2, 2, BRUTE);
    let f21 = small.gcd_sum(2, 2, 1, BRUTE);
    println!("gcd sums over odd pairs to {BRUTE}, the Jordan grouping against the gcd double sum");
    println!(
        "  S2 {f2:.10} against {b2:.10}  S4 {f4:.10} against {b4:.10}  S2(2,1) {f21:.10} against {b21:.10}  worst relative gap {:.2e}",
        [(f2, b2), (f4, b4), (f21, b21)]
            .iter()
            .map(|(a, b)| ((a - b) / b).abs())
            .fold(0.0, f64::max)
    );
    let big = tables(SPAN);
    let layers = ((SPAN + 1) / 2) as f64;
    let s2 = big.gcd_sum(2, 1, 1, SPAN);
    let s4 = big.gcd_sum(4, 2, 2, SPAN);
    let s21 = big.gcd_sum(2, 2, 1, SPAN);
    let (ho, h2) = (big.odd_one[SPAN], big.odd_two[SPAN]);
    let (odd_sum, pair, cube) = odd_euler();
    let closed_pair = PI * PI * LN_2 / 4.0;
    println!(
        "Euler sums with h_c the odd harmonic number, cut at {CUT} with the leading tail added"
    );
    println!(
        "  sum over odd c of h_c/c^2 {odd_sum:.13} against 7 zeta(3)/16 + pi^2 ln2/8 = {:.13}",
        7.0 * zeta3 / 16.0 + PI * PI * LN_2 / 8.0
    );
    println!(
        "  sum over all n of h_n/n^2 {:.13} against 7 zeta(3)/8 + pi^2 ln2/8 = {:.13}",
        all_euler(),
        7.0 * zeta3 / 8.0 + PI * PI * LN_2 / 8.0
    );
    println!(
        "  T2 = sum over odd a, b of 1/(a b max(a,b)) {pair:.13} against pi^2 ln2/4 = {closed_pair:.13}"
    );
    let limit_two = 4.0 * pair / (7.0 * zeta3);
    let closed_two = PI * PI * LN_2 / (7.0 * zeta3);
    let t4 = cube;
    let limit_four = 16.0 * t4 / (31.0 * zeta5);
    println!("  T4 = sum over odd k, l of 1/(k^2 l^2 max(k,l)) {t4:.10}");
    println!(
        "  S2/N limit 4 T2/(7 zeta(3)) = {limit_two:.10}, pi^2 ln2/(7 zeta(3)) = {closed_two:.10}; S4/N limit 16 T4/(31 zeta(5)) = {limit_four:.10}"
    );
    println!("flat variance constants, odd sides 1 to {SPAN}, {layers} layers, exact covariances summed by the Jordan grouping");
    let span = SPAN as f64;
    println!(
        "  S2/N {:.10}  S4/N {:.10}  at N = {SPAN}",
        s2 / span,
        s4 / span
    );
    let tree = (s2 - ho * ho) / (4.0 * layers);
    let void = (s4 - h2 * h2) / (4.0 * layers);
    let carpet = (2.0 * (s2 - ho * ho) - 4.0 * (s21 - h2 * ho) + (s4 - h2 * h2)) / (16.0 * layers);
    let net = (2.0 * (s2 - ho * ho) + 4.0 * (s21 - h2 * ho) + (s4 - h2 * h2)) / (16.0 * layers);
    let table = [
        ("tree", tree, closed_two / 2.0),
        ("void", void, limit_four / 2.0),
        ("carpet", carpet, closed_two / 4.0 + limit_four / 8.0),
        ("net", net, closed_two / 4.0 + limit_four / 8.0),
    ];
    let scale = span.ln().powi(2) / span;
    for (name, value, limit) in table {
        println!(
            "  {name:>6} L*Var {value:.10}  limit {limit:.10}  gap {:+.2e}  gap / (ln^2 N / N) {:+.3}",
            value - limit,
            (value - limit) / scale
        );
    }
    let count: u64 = (1..=SPAN)
        .step_by(2)
        .map(|b| u64::from(big.totient[b]))
        .sum();
    let predicted = 2.0 * span * span / (PI * PI);
    println!(
        "reduced fractions with odd denominator b <= {SPAN}: {count}, (2/pi^2) B^2 = {predicted:.1}, ratio {:.8}",
        count as f64 / predicted
    );
}
