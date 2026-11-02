use num_bigint::BigUint;
use std::env;
use std::time::Instant;

// DIGITS

fn ok_small(mut n: u128, base: u128, top: u128) -> bool {
    while n > 0 {
        if n % base > top {
            return false;
        }
        n /= base;
    }
    true
}

fn ok_big(n: &BigUint, base: u32, top: u8) -> bool {
    n.to_radix_le(base).iter().all(|&d| d <= top)
}

fn pow_big(base: u32, e: usize) -> BigUint {
    BigUint::from(base).pow(e as u32)
}

// SEARCH

struct Thin {
    cap: BigUint,
    top5: u8,
    top7: u8,
    p3: Vec<BigUint>,
    p5: Vec<BigUint>,
    p7: Vec<BigUint>,
    cut5: Vec<usize>,
    cut7: Vec<usize>,
    hits: Vec<BigUint>,
    nodes: u64,
}

impl Thin {
    fn new(cap: BigUint, top5: u8, top7: u8) -> Thin {
        let mut p3 = vec![BigUint::from(1u32)];
        while *p3.last().unwrap() < cap {
            let next = p3.last().unwrap() * 3u32;
            p3.push(next);
        }
        let l = p3.len() - 1;
        let one = BigUint::from(1u32);
        let rest: Vec<BigUint> = (0..=l).map(|k| (&p3[k] - &one) / 2u32).collect();
        let top = &rest[l];
        let mut p5 = vec![BigUint::from(1u32)];
        while p5.last().unwrap() <= top {
            let next = p5.last().unwrap() * 5u32;
            p5.push(next);
        }
        let mut p7 = vec![BigUint::from(1u32)];
        while p7.last().unwrap() <= top {
            let next = p7.last().unwrap() * 7u32;
            p7.push(next);
        }
        let cut5: Vec<usize> = rest
            .iter()
            .map(|r| p5.iter().position(|q| q > r).unwrap())
            .collect();
        let cut7: Vec<usize> = rest
            .iter()
            .map(|r| p7.iter().position(|q| q > r).unwrap())
            .collect();
        Thin {
            cap,
            top5,
            top7,
            p3,
            p5,
            p7,
            cut5,
            cut7,
            hits: Vec::new(),
            nodes: 0,
        }
    }

    fn walk(&mut self, k: usize, v: BigUint) {
        self.nodes += 1;
        if k == 0 {
            if ok_big(&v, 5, self.top5) && ok_big(&v, 7, self.top7) {
                self.hits.push(v);
            }
            return;
        }
        let one = BigUint::from(1u32);
        let h5 = &v / &self.p5[self.cut5[k]];
        if !ok_big(&h5, 5, self.top5) && !ok_big(&(&h5 + &one), 5, self.top5) {
            return;
        }
        let h7 = &v / &self.p7[self.cut7[k]];
        if !ok_big(&h7, 7, self.top7) && !ok_big(&(&h7 + &one), 7, self.top7) {
            return;
        }
        let w = &v + &self.p3[k - 1];
        self.walk(k - 1, v);
        if w < self.cap {
            self.walk(k - 1, w);
        }
    }

    fn run(mut self) -> (Vec<BigUint>, u64, usize) {
        let l = self.p3.len() - 1;
        self.walk(l, BigUint::from(0u32));
        (self.hits, self.nodes, l)
    }
}

fn members(cap: BigUint, top5: u8, top7: u8) -> Vec<BigUint> {
    Thin::new(cap, top5, top7).run().0
}

fn scan(cap: u128, top5: u128, top7: u128) -> Vec<BigUint> {
    let mut out = Vec::new();
    let mut n = 0u128;
    while n < cap {
        if ok_small(n, 3, 1) && ok_small(n, 5, top5) && ok_small(n, 7, top7) {
            out.push(BigUint::from(n));
        }
        n += 1;
    }
    out
}

fn shown(v: &[BigUint]) -> Vec<String> {
    v.iter().map(|x| x.to_string()).collect()
}

// OEIS CONTROL

const A030979: [u128; 23] = [
    0,
    1,
    10,
    756,
    757,
    3160,
    3186,
    3187,
    3250,
    7560,
    7561,
    7651,
    20007,
    59548377,
    59548401,
    45773612811,
    45775397187,
    237617431723407,
    24991943420078301,
    24991943420078302,
    24991943420078307,
    24991943715007536,
    24991943715007537,
];

fn a030979() -> Vec<BigUint> {
    A030979.iter().map(|&x| BigUint::from(x)).collect()
}

// BUDGETS

fn dim(base: f64, alphabet: f64) -> f64 {
    alphabet.ln() / base.ln()
}

fn up(x: f64) -> f64 {
    (x * 1e6).ceil() / 1e6
}

fn near(x: f64) -> f64 {
    (x * 1e6).round() / 1e6
}

fn down(x: f64) -> f64 {
    (x * 1e6).floor() / 1e6
}

fn budget() {
    let d3 = dim(3.0, 2.0);
    let d5 = dim(5.0, 3.0);
    let d7 = dim(7.0, 3.0);
    let d7w = dim(7.0, 4.0);
    println!("dim log_3 2 {:.6}", near(d3));
    println!("dim log_5 3 {:.6}", near(d5));
    println!("dim log_7 3 {:.6}", near(d7));
    println!("dim log_7 4 {:.6}", near(d7w));
    println!("pair 3 5 {:.6}", up(d3 + d5 - 1.0));
    println!("pair 3 7 {:.6}", up(d3 + d7 - 1.0));
    println!("pair 5 7 {:.6}", up(d5 + d7 - 1.0));
    println!("triple 3 5 7 {:.6}", up(d3 + d5 + d7 - 2.0));
    println!("triple 3 5 7 wide {:.6}", up(d3 + d5 + d7w - 2.0));
    println!("erdos 3 5 {:.6}", 1.0 / 2.0 + 2.0 / 4.0);
    println!("erdos 3 7 {:.6}", 1.0 / 2.0 + 2.0 / 6.0);
    println!("erdos 5 7 {:.6}", 2.0 / 4.0 + 2.0 / 6.0);
    println!("erdos 3 7 wide {:.6}", 1.0 / 2.0 + 3.0 / 6.0);
    println!("erdos 5 7 wide {:.6}", 2.0 / 4.0 + 3.0 / 6.0);
}

// VERBS

fn control() {
    let cap = BigUint::from(100_000_000u128);
    let thin = members(cap.clone(), 2, 2);
    let brute = scan(100_000_000u128, 2, 2);
    println!("thin below 10^8 pruned {:?}", shown(&thin));
    println!("thin agrees with the scan {}", thin == brute);
    let wide = members(cap, 2, 3);
    let wide_brute = scan(100_000_000u128, 2, 3);
    println!("wide below 10^8 pruned {:?}", shown(&wide));
    println!("wide agrees with the scan {}", wide == wide_brute);
    let far = BigUint::from(*A030979.last().unwrap() + 1);
    let wide_far = members(far, 2, 3);
    println!("wide terms below the last A030979 term {}", wide_far.len());
    println!("wide rebuilds A030979 {}", wide_far == a030979());
    println!(
        "thin sits inside wide {}",
        thin.iter().all(|x| wide.contains(x))
    );
}

fn report(base: u32, exp: usize, top7: u8, name: &str) {
    let cap = pow_big(base, exp);
    let t0 = Instant::now();
    let shownice = cap.to_string();
    let (hits, nodes, l) = Thin::new(cap, 2, top7).run();
    if shownice.len() <= 40 {
        println!("{} height {}^{} = {}", name, base, exp, shownice);
    } else {
        println!(
            "{} height {}^{} of {} decimal digits",
            name,
            base,
            exp,
            shownice.len()
        );
    }
    println!("{} powers of three {}", name, l);
    println!("{} nodes {}", name, nodes);
    println!("{} count {}", name, hits.len());
    if hits.len() <= 16 {
        println!("{} members {:?}", name, shown(&hits));
    } else {
        println!("{} largest {}", name, hits.last().unwrap());
    }
    let e = (hits.len() as f64).ln() / (exp as f64 * (base as f64).ln());
    println!("{} effective exponent {:.6}", name, down(e));
    println!("{} seconds {:.2}", name, t0.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let verb = args.get(1).map(|s| s.as_str()).unwrap_or("budget");
    let n: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(0);
    let arg = move || n;
    match verb {
        "budget" => budget(),
        "control" => control(),
        "reach" => big(move || report(3, arg(), 2, "thin")),
        "seven" => big(move || report(7, arg(), 2, "thin")),
        "wide" => big(move || report(3, arg(), 3, "wide")),
        "ten" => big(move || report(10, arg(), 3, "wide")),
        _ => println!("verbs budget control reach seven wide ten"),
    }
}

fn big<F: FnOnce() + Send + 'static>(f: F) {
    std::thread::Builder::new()
        .stack_size(1 << 30)
        .spawn(f)
        .unwrap()
        .join()
        .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(v: &[u128]) -> Vec<BigUint> {
        v.iter().map(|&x| BigUint::from(x)).collect()
    }

    #[test]
    fn the_thin_set_below_seven_to_the_seventeen() {
        let cap = pow_big(7, 17);
        assert_eq!(cap.to_string(), "232630513987207");
        assert_eq!(members(cap, 2, 2), list(&[0, 1, 3186, 3187, 20007]));
    }

    #[test]
    fn the_thin_set_gains_no_member_by_three_to_the_two_hundred() {
        assert_eq!(
            members(pow_big(3, 200), 2, 2),
            list(&[0, 1, 3186, 3187, 20007])
        );
    }

    #[test]
    fn the_pruned_walk_agrees_with_a_direct_scan() {
        let cap = 3u128.pow(15);
        assert_eq!(members(BigUint::from(cap), 2, 2), scan(cap, 2, 2));
        assert_eq!(members(BigUint::from(cap), 2, 3), scan(cap, 2, 3));
    }

    #[test]
    fn the_wide_variant_rebuilds_a030979() {
        let far = BigUint::from(*A030979.last().unwrap() + 1);
        assert_eq!(members(far, 2, 3), a030979());
    }

    #[test]
    fn every_member_passes_all_three_digit_rules() {
        for n in members(pow_big(3, 120), 2, 2) {
            assert!(ok_big(&n, 3, 1));
            assert!(ok_big(&n, 5, 2));
            assert!(ok_big(&n, 7, 2));
        }
    }

    #[test]
    fn the_pair_budgets_are_positive_and_the_triple_budget_is_not() {
        let d3 = dim(3.0, 2.0);
        let d5 = dim(5.0, 3.0);
        let d7 = dim(7.0, 3.0);
        assert_eq!(format!("{:.6}", up(d3 + d5 - 1.0)), "0.313536");
        assert_eq!(format!("{:.6}", up(d3 + d7 - 1.0)), "0.195505");
        assert_eq!(format!("{:.6}", up(d5 + d7 - 1.0)), "0.247182");
        assert_eq!(format!("{:.6}", up(d3 + d5 + d7 - 2.0)), "-0.121889");
        assert_eq!(
            format!("{:.6}", up(d3 + d5 + dim(7.0, 4.0) - 2.0)),
            "0.025951"
        );
    }

    #[test]
    fn the_dimensions_print_to_six_places() {
        assert_eq!(format!("{:.6}", near(dim(3.0, 2.0))), "0.630930");
        assert_eq!(format!("{:.6}", near(dim(5.0, 3.0))), "0.682606");
        assert_eq!(format!("{:.6}", near(dim(7.0, 3.0))), "0.564575");
        assert_eq!(format!("{:.6}", near(dim(7.0, 4.0))), "0.712414");
    }
}
