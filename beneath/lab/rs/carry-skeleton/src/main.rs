use mrly_math::memory::{allowed_windows, kappa, perron, transfer, Rule};
use std::collections::BTreeSet;
use std::time::Instant;

fn free(n: u64) -> u64 {
    n ^ (n << 1) ^ 1
}

fn defect(n: u64) -> u64 {
    (3 * n + 1) ^ free(n)
}

fn d_loc(n: u64) -> u32 {
    defect(n).count_ones()
}

fn carries(n: u64, bits: usize) -> u64 {
    let mut q = 1u64;
    let mut prev = 0u64;
    let mut out = 0u64;
    for i in 0..bits {
        let b = (n >> i) & 1;
        let next = u64::from(b + prev + q >= 2);
        if next == 1 {
            out |= 1 << (i + 1);
        }
        prev = b;
        q = next;
    }
    out
}

fn t_free(n: u64) -> u64 {
    if n.is_multiple_of(2) {
        n / 2
    } else {
        free(n) / 2
    }
}

fn support(q: u64) -> Vec<usize> {
    (0..64).filter(|j| (q >> j) & 1 == 1).collect()
}

fn differences(supp: &[usize]) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    for a in 0..supp.len() {
        for b in (a + 1)..supp.len() {
            out.insert(supp[b] - supp[a]);
        }
    }
    out
}

fn skeleton(n: u64, supp: &[usize]) -> u64 {
    let mut x = 1u64;
    for &j in supp {
        x ^= n << j;
    }
    x
}

fn zero_code(supp: &[usize]) -> (usize, u64) {
    let width = supp[supp.len() - 1] + 1;
    let gaps = differences(supp);
    let mut code = 0u64;
    for w in 0..(1usize << width) {
        let mut ok = true;
        for a in 0..width {
            for b in (a + 1)..width {
                let hit = (w >> (width - 1 - a)) & 1 == 1 && (w >> (width - 1 - b)) & 1 == 1;
                if hit && gaps.contains(&(b - a)) {
                    ok = false;
                }
            }
        }
        if ok {
            code |= 1 << w;
        }
    }
    (width, code)
}

fn word(n: u64, length: usize) -> Vec<usize> {
    (0..length).rev().map(|i| ((n >> i) & 1) as usize).collect()
}

fn longest_run(n: u64) -> u32 {
    let mut best = 0;
    let mut run = 0;
    for i in 0..64 {
        if (n >> i) & 1 == 1 {
            run += 1;
            if run > best {
                best = run;
            }
        } else {
            run = 0;
        }
    }
    best
}

fn v2(n: u64) -> u32 {
    if n == 0 {
        64
    } else {
        n.trailing_zeros()
    }
}

fn digits(n: u64) -> u32 {
    64 - n.leading_zeros()
}

fn depth_of(bits: &[u8]) -> usize {
    let mut q = 1u8;
    let mut prev = 0u8;
    let mut run = 0usize;
    let mut best = 0usize;
    let step = |b: u8, q: &mut u8, prev: &mut u8, run: &mut usize, best: &mut usize| {
        let next = u8::from(b + *prev + *q >= 2);
        if next == 1 {
            *run += 1;
            if *run > *best {
                *best = *run;
            }
        } else {
            *run = 0;
        }
        *prev = b;
        *q = next;
    };
    for &b in bits {
        step(b, &mut q, &mut prev, &mut run, &mut best);
    }
    for _ in 0..2 {
        step(0, &mut q, &mut prev, &mut run, &mut best);
    }
    best
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn report_identity() {
    let start = Instant::now();
    let limit = 1u64 << 18;
    let mut bad_carry = 0u64;
    let mut bad_conj = 0u64;
    let mut bad_rule = 0u64;
    for n in 0..limit {
        if carries(n, 22) != defect(n) {
            bad_carry += 1;
        }
        let m = 2 * n + 1;
        if m ^ (2 * m) != 2 * free(n) + 1 {
            bad_conj += 1;
        }
        let mut sixty = 0u64;
        for i in 0..24 {
            let hi = (m >> i) & 1;
            let lo = if i == 0 { 0 } else { (m >> (i - 1)) & 1 };
            sixty |= (hi ^ lo) << i;
        }
        if sixty != (m ^ (2 * m)) {
            bad_rule += 1;
        }
    }
    println!("identity limit 2^18");
    println!("identity carry_word_mismatches {bad_carry}");
    println!("identity conjugacy_mismatches {bad_conj}");
    println!("identity rule60_mismatches {bad_rule}");
    let mut worst = 0usize;
    for length in 2..=40usize {
        let mut u = 1u64;
        for j in 1..length {
            if j % 2 == 0 {
                u |= 1 << j;
            }
        }
        let v = u ^ 1;
        let diff = (3 * u + 1) ^ (3 * v + 1);
        let all = (2..=length).all(|i| (diff >> i) & 1 == 1);
        assert!(all, "radius witness fails at length {length}");
        worst = length;
    }
    println!("identity radius_witness_max_digit {worst}");
    println!("identity seconds {:.2}", start.elapsed().as_secs_f64());
}

fn report_rules() {
    let start = Instant::now();
    println!("rule q supp gaps width code windows rho kappa");
    for q in [3u64, 5, 7, 9, 11, 15] {
        let supp = support(q);
        let gaps = differences(&supp);
        let (width, code) = zero_code(&supp);
        let rule = Rule::new(1, width, code).expect("rule in range");
        let supp_text: Vec<String> = supp.iter().map(|j| j.to_string()).collect();
        let gap_text: Vec<String> = gaps.iter().map(|j| j.to_string()).collect();
        println!(
            "rule {q} {} {} {width} {code} {} {:.6} {:.6}",
            supp_text.join(","),
            gap_text.join(","),
            allowed_windows(&rule),
            perron(&rule),
            kappa(&rule)
        );
        let limit = 1u64 << 16;
        let mut mismatches = 0u64;
        for n in 0..limit {
            let free_carry = (q * n + 1) ^ skeleton(n, &supp) == 0;
            let padded = word(n, 16 + width);
            let accepted = n % 2 == 0 && rule.accepts(&padded);
            if free_carry != accepted {
                mismatches += 1;
            }
        }
        println!("rule {q} set_equality_mismatches_below_2^16 {mismatches}");
    }
    let golden = Rule::new(1, 2, 7).expect("golden in range");
    let supergolden = Rule::new(1, 3, 23).expect("supergolden in range");
    println!("kappa golden {:.6}", kappa(&golden));
    println!("kappa supergolden {:.6}", kappa(&supergolden));
    println!("kappa golden_transfer {:?}", transfer(&golden));
    println!("rule seconds {:.2}", start.elapsed().as_secs_f64());
}

fn report_density() {
    let start = Instant::now();
    let pi = [2i64, 1, 1, 2];
    let balance = [
        2 * pi[0] - (pi[0] + pi[1] + pi[2]),
        2 * pi[1] - pi[3],
        2 * pi[2] - pi[0],
        2 * pi[3] - (pi[1] + pi[2] + pi[3]),
    ];
    println!("density stationary_sixths {pi:?}");
    println!("density balance_residuals {balance:?}");
    println!(
        "density carry_on_mass {}/{}",
        pi[1] + pi[3],
        pi.iter().sum::<i64>()
    );
    println!("density digits mean mean_minus_half_L mean_over_L exact_residual");
    for l in 8..=22usize {
        let limit = 1u64 << l;
        let mut total = 0u128;
        for n in 0..limit {
            total += u128::from(d_loc(n));
        }
        let mean = total as f64 / limit as f64;
        let sign = if l % 2 == 0 { 1i128 } else { -1 };
        let closed = 3 * l as i128 * (limit / 2) as i128 + limit as i128 - sign;
        println!(
            "density {l} {mean:.6} {:.6} {:.6} {}",
            mean - l as f64 / 2.0,
            mean / l as f64,
            3 * total as i128 - closed
        );
    }
    println!("density seconds {:.2}", start.elapsed().as_secs_f64());
}

fn report_free() {
    let start = Instant::now();
    let limit = 1u64 << 20;
    let mut grew = 0u64;
    let mut unreached = 0u64;
    for n in 1..limit {
        if digits(t_free(n)) > digits(n) {
            grew += 1;
        }
        let mut x = n;
        let mut steps = 0;
        while x != 1 && steps < 4096 {
            x = t_free(x);
            steps += 1;
        }
        if x != 1 {
            unreached += 1;
        }
    }
    println!("free limit 2^20");
    println!("free digit_count_increases {grew}");
    println!("free values_not_reaching_one {unreached}");
    println!("free fixed_point_check {}", t_free(1));
    println!("free seconds {:.2}", start.elapsed().as_secs_f64());
}

type Stat = (&'static str, fn(u64) -> u64);

fn report_refutations() {
    let start = Instant::now();
    let limit = 1u64 << 16;
    let stats: [Stat; 6] = [
        ("popcount", |n| u64::from(n.count_ones())),
        ("longest_run", |n| u64::from(longest_run(n))),
        ("v2", |n| u64::from(v2(n))),
        ("digit_count", |n| u64::from(digits(n))),
        ("popcount_and_v2", |n| {
            u64::from(n.count_ones()) * 128 + u64::from(v2(n))
        }),
        ("all_four", |n| {
            ((u64::from(n.count_ones()) * 128 + u64::from(longest_run(n))) * 128 + u64::from(v2(n)))
                * 128
                + u64::from(digits(n))
        }),
    ];
    println!("refute statistic least_pair d_loc_pair disagreeing_pairs_below_2^16");
    for (name, key) in &stats {
        let mut first: std::collections::HashMap<u64, (u64, u32)> =
            std::collections::HashMap::new();
        let mut groups: std::collections::HashMap<u64, std::collections::HashMap<u32, u64>> =
            std::collections::HashMap::new();
        let mut witness: Option<(u64, u64, u32, u32)> = None;
        for n in 1..limit {
            let k = key(n);
            let d = d_loc(n);
            *groups.entry(k).or_default().entry(d).or_insert(0) += 1;
            match first.get(&k) {
                None => {
                    first.insert(k, (n, d));
                }
                Some(&(m, e)) => {
                    if e != d && witness.is_none() {
                        witness = Some((m, n, e, d));
                    }
                }
            }
        }
        let mut pairs = 0u64;
        for counts in groups.values() {
            let held: u64 = counts.values().sum();
            let agreeing: u64 = counts.values().map(|c| c * (c - 1) / 2).sum();
            pairs += held * (held - 1) / 2 - agreeing;
        }
        match witness {
            Some((m, n, e, d)) => println!("refute {name} {m},{n} {e},{d} {pairs}"),
            None => println!("refute {name} none none {pairs}"),
        }
    }
    println!("refute seconds {:.2}", start.elapsed().as_secs_f64());
}

fn report_depth() {
    let start = Instant::now();
    for l in [4usize, 8, 16, 32, 40] {
        let bits = vec![1u8; l];
        assert_eq!(depth_of(&bits), l + 1, "worst case fails at {l}");
    }
    println!("depth worst_case_is_L_plus_1_at_L 4,8,16,32,40");
    let golden = Rule::new(1, 2, 7).expect("golden in range");
    let phi = perron(&golden);
    let base = 2.0 / phi;
    println!("depth phi {phi:.9} two_over_phi {base:.9}");
    println!(
        "depth L samples mean_depth standard_error log_base_L offset increment_per_quadrupling"
    );
    let mut state = 0x5eed_1234_9abc_def0u64;
    let mut previous: Option<(usize, f64)> = None;
    for l in [16usize, 64, 256, 1024, 4096, 16384] {
        let samples = 200_000usize;
        let mut total = 0u64;
        let mut squares = 0u128;
        let mut bits = vec![0u8; l];
        for _ in 0..samples {
            let mut i = 0;
            while i < l {
                let chunk = splitmix(&mut state);
                let take = std::cmp::min(64, l - i);
                for j in 0..take {
                    bits[i + j] = ((chunk >> j) & 1) as u8;
                }
                i += take;
            }
            let depth = depth_of(&bits) as u64;
            total += depth;
            squares += u128::from(depth * depth);
        }
        let mean = total as f64 / samples as f64;
        let second = squares as f64 / samples as f64;
        let error = ((second - mean * mean) / (samples as f64 - 1.0)).sqrt();
        let predicted = (l as f64).ln() / base.ln();
        let increment = match previous {
            Some((pl, pm)) if l == 4 * pl => format!("{:.3}", mean - pm),
            _ => "-".to_string(),
        };
        println!(
            "depth {l} {samples} {mean:.3} {error:.3} {predicted:.3} {:.3} {increment}",
            mean - predicted
        );
        previous = Some((l, mean));
    }
    println!(
        "depth predicted_increment_per_quadrupling {:.6}",
        4.0f64.ln() / base.ln()
    );
    println!("depth seconds {:.2}", start.elapsed().as_secs_f64());
}

fn main() {
    let whole = Instant::now();
    report_identity();
    report_rules();
    report_density();
    report_free();
    report_refutations();
    report_depth();
    println!("run seconds {:.2}", whole.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    type Pinned = (&'static str, fn(u64) -> u64, u64, u64, u32, u32);

    #[test]
    fn the_carry_word_is_the_defect() {
        for n in 0..(1u64 << 16) {
            assert_eq!(carries(n, 20), defect(n));
        }
    }

    #[test]
    fn the_substitution_clears_the_constant() {
        for n in 0..(1u64 << 16) {
            let m = 2 * n + 1;
            assert_eq!(m ^ (2 * m), 2 * free(n) + 1);
        }
    }

    #[test]
    fn the_skeleton_is_rule_sixty() {
        for n in 0..(1u64 << 16) {
            let m = 2 * n + 1;
            let mut sixty = 0u64;
            for i in 0..20 {
                let hi = (m >> i) & 1;
                let lo = if i == 0 { 0 } else { (m >> (i - 1)) & 1 };
                sixty |= (hi ^ lo) << i;
            }
            assert_eq!(sixty, m ^ (2 * m));
        }
    }

    #[test]
    fn one_digit_reaches_every_digit() {
        for length in 2..=40usize {
            let mut u = 1u64;
            for j in 1..length {
                if j % 2 == 0 {
                    u |= 1 << j;
                }
            }
            let diff = (3 * u + 1) ^ (3 * (u ^ 1) + 1);
            assert!((2..=length).all(|i| (diff >> i) & 1 == 1));
        }
    }

    #[test]
    fn the_zero_carry_codes_are_pinned() {
        let pinned = [
            (3u64, 2usize, 7u64),
            (5, 3, 95),
            (7, 3, 23),
            (9, 4, 22015),
            (11, 4, 279),
            (15, 4, 279),
        ];
        for (q, width, code) in pinned {
            assert_eq!(zero_code(&support(q)), (width, code));
        }
    }

    #[test]
    fn the_zero_carry_set_is_the_rule() {
        for q in [3u64, 5, 7, 9, 11, 15] {
            let supp = support(q);
            let (width, code) = zero_code(&supp);
            let rule = Rule::new(1, width, code).unwrap();
            for n in 0..(1u64 << 14) {
                let carry_free = (q * n + 1) ^ skeleton(n, &supp) == 0;
                assert_eq!(carry_free, n % 2 == 0 && rule.accepts(&word(n, 14 + width)));
            }
        }
    }

    #[test]
    fn the_two_named_couplings_come_back() {
        let golden = Rule::new(1, 2, 7).unwrap();
        let supergolden = Rule::new(1, 3, 23).unwrap();
        assert!((kappa(&golden) - 0.098_239).abs() < 5e-7);
        assert!((kappa(&supergolden) - 0.115_204).abs() < 5e-7);
        assert!((perron(&golden) - 1.618_033_988_749_895).abs() < 1e-9);
        assert!((perron(&supergolden) - 1.465_571_231_876_768).abs() < 1e-9);
    }

    #[test]
    fn the_stationary_vector_balances() {
        let pi = [2i64, 1, 1, 2];
        assert_eq!(2 * pi[0], pi[0] + pi[1] + pi[2]);
        assert_eq!(2 * pi[1], pi[3]);
        assert_eq!(2 * pi[2], pi[0]);
        assert_eq!(2 * pi[3], pi[1] + pi[2] + pi[3]);
        assert_eq!(2 * (pi[1] + pi[3]), pi.iter().sum::<i64>());
    }

    #[test]
    fn the_carry_density_is_one_half() {
        for l in [16usize, 18, 20] {
            let limit = 1u64 << l;
            let mut total = 0u128;
            for n in 0..limit {
                total += u128::from(d_loc(n));
            }
            let sign = if l % 2 == 0 { 1i128 } else { -1 };
            let closed = 3 * l as i128 * (limit / 2) as i128 + limit as i128 - sign;
            assert_eq!(3 * total as i128, closed);
        }
    }

    #[test]
    fn the_carry_free_map_has_one_cycle() {
        for n in 1..(1u64 << 18) {
            assert!(digits(t_free(n)) <= digits(n));
            let mut x = n;
            let mut steps = 0;
            while x != 1 && steps < 4096 {
                x = t_free(x);
                steps += 1;
            }
            assert_eq!(x, 1);
        }
    }

    #[test]
    fn the_defect_is_no_function_of_its_statistics() {
        let pinned: [Pinned; 6] = [
            ("popcount", |n| u64::from(n.count_ones()), 1, 2, 2, 0),
            ("longest_run", |n| u64::from(longest_run(n)), 1, 2, 2, 0),
            ("v2", |n| u64::from(v2(n)), 1, 3, 2, 3),
            ("digit_count", |n| u64::from(digits(n)), 2, 3, 0, 3),
            (
                "popcount_and_v2",
                |n| u64::from(n.count_ones()) * 128 + u64::from(v2(n)),
                3,
                5,
                3,
                4,
            ),
            (
                "all_four",
                |n| {
                    ((u64::from(n.count_ones()) * 128 + u64::from(longest_run(n))) * 128
                        + u64::from(v2(n)))
                        * 128
                        + u64::from(digits(n))
                },
                19,
                25,
                3,
                4,
            ),
        ];
        for (_, key, m, n, dm, dn) in pinned {
            assert_eq!(key(m), key(n));
            assert_eq!(d_loc(m), dm);
            assert_eq!(d_loc(n), dn);
            assert_ne!(dm, dn);
        }
    }

    #[test]
    fn the_worst_depth_is_the_digit_count_plus_one() {
        for l in [1usize, 2, 4, 8, 16, 32, 40] {
            assert_eq!(depth_of(&vec![1u8; l]), l + 1);
        }
    }
}
