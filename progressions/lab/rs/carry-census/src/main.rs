use std::time::Instant;

// THE CRITERION

fn digits(mask: u64) -> Vec<u32> {
    (0..64).filter(|d| mask >> d & 1 == 1).collect()
}

fn sums(mask: u64) -> u128 {
    let ds = digits(mask);
    let mut ff = 0u128;
    for &f in &ds {
        for &g in &ds {
            ff |= 1 << (f + g);
        }
    }
    ff
}

fn signed(t: i64, mask: u64, ff: u128) -> bool {
    digits(mask).iter().any(|&h| {
        let v = t + 2 * h as i64;
        (0..128).contains(&v) && ff >> v & 1 == 1
    })
}

fn digit_free(mask: u64) -> bool {
    let ds = digits(mask);
    ds.iter().all(|&g| {
        ds.iter()
            .all(|&f| f >= g || 2 * g - f >= 64 || mask >> (2 * g - f) & 1 == 0)
    })
}

fn exits(base: u32, mask: u64) -> (bool, bool) {
    let ff = sums(mask);
    let b = base as i64;
    (signed(b, mask, ff), signed(-b, mask, ff))
}

fn closes(base: u32, mask: u64) -> bool {
    let ff = sums(mask);
    let b = base as i64;
    let s = |t: i64| signed(t, mask, ff);
    (s(b) && (s(-1) || (s(-b - 1) && s(1)))) || (s(-b) && (s(1) || (s(b + 1) && s(-1))))
}

fn free(base: u32, mask: u64) -> bool {
    digit_free(mask) && !closes(base, mask)
}

// THE DIGIT SETS

fn extend(base: u32, mask: u64, from: u32, out: &mut Vec<u64>) {
    out.push(mask);
    for d in from..base {
        let closes_digit = digits(mask)
            .iter()
            .any(|&g| 2 * g >= d && mask >> (2 * g - d) & 1 == 1);
        if !closes_digit {
            extend(base, mask | 1 << d, d + 1, out);
        }
    }
}

fn families(base: u32) -> Vec<u64> {
    let mut out = Vec::new();
    extend(base, 1, 1, &mut out);
    out
}

fn show(mask: u64) -> String {
    let words: Vec<String> = digits(mask).iter().map(|d| d.to_string()).collect();
    format!("{{{}}}", words.join(", "))
}

// BRUTE FORCE

fn level(base: u32, mask: u64, lvl: u32) -> Vec<u64> {
    let ds = digits(mask);
    let mut out = vec![0u64];
    for _ in 0..lvl {
        let mut next = Vec::with_capacity(out.len() * ds.len());
        for &x in &out {
            for &d in &ds {
                next.push(x * base as u64 + d as u64);
            }
        }
        out = next;
    }
    out.sort_unstable();
    out
}

fn brute(base: u32, mask: u64, lvl: u32) -> Option<(u64, u64, u64)> {
    let set = level(base, mask, lvl);
    let size = (base as u64).pow(lvl) as usize;
    let mut bits = vec![0u64; size.div_ceil(64)];
    for &x in &set {
        bits[x as usize / 64] |= 1 << (x % 64);
    }
    for (i, &x) in set.iter().enumerate() {
        for &z in &set[i + 1..] {
            if (x + z) % 2 == 0 {
                let y = (x + z) / 2;
                if bits[y as usize / 64] >> (y % 64) & 1 == 1 {
                    return Some((x, y, z));
                }
            }
        }
    }
    None
}

fn first_level(base: u32, mask: u64, top: u32) -> Option<u32> {
    (1..=top).find(|&lvl| brute(base, mask, lvl).is_some())
}

// THE CYCLIC GROUP

fn cyclic_extend(n: u32, set: &mut Vec<u32>, from: u32, best: &mut u32) {
    *best = (*best).max(set.len() as u32);
    for d in from..n {
        let closes_digit = set
            .iter()
            .any(|&g| set.iter().any(|&f| (f + g) % n == (2 * d) % n))
            || set.iter().any(|&h| (2 * h) % n == (2 * d) % n && h != d)
            || set
                .iter()
                .any(|&f| set.iter().any(|&h| (f + d) % n == (2 * h) % n));
        if !closes_digit {
            set.push(d);
            cyclic_extend(n, set, d + 1, best);
            set.pop();
        }
    }
}

fn cyclic(n: u32) -> u32 {
    let mut best = 0;
    cyclic_extend(n, &mut vec![0], 1, &mut best);
    best
}

// THE RUNS

fn control(bases: std::ops::RangeInclusive<u32>, top: u32) -> u32 {
    println!("control: brute force to level {top} against the criterion, every progression-free F with 0");
    println!(
        "| base | digit sets | with a progression | deepest first level | mismatches | seconds |"
    );
    println!("|---|---|---|---|---|---|");
    let mut total = 0;
    for base in bases {
        let start = Instant::now();
        let fams = families(base);
        let mut with = 0;
        let mut deepest = 0;
        let mut mismatches = 0;
        for &m in &fams {
            let first = first_level(base, m, top);
            if let Some(lvl) = first {
                with += 1;
                deepest = deepest.max(lvl);
            }
            if free(base, m) == first.is_some() {
                mismatches += 1;
                println!(
                    "mismatch at base {base}, F = {}: criterion free {}, brute first level {:?}",
                    show(m),
                    free(base, m),
                    first
                );
            }
        }
        total += mismatches;
        println!(
            "| {base} | {} | {with} | {deepest} | {mismatches} | {:.2} |",
            fams.len(),
            start.elapsed().as_secs_f64()
        );
    }
    total
}

fn verify(bases: std::ops::RangeInclusive<u32>) -> u32 {
    println!(
        "verify: brute force at level 3 against the criterion on every progression-free F with 0"
    );
    println!("| base | digit sets | with a progression | mismatches | seconds |");
    println!("|---|---|---|---|---|");
    let mut total = 0;
    let mut sets = 0;
    for base in bases {
        let start = Instant::now();
        let fams = families(base);
        sets += fams.len();
        let mut with = 0;
        let mut mismatches = 0;
        for &m in &fams {
            let found = brute(base, m, 3).is_some();
            with += u32::from(found);
            if free(base, m) == found {
                mismatches += 1;
                println!("mismatch at base {base}, F = {}", show(m));
            }
        }
        total += mismatches;
        println!(
            "| {base} | {} | {with} | {mismatches} | {:.2} |",
            fams.len(),
            start.elapsed().as_secs_f64()
        );
    }
    println!("verify: {sets} digit sets");
    total
}

fn census(bases: std::ops::RangeInclusive<u32>, checked: u32) -> Vec<u32> {
    let ratio = 2f64.ln() / 3f64.ln();
    println!("census: the largest progression-free design per base, every F with 0 by the criterion; maximisers above base {checked} brute-forced at level 3");
    println!("| base | half-base | cyclic | no exit | with an exit | best | maximisers | lift maximisers | least maximiser | base^(log 2/log 3) |");
    println!("|---|---|---|---|---|---|---|---|---|---|");
    let mut bests = Vec::new();
    let mut crossed = Vec::new();
    for base in bases {
        let fams = families(base);
        let half = base.div_ceil(2);
        let (mut half_best, mut noexit_best, mut exit_best, mut best) = (0u32, 0u32, 0u32, 0u32);
        let mut maximisers = Vec::new();
        for &m in &fams {
            let fill = m.count_ones();
            if m < 1 << half {
                assert!(
                    free(base, m),
                    "half-base F = {} at base {base} must be free",
                    show(m)
                );
                half_best = half_best.max(fill);
            }
            if !free(base, m) {
                continue;
            }
            let (up, down) = exits(base, m);
            if up || down {
                exit_best = exit_best.max(fill);
            } else {
                noexit_best = noexit_best.max(fill);
            }
            if fill > best {
                best = fill;
                maximisers = vec![m];
            } else if fill == best {
                maximisers.push(m);
            }
        }
        let bound = (base as f64).powf(ratio);
        if best as f64 > bound + 1e-9 {
            crossed.push(base);
        }
        if base > checked {
            for &m in &maximisers {
                assert_eq!(
                    brute(base, m, 3),
                    None,
                    "base {base}: the maximiser {} holds a progression at level 3",
                    show(m)
                );
            }
        }
        let cyc = cyclic(base);
        assert_eq!(
            cyc, noexit_best,
            "base {base}: the cyclic maximum and the no-exit maximum differ"
        );
        let lifts = maximisers
            .iter()
            .filter(|&&m| exits(base, m) == (false, false))
            .count();
        assert!(lifts > 0, "base {base}: no maximiser is a cyclic lift");
        println!("| {base} | {half_best} | {cyc} | {noexit_best} | {exit_best} | {best} | {} | {lifts} | {} | {bound:.3} |", maximisers.len(), show(maximisers[0]));
        if best > half_best {
            for &m in &maximisers {
                let (up, down) = exits(base, m);
                println!("  base {base} beats the half-base construction with {} ({} digits), exits up {up} down {down}", show(m), m.count_ones());
            }
        }
        bests.push(best);
    }
    let words: Vec<String> = crossed.iter().map(|b| b.to_string()).collect();
    println!(
        "bases where the best design beats dimension log 2/log 3: {}",
        if words.is_empty() {
            "none".to_string()
        } else {
            words.join(", ")
        }
    );
    bests
}

fn main() {
    let start = Instant::now();
    let mismatches = control(3..=12, 5);
    println!("control mismatches: {mismatches}");
    assert_eq!(mismatches, 0);
    let mismatches = verify(2..=27);
    println!("verify mismatches: {mismatches}");
    assert_eq!(mismatches, 0);
    let bests = census(2..=40, 27);
    let words: Vec<String> = bests.iter().map(|b| b.to_string()).collect();
    println!("best fill at bases 2..40: {}", words.join(", "));
    println!("total {:.2} seconds", start.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mask(ds: &[u32]) -> u64 {
        ds.iter().fold(0, |m, &d| m | 1 << d)
    }

    #[test]
    fn signed_set_is_not_symmetric() {
        let m = mask(&[0, 1, 3]);
        let ff = sums(m);
        assert!(signed(-5, m, ff));
        assert!(!signed(5, m, ff));
    }

    #[test]
    fn stanley_design_is_free() {
        assert!(free(3, mask(&[0, 1])));
        assert_eq!(brute(3, mask(&[0, 1]), 8), None);
    }

    #[test]
    fn full_digit_set_has_a_progression_in_f() {
        assert!(!digit_free(mask(&[0, 1, 2])));
        assert!(!free(3, mask(&[0, 1, 2])));
    }

    #[test]
    fn base_five_closes_by_a_carry() {
        let m = mask(&[0, 1, 3]);
        assert!(digit_free(m));
        assert!(!free(5, m));
        assert_eq!(brute(5, m, 2), Some((0, 3, 6)));
        assert!(level(5, m, 2).contains(&5) && level(5, m, 2).contains(&1));
    }

    #[test]
    fn even_digits_exit_without_closing() {
        let m = mask(&[0, 2]);
        assert_eq!(exits(4, m), (true, true));
        assert!(free(4, m));
        assert_eq!(brute(4, m, 6), None);
    }

    #[test]
    fn base_nineteen_record_is_free() {
        let m = mask(&[0, 1, 3, 12, 14, 15]);
        assert!(free(19, m));
        assert_eq!(exits(19, m), (false, false));
        assert_eq!(brute(19, m, 3), None);
    }

    #[test]
    fn levels_three_and_five_agree_with_the_criterion() {
        for base in 3..=9 {
            for m in families(base) {
                assert_eq!(
                    free(base, m),
                    brute(base, m, 3).is_none(),
                    "base {base} F {}",
                    show(m)
                );
            }
        }
        for base in 3..=6 {
            for m in families(base) {
                assert_eq!(
                    free(base, m),
                    brute(base, m, 5).is_none(),
                    "base {base} F {}",
                    show(m)
                );
            }
        }
    }

    #[test]
    fn sums_past_bit_63_are_kept() {
        let m = mask(&[0, 34]);
        assert_eq!(sums(m) >> 68 & 1, 1);
        assert!(signed(68, m, sums(m)));
        assert!(signed(-68, m, sums(m)));
        assert!(!signed(67, m, sums(m)));
    }

    #[test]
    fn base_thirty_five_record_is_nine() {
        let m = mask(&[0, 1, 3, 4, 10, 12, 22, 26, 28]);
        assert!(free(35, m));
        assert_eq!(exits(35, m), (false, false));
        assert_eq!(brute(35, m, 3), None);
        assert_eq!(cyclic(35), 9);
        let best = families(35)
            .into_iter()
            .filter(|&f| free(35, f))
            .map(|f| f.count_ones())
            .max();
        assert_eq!(best, Some(9));
    }
}
