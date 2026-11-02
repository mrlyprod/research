use crate::census::Entry;
use crate::sources::KMAX;

fn find<'a>(entries: &'a [Entry], name: &str) -> &'a Entry {
    entries.iter().find(|e| e.name == name).unwrap()
}

pub fn report(entries: &[Entry]) {
    println!("LAWS");
    let odds = find(entries, "odds");
    let mut want = vec![0u8; KMAX + 1];
    want[1] = 1;
    assert_eq!(odds.c, want, "odds is not the parity");
    println!("odds: c(w) = [w = 1] for w <= {KMAX}, the mess is the parity of the count");
    let evens = find(entries, "evens");
    want[0] = 1;
    assert_eq!(evens.c, want, "evens is not the complement of the parity");
    println!("evens: c(w) = [w <= 1] for w <= {KMAX}, the mess is the complement of the parity");
    let binary = find(entries, "binary");
    for w in 0..=KMAX {
        assert_eq!(
            binary.c[w] as u32,
            if w == 0 { 0 } else { w.count_ones() % 2 },
            "binary at {w}"
        );
    }
    println!("binary: c(w) = popcount(w) mod 2 for 1 <= w <= {KMAX}, c(0) = 0");
    let mut runs = 0;
    for k in 1..=KMAX {
        let d = binary.rows[k - 1].degree as usize;
        assert!(
            d == k || d == k - 1 || d == k - 2 && k >= 2,
            "binary degree at {k}"
        );
        runs += (d == k) as usize;
    }
    println!("binary: degree(k) is k, k-1 or k-2 for every k <= {KMAX}, and k at {runs} of them");
    let primes = find(entries, "primes");
    for w in (2..=KMAX).step_by(2) {
        assert_eq!(primes.c[w], (w % 4 == 2) as u8, "primes at {w}");
    }
    println!("primes: c(w) = [w = 2 mod 4] for even w <= {KMAX}");
    let rest = entries
        .iter()
        .filter(|e| !["odds", "evens"].contains(&e.name.as_str()))
        .filter(|e| e.rows.iter().all(|p| p.degree <= 1))
        .count();
    println!("sources of degree at most 1 at every k other than odds and evens: {rest}");
    let bent = entries
        .iter()
        .flat_map(|e| e.rows.iter().map(move |p| (e, p)))
        .filter(|(_, p)| p.k % 2 == 0 && p.nonlinearity == crate::mess::ceiling(p.k))
        .map(|(e, p)| format!("{}@{}", e.name, p.k))
        .collect::<Vec<_>>();
    println!(
        "sources at the ceiling at even k: {}",
        if bent.is_empty() {
            "none".to_string()
        } else {
            bent.join(" ")
        }
    );
    let odd_top = entries
        .iter()
        .flat_map(|e| e.rows.iter().map(move |p| (e, p)))
        .filter(|(_, p)| p.k % 2 == 1 && p.nonlinearity == crate::mess::ceiling(p.k))
        .map(|(e, p)| format!("{}@{}", e.name, p.k))
        .collect::<Vec<_>>();
    println!(
        "sources at the ceiling at odd k: {}",
        if odd_top.is_empty() {
            "none".to_string()
        } else {
            odd_top.join(" ")
        }
    );
}
