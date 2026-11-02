use std::collections::HashMap;

const LAWS: u32 = 6;
const SHOWN: u32 = 20;

pub fn masses(a: u64, b: u64, top: u32) -> Vec<u128> {
    let mut dist: HashMap<(u64, u64), u128> = HashMap::from([((0, 0), 1)]);
    let mut out = vec![0u128];
    for _ in 1..=top {
        let mut next: HashMap<(u64, u64), u128> = HashMap::new();
        for (&(ca, cb), &w) in &dist {
            for d in 0..3 {
                let (sa, sb) = (ca + a * d, cb + b * d);
                let (ra, rb) = (sa % 3, sb % 3);
                if ra > 1 || rb > 1 || (ra == 1 && rb == 1) {
                    continue;
                }
                let slot = next.entry((sa / 3, sb / 3)).or_insert(0);
                *slot = slot.checked_add(w).expect("the ray mass fits a u128");
            }
        }
        dist = next;
        out.push(dist.get(&(0, 0)).copied().unwrap_or(0) - 1);
    }
    out
}

fn law(s: &[u128], k: usize) -> Option<usize> {
    (k + 1..s.len())
        .rev()
        .take_while(|&m| s[m] == s[m - 1] + s[m - k])
        .last()
}

pub fn report(a: u64, b: u64, top: u32) {
    let m = masses(a, b, top);
    let s: Vec<u128> = m.iter().map(|v| v + 1).collect();
    println!(
        "ray ({a},{b}) M_n to n = {top}, carry states at most {}",
        a * b
    );
    let head: Vec<String> = m[1..=SHOWN.min(top) as usize]
        .iter()
        .map(|v| v.to_string())
        .collect();
    println!("M_1..M_{} {}", SHOWN.min(top), head.join(", "));
    println!("M_{top} = {}", m[top as usize]);
    for k in 1..=LAWS as usize {
        if let Some(from) = law(&s, k) {
            println!("M_n + 1 = (M_(n-1) + 1) + (M_(n-{k}) + 1) from n = {from} to {top}");
        }
    }
    let squares: Vec<String> = (1..=top as usize)
        .filter(|&n| {
            let r = (m[n] as f64).sqrt().round() as u128;
            m[n] > 0 && r * r == m[n]
        })
        .map(|n| format!("n {n}: {}", m[n]))
        .collect();
    println!("perfect squares: {}", squares.join(", "));
}
