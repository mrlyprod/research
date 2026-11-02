use crate::radii::{live, perron, v3, Live};
use mrly_num::factor::gcd;
use std::collections::HashSet;

const Y: u8 = 2;

pub struct Coded {
    pub live: Live,
    pub letters: Vec<Vec<u8>>,
    pub k: usize,
}

pub fn coded(a: u64, b: u64) -> Coded {
    let (a, b) = if a.is_multiple_of(3) { (b, a) } else { (a, b) };
    let live = live(a, b);
    let letters = live
        .edges
        .iter()
        .enumerate()
        .map(|(i, out)| {
            out.iter()
                .map(|&j| {
                    let eps = 3 * live.states[j] - live.states[i];
                    if eps == 0 {
                        0
                    } else if eps == b as i64 {
                        1
                    } else {
                        Y
                    }
                })
                .collect()
        })
        .collect();
    Coded {
        live,
        letters,
        k: v3(b) as usize,
    }
}

fn walks(c: &Coded, start: usize, len: usize, out: &mut Vec<Vec<u8>>) {
    let mut stack = vec![(start, Vec::new())];
    while let Some((s, word)) = stack.pop() {
        if word.len() == len {
            out.push(word);
            continue;
        }
        for (&j, &l) in c.live.edges[s].iter().zip(&c.letters[s]) {
            let mut w = word.clone();
            w.push(l);
            stack.push((j, w));
        }
    }
}

fn code(word: &[u8], k: usize, len: usize) -> u32 {
    let at = |i: usize| word.get(i).copied().unwrap_or(0);
    (1..=len)
        .filter(|&j| at(j + k - 1) == Y && at(j - 1) != Y)
        .fold(0u32, |s, j| s | 1 << (j - 1))
}

pub fn s_words(c: &Coded, len: usize) -> HashSet<u32> {
    let mut set = HashSet::new();
    let mut ws = Vec::new();
    for start in 0..c.live.states.len() {
        walks(c, start, len + c.k, &mut ws);
    }
    for w in &ws {
        set.insert(code(w, c.k, len));
    }
    set
}

pub fn free_words(k: usize, len: usize) -> Vec<u32> {
    (0u32..1 << len)
        .filter(|&w| len <= k || (w & (w >> k)) == 0)
        .collect()
}

fn contains(word: u32, len: usize, block: u32, blen: usize) -> bool {
    blen <= len && (0..=len - blen).any(|i| ((word >> i) & ((1 << blen) - 1)) == block)
}

pub fn missed(c: &Coded, max: usize) -> Vec<(u32, usize)> {
    let mut out = Vec::new();
    for len in 1..=max {
        let seen = s_words(c, len);
        for w in free_words(c.k, len) {
            if !seen.contains(&w) {
                out.push((w, len));
            }
        }
    }
    out
}

pub fn minimal(blocks: &[(u32, usize)]) -> Vec<(u32, usize)> {
    let mut out: Vec<(u32, usize)> = Vec::new();
    let mut sorted = blocks.to_vec();
    sorted.sort_by_key(|&(w, l)| (l, w));
    for (w, l) in sorted {
        if !out.iter().any(|&(b, bl)| contains(w, l, b, bl)) {
            out.push((w, l));
        }
    }
    out
}

pub fn show(w: u32, len: usize) -> String {
    (0..len)
        .map(|i| if (w >> i) & 1 == 1 { '1' } else { '0' })
        .collect()
}

pub fn closed(c: &Coded, n: usize) -> Vec<Vec<u8>> {
    let mut ws = Vec::new();
    walks(c, 0, n, &mut ws);
    let state_after = |w: &[u8]| {
        let mut s = 0usize;
        for &l in w {
            let pos = c.letters[s].iter().position(|&x| x == l).unwrap();
            s = c.live.edges[s][pos];
        }
        s
    };
    ws.into_iter().filter(|w| state_after(w) == 0).collect()
}

pub fn lemma(c: &Coded, top: usize) -> (bool, bool) {
    let mut injective = true;
    let mut free = true;
    for n in 1..=top {
        let walks = closed(c, n);
        let codes: HashSet<u32> = walks.iter().map(|w| code(w, c.k, n)).collect();
        injective &= codes.len() == walks.len();
        free &= codes.iter().all(|&w| n <= c.k || (w & (w >> c.k)) == 0);
    }
    (injective, free)
}

pub fn single(a: u64, b: u64, max: usize, probe: Option<Vec<u8>>) {
    let c = coded(a, b);
    let (mu, _, _) = perron(&c.live);
    println!(
        "ray {} {} k {} states {} mu {:.13}",
        a,
        b,
        c.k,
        c.live.states.len(),
        mu
    );
    for (i, out) in c.live.edges.iter().enumerate() {
        let arrows: Vec<String> = out
            .iter()
            .zip(&c.letters[i])
            .map(|(&j, &l)| format!("{}->{}", ['0', 'x', 'y'][l as usize], c.live.states[j]))
            .collect();
        println!("  {} : {}", c.live.states[i], arrows.join(" "));
    }
    let (inj, free) = lemma(&c, 12);
    println!(
        "coding lemma to level 12: injective on closed walks {}, distance-{}-free {}",
        inj, c.k, free
    );
    let mm = minimal(&missed(&c, max));
    let words: Vec<String> = mm.iter().map(|&(w, l)| show(w, l)).collect();
    println!(
        "minimal missed blocks to length {}: {}",
        max,
        if words.is_empty() {
            "none".to_string()
        } else {
            words.join(" ")
        }
    );
    match mm.first() {
        Some(&(_, l)) => println!(
            "every distance-{}-free block realised to length {}, the first miss at length {}",
            c.k,
            l - 1,
            l
        ),
        None => println!(
            "every distance-{}-free block realised to the cut {}, no miss",
            c.k, max
        ),
    }
    if let Some(block) = probe {
        println!(
            "probe {} length {}: {} by a path DP over the component",
            block
                .iter()
                .map(|&b| (b'0' + b) as char)
                .collect::<String>(),
            block.len(),
            if realised(&c, &block) {
                "realised"
            } else {
                "missed"
            }
        );
    }
}

pub fn census(height: u64, max: usize, probe: Option<Vec<u8>>) {
    let mut rows = Vec::new();
    let mut all_rays = 0usize;
    let mut lemma_ok = true;
    for b in 1..=height {
        for a in 1..=b {
            if gcd(a as u128, b as u128) != 1 {
                continue;
            }
            let c = coded(a, b);
            if c.live.states.len() == 1 {
                continue;
            }
            all_rays += 1;
            let (inj, free) = lemma(&c, 12);
            lemma_ok &= inj && free;
            if c.k != 2 {
                continue;
            }
            let (mu, _, _) = perron(&c.live);
            let blocks = missed(&c, max);
            rows.push((mu, a, b, c, blocks));
        }
    }
    rows.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    println!(
        "missed-block census, box 1 <= a <= b <= {}, S_j = [w_(j+k-1) = y][w_(j-1) != y] on the component of 0, blocks to length {}",
        height, max
    );
    println!(
        "coding lemma on all {} occupied rays of the box, every k, closed walks to level 12: injective and distance-k-free {}",
        all_rays, lemma_ok
    );
    println!("k = 2 occupied rays {}", rows.len());
    let mut shortest = vec![0usize; max + 2];
    let mut common: Option<HashSet<(u32, usize)>> = None;
    for (mu, a, b, c, blocks) in &rows {
        let mm = minimal(blocks);
        let shift = *a == 1 && *b == 9;
        if !shift {
            let set: HashSet<(u32, usize)> = blocks.iter().copied().collect();
            common = Some(match common {
                None => set,
                Some(c) => c.intersection(&set).copied().collect(),
            });
            shortest[mm.first().map(|&(_, l)| l).unwrap_or(max + 1)] += 1;
        }
        let words: Vec<String> = mm.iter().take(6).map(|&(w, l)| show(w, l)).collect();
        println!(
            "ray {} {} states {} mu {:.10} gap {:.6} minimal missed {}{}",
            a,
            b,
            c.live.states.len(),
            mu,
            phi - mu,
            if words.is_empty() {
                "none".to_string()
            } else {
                words.join(" ")
            },
            if shift { " SHIFT" } else { "" }
        );
    }
    println!(
        "off-shift rays by length of the shortest missed block: {}",
        (1..=max + 1)
            .filter(|&l| shortest[l] > 0)
            .map(|l| format!(
                "{} {}",
                if l > max {
                    "none".to_string()
                } else {
                    l.to_string()
                },
                shortest[l]
            ))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let common: Vec<(u32, usize)> = common.unwrap_or_default().into_iter().collect();
    let mm = minimal(&common);
    println!(
        "blocks missed by every off-shift k = 2 ray: {} (minimal: {})",
        common.len(),
        if mm.is_empty() {
            "none".to_string()
        } else {
            mm.iter()
                .map(|&(w, l)| show(w, l))
                .collect::<Vec<_>>()
                .join(" ")
        }
    );
    let off: Vec<&Coded> = rows
        .iter()
        .filter(|(_, a, b, _, _)| !(*a == 1 && *b == 9))
        .map(|(_, _, _, c, _)| c)
        .collect();
    let sets: Vec<HashSet<(u32, usize)>> = rows
        .iter()
        .filter(|(_, a, b, _, _)| !(*a == 1 && *b == 9))
        .map(|(_, _, _, _, blocks)| blocks.iter().copied().collect())
        .collect();
    let k = off.first().map(|c| c.k).unwrap_or(2);
    let chosen = cover(&sets);
    let word = join(&chosen, k);
    let missing = off.iter().filter(|c| !realised(c, &word)).count();
    println!(
        "greedy cover of the minimal missed blocks joined by 0^k: {} length {} from {}, missed by {} of {} off-shift rays by a path DP over each component",
        word.iter().map(|&b| (b'0' + b) as char).collect::<String>(),
        word.len(),
        chosen.iter().map(|&(w, l)| show(w, l)).collect::<Vec<_>>().join(" "),
        missing,
        off.len()
    );
    if let Some(block) = probe {
        let missing = off.iter().filter(|c| !realised(c, &block)).count();
        println!(
            "probe {} length {}: missed by {} of {} off-shift rays by a path DP over each component",
            block.iter().map(|&b| (b'0' + b) as char).collect::<String>(),
            block.len(),
            missing,
            off.len()
        );
    }
}

pub fn realised(c: &Coded, block: &[u8]) -> bool {
    let k = c.k;
    let mut front: HashSet<(usize, u32)> = (0..c.live.states.len()).map(|s| (s, 0)).collect();
    for step in 0..k + block.len() {
        let mut next = HashSet::new();
        for &(s, hist) in &front {
            for (&j, &l) in c.live.edges[s].iter().zip(&c.letters[s]) {
                let is_y = (l == Y) as u32;
                if step >= k && (is_y & (1 - (hist & 1))) != block[step - k] as u32 {
                    continue;
                }
                next.insert((j, (hist >> 1) | (is_y << (k - 1))));
            }
        }
        front = next;
    }
    !front.is_empty()
}

fn cover(sets: &[HashSet<(u32, usize)>]) -> Vec<(u32, usize)> {
    let mut uncovered: Vec<usize> = (0..sets.len()).collect();
    let mut chosen = Vec::new();
    while !uncovered.is_empty() {
        let mut candidates: Vec<(u32, usize)> = uncovered
            .iter()
            .flat_map(|&i| minimal(&sets[i].iter().copied().collect::<Vec<_>>()))
            .collect();
        candidates.sort_by_key(|&(w, l)| (l, w));
        candidates.dedup();
        let best = candidates
            .iter()
            .copied()
            .max_by_key(|&(w, l)| {
                (
                    uncovered
                        .iter()
                        .filter(|&&i| sets[i].contains(&(w, l)))
                        .count(),
                    std::cmp::Reverse((l, w)),
                )
            })
            .unwrap();
        uncovered.retain(|&i| !sets[i].contains(&best));
        chosen.push(best);
    }
    chosen
}

fn join(blocks: &[(u32, usize)], k: usize) -> Vec<u8> {
    let mut word = Vec::new();
    for (n, &(w, l)) in blocks.iter().enumerate() {
        if n > 0 {
            word.extend(std::iter::repeat_n(0u8, k));
        }
        word.extend((0..l).map(|i| ((w >> i) & 1) as u8));
    }
    word
}
