use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "research-erdos";
const BASES: [u64; 3] = [3, 4, 5];
const LEAST: u32 = 3;
const LIMIT: u64 = 10_000_000;
const TERMS: usize = 29;
const TOTAL: usize = 24_973_824;
const FROBENIUS: usize = 4_330_731;
const MISSING: usize = 45_704;
const WIDTH: usize = 860;

fn powers() -> Vec<usize> {
    let mut terms: Vec<usize> = BASES
        .iter()
        .flat_map(|&b| (LEAST..).map(move |j| b.pow(j)).take_while(|&p| p <= LIMIT))
        .map(|p| p as usize)
        .collect();
    terms.sort_unstable();
    terms
}

fn shift_or(bits: &mut [u64], a: usize, top: usize) {
    let (q, r) = (a / 64, a % 64);
    for i in (q..=top / 64).rev() {
        let mut word = bits[i - q] << r;
        if r > 0 && i > q {
            word |= bits[i - q - 1] >> (64 - r);
        }
        bits[i] |= word;
    }
}

fn has(bits: &[u64], k: usize) -> bool {
    bits[k / 64] >> (k % 64) & 1 == 1
}

fn count(bits: &[u64], lo: usize, hi: usize) -> usize {
    if lo >= hi {
        return 0;
    }
    let (a, b) = (lo / 64, (hi - 1) / 64);
    let head = !0u64 << (lo % 64);
    let tail = !0u64 >> (63 - (hi - 1) % 64);
    if a == b {
        return (bits[a] & head & tail).count_ones() as usize;
    }
    let inner: u32 = bits[a + 1..b].iter().map(|w| w.count_ones()).sum();
    (inner + (bits[a] & head).count_ones() + (bits[b] & tail).count_ones()) as usize
}

fn bins(bits: &[u64], total: usize) -> Vec<(usize, usize)> {
    let span = total + 1;
    (0..WIDTH / 2)
        .map(|c| {
            let (lo, hi) = if span >= WIDTH {
                (c * span / WIDTH, ((c + 1) * span).div_ceil(WIDTH))
            } else {
                let k = (2 * c + 1) * span / (2 * WIDTH);
                (k, k + 1)
            };
            (count(bits, lo, hi), hi - lo)
        })
        .collect()
}

fn main() -> Result<()> {
    let terms = powers();
    assert_eq!(terms.len(), TERMS);
    assert_eq!((terms[0], terms[TERMS - 1]), (27, 9_765_625));
    let mut bits = vec![0u64; TOTAL / 64 + 2];
    bits[0] = 1;
    let mut total = 0;
    let mut counts = Vec::new();
    let mut widths = Vec::new();
    for &a in &terms {
        total += a;
        shift_or(&mut bits, a, total);
        assert!(has(&bits, 0) && has(&bits, total));
        let (hit, wide): (Vec<usize>, Vec<usize>) = bins(&bits, total).into_iter().unzip();
        counts.push(hit);
        widths.push(wide);
    }
    assert_eq!(total, TOTAL);
    assert!((0..=TOTAL / 2).all(|k| has(&bits, k) == has(&bits, TOTAL - k)));
    let hole = (0..=TOTAL / 2).rev().find(|&k| !has(&bits, k));
    assert_eq!(hole, Some(FROBENIUS));
    assert_eq!(FROBENIUS - count(&bits, 1, FROBENIUS + 1), MISSING);
    assert_eq!(
        count(&bits, FROBENIUS + 1, TOTAL - FROBENIUS),
        TOTAL - 2 * FROBENIUS - 1
    );
    save(NAME, &json!({"counts": counts, "widths": widths}))
}
