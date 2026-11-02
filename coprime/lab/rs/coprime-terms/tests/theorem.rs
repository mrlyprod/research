fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let rest = a % b;
        a = b;
        b = rest;
    }
    a
}

fn points(q: u64, corners: &[Vec<u64>], level: u32) -> Vec<Vec<u64>> {
    let width = corners[0].len();
    let mut pts = vec![vec![0u64; width]];
    for _ in 0..level {
        let mut next = Vec::with_capacity(pts.len() * corners.len());
        for p in &pts {
            for v in corners {
                let mut x = Vec::with_capacity(width);
                for i in 0..width {
                    x.push(p[i] * q + v[i]);
                }
                next.push(x);
            }
        }
        pts = next;
    }
    pts
}

fn common(p: &[u64]) -> u64 {
    let mut g = 0;
    for c in p {
        g = gcd(g, *c);
    }
    g
}

fn grid(q: u64, width: usize, keep: impl Fn(&[u64]) -> bool) -> Vec<Vec<u64>> {
    let mut out = Vec::new();
    let count = (q as usize).pow(width as u32);
    for code in 0..count {
        let mut v = Vec::with_capacity(width);
        let mut rest = code as u64;
        for _ in 0..width {
            v.push(rest % q);
            rest /= q;
        }
        if keep(&v) {
            out.push(v);
        }
    }
    out
}

#[test]
fn bracket_identity() {
    let code: u64 = 34376528265;
    let mut base6 = Vec::new();
    let mut bit = 0;
    for a in 0..6u64 {
        for b in 0..6u64 {
            if (code >> bit) & 1 == 1 {
                base6.push(vec![b, a]);
            }
            bit += 1;
        }
    }
    assert_eq!(base6.len(), 8);
    for level in 1..=6u32 {
        let q6 = points(6, &base6, level)
            .iter()
            .filter(|p| gcd(common(p), 6) == 1)
            .count() as u64;
        assert_eq!(q6, 8u64.pow(level) / 2);
    }
    let menger = grid(3, 3, |v| v.iter().filter(|d| **d == 1).count() <= 1);
    for level in 1..=4u32 {
        let q3 = points(3, &menger, level)
            .iter()
            .filter(|p| !common(p).is_multiple_of(3))
            .count() as u64;
        assert_eq!(q3, 19 * 20u64.pow(level - 1));
    }
    let carpet = grid(3, 2, |v| v != [1, 1]);
    for level in 1..=6u32 {
        let q3 = points(3, &carpet, level)
            .iter()
            .filter(|p| !common(p).is_multiple_of(3))
            .count() as u64;
        assert_eq!(q3, 7 * 8u64.pow(level - 1));
    }
    let vicsek = grid(3, 2, |v| v.contains(&1));
    for level in 1..=7u32 {
        let q3 = points(3, &vicsek, level)
            .iter()
            .filter(|p| !common(p).is_multiple_of(3))
            .count() as u64;
        assert_eq!(q3, 5u64.pow(level));
    }
}

#[test]
fn box_bound() {
    let mut designs: Vec<(u64, Vec<Vec<u64>>)> = vec![
        (2, vec![vec![0, 0], vec![0, 1], vec![1, 0]]),
        (2, vec![vec![0, 1], vec![1, 0], vec![1, 1]]),
        (3, grid(3, 2, |v| v != [1, 1])),
        (3, grid(3, 2, |v| v.contains(&1))),
        (
            3,
            grid(3, 3, |v| v.iter().filter(|d| **d == 1).count() <= 1),
        ),
        (5, vec![vec![0, 0], vec![1, 0], vec![0, 1]]),
        (9, vec![vec![0, 0], vec![1, 0], vec![0, 1]]),
        (
            6,
            vec![
                vec![0, 0],
                vec![1, 1],
                vec![2, 3],
                vec![3, 2],
                vec![4, 5],
                vec![5, 4],
                vec![1, 3],
                vec![3, 1],
            ],
        ),
    ];
    let mut state: u64 = 7;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 33
    };
    while designs.len() < 48 {
        let q = 2 + next() % 9;
        let width = if next() % 3 == 2 { 3 } else { 2 };
        let cap = (q.pow(width as u32)).min(8);
        let k = 2 + next() % (cap - 1).max(1);
        let mut f: Vec<Vec<u64>> = Vec::new();
        while (f.len() as u64) < k {
            let v: Vec<u64> = (0..width).map(|_| next() % q).collect();
            if !f.contains(&v) {
                f.push(v);
            }
        }
        designs.push((q, f));
    }
    let mut checked = 0u64;
    for (q, f) in &designs {
        let k = f.len() as u64;
        let width = f[0].len() as u32;
        let mut level = 1u32;
        while k.pow(level + 1) <= 20000 && (*q as u128).pow(level + 1) < 1u128 << 60 {
            level += 1;
        }
        let pts = points(*q, f, level);
        let top = (*q as u128).pow(level).min(200) as u64;
        for m in 2..=top {
            let hits = pts.iter().filter(|p| p.iter().all(|c| c % m == 0)).count() as u128;
            let mut bound = u128::MAX;
            for h in 0..=level {
                let slab = (k as u128)
                    .saturating_pow(level - h)
                    .saturating_mul(((*q as u128).pow(h) / m as u128 + 1).saturating_pow(width));
                bound = bound.min(slab);
            }
            assert!(hits <= bound, "box bound violated at q={} m={}", q, m);
            checked += 1;
        }
    }
    assert!(checked > 5000);
}

#[test]
fn cantor_dust() {
    let dust = vec![vec![0, 0], vec![0, 2], vec![2, 0], vec![2, 2]];
    for level in 1..=8u32 {
        let pts = points(3, &dust, level);
        assert_eq!(pts.iter().filter(|p| common(p) == 1).count(), 0);
        if level == 8 {
            assert_eq!(pts.iter().filter(|p| common(p) == 2).count(), 33883);
        }
    }
}
