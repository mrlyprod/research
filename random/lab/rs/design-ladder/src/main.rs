use std::collections::{HashMap, HashSet};
use std::env;

mod reach;

// FRACTIONS

#[derive(Clone, Copy, PartialEq, Debug)]
struct Q {
    n: i128,
    d: i128,
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

impl Q {
    fn new(n: i128, d: i128) -> Q {
        let g = gcd(n, d).max(1);
        let s = if d < 0 { -1 } else { 1 };
        Q {
            n: s * n / g,
            d: s * d / g,
        }
    }
    fn int(n: i128) -> Q {
        Q { n, d: 1 }
    }
    fn add(self, o: Q) -> Q {
        Q::new(self.n * o.d + o.n * self.d, self.d * o.d)
    }
    fn mul(self, o: Q) -> Q {
        Q::new(self.n * o.n, self.d * o.d)
    }
    fn div(self, o: Q) -> Q {
        Q::new(self.n * o.d, self.d * o.n)
    }
    fn abs(self) -> Q {
        Q {
            n: self.n.abs(),
            d: self.d,
        }
    }
    fn lt(self, o: Q) -> bool {
        self.n * o.d < o.n * self.d
    }
    fn f(self) -> f64 {
        self.n as f64 / self.d as f64
    }
    fn show(self) -> String {
        if self.d == 1 {
            format!("{}", self.n)
        } else {
            format!("{}/{}", self.n, self.d)
        }
    }
}

// SEQUENCES

fn digit_count(mut n: u64, b: u64, d: u64) -> u32 {
    let mut c = 0;
    while n > 0 {
        if n % b == d {
            c += 1;
        }
        n /= b;
    }
    c
}

fn tm(n: u64) -> i8 {
    if n.count_ones() % 2 == 0 {
        1
    } else {
        -1
    }
}

fn rs(n: u64) -> i8 {
    if (n & (n >> 1)).count_ones() % 2 == 0 {
        1
    } else {
        -1
    }
}

fn t3(n: u64) -> i8 {
    if digit_count(n, 3, 2) % 2 == 0 {
        1
    } else {
        -1
    }
}

fn signed(v: &[i8], b: u64, mut n: u64, level: u32) -> i8 {
    let mut s = 1i8;
    for _ in 0..level {
        s *= v[(n % b) as usize];
        n /= b;
    }
    s
}

fn restricted(mask: u32, b: u64, mut n: u64, level: u32) -> bool {
    for _ in 0..level {
        if mask >> (n % b) & 1 == 0 {
            return false;
        }
        n /= b;
    }
    true
}

fn cantor(mut n: u64) -> bool {
    while n > 0 {
        if n % 3 == 1 {
            return false;
        }
        n /= 3;
    }
    true
}

fn carpet_morton(mut n: u64) -> bool {
    while n > 0 {
        if n % 9 == 4 {
            return false;
        }
        n /= 9;
    }
    true
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn sign(bit: bool) -> i8 {
    if bit {
        -1
    } else {
        1
    }
}

// CORRELATIONS

fn corr_exact(vp: &[i8], b: usize, h_max: usize) -> Vec<Q> {
    let al: Vec<i128> = (0..b)
        .map(|s| (0..b - s).map(|r| (vp[r] * vp[r + s]) as i128).sum())
        .collect();
    let bq = Q::int(b as i128);
    let mut c = vec![Q::int(0); h_max + 2];
    c[0] = Q::int(1);
    c[1] = Q::int(al[1]).div(Q::int(b as i128 - al[b - 1]));
    for s in 2..b {
        c[s] = Q::int(al[s]).add(Q::int(al[b - s]).mul(c[1])).div(bq);
    }
    for h in b..h_max + 2 {
        let (q, s) = (h / b, h % b);
        c[h] = if s == 0 {
            c[q]
        } else {
            Q::int(al[s])
                .mul(c[q])
                .add(Q::int(al[b - s]).mul(c[q + 1]))
                .div(bq)
        };
    }
    c
}

fn corr_float(vp: &[i8], b: usize, h_max: usize) -> Vec<f64> {
    let exact = corr_exact(vp, b, b + 1);
    let al: Vec<f64> = (0..b)
        .map(|s| (0..b - s).map(|r| (vp[r] * vp[r + s]) as f64).sum())
        .collect();
    let mut c = vec![0.0; h_max + 2];
    for h in 0..b.min(h_max + 2) {
        c[h] = exact[h].f();
    }
    for h in b..h_max + 2 {
        let (q, s) = (h / b, h % b);
        c[h] = if s == 0 {
            c[q]
        } else {
            (al[s] * c[q] + al[b - s] * c[q + 1]) / b as f64
        };
    }
    c
}

fn pearson(pairs: impl Iterator<Item = (i8, i8)>) -> Option<f64> {
    let (mut n, mut sx, mut sy, mut sxy) = (0f64, 0f64, 0f64, 0f64);
    for (x, y) in pairs {
        n += 1.0;
        sx += x as f64;
        sy += y as f64;
        sxy += (x as f64) * (y as f64);
    }
    let (mx, my) = (sx / n, sy / n);
    let vx = 1.0 - mx * mx;
    let vy = 1.0 - my * my;
    if vx < 1e-12 || vy < 1e-12 {
        return None;
    }
    Some((sxy / n - mx * my) / (vx * vy).sqrt())
}

fn mean(x: &[i8]) -> f64 {
    x.iter().map(|&v| v as f64).sum::<f64>() / x.len() as f64
}

// DILATION MATRIX

fn dilation_matrix(v: &[i8], b: usize, big: usize) -> Vec<Vec<i128>> {
    let mut m = vec![vec![0i128; big]; big];
    for r in 0..big {
        for d in 0..b {
            let t = big * d + r;
            let (q, s) = (t / b, t % b);
            m[r][q] += (v[d] * v[s]) as i128;
        }
    }
    m
}

fn mat_mul(a: &[Vec<i128>], b: &[Vec<i128>]) -> Vec<Vec<i128>> {
    let n = a.len();
    let mut c = vec![vec![0i128; n]; n];
    for i in 0..n {
        for k in 0..n {
            if a[i][k] == 0 {
                continue;
            }
            for j in 0..n {
                c[i][j] = c[i][j]
                    .checked_add(a[i][k].checked_mul(b[k][j]).expect("overflow"))
                    .expect("overflow");
            }
        }
    }
    c
}

fn radius_bound(m: &[Vec<i128>], squarings: u32) -> f64 {
    let mut p = m.to_vec();
    for _ in 0..squarings {
        p = mat_mul(&p, &p);
    }
    let norm = p
        .iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<i128>())
        .max()
        .unwrap();
    let k = 2f64.powi(squarings as i32);
    (norm as f64).powf(1.0 / k) * (1.0 + 1e-9)
}

fn ceil4(x: f64) -> f64 {
    (x * 1e4).ceil() / 1e4
}

// AUTOMATA

fn rule_word(rule: u8, l: u64, c: u64, r: u64) -> u64 {
    match rule {
        30 => l ^ (c | r),
        90 => l ^ r,
        150 => l ^ c ^ r,
        _ => {
            let mut out = 0u64;
            for p in 0..8 {
                if rule >> p & 1 == 1 {
                    let a = if p & 4 != 0 { l } else { !l };
                    let b = if p & 2 != 0 { c } else { !c };
                    let d = if p & 1 != 0 { r } else { !r };
                    out |= a & b & d;
                }
            }
            out
        }
    }
}

fn evolve(
    rule: u8,
    seed: &[bool],
    steps: usize,
    keep: &[usize],
    row_len: usize,
) -> (Vec<u8>, Vec<Vec<u8>>) {
    let width = (2 * steps + 2 * seed.len() + 512) / 64 * 64;
    let nw = width / 64;
    let centre = width / 2;
    let mut cur = vec![0u64; nw];
    for (n, &on) in seed.iter().enumerate() {
        if on {
            let p = centre + n;
            cur[p / 64] |= 1u64 << (p % 64);
        }
    }
    let mut lo = centre;
    let mut hi = centre + seed.len();
    let mut next = vec![0u64; nw];
    let mut column = Vec::with_capacity(steps);
    let mut rows = Vec::new();
    let bit = |row: &[u64], p: usize| (row[p / 64] >> (p % 64) & 1) as u8;
    let quiescent = rule & 1 == 0;
    for t in 0..steps {
        column.push(bit(&cur, centre));
        if keep.contains(&t) {
            rows.push((0..row_len).map(|n| bit(&cur, centre + n)).collect());
        }
        lo -= 1;
        hi += 1;
        let (a, z) = if quiescent {
            (
                (lo / 64).saturating_sub(1).max(1),
                (hi / 64 + 1).min(nw - 2),
            )
        } else {
            (1, nw - 2)
        };
        for i in a..=z {
            let c = cur[i];
            let l = (c << 1) | (cur[i - 1] >> 63);
            let r = (c >> 1) | (cur[i + 1] << 63);
            next[i] = rule_word(rule, l, c, r);
        }
        std::mem::swap(&mut cur, &mut next);
    }
    (column, rows)
}

fn kernel_counts(s: &[u8], b: usize, ell: usize) -> Vec<usize> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut out = Vec::new();
    let mut bk = 1usize;
    while bk * ell <= s.len() {
        for r in 0..bk {
            seen.insert((0..ell).map(|n| s[bk * n + r]).collect());
        }
        out.push(seen.len());
        bk *= b;
    }
    out
}

fn lz78(bits: &[u8]) -> usize {
    let mut dict: HashMap<(u32, u8), u32> = HashMap::new();
    let (mut node, mut next, mut phrases) = (0u32, 1u32, 0usize);
    for &b in bits {
        match dict.get(&(node, b)) {
            Some(&c) => node = c,
            None => {
                dict.insert((node, b), next);
                next += 1;
                phrases += 1;
                node = 0;
            }
        }
    }
    phrases + usize::from(node != 0)
}

// CODES

fn class_name(vp: &[i8], b: usize) -> &'static str {
    match (b, vp) {
        (_, v) if v.iter().all(|&x| x == 1) => "constant",
        (2, _) => "thue-morse",
        (3, [1, -1, 1]) => "(-1)^n",
        (3, [1, 1, -1]) => "parity of 2s",
        (3, [1, -1, -1]) => "parity of nonzero digits",
        _ => "?",
    }
}

fn codes() {
    let h_max = 4096usize;
    println!("signed designs: x_L(n) = prod over the L digits of v(digit), v = (-1)^code; class = v v(0), the infinite sequence up to a global sign");
    for b in [2usize, 3] {
        let level: u32 = if b == 2 { 20 } else { 13 };
        let top = (b as u64).pow(level - 1);
        for mask in 0u32..(1 << b) {
            let v: Vec<i8> = (0..b)
                .map(|d| if mask >> d & 1 == 1 { -1 } else { 1 })
                .collect();
            let vp: Vec<i8> = v.iter().map(|&x| x * v[0]).collect();
            let mut bad = 0u64;
            for n in 0..top {
                let x = signed(&v, b as u64, n, level);
                for r in 0..b {
                    if signed(&v, b as u64, b as u64 * n + r as u64, level) != v[0] * v[r] * x {
                        bad += 1;
                    }
                }
            }
            let bias: Vec<i8> = (0..b).map(|r| v[0] * v[r]).collect();
            let c = corr_exact(&vp, b, h_max);
            let (mut best, mut at) = (Q::int(0), 0usize);
            let mut period = 0usize;
            let mut first = 0usize;
            for h in 1..=h_max {
                if c[h].abs() == Q::int(1) && period == 0 {
                    period = h;
                }
                if c[h] != Q::int(0) && first == 0 {
                    first = h;
                }
                if best.lt(c[h].abs()) {
                    best = c[h].abs();
                    at = h;
                }
            }
            let mut dev = 0f64;
            for h in 1..=8u64 {
                let s: i64 = (0..top)
                    .map(|n| {
                        (signed(&v, b as u64, n, level) * signed(&v, b as u64, n + h, level)) as i64
                    })
                    .sum();
                dev = dev.max((s as f64 / top as f64 - c[h as usize].f()).abs());
            }
            println!(
                "b={b} code={mask:0w$b} v={v:?} class={} dilation violations={bad} of {} pairs, bias by r={bias:?}; C(1..b)=[{}] first nonzero lag {first} C={} ; max|C(h)| h<={h_max} = {} at h={at}; period {}; render C(h<=8) vs exact max dev {dev:.2e} at N={top}",
                class_name(&vp, b),
                top * b as u64,
                (1..=b).map(|h| c[h].show()).collect::<Vec<_>>().join(", "),
                c[first.max(1)].show(),
                best.show(),
                if period > 0 { period.to_string() } else { "none".into() },
                w = b
            );
        }
    }
    println!("digit-restricted designs: A_L(n) = all L digits in D; Pearson(A_L(n), A_L(bn + r)) over n < b^(L-1)");
    for b in [2usize, 3] {
        let level: u32 = if b == 2 { 20 } else { 13 };
        let top = (b as u64).pow(level - 1);
        for mask in 1u32..(1 << b) {
            let cells: Vec<String> = (0..b)
                .map(|r| {
                    let p = pearson((0..top).map(|n| {
                        (
                            sign(restricted(mask, b as u64, n, level)),
                            sign(restricted(mask, b as u64, b as u64 * n + r as u64, level)),
                        )
                    }));
                    match p {
                        Some(x) => format!("{x:.6}"),
                        None => "const".into(),
                    }
                })
                .collect();
            println!(
                "b={b} D={mask:0w$b} Pearson by r = [{}]",
                cells.join(", "),
                w = b
            );
        }
    }
}

// RUDIN-SHAPIRO

fn rudin() {
    let n_max = 1u64 << 20;
    let (mut bad2, mut bad4) = (0u64, 0u64);
    for n in 0..n_max {
        if rs(2 * n) != rs(n) {
            bad2 += 1;
        }
        if rs(4 * n + 1) != rs(n) {
            bad4 += 1;
        }
    }
    println!(
        "rudin-shapiro: x(2n) != x(n) at {bad2} of {n_max}; x(4n+1) != x(n) at {bad4} of {n_max}"
    );
    let h_max = 64usize;
    let x: Vec<i8> = (0..n_max + h_max as u64 + 1).map(rs).collect();
    let mut worst_at_power = 0i64;
    let mut worst_any = 0i64;
    let mut worst_h = 0usize;
    let mut table = Vec::new();
    for h in 1..=h_max {
        let mut s = 0i64;
        let mut mx = 0i64;
        let mut mp = 0i64;
        for n in 0..n_max as usize {
            s += (x[n] * x[n + h]) as i64;
            mx = mx.max(s.abs());
            if (n + 1).is_power_of_two() && n + 1 >= 1024 {
                mp = mp.max(s.abs());
            }
        }
        if h <= 8 {
            table.push(format!(
                "h={h}: max_N |S|={mx}, max over N=2^10..2^20 |S|={mp}, S(2^20)={s}"
            ));
        }
        if mx > worst_any {
            worst_any = mx;
            worst_h = h;
        }
        worst_at_power = worst_at_power.max(mp);
    }
    for line in table {
        println!("  {line}");
    }
    println!("  over 1 <= h <= {h_max} and N <= 2^20: max |S(N, h)| = {worst_any} at h = {worst_h}; at N = 2^k, k = 10..20: {worst_at_power}; sqrt(2^20) = 1024");
    let p = pearson((0..n_max).map(|n| (x[n as usize], rs(2 * n)))).unwrap();
    println!("  Pearson(x(n), x(2n)) over n < 2^20 = {p:.6}");
    let y = |n: u64| if n % 2 == 0 { tm(n) } else { -tm(n) };
    let mut dmax = 0f64;
    for d in 1..=3u32 {
        let big = 1u64 << d;
        for r in 0..big {
            if let Some(p) = pearson((0..n_max).map(|n| (y(n), y(big * n + r)))) {
                dmax = dmax.max(p.abs());
            }
        }
    }
    let k = pearson((0..n_max).map(|n| (y(2 * n), y(4 * n)))).unwrap();
    println!("(-1)^n x thue-morse: max |Pearson(y(n), y(2^d n + r))| over d <= 3, r < 2^d, n < 2^20 = {dmax:.6}; Pearson(y(2n), y(4n)) = {k:.6}");
}

// BATTERY

fn battery() {
    let n = 1usize << 15;
    let m = 9 * n + 16;
    let mut state = 1u64;
    let ctrl: Vec<i8> = (0..m)
        .map(|_| sign(splitmix(&mut state) >> 63 == 1))
        .collect();
    let (col, _) = evolve(30, &[true], m, &[], 0);
    let a051023 = [
        1u8, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1,
        1, 0, 0, 1, 1, 1,
    ];
    println!(
        "rule 30 centre column agrees with the first {} terms of A051023: {}",
        a051023.len(),
        col[..a051023.len()] == a051023
    );
    let seqs: Vec<(&str, Vec<i8>)> = vec![
        ("thue-morse", (0..m as u64).map(tm).collect()),
        ("rudin-shapiro", (0..m as u64).map(rs).collect()),
        ("cantor", (0..m as u64).map(|k| sign(cantor(k))).collect()),
        (
            "carpet-morton",
            (0..m as u64).map(|k| sign(carpet_morton(k))).collect(),
        ),
        (
            "rule-30-centre",
            col.iter().map(|&b| sign(b == 1)).collect(),
        ),
        ("moire-xor", (0..m as u64).map(|k| tm(k) * t3(k)).collect()),
        (
            "moire-rudin",
            (0..m as u64).map(|k| rs(k) * t3(k)).collect(),
        ),
        ("splitmix64-seed-1", ctrl),
    ];
    println!(
        "battery: N = {n} pairs per test; null band 4/sqrt(N) = {:.4}",
        4.0 / (n as f64).sqrt()
    );
    for (name, x) in &seqs {
        let mu = mean(&x[..n]);
        let (mut sh, mut sh_at) = (0f64, 0usize);
        for h in 1..=64 {
            if let Some(p) = pearson((0..n).map(|k| (x[k], x[k + h]))) {
                if p.abs() > sh {
                    sh = p.abs();
                    sh_at = h;
                }
            }
        }
        let (mut dl, mut dl_at) = (0f64, (0usize, 0usize));
        for big in 2..=9usize {
            for r in 0..big {
                if let Some(p) = pearson((0..n).map(|k| (x[k], x[big * k + r]))) {
                    if p.abs() > dl + 1e-12 {
                        dl = p.abs();
                        dl_at = (big, r);
                    }
                }
            }
        }
        let square = (0..n)
            .map(|k| (x[k] * x[2 * k] * x[3 * k] * x[6 * k]) as f64)
            .sum::<f64>()
            / n as f64;
        let mut ke = 0f64;
        for b in [2usize, 3] {
            if let Some(p) = pearson((0..n).map(|k| (x[b * k], x[b * b * k]))) {
                ke = ke.max(p.abs());
            }
        }
        let bits: Vec<u8> = x.iter().map(|&v| u8::from(v < 0)).collect();
        let lz: Vec<String> = [14u32, 16, 18]
            .iter()
            .map(|&e| {
                let len = 1usize << e;
                let c = lz78(&bits[..len]) as f64;
                format!("{:.3}", c * (c.log2() + 1.0) / len as f64)
            })
            .collect();
        println!(
            "  {name:<18} mean {mu:+.4} | shift max|r| {sh:.4} at h={sh_at} | dilation max|r| {dl:.4} at (B, r)=({}, {}) | kernel max|r| {ke:.4} | square mean {square:+.4} | LZ78 bits per bit at 2^14, 2^16, 2^18: {}",
            dl_at.0,
            dl_at.1,
            lz.join(", ")
        );
    }
}

// MOIRE

fn moire() {
    let v2 = [1i8, -1];
    let v3 = [1i8, 1, -1];
    let c2 = corr_exact(&v2, 2, 4096);
    let c3 = corr_exact(&v3, 3, 4096);
    let small: Vec<String> = (1..=12).map(|h| c2[h].mul(c3[h]).show()).collect();
    println!("C2(h) C3(h), h = 1..12: [{}]", small.join(", "));
    let h_top = 1usize << 21;
    let f2 = corr_float(&v2, 2, h_top);
    let f3 = corr_float(&v3, 3, h_top);
    let mut tail = vec![0f64; 22];
    let mut window = Vec::new();
    for k in 1..21usize {
        let (mut mx, mut at) = (0f64, 0usize);
        for h in (1usize << k)..(1usize << (k + 1)) {
            let p = (f2[h] * f3[h]).abs();
            if p > mx {
                mx = p;
                at = h;
            }
        }
        window.push((k, mx, at));
    }
    for k in (1..21usize).rev() {
        tail[k] = window[k - 1].1.max(tail[k + 1]);
    }
    for (k, mx, at) in &window {
        println!("  window [2^{k}, 2^{}): max |C2 C3| = {mx:.6e} at h = {at} = {}; sup over [2^{k}, 2^21) = {:.6e}", k + 1, base3(*at), tail[*k]);
    }
    let nz = 1usize << 22;
    let z: Vec<i8> = (0..(nz + 300) as u64).map(|k| tm(k) * t3(k)).collect();
    let mut dev = 0f64;
    for h in 1..=256usize {
        let s: i64 = (0..nz).map(|k| (z[k] * z[k + h]) as i64).sum();
        dev = dev.max((s as f64 / nz as f64 - f2[h] * f3[h]).abs());
    }
    println!("empirical C_z(h) at N = 2^22 against C2(h) C3(h), h <= 256: max |difference| = {dev:.3e}; 1/sqrt(N) = {:.3e}", 1.0 / (nz as f64).sqrt());
    println!("dilation sums of z: z(n) z(B n + r) = t2(r) t3(n) t3(B n + r) for B = 2^a, = t3(r) t2(n) t2(B n + r) for B = 3^c; bound on the spectral radius rho of the block matrix, theta = log_b rho");
    let cases: Vec<(&str, &[i8], usize, usize)> = vec![
        ("control t2 at B=2", &v2, 2, 2),
        ("z at B=2", &v3, 3, 2),
        ("z at B=4", &v3, 3, 4),
        ("z at B=8", &v3, 3, 8),
        ("z at B=16", &v3, 3, 16),
        ("z at B=3", &v2, 2, 3),
        ("z at B=9", &v2, 2, 9),
        ("z at B=27", &v2, 2, 27),
    ];
    for (name, v, b, big) in cases {
        let m = dilation_matrix(v, b, big);
        let rho = radius_bound(&m, 6);
        let theta = rho.ln() / (b as f64).ln();
        let k = if b == 2 { 60 } else { 38 };
        let series = pair_series(v, b, 1, big, k);
        let mean_k =
            series[k].iter().map(|x| x.abs()).max().unwrap() as f64 / (b as f64).powi(k as i32);
        println!("  {name:<18} rho <= {:.4}, theta <= {:.4}; exact max_r |S_r(b^{k})|/b^{k} = {mean_k:.3e}, growth exponent of S along b^k over the last 10 levels {:.4}", ceil4(rho), ceil4(theta), growth(&series, b, k, 10));
    }
    println!("kernel pairs of z: z(2n) z(4n) = t3(2n) t3(4n) and z(3n) z(9n) = t2(3n) t2(9n); sums S(b^k) exactly from the block recursion");
    let s24 = pair_series(&v3, 3, 2, 4, 38);
    let ones = s24.iter().all(|row| row[0] == 1);
    println!("  z(2n) z(4n): S(3^k) = 1 for every 0 <= k <= 38: {ones}");
    let s39 = pair_series(&v2, 2, 3, 9, 64);
    let law = |k: usize| {
        s39[k + 4][0] == 2 * s39[k + 3][0] - 2 * s39[k + 2][0] + s39[k + 1][0] + 2 * s39[k][0]
    };
    let holds = (0..=60).all(law);
    println!(
        "  z(3n) z(9n): S(2^k), k = 0..12: {:?}; S(2^(k+4)) = 2 S(2^(k+3)) - 2 S(2^(k+2)) + S(2^(k+1)) + 2 S(2^k) for every 0 <= k <= 60: {holds} (27 states); characteristic polynomial (x^2 - x - 1)(x^2 - x + 2); S(2^60)/2^60 = {:.3e}, S(2^60)/phi^60 = {:.6}, growth exponent {:.6}, log2 phi = {:.6}",
        (0..=12).map(|k| s39[k][0]).collect::<Vec<_>>(),
        s39[60][0] as f64 / 2f64.powi(60),
        s39[60][0] as f64 / ((1.0 + 5f64.sqrt()) / 2.0).powi(60),
        growth(&s39, 2, 60, 10),
        ((1.0 + 5f64.sqrt()) / 2.0).log2()
    );
    let naf = |mut h: i64, b: i64| {
        let mut w = 0;
        while h != 0 {
            let r = h.rem_euclid(b);
            let d = if b == 2 {
                if r == 0 {
                    0
                } else {
                    2 - (h.rem_euclid(4))
                }
            } else if r == 2 {
                -1
            } else {
                r
            };
            if d != 0 {
                w += 1;
            }
            h = (h - d) / b;
        }
        w
    };
    for (name, f, b) in [
        ("C2 by binary NAF weight", &f2, 2i64),
        ("C3 by balanced ternary weight", &f3, 3),
    ] {
        let mut by = vec![0f64; 24];
        for h in 1..(1usize << 20) {
            let w = naf(h as i64, b) as usize;
            by[w] = by[w].max(f[h].abs());
        }
        let cells: Vec<String> = (1..12).map(|w| format!("w{w} {:.2e}", by[w])).collect();
        println!("  max |{name}| over h < 2^20: {}", cells.join(", "));
    }
    complexity();
    moire_rudin();
}

fn moire_rudin() {
    let nz = 1usize << 22;
    let y: Vec<i8> = (0..(nz + 300) as u64).map(|k| rs(k) * t3(k)).collect();
    let mut worst = (0f64, 0usize);
    for h in 1..=256usize {
        let s: i64 = (0..nz).map(|k| (y[k] * y[k + h]) as i64).sum();
        let c = (s as f64 / nz as f64).abs();
        if c > worst.0 {
            worst = (c, h);
        }
    }
    println!("moire-rudin y = rs t3: max |C_y,N(h)| over h <= 256 at N = 2^22 = {:.3e} at h = {}; 1/sqrt(N) = {:.3e}", worst.0, worst.1, 1.0 / (nz as f64).sqrt());
    for (name, a, b, c) in [
        ("rs(n) rs(3n)", 1u64, 3u64, 0u64),
        ("rs(n) rs(3n+1)", 1, 3, 1),
        ("rs(n) rs(3n+2)", 1, 3, 2),
        ("rs(3n) rs(9n)", 3, 9, 0),
        ("rs(n) rs(2n+1)", 1, 2, 1),
    ] {
        let mut s = 0i64;
        let mut marks = Vec::new();
        for n in 0..(1u64 << 24) {
            s += (rs(a * n) * rs(b * n + c)) as i64;
            if (n + 1).is_power_of_two() && n + 1 >= 1 << 12 {
                marks.push(s);
            }
        }
        let e = ((marks[marks.len() - 1].abs().max(1) as f64).ln()
            - (marks[marks.len() - 5].abs().max(1) as f64).ln())
            / (4.0 * 2f64.ln());
        println!("  sum over n < 2^k of {name}, k = 12..24: {:?}; growth exponent over the last 4 levels {e:.4}", marks);
    }
    let big = 1usize << 22;
    let xr: Vec<i8> = (0..big as u64 + 64).map(rs).collect();
    let x3: Vec<i8> = (0..big as u64 + 64).map(t3).collect();
    for n in [8usize, 16, 32, 64] {
        let k2 = (n as f64).log2().ceil() as u32;
        let k3 = ((n as f64).ln() / 3f64.ln()).ceil() as u32;
        let fr_small = factors(&xr[..64 * 2usize.pow(k2) + n], n);
        let fr = factors(&xr[..big], n);
        let f3 = factors(&x3[..18 * 3usize.pow(k3) + n], n);
        let mut prod = HashSet::new();
        for a in &fr {
            for b in &f3 {
                prod.insert(a ^ b);
            }
        }
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let (mut a, mut b) = (0u64, 0u64);
        let mut pairs: HashSet<u128> = HashSet::new();
        for i in 0..big {
            a = ((a << 1) | u64::from(xr[i] < 0)) & mask;
            b = ((b << 1) | u64::from(x3[i] < 0)) & mask;
            if i + 1 >= n {
                pairs.insert(((a as u128) << 64) | b as u128);
            }
        }
        println!(
            "  n={n:<3} p_rs from prefix 64 2^k: {}, from 2^22: {} (8n - 8 = {}) | p3={} | joint pairs {} of {} | p_y = {}",
            fr_small.len(),
            fr.len(),
            8 * n - 8,
            f3.len(),
            pairs.len(),
            fr.len() * f3.len(),
            if pairs.len() == fr.len() * f3.len() { format!("{}", prod.len()) } else { "incomplete".into() }
        );
    }
}

fn pair_series(v: &[i8], b: usize, a_big: usize, b_big: usize, k_max: usize) -> Vec<Vec<i128>> {
    let x = |n: usize| signed(v, b as u64, n as u64, 64);
    let idx = |a: usize, r: usize| a * b_big + r;
    let size = a_big * b_big;
    let mut s: Vec<i128> = (0..size)
        .map(|i| (x(i / b_big) * x(i % b_big)) as i128)
        .collect();
    let mut out = vec![s.clone()];
    for _ in 0..k_max {
        let mut t = vec![0i128; size];
        for a in 0..a_big {
            for r in 0..b_big {
                for d in 0..b {
                    let (u, w) = (a_big * d + a, b_big * d + r);
                    t[idx(a, r)] += (v[u % b] * v[w % b]) as i128 * s[idx(u / b, w / b)];
                }
            }
        }
        s = t;
        out.push(s.clone());
    }
    out
}

fn growth(series: &[Vec<i128>], b: usize, k: usize, span: usize) -> f64 {
    let top = |k: usize| series[k].iter().map(|x| x.abs()).max().unwrap() as f64;
    ((top(k) / top(k - span)).ln() / span as f64) / (b as f64).ln()
}

fn base3(mut n: usize) -> String {
    let mut d = Vec::new();
    while n > 0 {
        d.push((b'0' + (n % 3) as u8) as char);
        n /= 3;
    }
    d.iter().rev().collect::<String>() + "_3"
}

fn factors(x: &[i8], n: usize) -> HashSet<u64> {
    let mut out = HashSet::new();
    for s in 0..=x.len() - n {
        let mut w = 0u64;
        for &v in &x[s..s + n] {
            w = (w << 1) | u64::from(v < 0);
        }
        out.insert(w);
    }
    out
}

fn complexity() {
    println!("factor complexity: p2, p3 exact from the prefix 18 b^k, b^k >= n; joint pairs read in the prefix of length 2^22; p_z = product set size when every pair occurs");
    let lens = [2usize, 4, 8, 12, 16, 24, 32, 48, 64];
    let big = 1usize << 22;
    let x2: Vec<i8> = (0..big as u64 + 64).map(tm).collect();
    let x3: Vec<i8> = (0..big as u64 + 64).map(t3).collect();
    let mut last: Option<(usize, usize)> = None;
    for &n in &lens {
        let k2 = (n as f64).log2().ceil() as u32;
        let k3 = ((n as f64).ln() / 3f64.ln()).ceil() as u32;
        let f2 = factors(&x2[..18 * 2usize.pow(k2) + n], n);
        let f3 = factors(&x3[..18 * 3usize.pow(k3) + n], n);
        let mut prod = HashSet::new();
        for a in &f2 {
            for b in &f3 {
                prod.insert(a ^ b);
            }
        }
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let (mut a, mut b) = (0u64, 0u64);
        let mut pairs: HashSet<u128> = HashSet::new();
        let mut zs: HashSet<u64> = HashSet::new();
        for i in 0..big {
            a = ((a << 1) | u64::from(x2[i] < 0)) & mask;
            b = ((b << 1) | u64::from(x3[i] < 0)) & mask;
            if i + 1 >= n {
                pairs.insert(((a as u128) << 64) | b as u128);
                zs.insert(a ^ b);
            }
        }
        let full = pairs.len() == f2.len() * f3.len();
        let pz = if full { prod.len() } else { zs.len() };
        let slope = last
            .map(|(m, p)| (pz as f64 / p as f64).ln() / (n as f64 / m as f64).ln())
            .map(|s| format!("{s:.3}"))
            .unwrap_or("-".into());
        println!(
            "  n={n:<3} p2={:<4} p3={:<4} joint pairs {} of {} ({}) | z factors in prefix {} | product set {} | p_z {} | local exponent {slope}",
            f2.len(),
            f3.len(),
            pairs.len(),
            f2.len() * f3.len(),
            if full { "complete" } else { "incomplete" },
            zs.len(),
            prod.len(),
            if full { format!("= {pz}") } else { format!("in [{}, {}]", zs.len(), prod.len()) }
        );
        last = Some((n, pz));
    }
}

// ROUNDS

fn rounds() {
    let steps = 1usize << 15;
    let ell = 24usize;
    let len = steps + 256;
    let single = vec![true];
    let tmb: Vec<bool> = (0..len as u64).map(|n| n.count_ones() % 2 == 1).collect();
    let cab: Vec<bool> = (0..len as u64).map(cantor).collect();
    println!("rounds: centre column c(t), t < 2^15; cumulative distinct kernel words of length {ell} over depth k = 0, 1, ...");
    let cases: Vec<(u8, &str, &Vec<bool>)> = vec![
        (90, "single cell", &single),
        (90, "thue-morse", &tmb),
        (150, "thue-morse", &tmb),
        (90, "cantor", &cab),
        (150, "cantor", &cab),
        (30, "single cell", &single),
        (30, "thue-morse", &tmb),
        (110, "single cell", &single),
        (254, "cantor", &cab),
        (254, "thue-morse", &tmb),
    ];
    for (rule, name, seed) in cases {
        let (col, _) = evolve(rule, seed, steps, &[], 0);
        println!(
            "  rule {rule:<3} seed {name:<12} base 2 {:?} | base 3 {:?}",
            kernel_counts(&col, 2, ell),
            kernel_counts(&col, 3, ell)
        );
    }
    let span = 1usize << 12;
    let keep: Vec<usize> = (0..span).collect();
    let (_, r90) = evolve(90, &single, span, &keep, span);
    for rule in [18u8, 22, 126, 146] {
        let (_, rr) = evolve(rule, &single, span, &keep, span);
        let off: usize = (0..span)
            .map(|t| (0..span).filter(|&n| rr[t][n] != r90[t][n]).count())
            .sum();
        let lucas: usize = (0..span)
            .map(|t| {
                (0..span)
                    .filter(|&n| {
                        let on = (t + n) % 2 == 0 && n <= t && ((t + n) / 2) & ((t - n) / 2) == 0;
                        (rr[t][n] == 1) != on
                    })
                    .count()
            })
            .sum();
        println!("  rule {rule:<3} single cell, t < {span}, n >= 0: cells differing from rule 90 {off}, from Pascal mod 2 {lucas}");
    }
    for (rule, name, seed) in [
        (18u8, "single cell", &single),
        (18, "thue-morse", &tmb),
        (18, "cantor", &cab),
    ] {
        let (col, _) = evolve(rule, seed, steps, &[], 0);
        println!(
            "  rule {rule:<3} seed {name:<12} base 2 {:?} | base 3 {:?}",
            kernel_counts(&col, 2, ell),
            kernel_counts(&col, 3, ell)
        );
    }
    println!("bounded rounds: row T of a rule on the thue-morse seed, n < 2^15");
    let (_, rows) = evolve(30, &tmb, 5, &[1, 2, 3, 4], steps);
    for (i, row) in rows.iter().enumerate() {
        println!(
            "  rule 30 row {}: base 2 {:?}",
            i + 1,
            kernel_counts(row, 2, ell)
        );
    }
    let (_, rows) = evolve(110, &tmb, 5, &[4], steps);
    println!(
        "  rule 110 row 4: base 2 {:?}",
        kernel_counts(&rows[0], 2, ell)
    );
    let tm_row: Vec<u8> = tmb[..steps].iter().map(|&b| u8::from(b)).collect();
    println!("  seed itself: base 2 {:?}", kernel_counts(&tm_row, 2, ell));
    let mut state = 1u64;
    let ctrl: Vec<u8> = (0..steps)
        .map(|_| (splitmix(&mut state) >> 63) as u8)
        .collect();
    println!(
        "  splitmix64 seed 1: base 2 {:?} | base 3 {:?}",
        kernel_counts(&ctrl, 2, ell),
        kernel_counts(&ctrl, 3, ell)
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("codes") => codes(),
        Some("rudin") => rudin(),
        Some("battery") => battery(),
        Some("moire") => moire(),
        Some("rounds") => rounds(),
        Some("shadow") => reach::shadow(),
        Some("escape") => reach::escape(),
        Some("lemma") => reach::lemma(),
        _ => eprintln!("verbs: codes | rudin | battery | moire | rounds | shadow | escape | lemma"),
    }
}
