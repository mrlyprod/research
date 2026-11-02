use std::cmp::Ordering;
use std::env;
use std::time::Instant;

// PACKING

fn binom(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for i in 0..k {
        r = r * (n - i) as u128 / (i + 1) as u128;
    }
    r as u64
}

struct Pack {
    off: Vec<u32>,
    wid: Vec<u32>,
    w: usize,
    bits: u32,
}

fn pack_for(d: usize) -> Pack {
    let mut off = Vec::with_capacity(d + 1);
    let mut wid = Vec::with_capacity(d + 1);
    let mut o = 0u32;
    for i in 0..=d {
        let b = (64 - binom(d, i).leading_zeros()).max(1);
        off.push(o);
        wid.push(b);
        o += b;
    }
    Pack {
        off,
        wid,
        w: (o as usize).div_ceil(64),
        bits: o,
    }
}

fn put(key: &mut [u64], i: usize, pk: &Pack, v: u32) {
    let bit = pk.off[i] as usize;
    let b = pk.wid[i] as usize;
    let w = bit >> 6;
    let off = bit & 63;
    key[w] |= (v as u64) << off;
    if off + b > 64 {
        key[w + 1] |= (v as u64) >> (64 - off);
    }
}

fn get(key: &[u64], i: usize, pk: &Pack) -> u32 {
    let bit = pk.off[i] as usize;
    let b = pk.wid[i] as usize;
    let w = bit >> 6;
    let off = bit & 63;
    let mut v = key[w] >> off;
    if off + b > 64 {
        v |= key[w + 1] << (64 - off);
    }
    (v & ((1u64 << b) - 1)) as u32
}

fn addk(dst: &mut [u64], src: &[u64]) {
    let mut c = 0u64;
    for i in 0..dst.len() {
        let (a, o1) = dst[i].overflowing_add(src[i]);
        let (a, o2) = a.overflowing_add(c);
        dst[i] = a;
        c = (o1 as u64) | (o2 as u64);
    }
}

fn subk(dst: &mut [u64], src: &[u64]) {
    let mut c = 0u64;
    for i in 0..dst.len() {
        let (a, o1) = dst[i].overflowing_sub(src[i]);
        let (a, o2) = a.overflowing_sub(c);
        dst[i] = a;
        c = (o1 as u64) | (o2 as u64);
    }
}

fn cmpk(a: &[u64], b: &[u64]) -> Ordering {
    for i in (0..a.len()).rev() {
        if a[i] != b[i] {
            return a[i].cmp(&b[i]);
        }
    }
    Ordering::Equal
}

fn hashk(k: &[u64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &w in k {
        h = (h ^ w).wrapping_mul(0x100_0000_01b3);
    }
    h ^ (h >> 31)
}

// TABLE

struct Table {
    w: usize,
    cap: usize,
    keys: Vec<u64>,
    vals: Vec<i64>,
    len: usize,
}

impl Table {
    fn new(w: usize, cap: usize) -> Table {
        Table {
            w,
            cap,
            keys: vec![0u64; w * cap],
            vals: vec![0i64; cap],
            len: 0,
        }
    }

    fn bytes(&self) -> usize {
        self.keys.len() * 8 + self.vals.len() * 8
    }

    fn empty(&self, s: usize) -> bool {
        self.keys[s * self.w..(s + 1) * self.w]
            .iter()
            .all(|&x| x == 0)
    }

    fn key(&self, s: usize) -> &[u64] {
        &self.keys[s * self.w..(s + 1) * self.w]
    }

    fn grow(&mut self) {
        let cap = self.cap * 2;
        let mut keys = vec![0u64; self.w * cap];
        let mut vals = vec![0i64; cap];
        for s in 0..self.cap {
            if self.empty(s) {
                continue;
            }
            let mut t = (hashk(self.key(s)) as usize) & (cap - 1);
            while keys[t * self.w..(t + 1) * self.w].iter().any(|&x| x != 0) {
                t = (t + 1) & (cap - 1);
            }
            keys[t * self.w..(t + 1) * self.w].copy_from_slice(self.key(s));
            vals[t] = self.vals[s];
        }
        self.cap = cap;
        self.keys = keys;
        self.vals = vals;
    }

    fn add(&mut self, k: &[u64], v: i64) {
        if (self.len + 1) * 10 > self.cap * 7 {
            self.grow();
        }
        let mut s = (hashk(k) as usize) & (self.cap - 1);
        loop {
            if self.empty(s) {
                self.keys[s * self.w..(s + 1) * self.w].copy_from_slice(k);
                self.vals[s] = v;
                self.len += 1;
                return;
            }
            if self.key(s) == k {
                self.vals[s] += v;
                return;
            }
            s = (s + 1) & (self.cap - 1);
        }
    }
}

// LAYERS

struct Layer {
    deg: usize,
    pk: Pack,
    w: usize,
    keys: Vec<u64>,
    nus: Vec<i32>,
    starts: Vec<u32>,
    count: u64,
}

impl Layer {
    fn new(deg: usize) -> Layer {
        let pk = pack_for(deg);
        let w = pk.w;
        let m2 = binom(deg, 2) as usize;
        Layer {
            deg,
            pk,
            w,
            keys: Vec::new(),
            nus: Vec::new(),
            starts: vec![0; (deg + 1) * (m2 + 1) + 1],
            count: 0,
        }
    }

    fn m2(&self) -> usize {
        binom(self.deg, 2) as usize
    }

    fn slot(&self, t1: usize, t2: usize) -> usize {
        t1 * (self.m2() + 1) + t2
    }

    fn bytes(&self) -> usize {
        self.keys.capacity() * 8 + self.nus.capacity() * 4 + self.starts.capacity() * 4
    }

    fn key(&self, i: usize) -> &[u64] {
        &self.keys[i * self.w..(i + 1) * self.w]
    }

    fn range(&self, t1: usize, t2: usize) -> (usize, usize) {
        if t1 > self.deg || t2 > self.m2() {
            return (0, 0);
        }
        let s = self.slot(t1, t2);
        (self.starts[s] as usize, self.starts[s + 1] as usize)
    }

    fn decode(&self, i: usize, out: &mut [u32]) {
        let k = self.key(i);
        for (j, o) in out.iter_mut().enumerate().take(self.deg + 1) {
            *o = get(k, j, &self.pk);
        }
    }
}

fn power(l: usize) -> i64 {
    if l == 0 {
        1
    } else if l == 1 {
        -1
    } else {
        0
    }
}

// LADDER

struct Stats {
    odd: Vec<u64>,
    cens: Vec<u64>,
    sums: Vec<i64>,
    mx: Vec<u64>,
    coef: Vec<u64>,
    bits: Vec<u32>,
    peak: usize,
    cut: usize,
    rev_bad: u64,
    rev_seen: u64,
}

fn ladder(top: usize, quiet: bool, budget: usize) -> Stats {
    let t0 = Instant::now();
    let mut layers: Vec<Layer> = Vec::new();
    let mut st = Stats {
        odd: Vec::new(),
        cens: Vec::new(),
        sums: Vec::new(),
        mx: Vec::new(),
        coef: Vec::new(),
        bits: Vec::new(),
        peak: 0,
        cut: top,
        rev_bad: 0,
        rev_seen: 0,
    };
    let mut run: i64 = 0;
    let mut best: u64 = 0;
    let mut total: u64 = 0;
    if !quiet {
        println!("   L  odd census    M* census M(q^L)       max    max/2^L  log2(max)/L    step  coef  bits heldMB     s");
    }
    let mut cut = false;
    for d in 0..=top {
        let mut l0 = Layer::new(d);
        let mut cmax: u64 = 0;
        if d == 0 {
            let mut k = vec![0u64; l0.w];
            put(&mut k, 0, &l0.pk, 1);
            l0.keys.extend_from_slice(&k);
            l0.nus.push(1);
            l0.starts[0] = 0;
            l0.starts[1] = 1;
            l0.count = 1;
            cmax = 1;
            run += 1;
            best = best.max(run.unsigned_abs());
        } else {
            let store = d < top;
            let pk = pack_for(d);
            let w = l0.w;
            let m2 = l0.m2();
            let mut pc = vec![0u32; top + 2];
            let mut shifts: Vec<Vec<u64>> = (0..=top + 1).map(|_| vec![0u64; w]).collect();
            let mut acc = vec![0u64; w];
            if store {
                let guess = layers[d - 1].count as usize * 5 / 2 + 16;
                l0.keys.reserve_exact(guess * w);
                l0.nus.reserve_exact(guess);
            }
            for v1 in 0..=d {
                for v2 in 0..=m2 {
                    let slot = v1 * (m2 + 1) + v2;
                    l0.starts[slot] = l0.count as u32;
                    let mut tab = Table::new(w, 1 << 6);
                    for (dp, layer) in layers.iter().enumerate().take(d) {
                        let j = d - dp;
                        for t1 in 0..=dp {
                            if v1 < t1 || v1 - t1 > 1 {
                                continue;
                            }
                            let n1 = v1 - t1;
                            if j == 1 && n1 != 1 {
                                continue;
                            }
                            for n2 in 0..=1usize {
                                if (j == 1 && n2 != 0) || (j == 2 && n2 != 1) {
                                    continue;
                                }
                                let off = t1 * n1 + n2;
                                if v2 < off {
                                    continue;
                                }
                                let (lo, hi) = layer.range(t1, v2 - off);
                                for idx in lo..hi {
                                    let nu = layer.nus[idx] as i64;
                                    layer.decode(idx, &mut pc);
                                    for (i, s) in shifts.iter_mut().enumerate().take(j + 1) {
                                        for x in s.iter_mut() {
                                            *x = 0;
                                        }
                                        for (k, &p) in pc.iter().enumerate().take(dp + 1) {
                                            if p != 0 {
                                                put(s, k + i, &pk, p);
                                            }
                                        }
                                    }
                                    for x in acc.iter_mut() {
                                        *x = 0;
                                    }
                                    addk(&mut acc, &shifts[0]);
                                    addk(&mut acc, &shifts[j]);
                                    if j >= 2 && n1 == 1 {
                                        addk(&mut acc, &shifts[j - 1]);
                                    }
                                    if j >= 3 && n2 == 1 {
                                        addk(&mut acc, &shifts[j - 2]);
                                    }
                                    tab.add(&acc, nu);
                                    let f = j.saturating_sub(3);
                                    for g in 1u64..(1u64 << f) {
                                        let cur = g ^ (g >> 1);
                                        let prev = (g - 1) ^ ((g - 1) >> 1);
                                        let dif = cur ^ prev;
                                        let i = 1 + dif.trailing_zeros() as usize;
                                        if cur & dif != 0 {
                                            addk(&mut acc, &shifts[i]);
                                        } else {
                                            subk(&mut acc, &shifts[i]);
                                        }
                                        tab.add(&acc, nu);
                                    }
                                }
                            }
                        }
                    }
                    let mut ord: Vec<u32> = (0..tab.cap as u32)
                        .filter(|&s| !tab.empty(s as usize))
                        .collect();
                    ord.sort_unstable_by(|&x, &y| cmpk(tab.key(x as usize), tab.key(y as usize)));
                    let (lo, hi) = layers[d - 1].range(v1, v2);
                    let mut sk = vec![0u64; (hi - lo) * w];
                    for (n, i) in (lo..hi).enumerate() {
                        layers[d - 1].decode(i, &mut pc);
                        for (c, &p) in pc.iter().enumerate().take(d) {
                            if p != 0 {
                                put(&mut sk[n * w..(n + 1) * w], c + 1, &pk, p);
                            }
                        }
                    }
                    let mut a = 0usize;
                    let mut c = 0usize;
                    while a < ord.len() || c < hi - lo {
                        let take_a = if a == ord.len() {
                            false
                        } else if c == hi - lo {
                            true
                        } else {
                            cmpk(tab.key(ord[a] as usize), &sk[c * w..(c + 1) * w])
                                == Ordering::Less
                        };
                        if take_a {
                            let s = ord[a] as usize;
                            let nu = -tab.vals[s];
                            run += nu;
                            best = best.max(run.unsigned_abs());
                            assert!(nu >= i32::MIN as i64 && nu <= i32::MAX as i64);
                            let k = tab.key(s);
                            assert_eq!(get(k, 0, &pk), 1);
                            assert_eq!(get(k, d, &pk), 1);
                            for i in 0..=d {
                                let a = get(k, i, &pk) as u64;
                                assert!(a <= binom(d, i));
                                cmax = cmax.max(a);
                            }
                            if store {
                                l0.keys.extend_from_slice(k);
                                l0.nus.push(nu as i32);
                            }
                            l0.count += 1;
                            a += 1;
                        } else {
                            run -= layers[d - 1].nus[lo + c] as i64;
                            best = best.max(run.unsigned_abs());
                            c += 1;
                        }
                    }
                    st.peak = st.peak.max(
                        layers.iter().map(|l| l.bytes()).sum::<usize>()
                            + l0.bytes()
                            + tab.bytes()
                            + ord.capacity() * 4
                            + sk.len() * 8,
                    );
                    if st.peak > budget {
                        cut = true;
                        break;
                    }
                }
                if cut {
                    break;
                }
            }
            if cut {
                if !quiet {
                    println!(
                        "  CUT inside level {} on the memory budget, {} MB",
                        d + 1,
                        st.peak / 1_048_576
                    );
                }
                break;
            }
            l0.starts[(d + 1) * (m2 + 1)] = l0.count as u32;
            assert!(cmax <= binom(d, d / 2));
            if store {
                let n = l0.count as usize;
                let mut rk = vec![0u64; w];
                for i in 0..n {
                    for x in rk.iter_mut() {
                        *x = 0;
                    }
                    for j in 0..=d {
                        let a = get(l0.key(i), j, &pk);
                        if a != 0 {
                            put(&mut rk, d - j, &pk, a);
                        }
                    }
                    let mut lo = 0usize;
                    let mut hi = n;
                    while lo < hi {
                        let mid = (lo + hi) / 2;
                        if cmpk(l0.key(mid), &rk) == Ordering::Less {
                            lo = mid + 1;
                        } else {
                            hi = mid;
                        }
                    }
                    if lo >= n
                        || cmpk(l0.key(lo), &rk) != Ordering::Equal
                        || l0.nus[lo] != l0.nus[i]
                    {
                        st.rev_bad += 1;
                    }
                    st.rev_seen += 1;
                }
            }
        }
        total += l0.count;
        let lv = d + 1;
        let sum = run + power(lv);
        let mx = best.max(sum.unsigned_abs());
        st.odd.push(l0.count);
        st.cens.push(total);
        st.sums.push(sum);
        st.mx.push(mx);
        st.coef.push(cmax);
        st.bits.push(l0.pk.bits);
        layers.push(l0);
        if !quiet {
            let step = if d > 0 {
                mx as f64 / st.mx[d - 1] as f64
            } else {
                0.0
            };
            println!(
                "  {:>2} {:>11} {:>12} {:>6} {:>9} {:>10.6} {:>12.6} {:>7.4} {:>5} {:>5} {:>6} {:>5.1}",
                lv,
                st.odd[d],
                total,
                sum,
                mx,
                mx as f64 / (1u64 << lv) as f64,
                (mx as f64).log2() / lv as f64,
                step,
                cmax,
                st.bits[d],
                st.peak / 1_048_576,
                t0.elapsed().as_secs_f64()
            );
        }
        st.cut = lv;
        if lv <= 18 {
            assert_eq!(st.mx[d], PIN_MX[d]);
            assert_eq!(st.cens[d], PIN_CENS[d]);
            assert_eq!(st.coef[d], PIN_COEF[d]);
            assert_eq!(st.sums[d], if lv == 1 { 0 } else { -1 });
            if lv == 18 && !quiet {
                println!("  the Python ladder to level 18 is matched term for term");
            }
        }
    }
    st
}

// REPORT

fn report(st: &Stats) {
    let n = st.mx.len();
    println!();
    println!("  running maxima  {:?}", st.mx);
    println!("  M(q^L)          {:?}", st.sums);
    println!("  monoid census   {:?}", st.cens);
    println!("  odd census      {:?}", st.odd);
    println!("  max coefficient {:?}", st.coef);
    println!();
    for r in 0..2usize {
        let cl: Vec<String> = (1..=n)
            .filter(|l| l % 2 == r)
            .map(|l| format!("{:.6}", st.mx[l - 1] as f64 / (1u64 << l) as f64))
            .collect();
        println!("  L = {} mod 2   max/2^L  {}", r, cl.join(", "));
    }
    println!();
    let mut first = 1;
    for l in 2..=n {
        let lo = st.mx[l - 2] as f64 / (1u64 << (l - 1)) as f64;
        let hi = st.mx[l - 1] as f64 / (1u64 << l) as f64;
        if hi <= lo {
            first = l;
        }
    }
    println!(
        "  max/2^L falls last at level {} and rises at every step to level {}",
        first, n
    );
    println!();
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    for w in [4usize, 6, 8, 10, 12, 14, 16, 18, 20] {
        if n > w {
            let g = (st.mx[n - 1] as f64 / st.mx[n - 1 - w] as f64).powf(1.0 / w as f64);
            lo = lo.min(g);
            hi = hi.max(g);
            println!(
                "  geometric mean step over the last {:>2} levels  {:.6}   log2 {:.6}",
                w,
                g,
                g.log2()
            );
        }
    }
    println!(
        "  the window hull of the step at depth {} is [{:.6}, {:.6}], the mass rate is 2",
        n, lo, hi
    );
    println!();
    for w in [4usize, 6, 8] {
        let mut row = Vec::new();
        let first = if n > w + 3 { n - 3 } else { w + 1 };
        for m in first..=n {
            let g = (st.mx[m - 1] as f64 / st.mx[m - 1 - w] as f64).powf(1.0 / w as f64);
            row.push(format!("{:.6}", g));
        }
        println!(
            "  the {:>2} level window read at depths {}..{}  {}",
            w,
            first,
            n,
            row.join(", ")
        );
    }
    println!();
    println!(
        "  the design mass rate is 2, the deepest level is {}, peak {} MB",
        st.cut,
        st.peak / 1_048_576
    );
    println!(
        "  every level to {} is complete, the rate below is a reading to that level",
        n
    );
    println!(
        "  the reciprocal Q -> x^(deg Q) Q(1/x) holds nu* fixed, {} mismatches over {} elements",
        st.rev_bad, st.rev_seen
    );
}

// PYTHON PIN

const PIN_MX: [u64; 18] = [
    1, 1, 2, 3, 4, 7, 15, 23, 45, 86, 162, 331, 741, 1665, 3173, 7508, 17753, 36147,
];

const PIN_CENS: [u64; 18] = [
    1, 2, 5, 11, 27, 61, 144, 331, 776, 1788, 4147, 9544, 22000, 50420, 115407, 263062, 598540,
    1357535,
];

const PIN_COEF: [u64; 18] = [
    1, 1, 2, 3, 6, 10, 20, 35, 70, 126, 252, 462, 924, 1716, 3432, 6435, 12870, 24310,
];

// VERBS

fn pin() {
    let st = ladder(17, false, usize::MAX);
    report(&st);
    println!();
    let mut bad = 0;
    for l in 1..=18usize {
        if st.mx[l - 1] != PIN_MX[l - 1] {
            bad += 1;
        }
        if st.cens[l - 1] != PIN_CENS[l - 1] {
            bad += 1;
        }
        if st.coef[l - 1] != PIN_COEF[l - 1] {
            bad += 1;
        }
        let want = if l == 1 { 0 } else { -1 };
        if st.sums[l - 1] != want {
            bad += 1;
        }
    }
    println!("  the Python ladder to level 18, mismatches {}", bad);
}

fn main() {
    let a: Vec<String> = env::args().collect();
    let verb = a.get(1).map(|s| s.as_str()).unwrap_or("pin");
    match verb {
        "pin" => pin(),
        "ladder" => {
            let top: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);
            let gb: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(2.7);
            let st = ladder(top, false, (gb * 1_073_741_824.0) as usize);
            report(&st);
        }
        _ => println!("verbs: pin, ladder <top degree> <GB>"),
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ladder_pins_to_the_python_at_level_eighteen() {
        let st = ladder(17, true, usize::MAX);
        assert_eq!(st.mx, PIN_MX.to_vec());
        assert_eq!(st.cens, PIN_CENS.to_vec());
        assert_eq!(st.coef, PIN_COEF.to_vec());
        assert_eq!(st.sums[0], 0);
        assert!(st.sums[1..].iter().all(|&s| s == -1));
    }
}
