use mrly_num::factor::gcd;

const KMAX: usize = 64;
const MATCH: f64 = 1e-9;
const PRIME_BITS: f64 = 61.0;

pub struct Live {
    pub states: Vec<i64>,
    pub edges: Vec<Vec<usize>>,
}

pub fn live(a: u64, b: u64) -> Live {
    let mut index = std::collections::HashMap::new();
    let mut order: Vec<i64> = vec![0];
    let mut edges: Vec<Vec<usize>> = Vec::new();
    index.insert(0i64, 0usize);
    let mut i = 0;
    while i < order.len() {
        let c = order[i];
        let mut out = Vec::new();
        for eps in [0i64, b as i64, -(a as i64)] {
            if (c + eps).rem_euclid(3) == 0 {
                let next = (c + eps).div_euclid(3);
                let j = *index.entry(next).or_insert_with(|| {
                    order.push(next);
                    order.len() - 1
                });
                out.push(j);
            }
        }
        edges.push(out);
        i += 1;
    }
    let mut back = vec![Vec::new(); order.len()];
    for (i, out) in edges.iter().enumerate() {
        for &j in out {
            back[j].push(i);
        }
    }
    let mut seen = vec![false; order.len()];
    seen[0] = true;
    let mut stack = vec![0usize];
    while let Some(x) = stack.pop() {
        for &y in &back[x] {
            if !seen[y] {
                seen[y] = true;
                stack.push(y);
            }
        }
    }
    let keep: Vec<usize> = (0..order.len()).filter(|&i| seen[i]).collect();
    let mut place = vec![usize::MAX; order.len()];
    for (k, &i) in keep.iter().enumerate() {
        place[i] = k;
    }
    Live {
        states: keep.iter().map(|&i| order[i]).collect(),
        edges: keep
            .iter()
            .map(|&i| {
                edges[i]
                    .iter()
                    .filter(|&&j| seen[j])
                    .map(|&j| place[j])
                    .collect()
            })
            .collect(),
    }
}

pub fn perron(live: &Live) -> (f64, f64, f64) {
    let n = live.edges.len();
    let mut v = vec![1.0f64; n];
    let (mut lo, mut hi) = (0.0, f64::INFINITY);
    for _ in 0..200_000 {
        let mut w = vec![0.0f64; n];
        for (i, out) in live.edges.iter().enumerate() {
            for &j in out {
                w[i] += v[j];
            }
        }
        lo = f64::INFINITY;
        hi = 0.0;
        for i in 0..n {
            let r = w[i] / v[i];
            lo = lo.min(r);
            hi = hi.max(r);
        }
        let norm = w.iter().sum::<f64>();
        for x in w.iter_mut() {
            *x /= norm;
        }
        v = w;
        if hi - lo < 1e-13 {
            break;
        }
    }
    ((lo + hi) / 2.0, lo, hi)
}

pub fn first_returns(live: &Live) -> Option<Vec<u128>> {
    let n = live.edges.len();
    let mut indeg = vec![0usize; n];
    for (i, out) in live.edges.iter().enumerate() {
        if i == 0 {
            continue;
        }
        for &j in out {
            if j != 0 {
                indeg[j] += 1;
            }
        }
    }
    let mut ready: Vec<usize> = (1..n).filter(|&i| indeg[i] == 0).collect();
    let mut topo = Vec::new();
    while let Some(x) = ready.pop() {
        topo.push(x);
        for &j in &live.edges[x] {
            if j != 0 {
                indeg[j] -= 1;
                if indeg[j] == 0 {
                    ready.push(j);
                }
            }
        }
    }
    if topo.len() != n - 1 {
        return None;
    }
    let mut depth = vec![0usize; n];
    for &x in topo.iter().rev() {
        depth[x] = 1 + live.edges[x]
            .iter()
            .map(|&j| if j == 0 { 0 } else { depth[j] })
            .max()
            .unwrap_or(0);
    }
    let longest = 1 + live.edges[0]
        .iter()
        .map(|&j| if j == 0 { 0 } else { depth[j] })
        .max()
        .unwrap_or(0);
    let mut f = vec![0u128; longest + 1];
    let mut cur = vec![0u128; n];
    cur[0] = 1;
    for slot in f.iter_mut().skip(1) {
        let mut next = vec![0u128; n];
        for (i, out) in live.edges.iter().enumerate() {
            if cur[i] == 0 {
                continue;
            }
            for &j in out {
                next[j] = next[j].checked_add(cur[i])?;
            }
        }
        *slot = next[0];
        next[0] = 0;
        cur = next;
    }
    Some(f)
}

pub fn trinomial(family: char, k: usize) -> Vec<i128> {
    let mut q = vec![0i128; k + 1];
    q[k] = 1;
    q[0] = -1;
    q[if family == 'A' { k - 1 } else { 1 }] -= 1;
    q
}

pub fn root(family: char, k: usize) -> f64 {
    let f = |x: f64| {
        if family == 'A' {
            x.powi(k as i32) - x.powi(k as i32 - 1) - 1.0
        } else {
            x.powi(k as i32) - x - 1.0
        }
    };
    let (mut lo, mut hi) = (1.0f64, 2.0f64);
    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

pub fn renewal_polynomial(f: &[u128]) -> Option<Vec<i128>> {
    let k = f.len() - 1;
    let mut p = vec![0i128; k + 1];
    p[k] = 1;
    for (j, &fj) in f.iter().enumerate().skip(1) {
        p[k - j] = -i128::try_from(fj).ok()?;
    }
    Some(p)
}

pub fn divides(q: &[i128], p: &[i128]) -> bool {
    let mut r: Vec<i128> = p.to_vec();
    let k = q.len() - 1;
    while r.len() > k {
        let lead = *r.last().unwrap();
        let shift = r.len() - 1 - k;
        for i in 0..=k {
            r[shift + i] -= lead * q[i];
        }
        r.pop();
    }
    r.iter().all(|&c| c == 0)
}

pub fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

pub fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(q) {
            return n == q;
        }
    }
    let d = (n - 1) >> (n - 1).trailing_zeros();
    'outer: for &base in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(base, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut dd = d;
        while dd != n - 1 {
            x = mulmod(x, x, n);
            dd <<= 1;
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

pub fn primes(count: usize) -> Vec<u64> {
    let mut out = Vec::new();
    let mut n = 1u64 << 61;
    while out.len() < count {
        if is_prime(n) {
            out.push(n);
        }
        n += 1;
    }
    out
}

pub fn sub(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}

pub fn charpoly_mod(live: &Live, p: u64) -> Vec<u64> {
    let n = live.edges.len();
    let mut h = vec![vec![0u64; n]; n];
    for (i, out) in live.edges.iter().enumerate() {
        for &j in out {
            h[i][j] = (h[i][j] + 1) % p;
        }
    }
    for m in 1..n.saturating_sub(1) {
        let Some(piv) = (m..n).find(|&i| h[i][m - 1] != 0) else {
            continue;
        };
        if piv != m {
            h.swap(piv, m);
            for row in h.iter_mut() {
                row.swap(piv, m);
            }
        }
        let tinv = powmod(h[m][m - 1], p - 2, p);
        for i in m + 1..n {
            if h[i][m - 1] == 0 {
                continue;
            }
            let u = mulmod(h[i][m - 1], tinv, p);
            let row_m = h[m].clone();
            for (cell, &r) in h[i].iter_mut().zip(&row_m) {
                *cell = sub(*cell, mulmod(u, r, p), p);
            }
            for row in h.iter_mut() {
                let add = mulmod(u, row[i], p);
                row[m] = (row[m] + add) % p;
            }
        }
    }
    let mut polys: Vec<Vec<u64>> = vec![vec![1u64]];
    for k in 0..n {
        let prev = &polys[k];
        let mut next = vec![0u64; prev.len() + 1];
        for (d, &c) in prev.iter().enumerate() {
            next[d + 1] = (next[d + 1] + c) % p;
            next[d] = sub(next[d], mulmod(h[k][k], c, p), p);
        }
        let mut prod = 1u64;
        for i in (0..k).rev() {
            prod = mulmod(prod, h[i + 1][i], p);
            if prod == 0 {
                break;
            }
            let coef = mulmod(h[i][k], prod, p);
            if coef == 0 {
                continue;
            }
            for (d, &c) in polys[i].iter().enumerate() {
                next[d] = sub(next[d], mulmod(coef, c, p), p);
            }
        }
        polys.push(next);
    }
    polys.pop().unwrap()
}

pub fn remainder_mod(c: &[u64], q: &[i128], p: u64) -> Vec<u64> {
    let k = q.len() - 1;
    let qm: Vec<u64> = q.iter().map(|&x| x.rem_euclid(p as i128) as u64).collect();
    let mut r: Vec<u64> = c.to_vec();
    while r.len() > k {
        let lead = *r.last().unwrap();
        let shift = r.len() - 1 - k;
        if lead != 0 {
            for i in 0..=k {
                r[shift + i] = sub(r[shift + i], mulmod(lead, qm[i], p), p);
            }
        }
        r.pop();
    }
    r
}

pub fn primes_needed(n: usize) -> usize {
    let bits = n as f64 * (2.0 + 2f64.sqrt()).log2() + 2.0;
    (bits / PRIME_BITS).ceil() as usize
}

pub fn exact_divisors(live: &Live, ps: &[u64]) -> Vec<(char, usize)> {
    let n = live.edges.len();
    let need = primes_needed(n);
    let polys: Vec<Vec<u64>> = ps[..need].iter().map(|&p| charpoly_mod(live, p)).collect();
    let mut out = Vec::new();
    for family in ['A', 'B'] {
        for k in 2..=KMAX.min(n) {
            if family == 'B' && k == 2 {
                continue;
            }
            let q = trinomial(family, k);
            let zero = polys
                .iter()
                .zip(ps)
                .all(|(c, &p)| remainder_mod(c, &q, p).iter().all(|&x| x == 0));
            if zero {
                out.push((family, k));
            }
        }
    }
    out
}

pub fn charpoly_exact(live: &Live, ps: &[u64]) -> Option<Vec<i128>> {
    let n = live.edges.len();
    if n > 70 {
        return None;
    }
    let (p0, p1) = (ps[0], ps[1]);
    let c0 = charpoly_mod(live, p0);
    let c1 = charpoly_mod(live, p1);
    let inv = powmod(p0 % p1, p1 - 2, p1);
    let modulus = p0 as i128 * p1 as i128;
    Some(
        c0.iter()
            .zip(&c1)
            .map(|(&a0, &a1)| {
                let t = mulmod(sub(a1, a0 % p1, p1), inv, p1);
                let v = a0 as i128 + p0 as i128 * t as i128;
                if v > modulus / 2 {
                    v - modulus
                } else {
                    v
                }
            })
            .collect(),
    )
}

pub fn v3(mut x: u64) -> u32 {
    let mut k = 0;
    while x > 0 && x.is_multiple_of(3) {
        x /= 3;
        k += 1;
    }
    k
}

pub struct Row {
    pub a: u64,
    pub b: u64,
    pub states: usize,
    pub rho: f64,
    pub width: f64,
    pub matched: Option<(char, usize, f64)>,
    pub exact: Vec<(char, usize)>,
    pub returns: Option<Vec<u128>>,
    pub renewal: Option<bool>,
    pub charpoly: Option<Vec<i128>>,
}

pub fn row(a: u64, b: u64, ps: &[u64]) -> Row {
    let live = live(a, b);
    let n = live.edges.len();
    let (rho, lo, hi) = perron(&live);
    let exact = if n >= 2 {
        exact_divisors(&live, ps)
    } else {
        Vec::new()
    };
    let mut near: Vec<(char, usize, f64)> = Vec::new();
    for family in ['A', 'B'] {
        for k in 2..=KMAX {
            if family == 'B' && k == 2 {
                continue;
            }
            let theta = root(family, k);
            if (rho - theta).abs() < MATCH {
                near.push((family, k, theta));
            }
        }
    }
    let matched = near
        .iter()
        .find(|m| exact.contains(&(m.0, m.1)))
        .or(near.first())
        .copied();
    let charpoly = if n >= 2 {
        charpoly_exact(&live, ps)
    } else {
        None
    };
    let returns = if n >= 2 { first_returns(&live) } else { None };
    let renewal = match (&returns, matched) {
        (Some(f), Some((family, k, _))) => {
            renewal_polynomial(f).map(|p| divides(&trinomial(family, k), &p))
        }
        _ => None,
    };
    Row {
        a,
        b,
        states: n,
        rho,
        width: hi - lo,
        matched,
        exact,
        returns,
        renewal,
        charpoly,
    }
}

pub fn components(edges: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = edges.len();
    let mut back = vec![Vec::new(); n];
    for (i, out) in edges.iter().enumerate() {
        for &j in out {
            back[j].push(i);
        }
    }
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for s in 0..n {
        if seen[s] {
            continue;
        }
        let mut stack = vec![(s, 0usize)];
        seen[s] = true;
        while let Some(&mut (v, ref mut k)) = stack.last_mut() {
            if *k < edges[v].len() {
                let w = edges[v][*k];
                *k += 1;
                if !seen[w] {
                    seen[w] = true;
                    stack.push((w, 0));
                }
            } else {
                order.push(v);
                stack.pop();
            }
        }
    }
    let mut part = vec![usize::MAX; n];
    let mut out = Vec::new();
    for &s in order.iter().rev() {
        if part[s] != usize::MAX {
            continue;
        }
        let id = out.len();
        let mut comp = vec![s];
        part[s] = id;
        let mut stack = vec![s];
        while let Some(v) = stack.pop() {
            for &w in &back[v] {
                if part[w] == usize::MAX {
                    part[w] = id;
                    comp.push(w);
                    stack.push(w);
                }
            }
        }
        out.push(comp);
    }
    out
}

pub fn box_rho(a: u64, b: u64) -> (f64, usize, usize) {
    let (a, b) = (a as usize, b as usize);
    let n = a * b;
    let mut edges: Vec<Vec<usize>> = vec![Vec::new(); n];
    for c1 in 0..a {
        for c2 in 0..b {
            for d in 0..3 {
                let (x, y) = (a * d + c1, b * d + c2);
                if matches!((x % 3, y % 3), (0, 0) | (1, 0) | (0, 1)) {
                    edges[c1 * b + c2].push((x / 3) * b + y / 3);
                }
            }
        }
    }
    let comps = components(&edges);
    let mut rho = 0.0f64;
    for comp in &comps {
        let place: std::collections::HashMap<usize, usize> =
            comp.iter().enumerate().map(|(k, &v)| (v, k)).collect();
        let local: Vec<Vec<usize>> = comp
            .iter()
            .map(|&v| {
                edges[v]
                    .iter()
                    .filter_map(|w| place.get(w).copied())
                    .collect()
            })
            .collect();
        if local.iter().all(|out| out.is_empty()) {
            continue;
        }
        let m = local.len();
        let mut v = vec![1.0f64; m];
        let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
        for _ in 0..200_000 {
            let mut w: Vec<f64> = v.clone();
            for (i, out) in local.iter().enumerate() {
                for &j in out {
                    w[i] += v[j];
                }
            }
            lo = f64::INFINITY;
            hi = 0.0;
            for i in 0..m {
                let r = w[i] / v[i];
                lo = lo.min(r);
                hi = hi.max(r);
            }
            let norm: f64 = w.iter().sum();
            for x in w.iter_mut() {
                *x /= norm;
            }
            v = w;
            if hi - lo < 1e-12 {
                break;
            }
        }
        rho = rho.max((lo + hi) / 2.0 - 1.0);
    }
    (rho, n, comps.len())
}

pub fn single(a: u64, b: u64) {
    let live = live(a, b);
    let ps = primes(primes_needed(live.edges.len()).max(2));
    println!("{}", describe(&row(a, b, &ps)));
    let (rho, states, comps) = box_rho(a, b);
    println!(
        "whole-box rho {:.10} on {} carry states in {} components, mu <= rho",
        rho, states, comps
    );
}

fn describe(row: &Row) -> String {
    let fam = match row.matched {
        Some((f, k, theta)) => format!(
            "{}{} theta {:.13} gap {:.1e}",
            f,
            k,
            theta,
            (row.rho - theta).abs()
        ),
        None => "neither".to_string(),
    };
    let exact: Vec<String> = row
        .exact
        .iter()
        .map(|(f, k)| format!("{}{}", f, k))
        .collect();
    let returns = match &row.returns {
        Some(f) => f
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(j, c)| format!("{}:{}", j, c))
            .collect::<Vec<_>>()
            .join(","),
        None => "cyclic".to_string(),
    };
    let renewal = match row.renewal {
        Some(true) => " renewal-divides",
        Some(false) => " renewal-fails",
        None => "",
    };
    let charpoly = match &row.charpoly {
        Some(c) => format!(
            " charpoly [{}]",
            c.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
        None => String::new(),
    };
    format!(
        "ray {} {} states {} mu {:.13} width {:.1e} match {} exact [{}] returns {}{}{}",
        row.a,
        row.b,
        row.states,
        row.rho,
        row.width,
        fam,
        exact.join(" "),
        returns,
        renewal,
        charpoly
    )
}

pub fn report(height: u64, threads: usize) {
    let ps = primes(primes_needed(height as usize + 2).max(2));
    let mut rays: Vec<(u64, u64)> = Vec::new();
    for b in 1..=height {
        for a in 1..=b {
            if gcd(a as u128, b as u128) == 1 {
                rays.push((a, b));
            }
        }
    }
    let total = rays.len();
    let rows: Vec<Row> = std::thread::scope(|scope| {
        let chunks: Vec<&[(u64, u64)]> = rays.chunks(total.div_ceil(threads)).collect();
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                let ps = &ps;
                scope.spawn(move || {
                    chunk
                        .iter()
                        .map(|&(a, b)| row(a, b, ps))
                        .collect::<Vec<Row>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let supergolden = root('A', 3);
    println!(
        "mass-rate census, mu the Perron root of the component of 0, box 1 <= a <= b <= {} primitive rays {} (swap de-duplicated, fibres excluded)",
        height, total
    );
    println!(
        "match tolerance {:.0e}, families A: x^k = x^(k-1) + 1, B: x^k = x + 1, k <= {}, exact: trinomial divides the live charpoly by modular remainder under the bound (2 + sqrt 2)^n with {} primes near 2^61",
        MATCH,
        KMAX,
        ps.len()
    );
    let mut no_three = 0usize;
    let mut congruence_fails = 0usize;
    let mut empty = 0usize;
    let mut occupied = Vec::new();
    for r in &rows {
        if r.a % 3 != 0 && r.b % 3 != 0 {
            no_three += 1;
        } else {
            let (p, q) = if r.a % 3 == 0 { (r.b, r.a) } else { (r.a, r.b) };
            let q1 = q / 3u64.pow(v3(q));
            if q1 % 3 != p % 3 {
                congruence_fails += 1;
            }
        }
        if r.states == 1 {
            empty += 1;
        } else {
            occupied.push(r);
        }
    }
    println!(
        "classes: 3 divides neither coordinate {}, congruence q1 != p mod 3 {}, mu = 1 (one live state) {}, occupied {}",
        no_three,
        congruence_fails,
        empty,
        occupied.len()
    );
    let mut shifts = 0usize;
    let mut top_off: Vec<&Row> = Vec::new();
    let mut certified = 0usize;
    let mut matched_only = 0usize;
    let mut neither = 0usize;
    let mut dag = 0usize;
    let mut worst_width = 0.0f64;
    for r in &occupied {
        let shift = r.a == 1 && r.b == 3u64.pow(v3(r.b)) && r.b > 1;
        if shift {
            shifts += 1;
        } else {
            top_off.push(r);
        }
        worst_width = worst_width.max(r.width);
        if r.returns.is_some() {
            dag += 1;
        }
        match r.matched {
            Some((f, k, _)) if r.exact.contains(&(f, k)) => certified += 1,
            Some(_) => matched_only += 1,
            None => neither += 1,
        }
    }
    println!(
        "occupied {}: certified (match and exact) {}, matched without exact {}, neither {}, acyclic off the start {}, shifts {}, worst Perron bracket {:.1e}",
        occupied.len(),
        certified,
        matched_only,
        neither,
        dag,
        shifts,
        worst_width
    );
    top_off.sort_by(|x, y| y.rho.partial_cmp(&x.rho).unwrap());
    let above: Vec<&&Row> = top_off
        .iter()
        .filter(|r| r.rho > supergolden + 1e-9)
        .collect();
    println!(
        "off-shift rays above the supergolden {:.13}: {}",
        supergolden,
        above.len()
    );
    for r in above {
        println!("WITNESS {}", describe(r));
    }
    println!("off-shift rays within 0.02 of phi {:.13}:", phi);
    for r in top_off.iter().filter(|r| (r.rho - phi).abs() < 0.02) {
        println!("  {}", describe(r));
    }
    println!("rays within 0.02 of the supergolden:");
    for r in top_off
        .iter()
        .filter(|r| (r.rho - supergolden).abs() < 0.02)
    {
        println!("  {}", describe(r));
    }
    println!("shift rays:");
    for r in occupied
        .iter()
        .filter(|r| r.a == 1 && r.b > 1 && r.b == 3u64.pow(v3(r.b)))
    {
        println!("  {}", describe(r));
    }
    for (a, b) in [(1u64, 12u64), (3, 10), (4, 9)] {
        let l = live(a, b);
        let poly = charpoly_exact(&l, &ps).unwrap();
        println!(
            "attainer ({},{}) states {:?} charpoly low to high {:?}",
            a, b, l.states, poly
        );
    }
    println!("rows:");
    for r in &occupied {
        println!("{}", describe(r));
    }
}
