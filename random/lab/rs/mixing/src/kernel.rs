use crate::catalog::Round;

// KERNEL

pub struct Kernel {
    pub m: usize,
    pub cb: Vec<u32>,
    pub la: Vec<u32>,
    pub start: Vec<usize>,
    pub ld: Vec<u32>,
    pub count: Vec<u32>,
    pub w: Vec<f64>,
}

pub fn wht(v: &mut [i64]) {
    let mut h = 1;
    while h < v.len() {
        for i in (0..v.len()).step_by(2 * h) {
            for j in i..i + h {
                let (a, b) = (v[j], v[j + h]);
                v[j] = a + b;
                v[j + h] = a - b;
            }
        }
        h *= 2;
    }
}

pub fn spots(mask: u32, bits: usize) -> Vec<usize> {
    (0..bits).filter(|&p| mask >> p & 1 == 1).collect()
}

pub fn expand(y: usize, cells: &[usize]) -> u32 {
    cells
        .iter()
        .enumerate()
        .fold(0, |m, (k, &p)| m | ((y >> k & 1) as u32) << p)
}

pub fn compress(x: u32, cells: &[usize]) -> usize {
    cells
        .iter()
        .enumerate()
        .fold(0, |m, (k, &p)| m | ((x >> p & 1) as usize) << k)
}

impl Kernel {
    pub fn new(round: &Round) -> Kernel {
        let bits = round.bits;
        let m = 1usize << bits;
        let acells = spots(round.fill, bits);
        let bcells = spots(round.void(), bits);
        let (ka, kb) = (acells.len(), bcells.len());
        let (sa, sb) = (1usize << ka, 1usize << kb);
        let mc: Vec<usize> = (0..sb)
            .map(|y| compress(round.mess(expand(y, &bcells)), &acells))
            .collect();
        let mut ac = vec![0i64; sa * sb];
        let mut f = vec![0i64; sb];
        for u in 0..sa {
            for y in 0..sb {
                f[y] = if (u & mc[y]).count_ones() % 2 == 0 {
                    1
                } else {
                    -1
                };
            }
            wht(&mut f);
            for x in f.iter_mut() {
                *x *= *x;
            }
            wht(&mut f);
            for a in 0..sb {
                ac[u * sb + a] = f[a] / sb as i64;
            }
        }
        let mut start = vec![0usize; sb + 1];
        let (mut ld, mut count) = (Vec::new(), Vec::new());
        let mut g = vec![0i64; sa];
        for a in 0..sb {
            for u in 0..sa {
                g[u] = ac[u * sb + a];
            }
            wht(&mut g);
            let mut total = 0i64;
            for (d, &c) in g.iter().enumerate() {
                assert!(
                    c >= 0 && c % sa as i64 == 0,
                    "a count that is not an integer"
                );
                let c = c / sa as i64;
                if c > 0 {
                    ld.push(round.spread(expand(d, &acells)));
                    count.push(c as u32);
                    total += c;
                }
            }
            assert_eq!(total, sb as i64, "a row that does not sum to the void");
            start[a + 1] = ld.len();
        }
        let w = count.iter().map(|&c| c as f64 / sb as f64).collect();
        Kernel {
            m,
            cb: (0..m as u32).map(|x| compress(x, &bcells) as u32).collect(),
            la: (0..m as u32).map(|x| round.spread(x)).collect(),
            start,
            ld,
            count,
            w,
        }
    }

    pub fn row(&self, alpha: usize) -> std::ops::Range<usize> {
        let a = self.cb[alpha] as usize;
        self.start[a]..self.start[a + 1]
    }

    pub fn edges(&self) -> usize {
        (1..self.m).map(|a| self.row(a).len()).sum()
    }

    pub fn forward(&self, v: &[f64], out: &mut [f64]) {
        out.fill(0.0);
        for alpha in 1..self.m {
            let x = v[alpha];
            if x == 0.0 {
                continue;
            }
            let base = self.la[alpha];
            for k in self.row(alpha) {
                out[(base ^ self.ld[k]) as usize] += x * self.w[k];
            }
        }
    }

    pub fn backward(&self, v: &[f64], out: &mut [f64]) {
        out[0] = 0.0;
        for alpha in 1..self.m {
            let base = self.la[alpha];
            out[alpha] = self
                .row(alpha)
                .map(|k| self.w[k] * v[(base ^ self.ld[k]) as usize])
                .sum();
        }
    }

    pub fn classes(&self) -> Vec<usize> {
        let mut label = vec![usize::MAX; self.m];
        let mut sizes = Vec::new();
        let mut stack = Vec::new();
        for s in 1..self.m {
            if label[s] != usize::MAX {
                continue;
            }
            let id = sizes.len();
            label[s] = id;
            stack.push(s);
            let mut size = 0;
            while let Some(a) = stack.pop() {
                size += 1;
                let base = self.la[a];
                for k in self.row(a) {
                    let b = (base ^ self.ld[k]) as usize;
                    if label[b] == usize::MAX {
                        label[b] = id;
                        stack.push(b);
                    }
                }
            }
            sizes.push(size);
        }
        sizes
    }

    pub fn irreducible(&self) -> bool {
        let mut seen = vec![false; self.m];
        let mut stack = vec![1usize];
        seen[1] = true;
        let mut reached = 1;
        while let Some(a) = stack.pop() {
            let base = self.la[a];
            for k in self.row(a) {
                let b = (base ^ self.ld[k]) as usize;
                if !seen[b] {
                    seen[b] = true;
                    reached += 1;
                    stack.push(b);
                }
            }
        }
        reached == self.m - 1
    }

    pub fn direct(&self, h: &[u32], alpha: usize) -> Vec<(u32, u32)> {
        let mut c = vec![0u32; self.m];
        for x in 0..self.m {
            c[(h[x] ^ h[x ^ alpha]) as usize] += 1;
        }
        let mut out = Vec::new();
        for (b, &v) in c.iter().enumerate() {
            if v > 0 {
                out.push((b as u32, v));
            }
        }
        out
    }

    pub fn sparse(&self, alpha: usize) -> Vec<(u32, u32)> {
        let scale = (self.m as u32) / (self.start.len() as u32 - 1);
        let base = self.la[alpha];
        let mut out: Vec<(u32, u32)> = self
            .row(alpha)
            .map(|k| (base ^ self.ld[k], self.count[k] * scale))
            .collect();
        out.sort();
        out
    }
}

// SPECTRUM

pub fn deflate(v: &mut [f64]) {
    let n = (v.len() - 1) as f64;
    let mean = v[1..].iter().sum::<f64>() / n;
    v[0] = 0.0;
    for x in v[1..].iter_mut() {
        *x -= mean;
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

pub fn seeded(m: usize, seed: u64) -> Vec<f64> {
    let mut s = seed | 1;
    let mut v: Vec<f64> = (0..m)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        })
        .collect();
    deflate(&mut v);
    let r = norm(&v);
    v.iter_mut().for_each(|x| *x /= r);
    v
}

fn sturm(alpha: &[f64], beta: &[f64], x: f64) -> usize {
    let mut count = 0;
    let mut q = 1.0;
    for i in 0..alpha.len() {
        let b2 = if i == 0 {
            0.0
        } else {
            beta[i - 1] * beta[i - 1]
        };
        q = alpha[i] - x - if i == 0 { 0.0 } else { b2 / q };
        if q == 0.0 {
            q = 1e-300;
        }
        if q < 0.0 {
            count += 1;
        }
    }
    count
}

fn top(alpha: &[f64], beta: &[f64]) -> f64 {
    let k = alpha.len();
    let bound = (0..k)
        .map(|i| {
            alpha[i].abs()
                + if i > 0 { beta[i - 1].abs() } else { 0.0 }
                + if i + 1 < k { beta[i].abs() } else { 0.0 }
        })
        .fold(0.0, f64::max);
    let (mut lo, mut hi) = (-bound, bound);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if sturm(alpha, beta, mid) < k {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn ritz(alpha: &[f64], beta: &[f64], theta: f64) -> Vec<f64> {
    let k = alpha.len();
    let shift = theta + 1e-10 * (1.0 + theta.abs());
    let mut x = vec![1.0; k];
    for _ in 0..3 {
        let mut c = vec![0.0; k];
        let mut d = x.clone();
        let mut diag: Vec<f64> = alpha.iter().map(|a| a - shift).collect();
        for i in 1..k {
            let f = beta[i - 1] / diag[i - 1];
            diag[i] -= f * beta[i - 1];
            d[i] -= f * d[i - 1];
            c[i - 1] = beta[i - 1];
        }
        let mut y = vec![0.0; k];
        y[k - 1] = d[k - 1] / diag[k - 1];
        for i in (0..k - 1).rev() {
            y[i] = (d[i] - c[i] * y[i + 1]) / diag[i];
        }
        let r = norm(&y);
        x = y.iter().map(|v| v / r).collect();
    }
    x
}

pub struct Spectrum {
    pub sigma: f64,
    pub residual: f64,
    pub steps: usize,
}

pub fn sigma2(apply: &dyn Fn(&[f64], &mut [f64]), m: usize, seed: u64, most: usize) -> Spectrum {
    let mut basis: Vec<Vec<f64>> = vec![seeded(m, seed)];
    let (mut alpha, mut beta) = (Vec::new(), Vec::new());
    let mut w = vec![0.0; m];
    let mut last = -1.0;
    let mut theta = 0.0;
    for step in 0..most {
        apply(&basis[step], &mut w);
        deflate(&mut w);
        let a = dot(&w, &basis[step]);
        alpha.push(a);
        for _ in 0..2 {
            for q in &basis {
                let c = dot(&w, q);
                w.iter_mut().zip(q).for_each(|(x, y)| *x -= c * y);
            }
        }
        theta = top(&alpha, &beta);
        let b = norm(&w);
        if (theta - last).abs() < 1e-14 * theta.abs().max(1e-300) && step >= 8 || b < 1e-13 {
            break;
        }
        last = theta;
        if step + 1 == most {
            break;
        }
        beta.push(b);
        basis.push(w.iter().map(|x| x / b).collect());
    }
    let k = alpha.len();
    let s = ritz(&alpha, &beta[..k - 1], theta);
    let mut y = vec![0.0; m];
    for (i, q) in basis.iter().take(k).enumerate() {
        y.iter_mut().zip(q).for_each(|(a, b)| *a += s[i] * b);
    }
    deflate(&mut y);
    let r = norm(&y);
    y.iter_mut().for_each(|x| *x /= r);
    apply(&y, &mut w);
    deflate(&mut w);
    let rayleigh = dot(&w, &y);
    let residual = w
        .iter()
        .zip(&y)
        .map(|(a, b)| (a - rayleigh * b).powi(2))
        .sum::<f64>()
        .sqrt();
    Spectrum {
        sigma: rayleigh.max(0.0).sqrt(),
        residual,
        steps: k,
    }
}

pub fn rate(
    apply: &dyn Fn(&[f64], &mut [f64]),
    m: usize,
    seed: u64,
    burn: usize,
    window: usize,
) -> f64 {
    let mut v = seeded(m, seed);
    let mut w = vec![0.0; m];
    let mut log = 0.0;
    for step in 0..burn + window {
        apply(&v, &mut w);
        deflate(&mut w);
        let r = norm(&w);
        if r == 0.0 {
            return 0.0;
        }
        if step >= burn {
            log += r.ln();
        }
        v.iter_mut().zip(&w).for_each(|(a, b)| *a = b / r);
    }
    (log / window as f64).exp()
}

// ARNOLDI

fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

fn hqr(h: &[Vec<f64>], n: usize) -> Vec<(f64, f64)> {
    let mut a = vec![vec![0.0; n + 1]; n + 1];
    for i in 1..=n {
        for j in 1..=n {
            a[i][j] = h[i - 1][j - 1];
        }
    }
    let (mut wr, mut wi) = (vec![0.0; n + 1], vec![0.0; n + 1]);
    let mut anorm = 0.0;
    for i in 1..=n {
        for j in i.max(2) - 1..=n {
            anorm += a[i][j].abs();
        }
    }
    let mut nn = n;
    let mut t = 0.0;
    while nn >= 1 {
        let mut its = 0;
        loop {
            let mut l = nn;
            while l >= 2 {
                let mut s = a[l - 1][l - 1].abs() + a[l][l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[l][l - 1].abs() + s == s {
                    a[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            let mut x = a[nn][nn];
            if l == nn {
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
            } else {
                let mut y = a[nn - 1][nn - 1];
                let mut w = a[nn][nn - 1] * a[nn - 1][nn];
                if l == nn - 1 {
                    let p = 0.5 * (y - x);
                    let q = p * p + w;
                    let mut z = q.abs().sqrt();
                    x += t;
                    if q >= 0.0 {
                        z = p + sign(z, p);
                        wr[nn - 1] = x + z;
                        wr[nn] = if z != 0.0 { x - w / z } else { x + z };
                        wi[nn - 1] = 0.0;
                        wi[nn] = 0.0;
                    } else {
                        wr[nn - 1] = x + p;
                        wr[nn] = x + p;
                        wi[nn - 1] = -z;
                        wi[nn] = z;
                    }
                    nn -= 2;
                } else {
                    assert!(its < 60, "hqr did not converge");
                    if its > 0 && its % 10 == 0 {
                        t += x;
                        for i in 1..=nn {
                            a[i][i] -= x;
                        }
                        let s = a[nn][nn - 1].abs() + a[nn - 1][nn - 2].abs();
                        x = 0.75 * s;
                        y = x;
                        w = -0.4375 * s * s;
                    }
                    its += 1;
                    let mut m = nn - 2;
                    let (mut p, mut q, mut r);
                    loop {
                        let z = a[m][m];
                        let (u, v) = (x - z, y - z);
                        p = (u * v - w) / a[m + 1][m] + a[m][m + 1];
                        q = a[m + 1][m + 1] - z - u - v;
                        r = a[m + 2][m + 1];
                        let s = p.abs() + q.abs() + r.abs();
                        p /= s;
                        q /= s;
                        r /= s;
                        if m == l {
                            break;
                        }
                        let u = a[m][m - 1].abs() * (q.abs() + r.abs());
                        let v = p.abs() * (a[m - 1][m - 1].abs() + z.abs() + a[m + 1][m + 1].abs());
                        if u + v == v {
                            break;
                        }
                        m -= 1;
                    }
                    for i in m + 2..=nn {
                        a[i][i - 2] = 0.0;
                        if i != m + 2 {
                            a[i][i - 3] = 0.0;
                        }
                    }
                    for k in m..nn {
                        if k != m {
                            p = a[k][k - 1];
                            q = a[k + 1][k - 1];
                            r = if k != nn - 1 { a[k + 2][k - 1] } else { 0.0 };
                            x = p.abs() + q.abs() + r.abs();
                            if x != 0.0 {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        let s = sign((p * p + q * q + r * r).sqrt(), p);
                        if s == 0.0 {
                            continue;
                        }
                        if k == m {
                            if l != m {
                                a[k][k - 1] = -a[k][k - 1];
                            }
                        } else {
                            a[k][k - 1] = -s * x;
                        }
                        p += s;
                        x = p / s;
                        y = q / s;
                        let z = r / s;
                        q /= p;
                        r /= p;
                        for j in k..=nn {
                            let mut p = a[k][j] + q * a[k + 1][j];
                            if k != nn - 1 {
                                p += r * a[k + 2][j];
                                a[k + 2][j] -= p * z;
                            }
                            a[k + 1][j] -= p * y;
                            a[k][j] -= p * x;
                        }
                        for i in l..=nn.min(k + 3) {
                            let mut p = x * a[i][k] + y * a[i][k + 1];
                            if k != nn - 1 {
                                p += z * a[i][k + 2];
                                a[i][k + 2] -= p * r;
                            }
                            a[i][k + 1] -= p * q;
                            a[i][k] -= p;
                        }
                    }
                }
            }
            if l + 1 >= nn {
                break;
            }
        }
    }
    (1..=n).map(|i| (wr[i], wi[i])).collect()
}

type C = (f64, f64);

fn cmul(a: C, b: C) -> C {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

fn cdiv(a: C, b: C) -> C {
    let d = b.0 * b.0 + b.1 * b.1;
    ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
}

fn eigvec(h: &[Vec<f64>], k: usize, mu: C) -> Vec<C> {
    let mut y: Vec<C> = vec![(1.0, 0.0); k];
    for _ in 0..3 {
        let mut a: Vec<Vec<C>> = (0..k)
            .map(|i| {
                (0..k)
                    .map(|j| {
                        (
                            h[i][j] - if i == j { mu.0 } else { 0.0 },
                            if i == j { -mu.1 } else { 0.0 },
                        )
                    })
                    .collect()
            })
            .collect();
        let mut b = y.clone();
        for c in 0..k {
            let piv = (c..k)
                .max_by(|&i, &j| {
                    (a[i][c].0.hypot(a[i][c].1)).total_cmp(&a[j][c].0.hypot(a[j][c].1))
                })
                .unwrap();
            a.swap(c, piv);
            b.swap(c, piv);
            if a[c][c] == (0.0, 0.0) {
                a[c][c] = (1e-300, 0.0);
            }
            for i in c + 1..k {
                let f = cdiv(a[i][c], a[c][c]);
                if f == (0.0, 0.0) {
                    continue;
                }
                for j in c..k {
                    let g = cmul(f, a[c][j]);
                    a[i][j] = (a[i][j].0 - g.0, a[i][j].1 - g.1);
                }
                let g = cmul(f, b[c]);
                b[i] = (b[i].0 - g.0, b[i].1 - g.1);
            }
        }
        for i in (0..k).rev() {
            let mut s = b[i];
            for j in i + 1..k {
                let g = cmul(a[i][j], y[j]);
                s = (s.0 - g.0, s.1 - g.1);
            }
            y[i] = cdiv(s, a[i][i]);
        }
        let r = y.iter().map(|c| c.0 * c.0 + c.1 * c.1).sum::<f64>().sqrt();
        y.iter_mut().for_each(|c| *c = (c.0 / r, c.1 / r));
    }
    y
}

pub struct Eigen {
    pub modulus: f64,
    pub value: C,
    pub residual: f64,
    pub matvecs: usize,
    pub restarts: usize,
}

pub fn radius(
    apply: &dyn Fn(&[f64], &mut [f64]),
    m: usize,
    seed: u64,
    dim: usize,
    tol: f64,
    most: usize,
) -> Eigen {
    let mut v = seeded(m, seed);
    let mut w = vec![0.0; m];
    let mut matvecs = 0;
    let mut best = Eigen {
        modulus: 0.0,
        value: (0.0, 0.0),
        residual: f64::INFINITY,
        matvecs: 0,
        restarts: 0,
    };
    for restart in 0..most {
        let mut basis = vec![v.clone()];
        let mut h = vec![vec![0.0; dim]; dim + 1];
        let mut k = dim;
        for j in 0..dim {
            apply(&basis[j], &mut w);
            matvecs += 1;
            deflate(&mut w);
            for _ in 0..2 {
                for (i, q) in basis.iter().enumerate() {
                    let c = dot(&w, q);
                    h[i][j] += c;
                    w.iter_mut().zip(q).for_each(|(x, y)| *x -= c * y);
                }
            }
            let b = norm(&w);
            h[j + 1][j] = b;
            if b < 1e-12 {
                k = j + 1;
                break;
            }
            if j + 1 < dim {
                basis.push(w.iter().map(|x| x / b).collect());
            }
        }
        let lambda = hqr(&h, k)
            .into_iter()
            .max_by(|a, b| {
                a.0.hypot(a.1)
                    .total_cmp(&b.0.hypot(b.1))
                    .then(a.1.total_cmp(&b.1))
            })
            .unwrap();
        let shift = 1e-12 * (1.0 + lambda.0.hypot(lambda.1));
        let y = eigvec(&h, k, (lambda.0 + shift, lambda.1));
        let (mut xr, mut xi) = (vec![0.0; m], vec![0.0; m]);
        for (i, q) in basis.iter().take(k).enumerate() {
            for t in 0..m {
                xr[t] += y[i].0 * q[t];
                xi[t] += y[i].1 * q[t];
            }
        }
        let r = (dot(&xr, &xr) + dot(&xi, &xi)).sqrt();
        xr.iter_mut().for_each(|x| *x /= r);
        xi.iter_mut().for_each(|x| *x /= r);
        let (mut pr, mut pi) = (vec![0.0; m], vec![0.0; m]);
        apply(&xr, &mut pr);
        apply(&xi, &mut pi);
        matvecs += 2;
        deflate(&mut pr);
        deflate(&mut pi);
        let (lr, li) = lambda;
        let residual = (0..m)
            .map(|t| {
                let a = pr[t] - (lr * xr[t] - li * xi[t]);
                let b = pi[t] - (lr * xi[t] + li * xr[t]);
                a * a + b * b
            })
            .sum::<f64>()
            .sqrt();
        best = Eigen {
            modulus: lr.hypot(li),
            value: lambda,
            residual,
            matvecs,
            restarts: restart,
        };
        if residual <= tol || k < dim {
            break;
        }
        v = xr.iter().zip(&xi).map(|(a, b)| a + b).collect();
        deflate(&mut v);
        let r = norm(&v);
        v.iter_mut().for_each(|x| *x /= r);
    }
    best
}

pub fn tv(v: &[f64]) -> f64 {
    let u = 1.0 / (v.len() - 1) as f64;
    0.5 * v[1..].iter().map(|x| (x - u).abs()).sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radius_finds_the_outer_complex_pair() {
        let op = |v: &[f64], out: &mut [f64]| {
            out[0] = 0.0;
            for i in 1..8 {
                out[i % 7 + 1] = 0.9 * v[i];
            }
            for i in 1..8 {
                out[i] += 0.05 * v[i];
            }
        };
        let angle = std::f64::consts::TAU / 7.0;
        let want = (0.05 + 0.9 * angle.cos()).hypot(0.9 * angle.sin());
        let e = radius(&op, 8, 1, 4, 1e-10, 200);
        assert!((e.modulus - want).abs() < 1e-9 && e.residual <= 1e-10);
    }
}
