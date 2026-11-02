const P: u128 = (1 << 61) - 1;

fn power(mut a: u128, mut e: u128) -> u128 {
    let mut r = 1;
    a %= P;
    while e > 0 {
        if e & 1 == 1 {
            r = r * a % P;
        }
        a = a * a % P;
        e >>= 1;
    }
    r
}

fn massey(s: &[u128]) -> Vec<u128> {
    let (mut c, mut b) = (vec![1u128], vec![1u128]);
    let (mut l, mut m, mut last) = (0usize, 1usize, 1u128);
    for n in 0..s.len() {
        let mut d = s[n] % P;
        for i in 1..=l.min(c.len() - 1) {
            d = (d + c[i] * (s[n - i] % P)) % P;
        }
        if d == 0 {
            m += 1;
            continue;
        }
        let t = c.clone();
        let coef = d * power(last, P - 2) % P;
        if c.len() < b.len() + m {
            c.resize(b.len() + m, 0);
        }
        for i in 0..b.len() {
            c[i + m] = (c[i + m] + P - coef * b[i] % P) % P;
        }
        if 2 * l <= n {
            l = n + 1 - l;
            b = t;
            last = d;
            m = 1;
        } else {
            m += 1;
        }
    }
    c.resize(l + 1, 0);
    c
}

pub struct Fit {
    pub poly: Vec<i128>,
    pub margin: i64,
}

pub fn fit(s: &[u128]) -> Option<Fit> {
    let c = massey(s);
    let poly: Vec<i128> = c
        .iter()
        .map(|&x| {
            if x > P / 2 {
                x as i128 - P as i128
            } else {
                x as i128
            }
        })
        .collect();
    let l = poly.len() - 1;
    for n in l..s.len() {
        let mut acc = 0i128;
        for (i, &k) in poly.iter().enumerate() {
            acc = acc.checked_add(k.checked_mul(s[n - i] as i128)?)?;
        }
        if acc != 0 {
            return None;
        }
    }
    Some(Fit {
        poly,
        margin: s.len() as i64 - 2 * l as i64,
    })
}

fn eval(poly: &[i128], r: i128) -> Option<i128> {
    poly.iter()
        .try_fold(0i128, |acc, &k| acc.checked_mul(r)?.checked_add(k))
}

pub fn roots(poly: &[i128]) -> (Vec<i128>, Vec<i128>) {
    let mut rest = poly.to_vec();
    let mut found = Vec::new();
    for r in (0..=64i128).flat_map(|r| if r == 0 { vec![0] } else { vec![r, -r] }) {
        while rest.len() > 1 && eval(&rest, r) == Some(0) {
            let mut q = Vec::with_capacity(rest.len() - 1);
            let mut acc = 0i128;
            for &k in &rest[..rest.len() - 1] {
                acc = acc * r + k;
                q.push(acc);
            }
            rest = q;
            found.push(r);
        }
    }
    found.sort_by(|a, b| b.cmp(a));
    (found, rest)
}

pub fn show(poly: &[i128]) -> String {
    let d = poly.len() - 1;
    let mut out = String::new();
    for (i, &k) in poly.iter().enumerate() {
        if k == 0 {
            continue;
        }
        let e = d - i;
        let sign = if k < 0 {
            " - "
        } else if out.is_empty() {
            ""
        } else {
            " + "
        };
        let a = k.abs();
        let lead = if out.is_empty() && k < 0 {
            "-".to_string()
        } else {
            sign.to_string()
        };
        let body = match (a, e) {
            (_, 0) => format!("{a}"),
            (1, 1) => "x".to_string(),
            (1, _) => format!("x^{e}"),
            (_, 1) => format!("{a} x"),
            _ => format!("{a} x^{e}"),
        };
        out.push_str(&lead);
        out.push_str(&body);
    }
    out
}
