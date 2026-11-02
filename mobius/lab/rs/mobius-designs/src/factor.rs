// SIEVES

pub fn small_primes(limit: usize) -> Vec<u64> {
    let mut comp = vec![false; limit + 1];
    let mut primes = Vec::new();
    for i in 2..=limit {
        if !comp[i] {
            primes.push(i as u64);
            let mut j = i * i;
            while j <= limit {
                comp[j] = true;
                j += i;
            }
        }
    }
    primes
}

pub fn mu_sieve(limit: usize) -> Vec<i8> {
    let mut mu = vec![0i8; limit + 1];
    let mut comp = vec![false; limit + 1];
    let mut primes: Vec<u32> = Vec::new();
    if limit >= 1 {
        mu[1] = 1;
    }
    for i in 2..=limit {
        if !comp[i] {
            primes.push(i as u32);
            mu[i] = -1;
        }
        for &p in &primes {
            let ip = i * p as usize;
            if ip > limit {
                break;
            }
            comp[ip] = true;
            if i % p as usize == 0 {
                mu[ip] = 0;
                break;
            }
            mu[ip] = -mu[i];
        }
    }
    mu
}

// MODULAR ARITHMETIC

fn mulmod(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1u64 % n;
    a %= n;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, n);
        }
        a = mulmod(a, a, n);
        e >>= 1;
    }
    r
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

// PRIMALITY

const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in WITNESSES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    'witness: for a in WITNESSES {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

// FACTOR SPLITTING

pub fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r > 0 && r as u128 * r as u128 > n as u128 {
        r -= 1;
    }
    while (r as u128 + 1) * (r as u128 + 1) <= n as u128 {
        r += 1;
    }
    r
}

fn brent(n: u64) -> u64 {
    if n & 1 == 0 {
        return 2;
    }
    let mut c = 1u64;
    loop {
        let step = |x: u64| {
            let y = mulmod(x, x, n) + c;
            if y >= n {
                y - n
            } else {
                y
            }
        };
        let mut y = 2u64;
        let mut x = y;
        let mut ys = y;
        let mut q = 1u64;
        let mut d = 1u64;
        let mut r = 1u64;
        let m = 128u64;
        while d == 1 {
            x = y;
            for _ in 0..r {
                y = step(y);
            }
            let mut k = 0u64;
            while k < r && d == 1 {
                ys = y;
                for _ in 0..m.min(r - k) {
                    y = step(y);
                    q = mulmod(q, x.abs_diff(y), n);
                }
                d = gcd(q, n);
                k += m;
            }
            r <<= 1;
        }
        if d == n {
            loop {
                ys = step(ys);
                d = gcd(x.abs_diff(ys), n);
                if d > 1 {
                    break;
                }
            }
        }
        if d < n {
            return d;
        }
        c += 1;
    }
}

// MOBIUS

pub fn mobius(n: u64, primes: &[u64]) -> i8 {
    if n == 1 {
        return 1;
    }
    let mut m = n;
    let mut parity = 0u32;
    let mut settled = false;
    for &p in primes {
        if m == 1 || p * p > m {
            settled = true;
            break;
        }
        if m.is_multiple_of(p) {
            m /= p;
            if m.is_multiple_of(p) {
                return 0;
            }
            parity ^= 1;
        }
    }
    if m > 1 {
        if settled {
            parity ^= 1;
        } else {
            let mut stack = vec![m];
            let mut found: Vec<u64> = Vec::new();
            while let Some(x) = stack.pop() {
                if is_prime(x) {
                    found.push(x);
                    continue;
                }
                let r = isqrt(x);
                if r * r == x {
                    return 0;
                }
                let d = brent(x);
                stack.push(d);
                stack.push(x / d);
            }
            found.sort_unstable();
            for w in found.windows(2) {
                if w[0] == w[1] {
                    return 0;
                }
            }
            parity ^= found.len() as u32 & 1;
        }
    }
    if parity & 1 == 0 {
        1
    } else {
        -1
    }
}
