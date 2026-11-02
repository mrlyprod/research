pub struct Tables {
    pub k: usize,
    pub binom: Vec<Vec<i128>>,
    pub kraw: Vec<Vec<i128>>,
}

pub struct Profile {
    pub k: usize,
    pub degree: i32,
    pub ones: u128,
    pub nonlinearity: u128,
    pub walsh_max: u128,
    pub autocorrelation_max: u128,
    pub differential: u128,
    pub immunity: usize,
}

// TABLES

pub fn pascal(n: usize) -> Vec<Vec<i128>> {
    let mut rows: Vec<Vec<i128>> = vec![vec![1]];
    for i in 1..=n {
        let mut row = vec![1i128; i + 1];
        for j in 1..i {
            row[j] = rows[i - 1][j - 1] + rows[i - 1][j];
        }
        rows.push(row);
    }
    rows
}

pub fn krawtchouk(k: usize) -> Vec<Vec<i128>> {
    (0..=k)
        .map(|v| {
            let mut poly = vec![0i128; k + 1];
            poly[0] = 1;
            for _ in 0..k - v {
                for i in (1..=k).rev() {
                    poly[i] += poly[i - 1];
                }
            }
            for _ in 0..v {
                for i in (1..=k).rev() {
                    poly[i] -= poly[i - 1];
                }
            }
            poly
        })
        .collect()
}

impl Tables {
    pub fn new(k: usize) -> Tables {
        Tables {
            k,
            binom: pascal(k),
            kraw: krawtchouk(k),
        }
    }
}

pub fn ceiling(k: usize) -> u128 {
    let half = 1u128 << (k - 1);
    if k % 2 == 0 {
        half - (1u128 << (k / 2 - 1))
    } else {
        half - (1u128 << ((k - 1) / 2))
    }
}

// ANF

pub fn anf(g: &[u8]) -> Vec<u8> {
    (0..g.len())
        .map(|w| {
            let mut acc = 0u8;
            let mut j = w;
            loop {
                acc ^= g[j];
                if j == 0 {
                    break;
                }
                j = (j - 1) & w;
            }
            acc
        })
        .collect()
}

pub fn degree(c: &[u8]) -> i32 {
    c.iter().rposition(|&bit| bit == 1).map_or(-1, |w| w as i32)
}

// SPECTRA

pub fn ones(t: &Tables, g: &[u8]) -> u128 {
    (0..=t.k)
        .filter(|&w| g[w] == 1)
        .map(|w| t.binom[t.k][w] as u128)
        .sum()
}

pub fn walsh(t: &Tables, g: &[u8]) -> Vec<i128> {
    (0..=t.k)
        .map(|v| {
            (0..=t.k)
                .map(|w| t.kraw[v][w] * (1 - 2 * g[w] as i128))
                .sum()
        })
        .collect()
}

pub fn autocorrelation(t: &Tables, g: &[u8]) -> Vec<i128> {
    let k = t.k;
    (0..=k)
        .map(|u| {
            let mut total = 0i128;
            for i in 0..=u {
                for j in 0..=k - u {
                    let sign = 1 - 2 * ((g[i + j] ^ g[u - i + j]) as i128);
                    total += t.binom[u][i] * t.binom[k - u][j] * sign;
                }
            }
            total
        })
        .collect()
}

pub fn profile(t: &Tables, g: &[u8]) -> Profile {
    let k = t.k;
    let w = walsh(t, g);
    let a = autocorrelation(t, g);
    let walsh_max = w.iter().map(|x| x.unsigned_abs()).max().unwrap_or(0);
    let autocorrelation_max = a[1..].iter().map(|x| x.unsigned_abs()).max().unwrap_or(0);
    Profile {
        k,
        degree: degree(&anf(g)),
        ones: ones(t, g),
        nonlinearity: (1u128 << (k - 1)) - walsh_max / 2,
        walsh_max,
        autocorrelation_max,
        differential: ((1u128 << k) + autocorrelation_max) / 2,
        immunity: (1..=k).take_while(|&v| w[v] == 0).count(),
    }
}

pub fn balanced(p: &Profile) -> bool {
    p.ones == 1u128 << (p.k - 1)
}

pub fn parseval(t: &Tables, g: &[u8]) -> bool {
    let w = walsh(t, g);
    let total: i128 = (0..=t.k).map(|v| t.binom[t.k][v] * w[v] * w[v]).sum();
    total == 1i128 << (2 * t.k)
}

// BRUTE FORCE

pub struct Table {
    pub degree: i32,
    pub ones: u128,
    pub nonlinearity: u128,
    pub autocorrelation_max: u128,
    pub differential: u128,
}

fn fwht(values: &mut [i64]) {
    let mut h = 1;
    while h < values.len() {
        for start in (0..values.len()).step_by(2 * h) {
            for i in start..start + h {
                let (a, b) = (values[i], values[i + h]);
                values[i] = a + b;
                values[i + h] = a - b;
            }
        }
        h *= 2;
    }
}

pub fn table(values: &[u8], k: usize) -> Table {
    let size = 1usize << k;
    assert_eq!(values.len(), size);
    let mut spectrum: Vec<i64> = values.iter().map(|&bit| 1 - 2 * bit as i64).collect();
    fwht(&mut spectrum);
    let walsh_max = spectrum.iter().map(|x| x.unsigned_abs()).max().unwrap() as u128;
    let mut squares: Vec<i64> = spectrum.iter().map(|x| x * x).collect();
    fwht(&mut squares);
    let autocorrelation_max = squares[1..]
        .iter()
        .map(|x| (x >> k).unsigned_abs())
        .max()
        .unwrap_or(0) as u128;
    let mut anf = values.to_vec();
    for bit in 0..k {
        for x in 0..size {
            if x >> bit & 1 == 1 {
                anf[x] ^= anf[x ^ (1 << bit)];
            }
        }
    }
    let degree = (0..size)
        .filter(|&x| anf[x] == 1)
        .map(|x| x.count_ones() as i32)
        .max()
        .unwrap_or(-1);
    Table {
        degree,
        ones: values.iter().map(|&bit| bit as u128).sum(),
        nonlinearity: (1u128 << (k - 1)) - walsh_max / 2,
        autocorrelation_max,
        differential: ((1u128 << k) + autocorrelation_max) / 2,
    }
}

pub fn symmetric_table(g: &[u8], k: usize) -> Vec<u8> {
    (0..1usize << k)
        .map(|x| g[x.count_ones() as usize])
        .collect()
}
