use crate::mess::krawtchouk;

pub const KFRONT: usize = 26;

#[derive(Clone)]
pub struct Row {
    pub count: u64,
    pub best_walsh: u64,
    pub best_count: u64,
    pub balanced_count: u64,
    pub balanced_best_walsh: u64,
    pub balanced_best_count: u64,
}

pub struct Front {
    pub k: usize,
    pub rows: Vec<Row>,
}

impl Row {
    fn empty() -> Row {
        Row {
            count: 0,
            best_walsh: u64::MAX,
            best_count: 0,
            balanced_count: 0,
            balanced_best_walsh: u64::MAX,
            balanced_best_count: 0,
        }
    }
    fn take(&mut self, walsh_max: u64, balanced: bool) {
        self.count += 1;
        if walsh_max < self.best_walsh {
            self.best_walsh = walsh_max;
            self.best_count = 0;
        }
        if walsh_max == self.best_walsh {
            self.best_count += 1;
        }
        if balanced {
            self.balanced_count += 1;
            if walsh_max < self.balanced_best_walsh {
                self.balanced_best_walsh = walsh_max;
                self.balanced_best_count = 0;
            }
            if walsh_max == self.balanced_best_walsh {
                self.balanced_best_count += 1;
            }
        }
    }
}

pub fn nonlinearity(k: usize, walsh_max: u64) -> u64 {
    (1u64 << (k - 1)) - walsh_max / 2
}

pub fn run(k: usize) -> Front {
    let kraw: Vec<Vec<i64>> = {
        let by_v = krawtchouk(k);
        (0..=k)
            .map(|w| (0..=k).map(|v| by_v[v][w] as i64).collect())
            .collect()
    };
    let sup: Vec<u64> = (0..=k)
        .map(|w| {
            (0..=k)
                .filter(|m| m & w == w)
                .fold(0u64, |acc, m| acc | 1 << m)
        })
        .collect();
    let mut sign = vec![1i64; k + 1];
    let mut spectrum = vec![0i64; k + 1];
    spectrum[0] = 1i64 << k;
    let mut cmask = 0u64;
    let mut rows = vec![Row::empty(); k + 1];
    let total = 1u64 << (k + 1);
    for step in 0..total {
        if step > 0 {
            let b = step.trailing_zeros() as usize;
            let old = sign[b];
            for v in 0..=k {
                spectrum[v] -= 2 * old * kraw[b][v];
            }
            sign[b] = -old;
            cmask ^= sup[b];
        }
        let walsh_max = spectrum.iter().map(|x| x.unsigned_abs()).max().unwrap();
        let degree = if cmask == 0 {
            0
        } else {
            63 - cmask.leading_zeros() as usize
        };
        rows[degree].take(walsh_max, spectrum[0] == 0);
    }
    Front { k, rows }
}

pub fn best(front: &Front) -> (u64, u64) {
    let walsh = front.rows.iter().map(|r| r.best_walsh).min().unwrap();
    let count = front
        .rows
        .iter()
        .filter(|r| r.best_walsh == walsh)
        .map(|r| r.best_count)
        .sum();
    (nonlinearity(front.k, walsh), count)
}
