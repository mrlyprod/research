use crate::catalog::Round;
use crate::kernel::{compress, expand, spots};

// TRIPLE

pub struct Triple {
    pub m: usize,
    pub size: usize,
    sb: usize,
    cb: Vec<usize>,
    la: Vec<u32>,
    ld: Vec<u32>,
}

impl Triple {
    pub fn new(round: &Round) -> Triple {
        let bits = round.bits;
        let m = 1usize << bits;
        let acells = spots(round.fill, bits);
        let bcells = spots(round.void(), bits);
        let sb = 1usize << bcells.len();
        let mess: Vec<u32> = (0..sb).map(|y| round.mess(expand(y, &bcells))).collect();
        let mut ld = vec![0u32; sb * sb];
        for a in 0..sb {
            for y in 0..sb {
                let d = mess[y] ^ mess[y ^ a];
                assert_eq!(d & !round.fill, 0);
                ld[a * sb + y] = round.spread(expand(compress(d, &acells), &acells));
            }
        }
        Triple {
            m,
            size: (m - 1) * (m - 2) + 1,
            sb,
            cb: (0..m as u32).map(|x| compress(x, &bcells)).collect(),
            la: (0..m as u32).map(|x| round.spread(x)).collect(),
            ld,
        }
    }

    pub fn index(&self, a: usize, b: usize) -> usize {
        1 + (a - 1) * (self.m - 2) + if b < a { b - 1 } else { b - 2 }
    }

    fn step(&self, a: usize, y: usize) -> usize {
        (self.la[a] ^ self.ld[self.cb[a] * self.sb + y]) as usize
    }

    pub fn forward(&self, v: &[f64], out: &mut [f64]) {
        out.fill(0.0);
        let w = 1.0 / self.sb as f64;
        for a in 1..self.m {
            for b in 1..self.m {
                if a == b {
                    continue;
                }
                let x = v[self.index(a, b)];
                if x == 0.0 {
                    continue;
                }
                for y in 0..self.sb {
                    out[self.index(self.step(a, y), self.step(b, y))] += x * w;
                }
            }
        }
    }

    pub fn backward(&self, v: &[f64], out: &mut [f64]) {
        out[0] = 0.0;
        let w = 1.0 / self.sb as f64;
        for a in 1..self.m {
            for b in 1..self.m {
                if a == b {
                    continue;
                }
                out[self.index(a, b)] = (0..self.sb)
                    .map(|y| v[self.index(self.step(a, y), self.step(b, y))])
                    .sum::<f64>()
                    * w;
            }
        }
    }

    pub fn marginal(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.m];
        for a in 1..self.m {
            for b in 1..self.m {
                if a != b {
                    out[a] += v[self.index(a, b)];
                }
            }
        }
        out
    }
}
