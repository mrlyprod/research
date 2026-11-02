use crate::word::{observe, series_value, CODES};
use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Frac {
    num: i128,
    den: i128,
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 {
        1
    } else {
        a
    }
}

impl Frac {
    pub fn new(num: i128, den: i128) -> Frac {
        assert!(den != 0, "a fraction needs a nonzero denominator");
        let sign = if den < 0 { -1 } else { 1 };
        let g = gcd(num, den);
        Frac {
            num: sign * num / g,
            den: sign * den / g,
        }
    }

    pub fn int(value: i64) -> Frac {
        Frac {
            num: value as i128,
            den: 1,
        }
    }

    pub fn zero() -> Frac {
        Frac { num: 0, den: 1 }
    }

    pub fn is_zero(&self) -> bool {
        self.num == 0
    }

    pub fn add(&self, other: &Frac) -> Frac {
        let left = self.num.checked_mul(other.den).expect("no overflow");
        let right = other.num.checked_mul(self.den).expect("no overflow");
        Frac::new(
            left.checked_add(right).expect("no overflow"),
            self.den.checked_mul(other.den).expect("no overflow"),
        )
    }

    pub fn mul(&self, other: &Frac) -> Frac {
        Frac::new(
            self.num.checked_mul(other.num).expect("no overflow"),
            self.den.checked_mul(other.den).expect("no overflow"),
        )
    }

    pub fn neg(&self) -> Frac {
        Frac {
            num: -self.num,
            den: self.den,
        }
    }

    pub fn div(&self, other: &Frac) -> Frac {
        assert!(!other.is_zero(), "no division by zero");
        Frac::new(
            self.num.checked_mul(other.den).expect("no overflow"),
            self.den.checked_mul(other.num).expect("no overflow"),
        )
    }

    pub fn integer(&self) -> i64 {
        assert!(self.den == 1, "the entry is not an integer");
        self.num as i64
    }

    pub fn parts(&self) -> (i128, i128) {
        (self.num, self.den)
    }

    pub fn below(&self, other: &Frac) -> bool {
        self.num * other.den < other.num * self.den
    }
}

pub struct Table {
    values: HashMap<Vec<u8>, [i64; 4]>,
}

impl Table {
    pub fn new() -> Table {
        Table {
            values: HashMap::new(),
        }
    }

    pub fn get(&mut self, word: &[u8], which: usize) -> i64 {
        if let Some(found) = self.values.get(word) {
            return found[which];
        }
        let obs = observe(word);
        let row = [
            series_value(&obs, 0),
            series_value(&obs, 1),
            series_value(&obs, 2),
            series_value(&obs, 3),
        ];
        self.values.insert(word.to_vec(), row);
        row[which]
    }
}

pub struct Rep {
    pub basis: Vec<Vec<u8>>,
    pub matrices: HashMap<u8, Vec<Vec<i64>>>,
    pub lambda: Vec<i64>,
    pub gamma: Vec<i64>,
}

fn suffixes() -> Vec<Vec<u8>> {
    let mut out = vec![Vec::new()];
    for &a in CODES.iter() {
        out.push(vec![a]);
    }
    for &a in CODES.iter() {
        for &b in CODES.iter() {
            out.push(vec![a, b]);
        }
    }
    out
}

struct Space {
    rows: Vec<Vec<Frac>>,
    echelon: Vec<(usize, Vec<Frac>, Vec<Frac>)>,
}

impl Space {
    fn new() -> Space {
        Space {
            rows: Vec::new(),
            echelon: Vec::new(),
        }
    }

    fn reduce(&self, target: &[Frac]) -> (Vec<Frac>, Vec<Frac>) {
        let mut rest: Vec<Frac> = target.to_vec();
        let mut coefficients = vec![Frac::zero(); self.rows.len()];
        for (pivot, row, combination) in self.echelon.iter() {
            if rest[*pivot].is_zero() {
                continue;
            }
            let factor = rest[*pivot].div(&row[*pivot]);
            for (slot, value) in rest.iter_mut().zip(row.iter()) {
                *slot = slot.add(&factor.mul(value).neg());
            }
            for (slot, value) in coefficients.iter_mut().zip(combination.iter()) {
                *slot = slot.add(&factor.mul(value));
            }
        }
        (coefficients, rest)
    }

    fn insert(&mut self, target: Vec<Frac>) -> bool {
        let (coefficients, rest) = self.reduce(&target);
        let pivot = rest.iter().position(|value| !value.is_zero());
        match pivot {
            None => false,
            Some(at) => {
                let mut combination: Vec<Frac> = coefficients.iter().map(|c| c.neg()).collect();
                combination.push(Frac::int(1));
                for (_, _, old) in self.echelon.iter_mut() {
                    old.push(Frac::zero());
                }
                self.rows.push(target);
                self.echelon.push((at, rest, combination));
                true
            }
        }
    }
}

fn row_of(word: &[u8], which: usize, suffix: &[Vec<u8>], table: &mut Table) -> Vec<Frac> {
    suffix
        .iter()
        .map(|tail| {
            let mut full = word.to_vec();
            full.extend_from_slice(tail);
            Frac::int(table.get(&full, which))
        })
        .collect()
}

pub fn build(which: usize, table: &mut Table) -> Rep {
    let suffix = suffixes();
    let mut space = Space::new();
    let mut basis: Vec<Vec<u8>> = Vec::new();
    let mut queue: Vec<Vec<u8>> = vec![Vec::new()];
    let mut head = 0usize;
    while head < queue.len() {
        let word = queue[head].clone();
        head += 1;
        let row = row_of(&word, which, &suffix, table);
        if !space.insert(row) {
            continue;
        }
        basis.push(word.clone());
        for &code in CODES.iter() {
            let mut next = word.clone();
            next.push(code);
            queue.push(next);
        }
    }
    let size = basis.len();
    let mut matrices: HashMap<u8, Vec<Vec<i64>>> = HashMap::new();
    for &code in CODES.iter() {
        let mut matrix = vec![vec![0i64; size]; size];
        for (index, word) in basis.iter().enumerate() {
            let mut next = word.clone();
            next.push(code);
            let row = row_of(&next, which, &suffix, table);
            let (coefficients, rest) = space.reduce(&row);
            assert!(
                rest.iter().all(|value| value.is_zero()),
                "the basis spans every extension"
            );
            for (slot, value) in matrix[index].iter_mut().zip(coefficients.iter()) {
                *slot = value.integer();
            }
        }
        matrices.insert(code, matrix);
    }
    let empty: Vec<u8> = Vec::new();
    let (coefficients, rest) = space.reduce(&row_of(&empty, which, &suffix, table));
    assert!(
        rest.iter().all(|value| value.is_zero()),
        "the empty word is in the span"
    );
    let lambda: Vec<i64> = coefficients.iter().map(|value| value.integer()).collect();
    let gamma: Vec<i64> = basis.iter().map(|word| table.get(word, which)).collect();
    Rep {
        basis,
        matrices,
        lambda,
        gamma,
    }
}

impl Rep {
    pub fn predict(&self, word: &[u8]) -> i64 {
        let mut state: Vec<i128> = self.lambda.iter().map(|value| *value as i128).collect();
        for code in word {
            let matrix = &self.matrices[code];
            let mut next = vec![0i128; state.len()];
            for (i, weight) in state.iter().enumerate() {
                if *weight == 0 {
                    continue;
                }
                for (j, slot) in next.iter_mut().enumerate() {
                    *slot += weight * matrix[i][j] as i128;
                }
            }
            state = next;
        }
        state
            .iter()
            .zip(self.gamma.iter())
            .map(|(a, b)| a * *b as i128)
            .sum::<i128>() as i64
    }

    pub fn classes(&self) -> Vec<(Vec<u8>, Vec<Vec<i64>>)> {
        let mut out: Vec<(Vec<u8>, Vec<Vec<i64>>)> = Vec::new();
        for &code in CODES.iter() {
            let matrix = self.matrices[&code].clone();
            match out.iter_mut().find(|(_, seen)| *seen == matrix) {
                Some((members, _)) => members.push(code),
                None => out.push((vec![code], matrix)),
            }
        }
        out
    }
}

pub fn product(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let size = a.len();
    let mut out = vec![vec![0i64; size]; size];
    for i in 0..size {
        for k in 0..size {
            if a[i][k] == 0 {
                continue;
            }
            for j in 0..size {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}
