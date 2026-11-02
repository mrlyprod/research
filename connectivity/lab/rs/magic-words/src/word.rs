pub const CODES: [u8; 15] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
pub const LIBRARY: [u8; 10] = [3, 5, 6, 7, 9, 10, 11, 12, 13, 14];

pub fn corners(code: u8) -> Vec<(usize, usize)> {
    (0..4usize)
        .filter(|i| (code >> i) & 1 == 1)
        .map(|i| (i / 2, i % 2))
        .collect()
}

pub struct Grid {
    pub side: usize,
    pub cells: Vec<bool>,
}

pub fn render(word: &[u8]) -> Grid {
    let mut side = 1usize;
    let mut cells = vec![true];
    for &code in word {
        let next = side * 2;
        let mut out = vec![false; next * next];
        let filled = corners(code);
        for r in 0..side {
            for c in 0..side {
                if cells[r * side + c] {
                    for &(a, b) in &filled {
                        out[(2 * r + a) * next + 2 * c + b] = true;
                    }
                }
            }
        }
        side = next;
        cells = out;
    }
    Grid { side, cells }
}

impl Grid {
    pub fn at(&self, r: usize, c: usize) -> bool {
        self.cells[r * self.side + c]
    }

    pub fn fill(&self) -> u64 {
        self.cells.iter().filter(|c| **c).count() as u64
    }

    pub fn diagonal(&self) -> u64 {
        (0..self.side).filter(|&i| self.at(i, i)).count() as u64
    }

    pub fn adjacent(&self) -> u64 {
        let mut count = 0u64;
        for r in 0..self.side {
            for c in 0..self.side {
                if !self.at(r, c) {
                    continue;
                }
                if r + 1 < self.side && self.at(r + 1, c) {
                    count += 1;
                }
                if c + 1 < self.side && self.at(r, c + 1) {
                    count += 1;
                }
            }
        }
        count
    }

    pub fn quads(&self) -> u64 {
        let mut count = 0u64;
        for r in 0..self.side.saturating_sub(1) {
            for c in 0..self.side.saturating_sub(1) {
                if self.at(r, c) && self.at(r + 1, c) && self.at(r, c + 1) && self.at(r + 1, c + 1)
                {
                    count += 1;
                }
            }
        }
        count
    }

    pub fn interior(&self) -> u64 {
        let mut count = 0u64;
        for r in 0..self.side {
            for c in 0..self.side {
                if !self.at(r, c) {
                    continue;
                }
                let inside = r > 0
                    && c > 0
                    && r + 1 < self.side
                    && c + 1 < self.side
                    && self.at(r - 1, c)
                    && self.at(r + 1, c)
                    && self.at(r, c - 1)
                    && self.at(r, c + 1);
                if inside {
                    count += 1;
                }
            }
        }
        count
    }

    pub fn boundary(&self) -> u64 {
        self.fill() - self.interior()
    }

    pub fn perimeter(&self) -> u64 {
        4 * self.fill() - 2 * self.adjacent()
    }

    pub fn euler(&self) -> i64 {
        self.fill() as i64 - self.adjacent() as i64 + self.quads() as i64
    }

    pub fn labels(&self) -> (Vec<i64>, u64) {
        let mut label = vec![-1i64; self.cells.len()];
        let mut count = 0i64;
        let mut stack: Vec<usize> = Vec::new();
        for start in 0..self.cells.len() {
            if !self.cells[start] || label[start] >= 0 {
                continue;
            }
            label[start] = count;
            stack.push(start);
            while let Some(at) = stack.pop() {
                let r = at / self.side;
                let c = at % self.side;
                let push = |rr: usize, cc: usize, label: &mut Vec<i64>, stack: &mut Vec<usize>| {
                    let next = rr * self.side + cc;
                    if self.cells[next] && label[next] < 0 {
                        label[next] = count;
                        stack.push(next);
                    }
                };
                if r > 0 {
                    push(r - 1, c, &mut label, &mut stack);
                }
                if r + 1 < self.side {
                    push(r + 1, c, &mut label, &mut stack);
                }
                if c > 0 {
                    push(r, c - 1, &mut label, &mut stack);
                }
                if c + 1 < self.side {
                    push(r, c + 1, &mut label, &mut stack);
                }
            }
            count += 1;
        }
        (label, count as u64)
    }

    pub fn components(&self) -> u64 {
        self.labels().1
    }

    pub fn holes(&self) -> u64 {
        let mut seen = vec![false; self.cells.len()];
        let mut stack: Vec<usize> = Vec::new();
        let open = |at: usize, seen: &mut Vec<bool>, stack: &mut Vec<usize>| {
            if !self.cells[at] && !seen[at] {
                seen[at] = true;
                stack.push(at);
            }
        };
        for i in 0..self.side {
            open(i, &mut seen, &mut stack);
            open((self.side - 1) * self.side + i, &mut seen, &mut stack);
            open(i * self.side, &mut seen, &mut stack);
            open(i * self.side + self.side - 1, &mut seen, &mut stack);
        }
        while let Some(at) = stack.pop() {
            let r = at / self.side;
            let c = at % self.side;
            let step = |rr: usize, cc: usize, seen: &mut Vec<bool>, stack: &mut Vec<usize>| {
                let next = rr * self.side + cc;
                if !self.cells[next] && !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            };
            if r > 0 {
                step(r - 1, c, &mut seen, &mut stack);
            }
            if r + 1 < self.side {
                step(r + 1, c, &mut seen, &mut stack);
            }
            if c > 0 {
                step(r, c - 1, &mut seen, &mut stack);
            }
            if c + 1 < self.side {
                step(r, c + 1, &mut seen, &mut stack);
            }
        }
        let mut count = 0u64;
        for start in 0..self.cells.len() {
            if self.cells[start] || seen[start] {
                continue;
            }
            count += 1;
            seen[start] = true;
            stack.push(start);
            while let Some(at) = stack.pop() {
                let r = at / self.side;
                let c = at % self.side;
                let step = |rr: usize, cc: usize, seen: &mut Vec<bool>, stack: &mut Vec<usize>| {
                    let next = rr * self.side + cc;
                    if !self.cells[next] && !seen[next] {
                        seen[next] = true;
                        stack.push(next);
                    }
                };
                if r > 0 {
                    step(r - 1, c, &mut seen, &mut stack);
                }
                if r + 1 < self.side {
                    step(r + 1, c, &mut seen, &mut stack);
                }
                if c > 0 {
                    step(r, c - 1, &mut seen, &mut stack);
                }
                if c + 1 < self.side {
                    step(r, c + 1, &mut seen, &mut stack);
                }
            }
        }
        count
    }

    pub fn profile(&self) -> Vec<u64> {
        let mut out = vec![0u64; 2 * self.side - 1];
        for r in 0..self.side {
            for c in 0..self.side {
                if self.at(r, c) {
                    out[r + c] += 1;
                }
            }
        }
        out
    }

    pub fn contacts(&self) -> (u64, u64) {
        let rows = (0..self.side)
            .filter(|&r| self.at(r, 0) && self.at(r, self.side - 1))
            .count() as u64;
        let cols = (0..self.side)
            .filter(|&c| self.at(0, c) && self.at(self.side - 1, c))
            .count() as u64;
        (rows, cols)
    }

    pub fn merging(&self) -> u64 {
        let (label, _) = self.labels();
        let mut touched: Vec<i64> = Vec::new();
        let mark = |at: usize, touched: &mut Vec<i64>| {
            if label[at] >= 0 && !touched.contains(&label[at]) {
                touched.push(label[at]);
            }
        };
        for r in 0..self.side {
            if self.at(r, 0) && self.at(r, self.side - 1) {
                mark(r * self.side, &mut touched);
                mark(r * self.side + self.side - 1, &mut touched);
            }
        }
        for c in 0..self.side {
            if self.at(0, c) && self.at(self.side - 1, c) {
                mark(c, &mut touched);
                mark((self.side - 1) * self.side + c, &mut touched);
            }
        }
        touched.len() as u64
    }
}

#[derive(Clone, PartialEq, Eq)]
pub struct Obs {
    pub fill: u64,
    pub diagonal: u64,
    pub boundary: u64,
    pub perimeter: u64,
    pub components: u64,
    pub euler: i64,
    pub holes: u64,
    pub profile: Vec<u64>,
}

pub fn observe(word: &[u8]) -> Obs {
    let grid = render(word);
    Obs {
        fill: grid.fill(),
        diagonal: grid.diagonal(),
        boundary: grid.boundary(),
        perimeter: grid.perimeter(),
        components: grid.components(),
        euler: grid.euler(),
        holes: grid.holes(),
        profile: grid.profile(),
    }
}

pub const SERIES: [&str; 4] = ["components", "Euler characteristic", "boundary", "holes"];

pub fn series_value(obs: &Obs, which: usize) -> i64 {
    match which {
        0 => obs.components as i64,
        1 => obs.euler,
        2 => obs.boundary as i64,
        _ => obs.holes as i64,
    }
}

pub fn spell(word: &[u8]) -> String {
    if word.is_empty() {
        return "e".to_string();
    }
    if word.len() == 1 {
        return format!("{}", word[0]);
    }
    let inner: Vec<String> = word.iter().map(|c| format!("{c}")).collect();
    format!("({})", inner.join(","))
}
