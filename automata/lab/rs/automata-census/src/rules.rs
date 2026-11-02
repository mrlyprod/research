pub const RULES: usize = 256;

pub fn output(rule: usize, l: u8, c: u8, r: u8) -> u8 {
    ((rule >> (4 * l as usize + 2 * c as usize + r as usize)) & 1) as u8
}

pub struct Diagram {
    pub steps: usize,
    pub width: usize,
    pub cells: Vec<u8>,
}

impl Diagram {
    pub fn at(&self, t: usize, i: usize) -> u8 {
        self.cells[t * self.width + i]
    }
    pub fn centre(&self) -> usize {
        self.width / 2
    }
    pub fn signed(&self, t: usize, offset: i64) -> u8 {
        let i = self.centre() as i64 + offset;
        if i < 0 || i >= self.width as i64 {
            0
        } else {
            self.at(t, i as usize)
        }
    }
}

pub fn evolve(rule: usize, steps: usize, pad: usize) -> (Diagram, u8) {
    let window = 2 * steps + 1;
    let width = window + 2 * pad;
    let mut cells = vec![0u8; (steps + 1) * width];
    cells[width / 2] = 1;
    let mut seen = 0u8;
    for t in 0..steps {
        for i in 0..width {
            let l = if i == 0 { 0 } else { cells[t * width + i - 1] };
            let c = cells[t * width + i];
            let r = if i + 1 == width {
                0
            } else {
                cells[t * width + i + 1]
            };
            seen |= 1 << (4 * l + 2 * c + r);
            cells[(t + 1) * width + i] = output(rule, l, c, r);
        }
    }
    for i in 0..width {
        let l = if i == 0 {
            0
        } else {
            cells[steps * width + i - 1]
        };
        let c = cells[steps * width + i];
        let r = if i + 1 == width {
            0
        } else {
            cells[steps * width + i + 1]
        };
        seen |= 1 << (4 * l + 2 * c + r);
    }
    let mut cropped = vec![0u8; (steps + 1) * window];
    for t in 0..=steps {
        cropped[t * window..(t + 1) * window]
            .copy_from_slice(&cells[t * width + pad..t * width + pad + window]);
    }
    (
        Diagram {
            steps,
            width: window,
            cells: cropped,
        },
        seen,
    )
}

pub fn single_seed(rule: usize, steps: usize) -> Diagram {
    evolve(rule, steps, steps).0
}

pub fn mirrored(diagram: &Diagram) -> Vec<u8> {
    let mut out = vec![0u8; diagram.cells.len()];
    for t in 0..=diagram.steps {
        for i in 0..diagram.width {
            out[t * diagram.width + i] = diagram.at(t, diagram.width - 1 - i);
        }
    }
    out
}
