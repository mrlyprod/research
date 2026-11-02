// DESIGN

#[derive(Clone, Copy, Debug)]
pub struct Design {
    pub name: &'static str,
    pub dimension: usize,
    pub invert: bool,
}

pub const MENGER: Design = Design {
    name: "menger",
    dimension: 3,
    invert: false,
};

pub const CARPET: Design = Design {
    name: "carpet",
    dimension: 2,
    invert: false,
};

pub const VICSEK: Design = Design {
    name: "vicsek",
    dimension: 2,
    invert: true,
};

pub const DESIGNS: [Design; 3] = [MENGER, CARPET, VICSEK];

impl Design {
    pub fn hit(&self, digit: u64) -> bool {
        (digit == 1) != self.invert
    }

    pub fn corners(&self) -> Vec<Vec<u64>> {
        let mut out = Vec::new();
        let total = 3u64.pow(self.dimension as u32);
        for code in 0..total {
            let mut vector = Vec::with_capacity(self.dimension);
            let mut rest = code;
            for _ in 0..self.dimension {
                vector.push(rest % 3);
                rest /= 3;
            }
            vector.reverse();
            if vector.iter().filter(|d| self.hit(**d)).count() <= 1 {
                out.push(vector);
            }
        }
        out
    }

    pub fn fill(&self) -> u64 {
        self.corners().len() as u64
    }

    pub fn zero_filled(&self) -> bool {
        !self.invert || self.dimension < 2
    }

    pub fn named(name: &str) -> Option<Design> {
        DESIGNS.iter().find(|d| d.name == name).copied()
    }
}
