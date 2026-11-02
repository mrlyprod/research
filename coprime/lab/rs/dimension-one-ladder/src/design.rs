use mrly_num::series::APERY;
use std::f64::consts::PI;

pub struct Design {
    pub name: &'static str,
    pub dimension: usize,
    pub invert: bool,
    pub fill: u64,
}

pub const CARPET: Design = Design {
    name: "carpet",
    dimension: 2,
    invert: false,
    fill: 8,
};
pub const MENGER: Design = Design {
    name: "menger",
    dimension: 3,
    invert: false,
    fill: 20,
};
pub const VICSEK: Design = Design {
    name: "vicsek",
    dimension: 2,
    invert: true,
    fill: 5,
};

impl Design {
    pub fn named(name: &str) -> Option<&'static Design> {
        [&CARPET, &MENGER, &VICSEK]
            .into_iter()
            .find(|d| d.name == name)
    }

    pub fn origin_filled(&self) -> bool {
        !self.invert
    }

    pub fn density(&self) -> f64 {
        match self.name {
            "carpet" => 189.0 / (32.0 * PI * PI),
            "menger" => (513.0 / 520.0) / APERY,
            _ => 27.0 / (4.0 * PI * PI),
        }
    }

    pub fn filled(&self, digits: &[u64]) -> bool {
        digits.iter().filter(|&&d| (d == 1) != self.invert).count() <= 1
    }

    pub fn corners(&self) -> Vec<Vec<u64>> {
        let mut out = Vec::new();
        for code in 0..3u64.pow(self.dimension as u32) {
            let mut rest = code;
            let mut digits = Vec::with_capacity(self.dimension);
            for _ in 0..self.dimension {
                digits.push(rest % 3);
                rest /= 3;
            }
            if self.filled(&digits) {
                out.push(digits);
            }
        }
        out
    }
}
