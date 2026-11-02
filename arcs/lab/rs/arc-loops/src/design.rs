use mrly_math::bang::factory::create;
use mrly_math::bang::Code;

#[derive(Clone)]
pub struct Design {
    pub number: usize,
    pub rule_base: usize,
    pub rule_code: u128,
    pub tile: Vec<bool>,
}

impl Design {
    pub fn new(rule_code: u128, number: usize, rule_base: usize) -> Design {
        let tensor = create(Code::from(rule_code), number, 2, rule_base, 1).expect("code in range");
        let tile = (0..number * number).map(|f| tensor.at(f) == 1).collect();
        Design {
            number,
            rule_base,
            rule_code,
            tile,
        }
    }

    pub fn full(code: u128, base: usize) -> Design {
        Design::new(code, base, base)
    }

    pub fn full_code(&self) -> u128 {
        self.tile
            .iter()
            .enumerate()
            .filter(|(_, &on)| on)
            .map(|(i, _)| 1u128 << i)
            .sum()
    }

    pub fn parity_code(&self) -> Option<u128> {
        let n = self.number;
        let mut code = 0u128;
        let mut seen = [None; 4];
        for i in 0..n {
            for j in 0..n {
                let bit = (i % 2) * 2 + (j % 2);
                let on = self.tile[i * n + j];
                match seen[bit] {
                    None => seen[bit] = Some(on),
                    Some(was) if was != on => return None,
                    _ => {}
                }
                if on {
                    code |= 1 << bit;
                }
            }
        }
        Some(code)
    }

    pub fn kept(&self) -> usize {
        self.tile.iter().filter(|&&on| on).count()
    }

    pub fn name(&self) -> String {
        let n = self.number;
        let full = if n == 2 {
            format!("bang dim 2, code {}", self.full_code())
        } else {
            format!("bang dim 2, base {n}, code {}", self.full_code())
        };
        match self.parity_code() {
            Some(p) if n > 2 => format!("{full} = bang dim 2, code {p} at side number {n}"),
            _ => full,
        }
    }

    pub fn cells(&self, level: usize) -> (usize, Vec<bool>) {
        let side = self.number.pow(level as u32);
        if level == 0 {
            return (1, vec![true]);
        }
        let tensor = create(
            Code::from(self.rule_code),
            self.number,
            2,
            self.rule_base,
            level,
        )
        .expect("code in range");
        (side, (0..side * side).map(|f| tensor.at(f) == 1).collect())
    }

    pub fn image(&self, g: usize) -> Design {
        let n = self.number;
        let mut tile = vec![false; n * n];
        for i in 0..n {
            for j in 0..n {
                let (a, b) = match g {
                    0 => (i, j),
                    1 => (n - 1 - i, n - 1 - j),
                    2 => (j, i),
                    _ => (n - 1 - j, n - 1 - i),
                };
                tile[a * n + b] = self.tile[i * n + j];
            }
        }
        let code: u128 = tile
            .iter()
            .enumerate()
            .filter(|(_, &on)| on)
            .map(|(i, _)| 1u128 << i)
            .sum();
        Design::full(code, n)
    }
}
