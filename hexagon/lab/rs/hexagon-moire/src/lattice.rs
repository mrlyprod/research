use mrly_math::three::Cell3d;
use mrly_tile::recipe::Design;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Carpet,
    Net,
    Tree,
    Void,
}

pub const FAMILIES: [Family; 4] = [Family::Carpet, Family::Net, Family::Tree, Family::Void];

impl Family {
    pub fn name(self) -> &'static str {
        match self {
            Family::Carpet => "carpet",
            Family::Net => "net",
            Family::Tree => "tree",
            Family::Void => "void",
        }
    }

    fn design(self) -> Design {
        match self {
            Family::Carpet => Design::Carpet,
            Family::Net => Design::Net,
            Family::Tree => Design::Ztree,
            Family::Void => Design::Void,
        }
    }

    pub fn cube(self, number: usize) -> Cell3d {
        mrly_math::three::named(self.design(), number, 1).expect("a named cube")
    }
}

pub struct Rule {
    corners: [bool; 8],
}

impl Rule {
    pub fn new(family: Family) -> Rule {
        let cube = family.cube(2);
        let mut corners = [false; 8];
        for (slot, corner) in corners.iter_mut().enumerate() {
            *corner = cube
                .types()
                .get(&[slot >> 2, (slot >> 1) & 1, slot & 1])
                .is_ok_and(|v| v == 1);
        }
        Rule { corners }
    }

    pub fn filled(&self, x: i64, y: i64, z: i64) -> bool {
        let bit = |c: i64| (c.div_euclid(4) & 1) as usize;
        self.corners[bit(x) << 2 | bit(y) << 1 | bit(z)]
    }
}

pub fn cell(rule: &Rule, n: i64, x: i64, z: i64) -> Option<bool> {
    let y = 6 * n - 2 - x - z;
    (0..4 * n).contains(&y).then(|| rule.filled(x, y, z))
}

pub fn column(point: f64, n: i64) -> i64 {
    ((point * 4.0 * n as f64).floor() as i64).min(4 * n - 1)
}

pub fn row(point: f64, n: i64) -> i64 {
    (2 * (point * 2.0 * n as f64).floor() as i64).min(4 * n - 2)
}
