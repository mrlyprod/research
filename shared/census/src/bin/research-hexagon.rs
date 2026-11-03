use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_tile::recipe::Design;

const NAME: &str = "research-hexagon";
const MESH: usize = 110;
const DEEPEST: i64 = 55;
const GRAIN: usize = 8;
const NONE: u16 = 65535;
const DIGITS: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

struct Rule {
    corners: [bool; 8],
}

impl Rule {
    fn new(design: Design) -> Result<Rule> {
        let cube = mrly_math::three::named(design, 2, 1)?;
        let mut corners = [false; 8];
        for (slot, corner) in corners.iter_mut().enumerate() {
            *corner = cube.types().get(&[slot >> 2, (slot >> 1) & 1, slot & 1])? == 1;
        }
        Ok(Rule { corners })
    }

    fn filled(&self, x: i64, y: i64, z: i64) -> bool {
        let bit = |c: i64| (c.div_euclid(4) & 1) as usize;
        self.corners[bit(x) << 2 | bit(y) << 1 | bit(z)]
    }

    fn at(&self, n: i64, x: i64, z: i64) -> Option<bool> {
        let y = 6 * n - 2 - x - z;
        (0..4 * n).contains(&y).then(|| self.filled(x, y, z))
    }

    fn sample(&self, n: i64, x: f64, z: f64) -> Option<bool> {
        if !(0.0..1.0).contains(&x) || !(0.0..1.0).contains(&z) {
            return None;
        }
        let column = ((x * 4.0 * n as f64).floor() as i64).min(4 * n - 1);
        let row = (2 * (z * 2.0 * n as f64).floor() as i64).min(4 * n - 2);
        self.at(n, column, row)
    }
}

fn layers() -> impl Iterator<Item = i64> {
    (3..=DEEPEST).step_by(2)
}

fn stack(rule: &Rule, x: f64, z: f64) -> Option<f64> {
    let mut filled = 0usize;
    let mut seen = 0usize;
    for n in layers() {
        if let Some(cell) = rule.sample(n, x, z) {
            seen += 1;
            filled += usize::from(cell);
        }
    }
    (seen > 0).then(|| filled as f64 / seen as f64)
}

fn unit(across: f64, down: f64) -> (f64, f64) {
    let span = 2.0 * MESH as f64;
    let z = 1.0 - down / span;
    (0.25 + across / span - 0.5 * z, z)
}

fn row_len(row: usize) -> usize {
    let reach = if row < MESH { row } else { 2 * MESH - 1 - row };
    2 * (MESH + reach) + 1
}

fn corners(row: usize, col: usize) -> [(f64, f64); 3] {
    let up = row < MESH;
    let reach = if up { row } else { 2 * MESH - 1 - row };
    let (long, short) = (MESH + reach + 1, MESH + reach);
    let (top_len, bot_len) = if up { (short, long) } else { (long, short) };
    let step = (col / 2) as f64;
    let (down, top) = (row as f64, row as f64 + 1.0);
    if col.is_multiple_of(2) == up {
        let edge = (2 * MESH - bot_len) as f64 / 2.0 + step;
        [
            unit(edge, top),
            unit(edge + 1.0, top),
            unit(edge + 0.5, down),
        ]
    } else {
        let edge = (2 * MESH - top_len) as f64 / 2.0 + step;
        [
            unit(edge, down),
            unit(edge + 1.0, down),
            unit(edge + 0.5, top),
        ]
    }
}

fn shade(rule: &Rule, tri: [(f64, f64); 3]) -> Option<f64> {
    let mut total = 0.0;
    let mut taken = 0usize;
    for i in 0..GRAIN {
        for j in 0..GRAIN - i {
            let a = (i as f64 + 1.0 / 3.0) / GRAIN as f64;
            let b = (j as f64 + 1.0 / 3.0) / GRAIN as f64;
            let c = 1.0 - a - b;
            let x = a * tri[0].0 + b * tri[1].0 + c * tri[2].0;
            let z = a * tri[0].1 + b * tri[1].1 + c * tri[2].1;
            if let Some(value) = stack(rule, x, z) {
                total += value;
                taken += 1;
            }
        }
    }
    (taken > 0).then(|| total / taken as f64)
}

fn ink_law(rule: &Rule, n: i64) -> bool {
    let chi = if (3 * n - 1) / 2 % 2 == 0 { 1 } else { -1 };
    let (mut filled, mut inside) = (0i64, 0i64);
    for x in 0..4 * n {
        for step in 0..2 * n {
            if let Some(cell) = rule.at(n, x, 2 * step) {
                inside += 1;
                filled += i64::from(cell);
            }
        }
    }
    filled * 8 * n * n == inside * (4 * n * n + chi * n * n + 4 * n - chi)
}

fn base64(bytes: &[u8]) -> String {
    let mut text = String::new();
    for chunk in bytes.chunks(3) {
        let word = chunk.iter().enumerate().fold(0u32, |word, (i, &byte)| {
            word | u32::from(byte) << (16 - 8 * i)
        });
        for k in 0..=chunk.len() {
            text.push(char::from(DIGITS[(word >> (18 - 6 * k) & 63) as usize]));
        }
    }
    text
}

fn main() -> Result<()> {
    let carpet = Rule::new(Design::Carpet)?;
    assert_eq!(layers().count(), 27);
    assert!(layers().all(|n| ink_law(&carpet, n)));

    let mut field: Vec<Vec<Option<f64>>> = Vec::with_capacity(2 * MESH);
    for row in 0..2 * MESH {
        let mut line = Vec::with_capacity(row_len(row));
        for col in 0..row_len(row) {
            line.push(shade(&carpet, corners(row, col)));
        }
        field.push(line);
    }
    assert_eq!(field.iter().map(Vec::len).sum::<usize>(), 6 * MESH * MESH);
    assert!(layers().all(|n| carpet.sample(n, 0.5, 0.5).is_some()));

    let drawn: Vec<f64> = field.iter().flatten().flatten().copied().collect();
    let background = drawn.iter().sum::<f64>() / drawn.len() as f64;
    let mut spread: Vec<f64> = drawn
        .iter()
        .map(|value| (value - background).abs())
        .collect();
    spread.sort_by(|a, b| a.total_cmp(b));
    let reach = spread[spread.len() * 9 / 10];

    let mut shades = drawn.clone();
    shades.sort_by(|a, b| a.total_cmp(b));
    shades.dedup_by(|a, b| a.to_bits() == b.to_bits());
    assert!(shades.len() < usize::from(NONE));
    let mut codes = Vec::with_capacity(12 * MESH * MESH);
    for cell in field.iter().flatten() {
        let code = cell.map_or(NONE, |value| {
            shades
                .binary_search_by(|shade| shade.total_cmp(&value))
                .expect("a shade of the table") as u16
        });
        codes.extend(code.to_le_bytes());
    }
    let table: Vec<u8> = shades
        .iter()
        .flat_map(|shade| shade.to_le_bytes())
        .collect();
    save(
        NAME,
        &json!({
            "mesh": MESH,
            "background": background,
            "reach": reach,
            "shades": base64(&table),
            "field": base64(&codes),
        }),
    )
}
