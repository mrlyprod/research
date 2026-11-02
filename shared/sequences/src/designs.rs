use crate::row::{Row, SECONDS, TERMS};
use mrly_core::error::Result;
use mrly_math::bang::{baseq, code_to_corners, factory, Code};
use mrly_math::counts::{self, Exposure};
use mrly_math::{six, three, two};
use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

pub const SPACES: [(usize, usize); 9] = [
    (1, 2),
    (2, 2),
    (3, 2),
    (4, 2),
    (1, 3),
    (2, 3),
    (1, 4),
    (2, 4),
    (1, 5),
];
pub const CELLS: u128 = 600_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Axis {
    Level,
    Side,
}

impl Axis {
    pub fn slug(self) -> &'static str {
        match self {
            Axis::Level => "level",
            Axis::Side => "side",
        }
    }

    fn place(self, index: usize, number: usize) -> (usize, u32) {
        match self {
            Axis::Level => (number, index as u32 + 1),
            Axis::Side => (2 * index + 3, 1),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Group {
    Closed,
    Cut,
    Profile,
    Render,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Measure {
    Fills,
    Voids,
    Surface,
    Peak,
    Heights,
    Vertices,
    Edges,
    Faces,
    Euler,
    Perimeter,
    Hidden,
    Triangles,
    Gaps,
    ProFills,
    ProVoids,
    Holes,
    Pieces,
    SliceEuler,
}

impl Measure {
    pub const ALL: [Measure; 18] = [
        Measure::Fills,
        Measure::Voids,
        Measure::Surface,
        Measure::Peak,
        Measure::Heights,
        Measure::Vertices,
        Measure::Edges,
        Measure::Faces,
        Measure::Euler,
        Measure::Perimeter,
        Measure::Hidden,
        Measure::Triangles,
        Measure::Gaps,
        Measure::ProFills,
        Measure::ProVoids,
        Measure::Holes,
        Measure::Pieces,
        Measure::SliceEuler,
    ];

    pub fn slug(self) -> &'static str {
        match self {
            Measure::Fills => "fills",
            Measure::Voids => "voids",
            Measure::Surface => "surface",
            Measure::Peak => "peak",
            Measure::Heights => "heights",
            Measure::Vertices => "vertices",
            Measure::Edges => "edges",
            Measure::Faces => "faces",
            Measure::Euler => "euler",
            Measure::Perimeter => "perimeter",
            Measure::Hidden => "hidden",
            Measure::Triangles => "triangles",
            Measure::Gaps => "gaps",
            Measure::ProFills => "pro-fills",
            Measure::ProVoids => "pro-voids",
            Measure::Holes => "holes",
            Measure::Pieces => "pieces",
            Measure::SliceEuler => "slice-euler",
        }
    }

    fn group(self) -> Group {
        match self {
            Measure::Fills | Measure::Voids | Measure::Surface => Group::Closed,
            Measure::Triangles | Measure::Gaps | Measure::ProFills | Measure::ProVoids => {
                Group::Cut
            }
            Measure::Peak | Measure::Heights => Group::Profile,
            _ => Group::Render,
        }
    }

    pub fn applies(self, dimension: usize, base: usize) -> bool {
        match self {
            Measure::Fills | Measure::Voids | Measure::Surface => true,
            Measure::Peak | Measure::Heights => dimension >= 2,
            Measure::Vertices | Measure::Edges | Measure::Euler => dimension == 2 || dimension == 3,
            Measure::Perimeter => dimension == 2,
            Measure::Faces | Measure::Hidden => dimension == 3,
            _ => dimension == 3 && base == 2,
        }
    }

    pub fn generator(self, dimension: usize) -> &'static str {
        let plane = dimension == 2;
        match self {
            Measure::Fills | Measure::Voids => "mrly_math::counts::fill",
            Measure::Surface => "mrly_math::counts::Exposure::at",
            Measure::Peak | Measure::Heights => "mrly_math::counts::profile_of_tile",
            Measure::Vertices | Measure::Edges | Measure::Euler if plane => {
                "mrly_math::two::census"
            }
            Measure::Perimeter => "mrly_math::two::census",
            Measure::Vertices | Measure::Edges | Measure::Euler | Measure::Faces => {
                "mrly_math::three::census"
            }
            Measure::Hidden => "mrly_math::three::hidden",
            Measure::Triangles => "mrly_math::counts::cut_fills",
            Measure::Gaps => "mrly_math::counts::cut_voids",
            Measure::ProFills => "mrly_math::counts::pro_fills",
            Measure::ProVoids => "mrly_math::counts::pro_voids",
            Measure::Holes => "mrly_math::six::holes",
            Measure::Pieces => "mrly_math::six::components",
            Measure::SliceEuler => "mrly_math::six::euler",
        }
    }
}

pub fn name(code: u128, dimension: usize, base: usize, measure: Measure, axis: Axis) -> String {
    let base = if base == 2 {
        String::new()
    } else {
        format!("_base={base}")
    };
    format!(
        "sequence_dim={dimension}{base}_code={code}_measure={}_axis={}",
        measure.slug(),
        axis.slug()
    )
}

fn within(cost: Option<u128>) -> bool {
    cost.is_some_and(|cost| cost <= CELLS)
}

fn values(
    group: Group,
    code: u128,
    dimension: usize,
    base: usize,
    number: usize,
    level: u32,
    exposure: &Exposure,
) -> Option<Vec<(Measure, Option<i128>)>> {
    let id = Code::from(code);
    let grid = (number as u128).checked_pow(dimension as u32 * level);
    let out = match group {
        Group::Closed => {
            let fill = counts::fill(id, number, dimension, 1, base)
                .ok()?
                .checked_pow(level);
            vec![
                (Measure::Fills, fill.and_then(|f| i128::try_from(f).ok())),
                (
                    Measure::Voids,
                    grid.zip(fill).and_then(|(g, f)| i128::try_from(g - f).ok()),
                ),
                (
                    Measure::Surface,
                    exposure.at(level).and_then(|s| i128::try_from(s).ok()),
                ),
            ]
        }
        Group::Cut => {
            let cost = (number as u128)
                .checked_pow(2 * level)
                .and_then(|side| side.checked_mul(16));
            if !within(cost) {
                return None;
            }
            vec![
                (
                    Measure::Triangles,
                    Some(counts::cut_fills(id, number, level).ok()? as i128),
                ),
                (
                    Measure::Gaps,
                    Some(counts::cut_voids(id, number, level).ok()? as i128),
                ),
                (
                    Measure::ProFills,
                    Some(counts::pro_fills(id, number, level).ok()? as i128),
                ),
                (
                    Measure::ProVoids,
                    Some(counts::pro_voids(id, number, level).ok()? as i128),
                ),
            ]
        }
        Group::Profile => {
            let span = (number as u128)
                .checked_pow(level)
                .map(|side| dimension as u128 * (side - 1) + 1);
            if !within(span) {
                return None;
            }
            let tile = factory::create(id, number, dimension, base, 1).ok()?;
            let profile = counts::profile_of_tile(&tile, level).ok()?;
            vec![
                (
                    Measure::Peak,
                    Some(profile.iter().max().copied().unwrap_or(0) as i128),
                ),
                (
                    Measure::Heights,
                    Some(profile.iter().filter(|&&c| c > 0).count() as i128),
                ),
            ]
        }
        Group::Render => {
            if !within(grid) {
                return None;
            }
            if dimension == 2 {
                let tally =
                    two::census(&two::create(id, number, level as usize, 0, base).ok()?).ok()?;
                vec![
                    (Measure::Vertices, Some(tally.vertices as i128)),
                    (Measure::Edges, Some(tally.edges as i128)),
                    (Measure::Euler, Some(tally.euler as i128)),
                    (Measure::Perimeter, Some(tally.perimeter as i128)),
                ]
            } else {
                let cell = three::create(id, number, level as usize, base).ok()?;
                let tally = three::census(&cell).ok()?;
                let mut out = vec![
                    (Measure::Vertices, Some(tally.vertices as i128)),
                    (Measure::Edges, Some(tally.edges as i128)),
                    (Measure::Faces, Some(tally.faces as i128)),
                    (Measure::Euler, Some(tally.euler as i128)),
                    (Measure::Hidden, Some(three::hidden(&cell) as i128)),
                ];
                if base == 2 {
                    let slice = six::cut(&cell).ok()?;
                    out.push((Measure::Holes, Some(six::holes(&slice).ok()? as i128)));
                    out.push((Measure::Pieces, Some(six::components(&slice).ok()? as i128)));
                    out.push((Measure::SliceEuler, Some(six::euler(&slice, false) as i128)));
                }
                out
            }
        }
    };
    Some(out)
}

pub fn sequences(
    code: u128,
    dimension: usize,
    base: usize,
    axis: Axis,
    full: bool,
) -> Result<Vec<(Measure, Vec<i128>)>> {
    let corners = code_to_corners(Code::from(code), dimension, base)?;
    let measures: Vec<Measure> = Measure::ALL
        .into_iter()
        .filter(|m| m.applies(dimension, base) && (full || m.group() == Group::Closed))
        .collect();
    let mut terms: HashMap<Measure, Vec<i128>> = measures.iter().map(|&m| (m, vec![])).collect();
    let mut dead: HashSet<Measure> = HashSet::new();
    let mut spent: HashMap<Group, Duration> = HashMap::new();
    let groups: Vec<Group> = [Group::Closed, Group::Cut, Group::Profile, Group::Render]
        .into_iter()
        .filter(|g| measures.iter().any(|m| m.group() == *g))
        .collect();
    let level_number = base.max(3);
    let level_exposure = Exposure::from_corners(&corners, level_number, dimension, base);
    for index in 0..TERMS {
        let (number, level) = axis.place(index, level_number);
        let exposure = match axis {
            Axis::Level => level_exposure.clone(),
            Axis::Side => Exposure::from_corners(&corners, number, dimension, base),
        };
        for &group in &groups {
            let members: Vec<Measure> = measures
                .iter()
                .copied()
                .filter(|m| m.group() == group && !dead.contains(m))
                .collect();
            if members.is_empty() {
                continue;
            }
            let start = Instant::now();
            let found = values(group, code, dimension, base, number, level, &exposure);
            let total = spent.entry(group).or_default();
            *total += start.elapsed();
            let over = *total > SECONDS;
            let found: HashMap<Measure, Option<i128>> =
                found.unwrap_or_default().into_iter().collect();
            for measure in members {
                match found.get(&measure).copied().flatten() {
                    Some(value) => terms.get_mut(&measure).unwrap().push(value),
                    None => {
                        dead.insert(measure);
                    }
                }
                if over {
                    dead.insert(measure);
                }
            }
        }
        if measures.iter().all(|m| dead.contains(m)) {
            break;
        }
    }
    Ok(measures
        .into_iter()
        .map(|m| (m, terms.remove(&m).unwrap()))
        .collect())
}

fn binomial(n: usize, k: usize) -> i128 {
    (0..k).fold(1i128, |acc, i| acc * (n - i) as i128 / (i as i128 + 1))
}

pub fn fill_polynomial(corners: &[Vec<u8>], dimension: usize) -> Vec<i128> {
    let mut out = vec![0i128; dimension + 1];
    for corner in corners {
        let ones = corner.iter().filter(|&&r| r != 0).count();
        for j in 0..=ones {
            let sign = if (ones - j) % 2 == 0 { 1 } else { -1 };
            out[dimension - ones + j] += sign * binomial(ones, j);
        }
    }
    out
}

fn odd_power(dimension: usize) -> Vec<i128> {
    (0..=dimension)
        .map(|j| {
            let sign = if (dimension - j) % 2 == 0 { 1 } else { -1 };
            sign * binomial(dimension, j) * (1i128 << j)
        })
        .collect()
}

pub fn polynomial_text(rising: &[i128]) -> String {
    let mut out = String::new();
    for (power, &c) in rising.iter().enumerate().rev() {
        if c == 0 {
            continue;
        }
        let magnitude = c.abs();
        let body = match (power, magnitude) {
            (0, m) => m.to_string(),
            (1, 1) => "k".to_string(),
            (1, m) => format!("{m}k"),
            (p, 1) => format!("k^{p}"),
            (p, m) => format!("{m}k^{p}"),
        };
        if out.is_empty() {
            out = if c < 0 { format!("-{body}") } else { body };
        } else {
            out += &format!(" {} {body}", if c < 0 { "-" } else { "+" });
        }
    }
    if out.is_empty() {
        "0".to_string()
    } else {
        out
    }
}

pub fn recurrence_text(coefficients: &[i128]) -> String {
    let mut out = String::from("a(level) =");
    for (i, &c) in coefficients.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = format!("a(level-{})", i + 1);
        let magnitude = c.abs();
        let body = if magnitude == 1 {
            term
        } else {
            format!("{magnitude} {term}")
        };
        if out.ends_with('=') {
            out += &format!(" {}{body}", if c < 0 { "-" } else { "" });
        } else {
            out += &format!(" {} {body}", if c < 0 { "-" } else { "+" });
        }
    }
    out
}

pub fn closed(
    code: u128,
    dimension: usize,
    base: usize,
    measure: Measure,
    axis: Axis,
) -> Result<Option<String>> {
    let corners = code_to_corners(Code::from(code), dimension, base)?;
    let number = base.max(3);
    let fill = counts::fill(Code::from(code), number, dimension, 1, base)?;
    let grid = (number as u128).pow(dimension as u32);
    Ok(match (measure, axis) {
        (Measure::Fills, Axis::Level) => Some(format!("{fill}^level")),
        (Measure::Voids, Axis::Level) => Some(format!("{grid}^level - {fill}^level")),
        (Measure::Fills, Axis::Side) if base == 2 => {
            Some(polynomial_text(&fill_polynomial(&corners, dimension)))
        }
        (Measure::Voids, Axis::Side) if base == 2 => {
            let fills = fill_polynomial(&corners, dimension);
            let voids: Vec<i128> = odd_power(dimension)
                .iter()
                .zip(&fills)
                .map(|(g, f)| g - f)
                .collect();
            Some(polynomial_text(&voids))
        }
        (Measure::Surface, Axis::Level) => Some(recurrence_text(
            &Exposure::from_corners(&corners, number, dimension, base).recurrence(),
        )),
        _ => None,
    })
}

pub fn orbit_rows(
    dimension: usize,
    base: usize,
    rep: u128,
    group: &[Vec<usize>],
) -> Result<Vec<Row>> {
    let mates: Vec<u128> = baseq::orbit(group, Code::from(rep))
        .into_iter()
        .map(|c| c.get())
        .filter(|&c| c != rep)
        .collect();
    let full_mates = base == 2 && dimension <= 3;
    let mut rows = Vec::new();
    for axis in [Axis::Level, Axis::Side] {
        let mut seen: HashMap<Measure, HashSet<Vec<i128>>> = HashMap::new();
        for (code, full) in
            std::iter::once((rep, true)).chain(mates.iter().map(|&m| (m, full_mates)))
        {
            for (measure, terms) in sequences(code, dimension, base, axis, full)? {
                if !seen.entry(measure).or_default().insert(terms.clone()) {
                    continue;
                }
                rows.push(
                    Row::new(
                        name(code, dimension, base, measure, axis),
                        measure.generator(dimension),
                        &format!("{}/{}", measure.slug(), axis.slug()),
                        terms,
                    )
                    .closed(closed(code, dimension, base, measure, axis)?),
                );
            }
        }
    }
    Ok(rows)
}

pub fn jobs() -> Result<Vec<(usize, usize, u128)>> {
    let mut out = Vec::new();
    for (dimension, base) in SPACES {
        for (code, _) in baseq::representatives(base, dimension)? {
            out.push((dimension, base, code.get()));
        }
    }
    Ok(out)
}
