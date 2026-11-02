use crate::design::Graph;
use std::collections::{BTreeSet, HashMap};

pub type Corner = (i64, i64);

pub fn edges(level: u32) -> BTreeSet<(Corner, Corner)> {
    let mut out = BTreeSet::new();
    let side = 1i64 << level;
    split(&mut out, (0, 0), (side, 0), (0, side), level);
    out
}

fn split(out: &mut BTreeSet<(Corner, Corner)>, a: Corner, b: Corner, c: Corner, level: u32) {
    if level == 0 {
        for (one, two) in [(a, b), (b, c), (a, c)] {
            out.insert(if one < two { (one, two) } else { (two, one) });
        }
        return;
    }
    let ab = ((a.0 + b.0) / 2, (a.1 + b.1) / 2);
    let bc = ((b.0 + c.0) / 2, (b.1 + c.1) / 2);
    let ac = ((a.0 + c.0) / 2, (a.1 + c.1) / 2);
    split(out, a, ab, ac, level - 1);
    split(out, ab, b, bc, level - 1);
    split(out, ac, bc, c, level - 1);
}

pub struct Gasket {
    pub points: Vec<Corner>,
    pub graph: Graph,
}

pub fn build(level: u32) -> Gasket {
    let edges = edges(level);
    let points: Vec<Corner> = edges
        .iter()
        .flat_map(|edge| [edge.0, edge.1])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let seat: HashMap<Corner, u32> = points
        .iter()
        .enumerate()
        .map(|(index, point)| (*point, index as u32))
        .collect();
    let mut adjacency: Vec<Vec<u32>> = vec![Vec::new(); points.len()];
    for edge in &edges {
        let (here, there) = (seat[&edge.0], seat[&edge.1]);
        adjacency[here as usize].push(there);
        adjacency[there as usize].push(here);
    }
    for row in adjacency.iter_mut() {
        row.sort_unstable();
    }
    Gasket {
        points,
        graph: Graph {
            shape: Vec::new(),
            cells: Vec::new(),
            adjacency,
        },
    }
}
