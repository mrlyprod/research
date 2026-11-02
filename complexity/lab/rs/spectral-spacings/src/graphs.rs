use mrly_core::{Rng, Tensor};
use mrly_math::bang::factory::{corners_to_code, residue_corners};
use mrly_math::bang::Code;
use mrly_math::graph::core_graph;
use std::collections::HashMap;

pub struct Graph {
    pub adjacency: Vec<Vec<u32>>,
}

impl Graph {
    pub fn nodes(&self) -> usize {
        self.adjacency.len()
    }

    pub fn edges(&self) -> usize {
        self.adjacency.iter().map(|row| row.len()).sum::<usize>() / 2
    }

    pub fn from_pairs(nodes: usize, pairs: &[(usize, usize)]) -> Graph {
        let mut adjacency = vec![Vec::new(); nodes];
        for (a, b) in pairs {
            adjacency[*a].push(*b as u32);
            adjacency[*b].push(*a as u32);
        }
        Graph { adjacency }
    }

    pub fn of(grid: &Tensor) -> Graph {
        let network = core_graph(grid).expect("a grid has a core graph");
        let pairs: Vec<(usize, usize)> = network
            .branches
            .iter()
            .map(|b| (b.parent, b.child))
            .collect();
        Graph::from_pairs(network.nodes.len(), &pairs)
    }

    pub fn components(&self) -> usize {
        let mut parent: Vec<usize> = (0..self.nodes()).collect();
        for (node, row) in self.adjacency.iter().enumerate() {
            for other in row {
                let (a, b) = (root(&mut parent, node), root(&mut parent, *other as usize));
                if a != b {
                    parent[a] = b;
                }
            }
        }
        (0..self.nodes())
            .filter(|node| root(&mut parent, *node) == *node)
            .count()
    }
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

pub fn carpet(level: usize) -> Graph {
    Graph::of(
        mrly_math::two::create(Code::from(495u64), 3, level, 0, 3)
            .expect("the carpet renders")
            .types(),
    )
}

pub fn sierpinski(level: usize) -> Graph {
    Graph::of(
        mrly_math::two::create(Code::from(7u64), 2, level, 0, 2)
            .expect("the gasket renders")
            .types(),
    )
}

fn sponge_code() -> u128 {
    let filled: Vec<Vec<u8>> = residue_corners(3, 3)
        .into_iter()
        .filter(|corner| corner.iter().filter(|digit| **digit == 1).count() <= 1)
        .collect();
    corners_to_code(&filled, 3, 3).get()
}

pub fn sponge(level: usize) -> Graph {
    Graph::of(
        mrly_math::three::create(Code::from(sponge_code()), 3, level, 3)
            .expect("the sponge renders")
            .types(),
    )
}

pub fn square(side: usize) -> Graph {
    Graph::of(&Tensor::full(vec![side, side], 1))
}

pub fn random(nodes: usize, p: f64, seed: u64) -> Graph {
    let mut rng = Rng::new(seed);
    let mut pairs = Vec::new();
    for a in 0..nodes {
        for b in a + 1..nodes {
            if rng.chance(p) {
                pairs.push((a, b));
            }
        }
    }
    Graph::from_pairs(nodes, &pairs)
}

type Point = [i32; 3];
type Triangle = [Point; 3];

fn sponge_cell(mut cell: [usize; 3], level: usize) -> bool {
    for _ in 0..level {
        let ones = cell.iter().filter(|c| **c % 3 == 1).count();
        if ones > 1 {
            return false;
        }
        cell = [cell[0] / 3, cell[1] / 3, cell[2] / 3];
    }
    true
}

fn cut_points(cell: [usize; 3], d: i32) -> Vec<Point> {
    let base = [2 * cell[0] as i32, 2 * cell[1] as i32, 2 * cell[2] as i32];
    let mut out = Vec::new();
    for axis in 0..3 {
        let others: Vec<usize> = (0..3).filter(|b| *b != axis).collect();
        for b1 in [0, 2] {
            for b2 in [0, 2] {
                if d - b1 - b2 == 1 {
                    let mut p = base;
                    p[axis] += 1;
                    p[others[0]] += b1;
                    p[others[1]] += b2;
                    out.push(p);
                }
            }
        }
    }
    out
}

fn sorted(mut t: Triangle) -> Triangle {
    t.sort();
    t
}

fn dist2(u: &Point, v: &Point) -> i32 {
    (0..3).map(|k| (u[k] - v[k]) * (u[k] - v[k])).sum()
}

fn cell_pieces(cell: [usize; 3], d: i32) -> Vec<Triangle> {
    let pts = cut_points(cell, d);
    if pts.len() == 3 {
        return vec![sorted([pts[0], pts[1], pts[2]])];
    }
    let mid = [
        2 * cell[0] as i32 + 1,
        2 * cell[1] as i32 + 1,
        2 * cell[2] as i32 + 1,
    ];
    let mut out = Vec::new();
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            if dist2(&pts[i], &pts[j]) == 2 {
                out.push(sorted([mid, pts[i], pts[j]]));
            }
        }
    }
    out
}

pub fn slice(level: usize) -> Graph {
    let n = 3usize.pow(level as u32);
    let sigma = (3 * n - 1) / 2;
    let mut triangles: Vec<Triangle> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for s in [sigma - 2, sigma - 1, sigma] {
                if s < i + j || s - i - j >= n {
                    continue;
                }
                let cell = [i, j, s - i - j];
                if !sponge_cell(cell, level) {
                    continue;
                }
                triangles.extend(cell_pieces(cell, 3 * n as i32 - 2 * s as i32));
            }
        }
    }
    let mut owners: HashMap<(Point, Point), Vec<usize>> = HashMap::new();
    for (index, t) in triangles.iter().enumerate() {
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            owners.entry((t[a], t[b])).or_default().push(index);
        }
    }
    let pairs: Vec<(usize, usize)> = owners
        .values()
        .filter(|o| o.len() == 2)
        .map(|o| (o[0], o[1]))
        .collect();
    Graph::from_pairs(triangles.len(), &pairs)
}
