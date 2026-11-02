use mrly_core::Tensor;
use mrly_math::bang::baseq::axis_maps;
use mrly_math::bang::factory::{corners_to_code, residue_corners};
use mrly_math::bang::universe::permutations;
use mrly_math::bang::Code;
use mrly_math::graph::core_graph;
use std::collections::BTreeSet;

pub const BASE: usize = 3;

pub fn plane(code: u128, base: usize, level: usize) -> Tensor {
    mrly_math::two::create(Code::from(code), base, level, 0, base)
        .expect("a plane design renders")
        .types()
        .clone()
}

pub fn sponge_code() -> u128 {
    let filled: Vec<Vec<u8>> = residue_corners(3, BASE)
        .into_iter()
        .filter(|corner| corner.iter().filter(|digit| **digit == 1).count() <= 1)
        .collect();
    corners_to_code(&filled, 3, BASE).get()
}

pub fn sponge(level: usize) -> Tensor {
    mrly_math::three::create(Code::from(sponge_code()), BASE, level, BASE)
        .expect("the sponge renders")
        .types()
        .clone()
}

pub fn filled(grid: &Tensor) -> Vec<usize> {
    (0..grid.size())
        .filter(|&flat| grid.at(flat) != 0)
        .collect()
}

pub fn strides(shape: &[usize]) -> Vec<usize> {
    let mut out = vec![1; shape.len()];
    for axis in (0..shape.len().saturating_sub(1)).rev() {
        out[axis] = out[axis + 1] * shape[axis + 1];
    }
    out
}

pub fn coords(flat: usize, shape: &[usize]) -> Vec<usize> {
    let mut left = flat;
    strides(shape)
        .iter()
        .map(|stride| {
            let index = left / stride;
            left %= stride;
            index
        })
        .collect()
}

pub struct Graph {
    pub shape: Vec<usize>,
    pub cells: Vec<usize>,
    pub adjacency: Vec<Vec<u32>>,
}

impl Graph {
    pub fn of(grid: &Tensor) -> Graph {
        let network = core_graph(grid).expect("a grid has a core graph");
        let mut adjacency: Vec<Vec<u32>> = vec![Vec::new(); network.nodes.len()];
        for branch in &network.branches {
            adjacency[branch.parent].push(branch.child as u32);
            adjacency[branch.child].push(branch.parent as u32);
        }
        for row in adjacency.iter_mut() {
            row.sort_unstable();
        }
        Graph {
            shape: grid.shape.clone(),
            cells: filled(grid),
            adjacency,
        }
    }

    pub fn nodes(&self) -> usize {
        self.adjacency.len()
    }

    pub fn labels(&self) -> Vec<usize> {
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
            .map(|node| root(&mut parent, node))
            .collect()
    }
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

pub struct Components {
    pub count: usize,
    pub giant: Tensor,
    pub share: f64,
    pub spanning: bool,
}

pub fn components(grid: &Tensor) -> Components {
    let graph = Graph::of(grid);
    let labels = graph.labels();
    let mut sizes = vec![0usize; graph.nodes()];
    for label in &labels {
        sizes[*label] += 1;
    }
    let count = sizes.iter().filter(|size| **size > 0).count();
    let best = (0..graph.nodes())
        .max_by_key(|node| (sizes[*node], usize::MAX - node))
        .unwrap_or(0);
    let mut giant = Tensor::new(grid.shape.clone());
    let mut low = vec![usize::MAX; grid.shape.len()];
    let mut high = vec![0usize; grid.shape.len()];
    for (node, flat) in graph.cells.iter().enumerate() {
        if labels[node] == best {
            giant.put(*flat, 1);
            for (axis, at) in coords(*flat, &grid.shape).iter().enumerate() {
                low[axis] = low[axis].min(*at);
                high[axis] = high[axis].max(*at);
            }
        }
    }
    let walls =
        (0..grid.shape.len()).all(|axis| low[axis] == 0 && high[axis] + 1 == grid.shape[axis]);
    let held = if graph.nodes() == 0 { 0 } else { sizes[best] };
    Components {
        count,
        giant,
        share: if graph.nodes() == 0 {
            0.0
        } else {
            held as f64 / graph.nodes() as f64
        },
        spanning: walls,
    }
}

pub fn group() -> Vec<Vec<usize>> {
    let cells = residue_corners(2, BASE);
    let maps = axis_maps(BASE);
    let mut out = Vec::new();
    for order in permutations(2) {
        for first in &maps {
            for second in &maps {
                out.push(
                    cells
                        .iter()
                        .map(|cell| {
                            first[cell[order[0]] as usize] * BASE + second[cell[order[1]] as usize]
                        })
                        .collect(),
                );
            }
        }
    }
    out
}

pub fn carry(element: &[usize], code: u128) -> u128 {
    element
        .iter()
        .enumerate()
        .filter(|(index, _)| code >> index & 1 == 1)
        .map(|(_, image)| 1u128 << image)
        .sum()
}

pub fn orbit(group: &[Vec<usize>], code: u128) -> BTreeSet<u128> {
    group.iter().map(|element| carry(element, code)).collect()
}

pub fn classes(group: &[Vec<usize>]) -> Vec<(u128, usize)> {
    let mut seen = vec![false; 512];
    let mut out = Vec::new();
    for code in 0..512u128 {
        if seen[code as usize] {
            continue;
        }
        let members = orbit(group, code);
        out.push((code, members.len()));
        for member in members {
            seen[member as usize] = true;
        }
    }
    out
}
