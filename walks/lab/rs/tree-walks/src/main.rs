use mrly_core::Tensor;
use mrly_math::bang::Code;
use mrly_math::graph::core_graph;

const SUBJECTS: [(u128, usize, usize); 2] = [(79, 3, 6), (7, 2, 9)];

fn plane(code: u128, base: usize, level: usize) -> Tensor {
    mrly_math::two::create(Code::from(code), base, level, 0, base)
        .expect("a plane design renders")
        .types()
        .clone()
}

fn fold(hash: u64, word: u64) -> u64 {
    (hash ^ word).wrapping_mul(0x100000001b3)
}

fn main() {
    for (code, base, top) in SUBJECTS {
        for level in 1..=top {
            let grid = plane(code, base, level);
            let cells: Vec<usize> = (0..grid.size()).filter(|&f| grid.at(f) != 0).collect();
            let network = core_graph(&grid).expect("a grid has a core graph");
            let mut edges: Vec<(usize, usize)> = network
                .branches
                .iter()
                .map(|b| {
                    let (u, v) = (cells[b.parent], cells[b.child]);
                    (u.min(v), u.max(v))
                })
                .collect();
            edges.sort_unstable();
            let cell_hash = cells
                .iter()
                .fold(0xcbf29ce484222325, |h, &c| fold(h, c as u64));
            let edge_hash = edges.iter().fold(0xcbf29ce484222325, |h, &(u, v)| {
                fold(fold(h, u as u64), v as u64)
            });
            println!(
                "base {base} code {code} level {level} nodes {} edges {} cells {cell_hash:016x} edge-list {edge_hash:016x}",
                cells.len(),
                edges.len()
            );
        }
    }
}
