use crate::design::Graph;
use faer::{Mat, Side};
use mrly_core::Rng;

pub const MODES: usize = 4;
pub const DENSE_LIMIT: usize = 2000;
const BLOCK: usize = 6;
const ROUNDS: usize = 14;
const SOLVE_TOL: f64 = 1e-9;
const SOLVE_CAP: usize = 200000;
const RITZ_TOL: f64 = 1e-11;

pub fn apply(graph: &Graph, x: &[f64], out: &mut [f64]) {
    for (node, row) in graph.adjacency.iter().enumerate() {
        let mut total = row.len() as f64 * x[node];
        for other in row {
            total -= x[*other as usize];
        }
        out[node] = total;
    }
}

pub fn dense(graph: &Graph) -> Mat<f64> {
    let n = graph.nodes();
    let mut out = Mat::<f64>::zeros(n, n);
    for (node, row) in graph.adjacency.iter().enumerate() {
        out[(node, node)] = row.len() as f64;
        for other in row {
            out[(node, *other as usize)] = -1.0;
        }
    }
    out
}

pub fn low(graph: &Graph) -> Vec<f64> {
    let n = graph.nodes();
    let want = MODES.min(n.saturating_sub(1));
    if want == 0 {
        return vec![0.0];
    }
    if n <= DENSE_LIMIT {
        let mut values = dense(graph)
            .as_ref()
            .self_adjoint_eigenvalues(Side::Lower)
            .expect("the dense eigensolver converges");
        values.truncate(want + 1);
        return values;
    }
    let mut values = vec![0.0];
    values.extend(krylov(graph, want));
    values
}

pub fn slow_mode(graph: &Graph) -> Vec<f64> {
    let eigen = dense(graph)
        .as_ref()
        .self_adjoint_eigen(Side::Lower)
        .expect("the dense eigensolver converges");
    let u = eigen.U();
    (0..graph.nodes()).map(|node| u[(node, 1)]).collect()
}

fn krylov(graph: &Graph, want: usize) -> Vec<f64> {
    let n = graph.nodes();
    let mut rng = Rng::new(12345);
    let mut block: Vec<Vec<f64>> = (0..BLOCK)
        .map(|_| (0..n).map(|_| rng.unit() * 2.0 - 1.0).collect())
        .collect();
    orthonormalise(&mut block, &[]);
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut gram: Vec<Vec<f64>> = Vec::new();
    let mut work = vec![0.0; n];
    let mut latest: Vec<f64> = Vec::new();
    for _ in 0..ROUNDS {
        let mut grown: Vec<Vec<f64>> = block
            .iter()
            .map(|column| solve(graph, column, &mut work))
            .collect();
        orthonormalise(&mut grown, &basis);
        if grown.is_empty() {
            break;
        }
        for column in &grown {
            apply(graph, column, &mut work);
            let mut row: Vec<f64> = basis.iter().map(|seat| dot(seat, &work)).collect();
            row.push(dot(column, &work));
            for (seat, entry) in row.iter().enumerate().take(gram.len()) {
                gram[seat].push(*entry);
            }
            gram.push(row);
            basis.push(column.clone());
        }
        block = grown;
        if basis.len() < want + 2 {
            continue;
        }
        let candidate = ritz(&gram, want);
        let settled = latest.len() == candidate.len()
            && latest
                .iter()
                .zip(&candidate)
                .all(|(old, new)| (old - new).abs() <= RITZ_TOL * new.abs());
        latest = candidate;
        if settled {
            break;
        }
    }
    latest
}

fn ritz(gram: &[Vec<f64>], want: usize) -> Vec<f64> {
    let size = gram.len();
    let small = Mat::<f64>::from_fn(size, size, |row, column| {
        0.5 * (gram[row][column] + gram[column][row])
    });
    let mut values = small
        .as_ref()
        .self_adjoint_eigenvalues(Side::Lower)
        .expect("the Ritz eigensolver converges");
    values.truncate(want);
    values
}

fn solve(graph: &Graph, rhs: &[f64], work: &mut [f64]) -> Vec<f64> {
    let n = graph.nodes();
    let mut residual = rhs.to_vec();
    centre(&mut residual);
    let target = dot(&residual, &residual).sqrt();
    let mut solution = vec![0.0; n];
    if target == 0.0 {
        return solution;
    }
    let mut direction = residual.clone();
    let mut square = dot(&residual, &residual);
    for _ in 0..SOLVE_CAP {
        apply(graph, &direction, work);
        let curvature = dot(&direction, work);
        if curvature <= 0.0 {
            break;
        }
        let stride = square / curvature;
        for seat in 0..n {
            solution[seat] += stride * direction[seat];
            residual[seat] -= stride * work[seat];
        }
        centre(&mut residual);
        let next = dot(&residual, &residual);
        if next.sqrt() <= SOLVE_TOL * target {
            break;
        }
        let blend = next / square;
        for seat in 0..n {
            direction[seat] = residual[seat] + blend * direction[seat];
        }
        square = next;
    }
    centre(&mut solution);
    solution
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn centre(column: &mut [f64]) {
    let mean = column.iter().sum::<f64>() / column.len() as f64;
    for entry in column.iter_mut() {
        *entry -= mean;
    }
}

fn orthonormalise(block: &mut Vec<Vec<f64>>, against: &[Vec<f64>]) {
    let mut kept: Vec<Vec<f64>> = Vec::new();
    for column in block.iter_mut() {
        for _ in 0..2 {
            centre(column);
            for seat in against.iter().chain(kept.iter()) {
                let overlap = dot(seat, column);
                for (entry, base) in column.iter_mut().zip(seat) {
                    *entry -= overlap * base;
                }
            }
        }
        let norm = dot(column, column).sqrt();
        if norm > 1e-8 {
            for entry in column.iter_mut() {
                *entry /= norm;
            }
            kept.push(column.clone());
        }
    }
    *block = kept;
}

pub fn exponents(coarse: &[f64], fine: &[f64], scale: f64) -> Vec<f64> {
    let modes = MODES.min(coarse.len() - 1).min(fine.len() - 1);
    (1..=modes)
        .map(|mode| (coarse[mode] / fine[mode]).ln() / scale.ln())
        .collect()
}
