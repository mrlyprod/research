use crate::graphs::Graph;
use faer::{Mat, Side};

pub fn eigenvalues(graph: &Graph, normalised: bool) -> Vec<f64> {
    let n = graph.nodes();
    let mut matrix = Mat::<f64>::zeros(n, n);
    for (node, row) in graph.adjacency.iter().enumerate() {
        let degree = row.len() as f64;
        if normalised {
            matrix[(node, node)] = if degree > 0.0 { 1.0 } else { 0.0 };
        } else {
            matrix[(node, node)] = degree;
        }
        for other in row {
            let weight = if normalised {
                -1.0 / (degree * graph.adjacency[*other as usize].len() as f64).sqrt()
            } else {
                -1.0
            };
            matrix[(node, *other as usize)] = weight;
        }
    }
    let mut values = matrix
        .as_ref()
        .self_adjoint_eigenvalues(Side::Lower)
        .expect("the dense eigensolver converges");
    values.sort_by(|a, b| a.partial_cmp(b).expect("finite eigenvalues"));
    values
}
