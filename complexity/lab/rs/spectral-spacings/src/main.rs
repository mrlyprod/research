mod graphs;
mod multiplicity;
mod numerics;
mod spectrum;
mod unfolding;

use graphs::Graph;
use numerics::{fraction_below, goe, ks_distance, ks_pvalue, poisson};

const DEGREE: usize = 12;
const HALF: usize = 10;
const TOLERANCE: f64 = 1e-9;

struct Row {
    name: String,
    graph: Graph,
}

fn subjects() -> Vec<Row> {
    let mut out = Vec::new();
    let mut push = |name: String, graph: Graph| out.push(Row { name, graph });
    push(
        "random graph n=400 p=0.1".into(),
        graphs::random(400, 0.1, 7),
    );
    push("square lattice 20x20".into(), graphs::square(20));
    push("carpet L=3".into(), graphs::carpet(3));
    push("carpet L=4".into(), graphs::carpet(4));
    push("sponge L=2".into(), graphs::sponge(2));
    push("slice L=2".into(), graphs::slice(2));
    push("slice L=3".into(), graphs::slice(3));
    push("sierpinski L=5".into(), graphs::sierpinski(5));
    push("sierpinski L=6".into(), graphs::sierpinski(6));
    out
}

fn report(label: &str, spacings: &unfolding::Spacings) {
    let m = spacings.values.len();
    let d_goe = ks_distance(&spacings.values, goe);
    let d_poi = ks_distance(&spacings.values, poisson);
    println!(
        "    {label}: P(s<0.5) = {:.4}  KS goe = {:.4} (p = {:.3})  KS poisson = {:.4} (p = {:.3})  negative steps = {}",
        fraction_below(&spacings.values, 0.5),
        d_goe,
        ks_pvalue(d_goe, m),
        d_poi,
        ks_pvalue(d_poi, m),
        spacings.negatives
    );
}

fn main() {
    println!("GOE P(s<0.5) = 1 - exp(-pi/16) = {:.5}", goe(0.5));
    println!("POISSON P(s<0.5) = 1 - exp(-1/2) = {:.5}", poisson(0.5));
    println!("polynomial unfolder degree {DEGREE}, window unfolder half width {HALF}, clustering tolerance {TOLERANCE:e}");
    for row in subjects() {
        let graph = &row.graph;
        println!(
            "{}: {} nodes, {} edges, {} components",
            row.name,
            graph.nodes(),
            graph.edges(),
            graph.components()
        );
        for (kind, normalised) in [("combinatorial", false), ("normalised", true)] {
            let values = spectrum::eigenvalues(graph, normalised);
            let classes = multiplicity::classes(&values, TOLERANCE);
            println!(
                "  {kind} laplacian: distinct {} ({:.2}%), zero spacings {} ({:.2}%), in repeated classes {} ({:.2}%), largest multiplicity {}",
                classes.distinct,
                100.0 * classes.distinct as f64 / values.len() as f64,
                values.len() - classes.distinct,
                100.0 * (values.len() - classes.distinct) as f64 / values.len() as f64,
                classes.repeated,
                100.0 * classes.repeated as f64 / values.len() as f64,
                classes.largest
            );
            report("polynomial", &unfolding::polynomial(&values, DEGREE));
            report("window", &unfolding::window(&values, HALF));
        }
    }
    let big = graphs::slice(4);
    println!(
        "slice L=4: {} nodes, {} edges, {} components, spectrum not computed",
        big.nodes(),
        big.edges(),
        big.components()
    );
}
