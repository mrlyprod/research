pub fn pairwise(values: &[f64]) -> f64 {
    if values.len() <= 128 {
        return values.iter().sum();
    }
    let half = values.len() / 2;
    pairwise(&values[..half]) + pairwise(&values[half..])
}

pub fn mean(values: &[f64]) -> f64 {
    pairwise(values) / values.len() as f64
}

pub fn odds(limit: usize) -> impl Iterator<Item = usize> {
    (1..=limit).step_by(2)
}
