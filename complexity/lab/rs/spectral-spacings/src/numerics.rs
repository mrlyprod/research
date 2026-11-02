pub fn sum(values: &[f64]) -> f64 {
    if values.len() <= 8 {
        return values.iter().fold(0.0, |acc, v| acc + v);
    }
    let (left, right) = values.split_at(values.len() / 2);
    sum(left) + sum(right)
}

pub fn mean(values: &[f64]) -> f64 {
    sum(values) / values.len() as f64
}

pub fn fraction_below(values: &[f64], bound: f64) -> f64 {
    values.iter().filter(|v| **v < bound).count() as f64 / values.len() as f64
}

pub fn goe(s: f64) -> f64 {
    1.0 - (-std::f64::consts::PI * s * s / 4.0).exp()
}

pub fn poisson(s: f64) -> f64 {
    1.0 - (-s).exp()
}

pub fn ks_distance(values: &[f64], law: fn(f64) -> f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite spacings"));
    let m = sorted.len() as f64;
    let mut worst: f64 = 0.0;
    for (i, s) in sorted.iter().enumerate() {
        let f = law(*s);
        worst = worst.max((i + 1) as f64 / m - f).max(f - i as f64 / m);
    }
    worst
}

pub fn ks_pvalue(distance: f64, count: usize) -> f64 {
    let z = distance * (count as f64).sqrt();
    if z < 1e-6 {
        return 1.0;
    }
    let mut total = 0.0;
    for k in 1..200 {
        let sign = if k % 2 == 1 { 2.0 } else { -2.0 };
        total += sign * (-2.0 * (k * k) as f64 * z * z).exp();
    }
    total.clamp(0.0, 1.0)
}
