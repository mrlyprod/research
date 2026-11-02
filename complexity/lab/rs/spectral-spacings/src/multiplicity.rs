pub struct Classes {
    pub distinct: usize,
    pub repeated: usize,
    pub largest: usize,
}

pub fn classes(values: &[f64], tolerance: f64) -> Classes {
    let mut sizes = Vec::new();
    let mut run = 1;
    for pair in values.windows(2) {
        if pair[1] - pair[0] > tolerance {
            sizes.push(run);
            run = 1;
        } else {
            run += 1;
        }
    }
    sizes.push(run);
    Classes {
        distinct: sizes.len(),
        repeated: sizes.iter().filter(|s| **s > 1).sum(),
        largest: *sizes.iter().max().unwrap_or(&0),
    }
}
