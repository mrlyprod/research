use mrly_core::tensor::Tensor;
use mrly_math::bang::factory::create;
use mrly_math::bang::Code;
use mrly_math::shape::{census, crop, named, Frac, Shape};

// DESIGNS

fn design(dimension: usize, level: usize) -> Tensor {
    if level == 0 {
        return Tensor::full(vec![1; dimension], 1);
    }
    match dimension {
        2 => create(Code::from(7u128), 3, 2, 2, level).unwrap(),
        _ => create(Code::from(23u128), 3, 3, 2, level).unwrap(),
    }
}

// EXPOSURE

fn exposed(types: &Tensor) -> usize {
    let dims = types.shape.clone();
    let rank = dims.len();
    let mut index = vec![0usize; rank];
    let mut probe = vec![0usize; rank];
    let mut count = 0;
    for flat in 0..types.size() {
        let mut rem = flat;
        for axis in (0..rank).rev() {
            index[axis] = rem % dims[axis];
            rem /= dims[axis];
        }
        if types.at(flat) == 0 {
            continue;
        }
        for axis in 0..rank {
            for step in [-1i64, 1] {
                let next = index[axis] as i64 + step;
                if next < 0 || next as usize >= dims[axis] {
                    count += 1;
                    continue;
                }
                probe.copy_from_slice(&index);
                probe[axis] = next as usize;
                if types.at(types.index(&probe)) == 0 {
                    count += 1;
                }
            }
        }
    }
    count
}

// LINES

fn line(label: &str, name: &str, dimension: usize, radius: Frac, types: &Tensor) {
    let shape = named(name, dimension, radius).unwrap();
    let tally = census(&shape, types).unwrap();
    let kept = crop(types, &shape, true).unwrap();
    let inner = crop(types, &shape, false).unwrap();
    let anti_kept = crop(types, &Shape::Anti(Box::new(shape.clone())), true).unwrap();
    let anti_inner = crop(types, &Shape::Anti(Box::new(shape)), false).unwrap();
    let filled = types.sum();
    assert_eq!(tally.filled.iter().sum::<usize>() as u64, filled);
    assert_eq!(kept.sum() as usize, tally.filled[1] + tally.filled[2]);
    assert_eq!(inner.sum() as usize, tally.filled[2]);
    assert_eq!(anti_kept.sum() as usize, tally.filled[0] + tally.filled[1]);
    assert_eq!(anti_inner.sum() as usize, tally.filled[0]);
    assert_eq!(inner.sum() + anti_kept.sum(), filled);
    assert_eq!(kept.sum() + anti_inner.sum(), filled);
    println!(
        "crop-counts {label} {name} r={}/{} filled_in={} filled_cut={} exposed_after={}",
        radius.num,
        radius.den,
        tally.filled[2],
        tally.filled[1],
        exposed(&kept)
    );
}

fn main() {
    println!("crop-counts generator: cargo run --release --manifest-path research/crop/lab/rs/crop-counts/Cargo.toml");
    let half = Frac::new(1, 2).unwrap();
    for level in 0..=5usize {
        let carpet = design(2, level);
        for name in ["ball", "diamond"] {
            line(&format!("levels carpet L={level}"), name, 2, half, &carpet);
        }
    }
    for level in 0..=4usize {
        let sponge = design(3, level);
        for name in ["ball", "diamond"] {
            line(&format!("levels sponge L={level}"), name, 3, half, &sponge);
        }
    }
    let carpet = design(2, 4);
    let sponge = design(3, 3);
    for rnum in 1..=24i64 {
        let radius = Frac::new(rnum, 24).unwrap();
        for name in ["ball", "diamond"] {
            line("sweep carpet L=4", name, 2, radius, &carpet);
            line("sweep sponge L=3", name, 3, radius, &sponge);
        }
    }
}
