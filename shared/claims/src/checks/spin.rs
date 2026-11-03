use mrly_math::bang::corners_to_code;
use mrly_math::moire::{stack, Combine, Field, Lattice, Spec};
use mrly_math::spin::{profile, reach};

const RINGS: [(usize, usize, f64); 3] = [(3, 5, 0.38), (5, 7, -0.33), (9, 13, 0.38)];

const SIZE: usize = 135;

const STEPS: usize = 400;

pub fn the_coprime_law_survives_the_spin() -> Result<(), String> {
    let far = reach(SIZE);
    let disc: Vec<f64> = (0..STEPS)
        .map(|k| {
            let radius = k as f64 * far / (STEPS - 1) as f64;
            if radius <= SIZE as f64 / 2.0 {
                radius
            } else {
                0.0
            }
        })
        .collect();
    let ones = vec![1.0; SIZE * SIZE];
    for (m, n, ring) in RINGS {
        let (a, b) = (layer(m)?, layer(n)?);
        let flat = pearson(&a.as_f64(), &b.as_f64(), &ones);
        if flat.abs() >= 0.03 {
            return Err(format!("the flat layers {m} and {n} correlate at {flat}"));
        }
        let spun = pearson(&profiled(&a), &profiled(&b), &disc);
        if (spun - ring).abs() >= 0.02 {
            return Err(format!("the ring profiles {m} and {n} correlate at {spun}"));
        }
    }
    Ok(())
}

fn layer(scale: usize) -> Result<Field, String> {
    let spec = Spec::new(corners_to_code(&[vec![0, 0]], 2, 2).get(), 2, 2);
    stack(spec, &[scale], Combine::Sum, 1, Lattice::Square, SIZE, &[])
        .map_err(|_| format!("the layer at {scale} does not build"))
}

fn profiled(field: &Field) -> Vec<f64> {
    profile(&field.data, SIZE, STEPS)
        .unwrap_or_default()
        .iter()
        .map(|&value| value as f64)
        .collect()
}

fn pearson(a: &[f64], b: &[f64], weight: &[f64]) -> f64 {
    let total: f64 = weight.iter().sum();
    let mean = |x: &[f64]| x.iter().zip(weight).map(|(v, w)| v * w).sum::<f64>() / total;
    let (ma, mb) = (mean(a), mean(b));
    let covariance: f64 = a
        .iter()
        .zip(b)
        .zip(weight)
        .map(|((x, y), w)| w * (x - ma) * (y - mb))
        .sum();
    let spread = |x: &[f64], m: f64| {
        x.iter()
            .zip(weight)
            .map(|(v, w)| w * (v - m) * (v - m))
            .sum::<f64>()
    };
    covariance / (spread(a, ma) * spread(b, mb)).sqrt()
}
