mod derive;
mod transfer;

use mrly_core::Tensor;
use mrly_math::bang::factory::create;
use mrly_math::bang::Code;
use mrly_math::shape::{census, named, Frac, Shape};

fn ceil_sqrt(value: u64) -> u64 {
    if value == 0 {
        return 0;
    }
    let mut root = (value as f64).sqrt() as u64;
    while root.saturating_mul(root) < value {
        root += 1;
    }
    while root > 0 && (root - 1) * (root - 1) >= value {
        root -= 1;
    }
    root
}

fn bucket(quad: u64) -> u64 {
    ceil_sqrt(quad).div_ceil(2)
}

fn design(dimension: usize, level: usize) -> Tensor {
    match dimension {
        2 => create(Code::from(7u128), 3, 2, 2, level).unwrap(),
        _ => create(Code::from(23u128), 3, 3, 2, level).unwrap(),
    }
}

struct Sweep {
    r_max: u64,
    seen: Vec<u64>,
    inside: Vec<u64>,
    touch: Vec<u64>,
    all_seen: Vec<u64>,
    all_inside: Vec<u64>,
    all_touch: Vec<u64>,
}

fn sweep(types: &Tensor, shift: i64, r_max: u64) -> Sweep {
    let dims = types.shape.clone();
    let rank = dims.len();
    let width = (r_max + 2) as usize;
    let mut out = Sweep {
        r_max,
        seen: vec![0; width],
        inside: vec![0; width],
        touch: vec![0; width],
        all_seen: vec![0; width],
        all_inside: vec![0; width],
        all_touch: vec![0; width],
    };
    let mut index = vec![0i64; rank];
    for flat in 0..types.size() {
        let mut centre = 0u64;
        let mut far = 0u64;
        let mut near = 0u64;
        for coordinate in &index {
            let offset = (2 * coordinate + 1 - shift).unsigned_abs();
            centre += offset * offset;
            far += (offset + 1) * (offset + 1);
            let low = offset.saturating_sub(1);
            near += low * low;
        }
        let slots = [bucket(centre), bucket(far), bucket(near)];
        let filled = types.at(flat) != 0;
        for (which, slot) in slots.iter().enumerate() {
            if *slot > r_max {
                continue;
            }
            let at = *slot as usize;
            match which {
                0 => {
                    out.all_seen[at] += 1;
                    if filled {
                        out.seen[at] += 1;
                    }
                }
                1 => {
                    out.all_inside[at] += 1;
                    if filled {
                        out.inside[at] += 1;
                    }
                }
                _ => {
                    out.all_touch[at] += 1;
                    if filled {
                        out.touch[at] += 1;
                    }
                }
            }
        }
        for axis in (0..rank).rev() {
            index[axis] += 1;
            if (index[axis] as usize) < dims[axis] {
                break;
            }
            index[axis] = 0;
        }
    }
    for column in [
        &mut out.seen,
        &mut out.inside,
        &mut out.touch,
        &mut out.all_seen,
        &mut out.all_inside,
        &mut out.all_touch,
    ] {
        for r in 1..column.len() {
            column[r] += column[r - 1];
        }
    }
    out
}

impl Sweep {
    fn cut(&self, r: u64) -> u64 {
        self.touch[r as usize] - self.inside[r as usize]
    }
    fn all_cut(&self, r: u64) -> u64 {
        self.all_touch[r as usize] - self.all_inside[r as usize]
    }
}

fn ball(dimension: usize, shift: i64, radius: Frac) -> Shape {
    if shift == 0 {
        Shape::Ball {
            center: vec![Frac::whole(0); dimension],
            radius,
        }
    } else {
        named("ball", dimension, radius).unwrap()
    }
}

fn oracle(label: &str, types: &Tensor, shift: i64, side: usize, table: &Sweep, radii: &[u64]) {
    for r in radii {
        let shape = ball(
            types.shape.len(),
            shift,
            Frac::new(*r as i64, side as i64).unwrap(),
        );
        let tally = census(&shape, types).unwrap();
        assert_eq!(tally.filled[2] as u64, table.inside[*r as usize]);
        assert_eq!(tally.filled[1] as u64, table.cut(*r));
        assert_eq!(tally.cells[2] as u64, table.all_inside[*r as usize]);
        assert_eq!(tally.cells[1] as u64, table.all_cut(*r));
        println!(
            "circle-crop oracle {label} r={r} census_in={} census_cut={} sweep_in={} sweep_cut={}",
            tally.filled[2],
            tally.filled[1],
            table.inside[*r as usize],
            table.cut(*r)
        );
    }
}

fn trend(label: &str, name: &str, running: &[f64], from: u64, to: u64) {
    if to < from * 3 {
        return;
    }
    let mut least = f64::INFINITY;
    let mut most = f64::NEG_INFINITY;
    let mut total = 0.0f64;
    let mut count = 0u64;
    for r in from..=(to / 3) {
        let here = running[r as usize];
        let there = running[(r * 3) as usize];
        if here <= 0.0 || there <= 0.0 {
            continue;
        }
        let step = (there / here).ln() / 3.0f64.ln();
        least = least.min(step);
        most = most.max(step);
        total += step;
        count += 1;
    }
    if count == 0 {
        return;
    }
    let mut sum_x = 0.0f64;
    let mut sum_y = 0.0f64;
    let mut sum_xx = 0.0f64;
    let mut sum_xy = 0.0f64;
    let mut points = 0.0f64;
    for r in from..=to {
        if running[r as usize] <= 0.0 {
            continue;
        }
        let x = (r as f64).ln();
        let y = running[r as usize].ln();
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        points += 1.0;
    }
    let fitted = (points * sum_xy - sum_x * sum_y) / (points * sum_xx - sum_x * sum_x);
    let ends =
        (running[to as usize] / running[from as usize]).ln() / (to as f64 / from as f64).ln();
    println!(
        "circle-crop trend {label} {name} r={from}..{} steps={count} min={:.6} mean={:.6} max={:.6} fitted={fitted:.6} ends={ends:.6}",
        to,
        least,
        total / count as f64,
        most
    );
}

fn windows(label: &str, series: &[(&str, Vec<f64>, u64)]) {
    for (name, values, limit) in series {
        let mut lows: Vec<(u64, f64, u64)> = Vec::new();
        let mut power = 1u64;
        while power <= *limit {
            let top = (power * 3 - 1).min(*limit);
            let mut best = 0.0f64;
            let mut spot = power;
            for r in power..=top {
                if values[r as usize] > best {
                    best = values[r as usize];
                    spot = r;
                }
            }
            lows.push((power, best, spot));
            power *= 3;
        }
        for step in 0..lows.len() {
            let (start, best, spot) = lows[step];
            let raw = if start > 1 && best > 0.0 {
                best.ln() / (start as f64).ln()
            } else {
                f64::NAN
            };
            let slope = if step + 1 < lows.len() && best > 0.0 && lows[step + 1].1 > 0.0 {
                (lows[step + 1].1 / best).ln() / 3.0f64.ln()
            } else {
                f64::NAN
            };
            println!(
                "circle-crop window {label} {name} k={step} r={start}..{} max={best:.3} at={spot} at_over_start={:.4} raw={raw:.6} slope={slope:.6}",
                (start * 3 - 1).min(*limit),
                spot as f64 / start as f64
            );
        }
    }
}

fn cap(dimension: usize, r: u64) -> f64 {
    match dimension {
        2 => (3 * r + 5) as f64,
        _ => std::f64::consts::PI * 3.0f64.sqrt() * (r * r + 1) as f64,
    }
}

fn means(label: &str, table: &Sweep, fill: u64, dimension: usize) {
    let r_max = table.r_max;
    let mass = fill as f64;
    let share = 3.0f64.powi(dimension as i32) / mass;
    let edge = dimension as i32 - 1;
    let mut prior_least = f64::NAN;
    let mut prior_most = f64::NAN;
    let mut power = 1u64;
    while power * 3 <= r_max {
        let start = power;
        let stop = power * 3 - 1;
        let span = (stop - start + 1) as f64;
        let level = power.ilog(3) + 1;
        let weight = share.powi(level as i32);
        let mut total = 0u64;
        let mut least = f64::INFINITY;
        let mut most = f64::NEG_INFINITY;
        let mut full_least = f64::INFINITY;
        let mut full_most = f64::NEG_INFINITY;
        let mut carried = 0.0f64;
        for r in start..=stop {
            assert!(table.cut(r) <= table.all_cut(r));
            total += table.cut(r);
            let phi = table.cut(r) as f64 * weight / table.all_cut(r) as f64;
            least = least.min(phi);
            most = most.max(phi);
            carried += phi;
            let density = table.all_cut(r) as f64 / (r as f64).powi(edge);
            full_least = full_least.min(density);
            full_most = full_most.max(density);
        }
        let exact_low =
            table.inside[(power * 3) as usize] as i64 - table.touch[start as usize] as i64;
        let exact_high =
            2 * (table.touch[(power * 3) as usize] as i64 - table.inside[start as usize] as i64);
        assert!(exact_low <= total as i64);
        assert!(total as i64 <= exact_high);
        let mut depth = 0u32;
        let mut reach = power;
        while reach * 3 <= r_max {
            reach *= 3;
            depth += 1;
        }
        let scale = mass.powi(depth as i32);
        let main_low = table.inside[reach as usize] as f64 / scale;
        let main_high = table.touch[reach as usize] as f64 / scale;
        let slack = (table.cut(power * 3) + table.cut(start)) as f64;
        let form_low = (mass - 1.0) * main_low - slack;
        let form_high = 2.0 * ((mass - 1.0) * main_high + slack);
        assert!(form_low <= total as f64);
        assert!(total as f64 <= form_high);
        let kappa_low = total as f64 / ((mass - 1.0) * main_high);
        let kappa_high = total as f64 / ((mass - 1.0) * main_low);
        println!(
            "circle-crop mean {label} k={} r={start}..{stop} sum={total} mean={:.6} exact_low={exact_low} exact_high={exact_high} form_low={:.6} form_high={:.6} kappa_low={:.6} kappa_high={:.6}",
            level - 1,
            total as f64 / span,
            (form_low * 1e6).floor() / 1e6,
            (form_high * 1e6).ceil() / 1e6,
            (kappa_low * 1e6).floor() / 1e6,
            (kappa_high * 1e6).ceil() / 1e6
        );
        let mean_phi = carried / span;
        let ground = weight * form_low.max(0.0) / span / cap(dimension, stop);
        assert!(mean_phi >= ground);
        println!(
            "circle-crop factor {label} k={} r={start}..{stop} level={level} min={:.6} mean={:.6} max={:.6} ground={:.6} min_step={:.6} max_step={:.6} full_low={:.6} full_high={:.6}",
            level - 1,
            (least * 1e6).floor() / 1e6,
            mean_phi,
            (most * 1e6).ceil() / 1e6,
            (ground * 1e6).floor() / 1e6,
            least - prior_least,
            most - prior_most,
            (full_least * 1e6).floor() / 1e6,
            (full_most * 1e6).ceil() / 1e6
        );
        prior_least = least;
        prior_most = most;
        power *= 3;
    }
}

fn root_floor(value: u64) -> u64 {
    if value == 0 {
        return 0;
    }
    let mut root = (value as f64).sqrt() as u64;
    while root * root > value {
        root -= 1;
    }
    while (root + 1) * (root + 1) <= value {
        root += 1;
    }
    root
}

fn ones_table(level: usize) -> Vec<u32> {
    let size = 3usize.pow(level as u32);
    let mut out = vec![0u32; size];
    for value in 1..size {
        out[value] = (out[value / 3] << 1) | u32::from(value % 3 == 1);
    }
    out
}

fn trim(value: f64, up: bool) -> f64 {
    if up {
        (value * 1e6).ceil() / 1e6
    } else {
        (value * 1e6).floor() / 1e6
    }
}

fn list(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:.6}"))
        .collect::<Vec<String>>()
        .join(",")
}

struct Shell {
    total: u64,
    filled: u64,
    bad: Vec<u64>,
    pair: Vec<u64>,
}

fn record(mask: u32, level: usize, out: &mut Shell) {
    out.total += 1;
    if mask == 0 {
        out.filled += 1;
        return;
    }
    let mut rest = mask;
    while rest != 0 {
        let low = rest.trailing_zeros() as usize;
        out.bad[low] += 1;
        let mut other = rest & (rest - 1);
        while other != 0 {
            out.pair[low * level + other.trailing_zeros() as usize] += 1;
            other &= other - 1;
        }
        rest &= rest - 1;
    }
}

fn shell(ones: &[u32], dimension: usize, radius: u64, level: usize, out: &mut Shell) {
    out.total = 0;
    out.filled = 0;
    for slot in out.bad.iter_mut() {
        *slot = 0;
    }
    for slot in out.pair.iter_mut() {
        *slot = 0;
    }
    let square = radius * radius;
    if dimension == 2 {
        for i in 0..=radius {
            let high = root_floor(square - i * i);
            let step = (i + 1) * (i + 1);
            let low = if step >= square {
                0
            } else {
                root_floor(square - step)
            };
            let first = ones[i as usize];
            for y in low..=high {
                record(first & ones[y as usize], level, out);
            }
        }
    } else {
        for i in 0..=radius {
            let first = ones[i as usize];
            let span = root_floor(square - i * i);
            for j in 0..=span {
                let base = i * i + j * j;
                let high = root_floor(square - base);
                let step = (i + 1) * (i + 1) + (j + 1) * (j + 1);
                let low = if step >= square {
                    0
                } else {
                    root_floor(square - step)
                };
                let second = ones[j as usize];
                let both = first & second;
                let either = first | second;
                for z in low..=high {
                    record(both | (either & ones[z as usize]), level, out);
                }
            }
        }
    }
}

fn digit_census(label: &str, table: &Sweep, fill: u64, dimension: usize) {
    let r_max = table.r_max;
    let mass = fill as f64;
    let room = 3.0f64.powi(dimension as i32);
    let share = room / mass;
    let null = 1.0 - mass / room;
    let ones = ones_table((r_max + 1).ilog(3) as usize);
    let mut fine_low = f64::INFINITY;
    let mut fine_high = f64::NEG_INFINITY;
    let mut scale_low = f64::INFINITY;
    let mut scale_high = f64::NEG_INFINITY;
    let mut near_low = f64::INFINITY;
    let mut near_high = f64::NEG_INFINITY;
    let mut far_low = f64::INFINITY;
    let mut far_high = f64::NEG_INFINITY;
    let mut load_high = f64::NEG_INFINITY;
    let mut ind_least = f64::INFINITY;
    let mut ind_most = f64::NEG_INFINITY;
    let mut psi_least = f64::INFINITY;
    let mut psi_most = f64::NEG_INFINITY;
    let mut count = 0u64;
    let mut power = 1u64;
    while power * 3 <= r_max {
        let start = power;
        let stop = power * 3 - 1;
        let span = (stop - start + 1) as f64;
        let level = power.ilog(3) as usize + 1;
        let mut state = Shell {
            total: 0,
            filled: 0,
            bad: vec![0; level],
            pair: vec![0; level * level],
        };
        let mut sum_total = 0u64;
        let mut sum_bad = vec![0u64; level];
        let mut sum_pair = vec![0u64; level * level];
        let mut rate_low = vec![f64::INFINITY; level];
        let mut rate_high = vec![f64::NEG_INFINITY; level];
        let mut rate_sum = vec![0.0f64; level];
        let mut ind_low = f64::INFINITY;
        let mut ind_high = f64::NEG_INFINITY;
        let mut ind_sum = 0.0f64;
        let mut psi_low = f64::INFINITY;
        let mut psi_high = f64::NEG_INFINITY;
        let mut psi_sum = 0.0f64;
        let mut load_sum = 0.0f64;
        let mut load_most = f64::NEG_INFINITY;
        for r in start..=stop {
            shell(&ones, dimension, r, level, &mut state);
            assert_eq!(state.total, table.all_cut(r));
            assert_eq!(state.filled, table.cut(r));
            if dimension == 2 {
                assert_eq!(state.total, 2 * r + 1);
            }
            let mut union = 0u64;
            let mut product = 1.0f64;
            for j in 0..level {
                if dimension == 2 {
                    let parents = 2 * (r / 3u64.pow(j as u32 + 1)) + 1;
                    assert!(state.bad[j] <= 2 * 3u64.pow(j as u32) * parents);
                }
                union += state.bad[j];
                sum_bad[j] += state.bad[j];
                let rate = state.bad[j] as f64 / state.total as f64;
                rate_low[j] = rate_low[j].min(rate);
                rate_high[j] = rate_high[j].max(rate);
                rate_sum[j] += rate;
                product *= 1.0 - rate;
            }
            assert!(state.filled + union >= state.total);
            assert!(product > 0.0);
            sum_total += state.total;
            for (sum, &pair) in sum_pair.iter_mut().zip(&state.pair) {
                *sum += pair;
            }
            let load = union as f64 / state.total as f64;
            load_sum += load;
            load_most = load_most.max(load);
            let independent = product * share.powi(level as i32);
            let correction = state.filled as f64 / state.total as f64 / product;
            ind_low = ind_low.min(independent);
            ind_high = ind_high.max(independent);
            ind_sum += independent;
            psi_low = psi_low.min(correction);
            psi_high = psi_high.max(correction);
            psi_sum += correction;
        }
        let scale = sum_total as f64;
        let mean: Vec<f64> = rate_sum.iter().map(|value| value / span).collect();
        let scaled: Vec<f64> = mean
            .iter()
            .enumerate()
            .map(|(j, value)| (value - null) * start as f64 / 3.0f64.powi(j as i32))
            .collect();
        let least: Vec<f64> = rate_low.iter().map(|value| trim(*value, false)).collect();
        let most: Vec<f64> = rate_high.iter().map(|value| trim(*value, true)).collect();
        let mut near: Vec<f64> = Vec::new();
        let mut far: Vec<f64> = Vec::new();
        for j in 0..level {
            for gap in 1..=2usize {
                if j + gap >= level {
                    continue;
                }
                let joint = sum_pair[j * level + j + gap] as f64 / scale;
                let apart = (sum_bad[j] as f64 / scale) * (sum_bad[j + gap] as f64 / scale);
                let value = if apart > 0.0 { joint / apart } else { f64::NAN };
                if gap == 1 {
                    near_low = near_low.min(value);
                    near_high = near_high.max(value);
                    near.push(value);
                } else {
                    far_low = far_low.min(value);
                    far_high = far_high.max(value);
                    far.push(value);
                }
            }
        }
        for j in 0..level {
            scale_low = scale_low.min(scaled[j]);
            scale_high = scale_high.max(scaled[j]);
            if j + 3 < level {
                fine_low = fine_low.min(mean[j]);
                fine_high = fine_high.max(mean[j]);
            }
        }
        load_high = load_high.max(load_most);
        ind_least = ind_least.min(ind_low);
        ind_most = ind_most.max(ind_high);
        psi_least = psi_least.min(psi_low);
        psi_most = psi_most.max(psi_high);
        count += 1;
        println!(
            "circle-crop digits {label} k={} r={start}..{stop} level={level} null={null:.6} ind_low={:.6} ind_mean={:.6} ind_high={:.6} psi_low={:.6} psi_mean={:.6} psi_high={:.6} union_mean={:.6} union_high={:.6}",
            level - 1,
            trim(ind_low, false),
            ind_sum / span,
            trim(ind_high, true),
            trim(psi_low, false),
            psi_sum / span,
            trim(psi_high, true),
            load_sum / span,
            trim(load_most, true)
        );
        println!(
            "circle-crop digitrate {label} k={} r={start}..{stop} p_low={} p_mean={} p_high={} p_scaled={}",
            level - 1,
            list(&least),
            list(&mean),
            list(&most),
            list(&scaled)
        );
        println!(
            "circle-crop digitpair {label} k={} r={start}..{stop} rho1={} rho2={}",
            level - 1,
            list(&near),
            list(&far)
        );
        power *= 3;
    }
    println!(
        "circle-crop digittotal {label} windows={count} fine_low={:.6} fine_high={:.6} scaled_low={:.6} scaled_high={:.6} rho1_low={:.6} rho1_high={:.6} rho2_low={:.6} rho2_high={:.6} union_high={:.6} ind_low={:.6} ind_high={:.6} psi_low={:.6} psi_high={:.6}",
        trim(fine_low, false),
        trim(fine_high, true),
        trim(scale_low, false),
        trim(scale_high, true),
        trim(near_low, false),
        trim(near_high, true),
        trim(far_low, false),
        trim(far_high, true),
        trim(load_high, true),
        trim(ind_least, false),
        trim(ind_most, true),
        trim(psi_least, false),
        trim(psi_most, true)
    );
}

fn corner(tag: &str, dimension: usize, levels: &[usize], fill: u64) -> (u64, u64) {
    let mut prior: Option<Sweep> = None;
    let mut tally = (0u64, 0u64);
    for (step, level) in levels.iter().enumerate() {
        let types = design(dimension, *level);
        let side = types.shape[0];
        let r_max = (side - 1) as u64;
        let table = sweep(&types, 0, r_max);
        let label = format!("{tag} corner L={level}");
        for r in 1..=r_max {
            let at = r as usize;
            assert!(table.inside[at] <= table.seen[at]);
            assert!(table.seen[at] <= table.touch[at]);
            assert!(table.all_inside[at] <= table.all_seen[at]);
            assert!(table.all_seen[at] <= table.all_touch[at]);
            assert!(table.all_cut(r) as f64 <= cap(dimension, r));
            if dimension == 2 {
                assert_eq!(table.all_cut(r), 2 * r + 1);
            }
            let column = (r as f64 / ((dimension - 1) as f64).sqrt()).powi(dimension as i32 - 1);
            assert!(table.all_cut(r) as f64 >= column);
        }
        if let Some(past) = &prior {
            for r in 1..=past.r_max {
                let at = r as usize;
                assert_eq!(table.seen[at], past.seen[at]);
                assert_eq!(table.inside[at], past.inside[at]);
                assert_eq!(table.touch[at], past.touch[at]);
            }
            println!(
                "circle-crop levelfree {label} matches L={} on r=1..{}",
                levels[step - 1],
                past.r_max
            );
        }
        let width = match (dimension, side) {
            (_, 0..=243) => 5,
            (2, 244..=729) => 4,
            (2, 730..=2187) => 3,
            (2, 2188..=6561) => 2,
            (2, _) => 0,
            (_, _) => 0,
        };
        let sample: Vec<u64> = [13u64, 7, 5, 3, 2][..width]
            .iter()
            .map(|d| (r_max / d).max(1))
            .collect();
        oracle(&label, &types, 0, side, &table, &sample);
        if step + 1 == levels.len() {
            tally = report(&label, &table, fill, dimension);
            digit_census(&label, &table, fill, dimension);
        }
        prior = Some(table);
    }
    tally
}

fn report(label: &str, table: &Sweep, fill: u64, dimension: usize) -> (u64, u64) {
    let r_max = table.r_max;
    let reach_max = r_max / 3;
    let width = (r_max + 2) as usize;
    let mut low = vec![0.0f64; width];
    let mut high = vec![0.0f64; width];
    let mut swing = vec![0.0f64; width];
    let mut delta = vec![0.0f64; width];
    let mut rows: Vec<String> = Vec::new();
    let mut banded = 0u64;
    for r in 1..=r_max {
        let at = r as usize;
        let mut depth = 0u32;
        let mut reach = r;
        while reach * 3 <= r_max {
            reach *= 3;
            depth += 1;
        }
        if depth > 0 {
            banded += 1;
        }
        let scale = (fill as f64).powi(depth as i32);
        let centre = table.seen[at] as f64 - table.seen[reach as usize] as f64 / scale;
        let band = table.cut(reach) as f64 / scale;
        swing[at] = centre;
        low[at] = (centre.abs() - band).max(0.0);
        high[at] = centre.abs() + band;
        assert!(low[at] <= table.cut(r) as f64);
        let mut line = format!(
            "circle-crop row {label} r={r} res={} N={} in={} cut={} Nfull={} Nfull_in={} Nfull_cut={}",
            u64::from(r.is_power_of_three()),
            table.seen[at],
            table.inside[at],
            table.cut(r),
            table.all_seen[at],
            table.all_inside[at],
            table.all_cut(r)
        );
        if r <= reach_max {
            let jump = table.seen[(r * 3) as usize] as i64 - fill as i64 * table.seen[at] as i64;
            assert!(jump.unsigned_abs() <= table.cut(r * 3) + fill * table.cut(r));
            delta[at] = jump.abs() as f64;
            line.push_str(&format!(" delta={jump}"));
        } else {
            line.push_str(" delta=na");
        }
        line.push_str(&format!(
            " E={centre:.6} Elow={:.6} Ehigh={:.6} depth={depth}",
            (low[at] * 1e6).floor() / 1e6,
            (high[at] * 1e6).ceil() / 1e6
        ));
        rows.push(line);
    }
    let mut running_delta = vec![0.0f64; width];
    let mut running_cut = vec![0.0f64; width];
    for r in 1..=r_max {
        let at = r as usize;
        running_delta[at] = running_delta[at - 1].max(delta[at]);
        running_cut[at] = running_cut[at - 1].max(table.cut(r) as f64);
        println!(
            "{} run_delta={:.0} run_cut={:.0}",
            rows[at - 1],
            running_delta[at],
            running_cut[at]
        );
    }
    trend(label, "delta", &running_delta, 27, reach_max);
    trend(label, "cut", &running_cut, 27, r_max);
    let mut power = 1u64;
    while power <= r_max {
        let at = power as usize;
        let over = if power >= 3 {
            table.seen[at] as f64 / table.seen[(power / 3) as usize] as f64
        } else {
            f64::NAN
        };
        let mut best = 0.0f64;
        let mut under = 0u64;
        let mut span = 0u64;
        let top = (power * 3 - 1).min(reach_max);
        if power <= reach_max {
            for r in power..=top {
                best = best.max(delta[r as usize]);
                span += 1;
                if delta[r as usize] <= delta[at] {
                    under += 1;
                }
            }
        }
        let rank = if span > 0 {
            under as f64 / span as f64
        } else {
            f64::NAN
        };
        let mut quiet = 0u64;
        let mut reach_span = 0u64;
        for r in power..=(power * 3 - 1).min(r_max) {
            reach_span += 1;
            if table.cut(r) <= table.cut(power) {
                quiet += 1;
            }
        }
        let cut_rank = quiet as f64 / reach_span as f64;
        let defect = if power <= reach_max {
            format!(
                "delta={:.0} window_max_delta={best:.0} delta_rank={rank:.4} E={:.6}",
                delta[at], swing[at]
            )
        } else {
            String::from(
                "delta=beyond_reach window_max_delta=beyond_reach delta_rank=beyond_reach E=beyond_reach",
            )
        };
        let mass = (fill as f64).powi(power.ilog(3) as i32);
        let mu_low = (table.inside[at] as f64 / mass * 1e9).floor() / 1e9;
        let mu_high = ((table.inside[at] + table.cut(power)) as f64 / mass * 1e9).ceil() / 1e9;
        println!(
            "circle-crop resonance {label} r={power} N={} ratio={over:.6} in={} cut={} {defect} cut_rank={cut_rank:.4} mu_low={mu_low:.9} mu_high={mu_high:.9}",
            table.seen[at],
            table.inside[at],
            table.cut(power)
        );
        power *= 3;
    }
    let cuts: Vec<f64> = (0..width)
        .map(|r| table.cut((r as u64).min(r_max)) as f64)
        .collect();
    windows(
        label,
        &[
            ("cut", cuts, r_max),
            ("delta", delta, reach_max),
            ("errorlow", low, reach_max),
            ("errorhigh", high, reach_max),
        ],
    );
    means(label, table, fill, dimension);
    println!("circle-crop bands {label} rows={r_max} banded={banded}");
    (r_max, banded)
}

trait Triadic {
    fn is_power_of_three(&self) -> bool;
}

impl Triadic for u64 {
    fn is_power_of_three(&self) -> bool {
        let mut value = *self;
        if value == 0 {
            return false;
        }
        while value.is_multiple_of(3) {
            value /= 3;
        }
        value == 1
    }
}

fn centre(tag: &str, dimension: usize, level: usize) -> u64 {
    let types = design(dimension, level);
    let side = types.shape[0];
    let r_max = ((side - 1) / 2) as u64;
    let table = sweep(&types, side as i64, r_max);
    let label = format!("{tag} centre L={level}");
    let cells = (side as u64).pow(dimension as u32);
    let filled = types.sum();
    println!("circle-crop density {label} side={side} cells={cells} fill={filled}");
    let sample: Vec<u64> = [11u64, 5, 3, 2]
        .iter()
        .map(|d| (r_max / d).max(1))
        .collect();
    let sample = if level >= 7 || (dimension == 3 && level >= 5) {
        sample[..2].to_vec()
    } else {
        sample
    };
    oracle(&label, &types, side as i64, side, &table, &sample);
    let mut first = 0u64;
    let mut relative = vec![0.0f64; (r_max + 2) as usize];
    for r in 1..=r_max {
        let at = r as usize;
        assert!(table.inside[at] <= table.seen[at]);
        assert!(table.seen[at] <= table.touch[at]);
        if first == 0 && table.seen[at] > 0 {
            first = r;
        }
        let exact =
            cells as i64 * table.seen[at] as i64 - filled as i64 * table.all_seen[at] as i64;
        let error = exact as f64 / cells as f64;
        let main = filled as f64 / cells as f64 * table.all_seen[at] as f64;
        relative[at] = if main > 0.0 {
            (error / main).abs()
        } else {
            0.0
        };
        println!(
            "circle-crop row {label} r={r} res={} N={} in={} cut={} Nfull={} err_num={exact} err={error:.6} rel={:.6}",
            u64::from(r.is_power_of_three()),
            table.seen[at],
            table.inside[at],
            table.cut(r),
            table.all_seen[at],
            relative[at]
        );
    }
    let hole = (side / 3 - 1) as u64 / 2;
    assert!(first > hole);
    if dimension == 2 {
        assert_eq!(first, hole + 1);
    }
    println!("circle-crop hole {label} block_inradius={hole} first_hit={first}");
    let cuts: Vec<f64> = (0..=(r_max + 1))
        .map(|r| table.cut(r.min(r_max)) as f64)
        .collect();
    windows(
        &label,
        &[("cut", cuts, r_max), ("relative", relative, r_max)],
    );
    r_max
}

fn gauss(dimension: usize, level: usize) {
    let types = design(dimension, level);
    let side = types.shape[0];
    let r_max = (side - 1) as u64;
    let table = sweep(&types, 0, r_max);
    let pi = std::f64::consts::PI;
    let mut worst = 0.0f64;
    for r in 1..=r_max {
        let volume = if dimension == 2 {
            pi * (r as f64) * (r as f64) / 4.0
        } else {
            pi * (r as f64).powi(3) / 6.0
        };
        let gap = (table.all_seen[r as usize] as f64 - volume).abs()
            / (r as f64).powi(dimension as i32 - 1);
        worst = worst.max(gap);
    }
    println!(
        "circle-crop gauss D={dimension} L={level} r_max={r_max} max|Nfull-vol|/r^(D-1)={worst:.6}"
    );
}

fn digits(dimension: usize) -> Vec<Vec<f64>> {
    let types = design(dimension, 1);
    let mut out: Vec<Vec<f64>> = Vec::new();
    let mut index = vec![0usize; dimension];
    for flat in 0..types.size() {
        if types.at(flat) != 0 {
            out.push(index.iter().map(|value| *value as f64).collect());
        }
        for axis in (0..dimension).rev() {
            index[axis] += 1;
            if index[axis] < 3 {
                break;
            }
            index[axis] = 0;
        }
    }
    out
}

fn factor(set: &[Vec<f64>], point: &[f64]) -> (f64, f64) {
    let mut real = 0.0f64;
    let mut imaginary = 0.0f64;
    for digit in set {
        let dot: f64 = digit.iter().zip(point).map(|(a, b)| a * b).sum();
        let angle = -2.0 * std::f64::consts::PI * dot;
        real += angle.cos();
        imaginary += angle.sin();
    }
    let mass = set.len() as f64;
    (real / mass, imaginary / mass)
}

fn hat(set: &[Vec<f64>], point: &[f64], terms: u32) -> (f64, f64) {
    let mut real = 1.0f64;
    let mut imaginary = 0.0f64;
    for term in 1..=terms {
        let scale = 3.0f64.powi(term as i32);
        let scaled: Vec<f64> = point.iter().map(|value| value / scale).collect();
        let (pr, pi) = factor(set, &scaled);
        let next = (real * pr - imaginary * pi, real * pi + imaginary * pr);
        real = next.0;
        imaginary = next.1;
    }
    (real, imaginary)
}

fn step(set: &[Vec<f64>], point: &[f64]) -> f64 {
    let dimension = point.len();
    let mut total = 0.0f64;
    let count = 3usize.pow(dimension as u32);
    for code in 0..count {
        let mut shifted = point.to_vec();
        let mut rest = code;
        for value in shifted.iter_mut() {
            *value += (rest % 3) as f64 / 3.0;
            rest /= 3;
        }
        let (re, im) = factor(set, &shifted);
        total += (re * re + im * im).sqrt();
    }
    total
}

fn mass(tag: &str, dimension: usize, top: u32, sharp: f64) {
    let set = digits(dimension);
    let zero = vec![0.0f64; dimension];
    let ground = step(&set, &zero);
    let witness: Vec<f64> = vec![155.0 / 243.0; dimension];
    let broken = step(&set, &witness);
    assert!((ground - 2.0).abs() < 1e-9);
    assert!(broken < 2.0);
    println!(
        "circle-crop transform {tag} step h_lattice={ground:.6} h_witness={:.6} at=155/243",
        (broken * 1e6).ceil() / 1e6
    );
    let mut prior = f64::NAN;
    for level in 1..=top {
        let side = 3usize.pow(level);
        let count = side.pow(dimension as u32);
        let mut total = 0.0f64;
        let mut index = vec![0usize; dimension];
        for _ in 0..count {
            let point: Vec<f64> = index.iter().map(|value| *value as f64).collect();
            let mut piece = 1.0f64;
            for term in 1..=level {
                let scale = 3.0f64.powi(term as i32);
                let scaled: Vec<f64> = point.iter().map(|value| value / scale).collect();
                let (re, im) = factor(&set, &scaled);
                piece *= (re * re + im * im).sqrt();
            }
            total += piece;
            for axis in (0..dimension).rev() {
                index[axis] += 1;
                if index[axis] < side {
                    break;
                }
                index[axis] = 0;
            }
        }
        assert!(total >= 2.0f64.powi(level as i32) - 1e-9);
        let bulk = total - 1.0;
        let ratio = bulk / prior;
        let stride = ratio.ln() / 3.0f64.ln();
        println!(
            "circle-crop transform {tag} mass L={level} lambda={bulk:.6} step={ratio:.6} log3={stride:.6} cost={:.6} error={:.6}",
            stride + 0.5,
            stride + 0.5 + sharp - 2.0
        );
        prior = bulk;
    }
}

fn transform(tag: &str, dimension: usize, point: &[f64]) {
    let set = digits(dimension);
    let tripled: Vec<f64> = point.iter().map(|value| value * 3.0).collect();
    let (ar, ai) = hat(&set, point, 60);
    let (br, bi) = hat(&set, &tripled, 60);
    let (pr, pi) = factor(&set, point);
    assert!((br - (pr * ar - pi * ai)).abs() < 1e-12);
    assert!((bi - (pr * ai + pi * ar)).abs() < 1e-12);
    let integral = point
        .iter()
        .all(|value| (value - value.round()).abs() < 1e-12);
    let size = (pr * pr + pi * pi).sqrt();
    assert_eq!(integral, (size - 1.0).abs() < 1e-12);
    let name: Vec<String> = point.iter().map(|value| format!("{value}")).collect();
    println!(
        "circle-crop transform {tag} t=({}) integral={integral} hat_mu={:.8} hat_mu_3t={:.8} P={size:.8}",
        name.join(","),
        (ar * ar + ai * ai).sqrt(),
        (br * br + bi * bi).sqrt()
    );
}

fn main() {
    println!("circle-crop generator: cargo run --release --manifest-path research/crop/lab/rs/circle-crop/Cargo.toml");
    println!("circle-crop convention: cells are indexed x in [0,3^L)^D, a cell counts when its centre x+1/2 lies in the closed ball, corner balls sit at the lattice corner 0 and centre balls at the grid centre 3^L/2");
    gauss(2, 5);
    gauss(3, 3);
    transform("carpet", 2, &[0.5, 0.0]);
    transform("carpet", 2, &[1.0, 0.0]);
    transform("sponge", 3, &[0.5, 0.0, 0.0]);
    transform("sponge", 3, &[1.0, 0.0, 0.0]);
    mass("carpet", 2, 6, 8.0f64.ln() / 3.0f64.ln());
    let carpet = corner("carpet", 2, &[5, 6, 7, 8, 9], 8);
    let sponge = corner("sponge", 3, &[3, 4, 5, 6], 20);
    let mut rows = carpet.0 + sponge.0;
    let banded = carpet.1 + sponge.1;
    rows += centre("carpet", 2, 6);
    rows += centre("carpet", 2, 7);
    rows += centre("sponge", 3, 4);
    rows += centre("sponge", 3, 5);
    transfer::transfer();
    derive::derive();
    println!("circle-crop totals rows={rows} banded={banded}");
}
