use crate::root_floor;
use crate::transfer::{levels, masks_seen, rung, show};

// PATTERN LAW

fn mask_of(steps: [i64; 4], reach: usize) -> u16 {
    let mut mask = 0u16;
    for column in 0..reach.min(3) {
        let low = steps[column + 1].max(0);
        let high = steps[column].min(2);
        let mut row = low;
        while row <= high {
            mask |= 1 << (3 * column as u16 + row as u16);
            row += 1;
        }
    }
    mask
}

fn window(steps: [i64; 4]) -> Option<(i64, i64)> {
    let mut low = i64::MIN;
    let mut high = i64::MAX;
    for near in 0..4usize {
        for far in 0..4usize {
            if far > near {
                let gap = (far - near) as i64;
                low = low.max(6 * (steps[near] - steps[far] - 1) / gap);
            }
            if far < near {
                let gap = (near - far) as i64;
                high = high.min(6 * (steps[far] - steps[near] + 1) / gap);
            }
        }
    }
    let floor = low.max(0);
    if high > floor {
        Some((floor, high))
    } else {
        None
    }
}

// LINE ALPHABET

fn line_alphabet(span: i64) -> Vec<(u16, i64, i64)> {
    let mut found: Vec<Option<(i64, i64)>> = vec![None; 512];
    for a in -span..=span {
        for b in -span..=span {
            for c in -span..=span {
                for d in -span..=span {
                    let steps = [a, b, c, d];
                    let Some((low, high)) = window(steps) else {
                        continue;
                    };
                    let mask = mask_of(steps, 3);
                    if mask == 0 {
                        continue;
                    }
                    let slot = &mut found[mask as usize];
                    match slot {
                        None => *slot = Some((low, high)),
                        Some(range) => {
                            range.0 = range.0.min(low);
                            range.1 = range.1.max(high);
                        }
                    }
                }
            }
        }
    }
    (0..512)
        .filter_map(|mask| found[mask].map(|range| (mask as u16, range.0, range.1)))
        .collect()
}

// SHELL AGAINST THE LINE

fn compare(label: &str, radius: u64) {
    let deep = levels(radius);
    let seen = masks_seen(radius, 1, deep - 3);
    let all = masks_seen(radius, 1, deep);
    let alphabet = line_alphabet(9);
    let line: Vec<u16> = alphabet.iter().map(|entry| entry.0).collect();
    let outside: Vec<u16> = seen.iter().cloned().filter(|m| !line.contains(m)).collect();
    let missing: Vec<u16> = line.iter().cloned().filter(|m| !seen.contains(m)).collect();
    let curved: Vec<u16> = all.iter().cloned().filter(|m| !line.contains(m)).collect();
    assert!(outside.is_empty() && missing.is_empty());
    assert_eq!(curved.len(), 1);
    println!(
        "alphabet {} r={} line={} shell={} shell_outside_line=[{}] line_missing_from_shell=[{}] every_level_outside_line=[{}]",
        label,
        radius,
        line.len(),
        seen.len(),
        show(&outside),
        show(&missing),
        show(&curved)
    );
}

// FROZEN SLOPE

fn frozen(label: &str, radius: u64) {
    let deep = levels(radius);
    for level in 0..deep.saturating_sub(1) {
        let below = rung(radius, level);
        let here = rung(radius, level + 1);
        let scale = 3u64.pow(level) as f64;
        let far = (radius * radius) as f64;
        let top_here = here.high.len() - 2;
        let top_below = below.high.len() - 2;
        let mut columns = 0u64;
        let mut column_hits = 0u64;
        let mut boxes = 0u64;
        let mut faults = 0u64;
        let mut chord_faults = 0u64;
        let mut shallow = 0u64;
        let mut shallow_faults = 0u64;
        let mut worst = 0usize;
        for i in 0..=top_here {
            columns += 1;
            let edge = (3 * i) as f64 * scale;
            let height = (far - edge * edge).max(0.0).sqrt() / scale;
            if height <= 0.0 {
                continue;
            }
            let slope = (3 * i) as f64 / height;
            let far_edge = (3 * (i + 1)) as f64 * scale;
            let end = (far - far_edge * far_edge).max(0.0).sqrt() / scale;
            let chord = (height - end) / 3.0;
            let mut guess = [0i64; 4];
            let mut rope = [0i64; 4];
            guess[0] = below.high[3 * i] as i64;
            rope[0] = guess[0];
            for k in 1..4usize {
                guess[k] = (height - k as f64 * slope).floor() as i64;
                rope[k] = (height - k as f64 * chord).floor() as i64;
            }
            let reach = (top_below + 1).saturating_sub(3 * i).min(3);
            let mut exact = true;
            for (k, &g) in guess.iter().enumerate().skip(1) {
                if 3 * i + k <= top_below + 1 && g != below.high[3 * i + k] as i64 {
                    exact = false;
                }
            }
            if exact {
                column_hits += 1;
            }
            let low = here.high[i + 1];
            let high = here.high[i];
            for y in low..=high {
                boxes += 1;
                let mut truth = [0i64; 4];
                for (k, t) in truth.iter_mut().enumerate() {
                    let column = 3 * i + k;
                    *t = if column <= top_below + 1 {
                        below.high[column] as i64 - 3 * y as i64
                    } else {
                        -1
                    };
                }
                let shifted = [
                    guess[0] - 3 * y as i64,
                    guess[1] - 3 * y as i64,
                    guess[2] - 3 * y as i64,
                    guess[3] - 3 * y as i64,
                ];
                let roped = [
                    rope[0] - 3 * y as i64,
                    rope[1] - 3 * y as i64,
                    rope[2] - 3 * y as i64,
                    rope[3] - 3 * y as i64,
                ];
                let seen = mask_of(truth, reach);
                if seen != mask_of(shifted, reach) {
                    faults += 1;
                    worst = worst.max(i);
                }
                let missed = seen != mask_of(roped, reach);
                if missed {
                    chord_faults += 1;
                }
                if slope <= 1.0 {
                    shallow += 1;
                    if missed {
                        shallow_faults += 1;
                    }
                }
            }
        }
        println!(
            "frozen {} r={} level={} columns={} column_exact={} boxes={} tangent_faults={} tangent_rate={:.6} chord_faults={} chord_rate={:.6} shallow_boxes={} shallow_chord_faults={} shallow_rate={:.6} last_fault_column={} of {}",
            label,
            radius,
            level,
            columns,
            column_hits,
            boxes,
            faults,
            faults as f64 / boxes.max(1) as f64,
            chord_faults,
            chord_faults as f64 / boxes.max(1) as f64,
            shallow,
            shallow_faults,
            shallow_faults as f64 / shallow.max(1) as f64,
            worst,
            top_here
        );
    }
}

// CENTRE SEAT

fn centre(label: &str, radius: u64) {
    let deep = levels(radius);
    for level in 0..deep.saturating_sub(1) {
        let below = rung(radius, level);
        let here = rung(radius, level + 1);
        let top_here = here.high.len() - 2;
        let top_below = below.high.len() - 2;
        let mut boxes = 0u64;
        let mut seated = 0u64;
        for i in 0..=top_here {
            let low = here.high[i + 1];
            let high = here.high[i];
            for y in low..=high {
                boxes += 1;
                let column = 3 * i + 1;
                if column <= top_below
                    && below.high[column + 1] <= 3 * y + 1
                    && 3 * y < below.high[column]
                {
                    seated += 1;
                }
            }
        }
        println!(
            "centre {} r={} level={} boxes={} centre_crossed={} rate={:.6} derived={:.6} ratio={:.6}",
            label,
            radius,
            level,
            boxes,
            seated,
            seated as f64 / boxes as f64,
            1.0 / 3.0,
            3.0 * seated as f64 / boxes as f64
        );
    }
}

// STRAIGHT LINE LADDER

fn ones(value: u64, depth: u32) -> u32 {
    let mut mask = 0u32;
    let mut left = value;
    for slot in 0..depth {
        if left % 3 == 1 {
            mask |= 1 << slot;
        }
        left /= 3;
    }
    mask
}

fn ladder(name: &str, num: i128, den: i128, off: i128, depth: u32) {
    let width = 3i128.pow(depth);
    let height = |x: i128| -> i128 { (num * (width - x) + off) / den };
    let mut total = 0u64;
    let mut alive = 0u64;
    let mut hits = vec![0u64; depth as usize];
    for x in 0..width {
        let column = ones(x as u64, depth);
        let low = height(x + 1);
        let high = height(x);
        for z in low..=high {
            total += 1;
            let both = column & ones(z as u64, depth);
            if both == 0 {
                alive += 1;
            }
            let mut left = both;
            while left != 0 {
                let slot = left.trailing_zeros() as usize;
                hits[slot] += 1;
                left &= left - 1;
            }
        }
    }
    let mut product = 1.0f64;
    for &h in &hits {
        product *= 1.0 - h as f64 / total as f64;
    }
    let survival = alive as f64 / total as f64;
    let psi = survival / product;
    let independent = (8.0f64 / 9.0).powi(depth as i32);
    println!(
        "line {} slope={:.12} shift={:.6} depth={} cells={} alive={} survival={:.9} marginal={:.9} psi={:.9} logpsi={:.6} rate={:.9} eight_ninths={:.9}",
        name,
        num as f64 / den as f64,
        off as f64 / den as f64,
        depth,
        total,
        alive,
        survival,
        product,
        psi,
        psi.ln(),
        survival.powf(1.0 / depth as f64),
        independent.powf(1.0 / depth as f64)
    );
}

// RESONANCE TRACKING

fn coprime(a: i128, b: i128) -> bool {
    let mut x = a;
    let mut y = b;
    while y != 0 {
        let t = x % y;
        x = y;
        y = t;
    }
    x == 1
}

fn steeper(u: i128, r: i128, p: i128, q: i128) -> bool {
    if p < 0 {
        return true;
    }
    u * u * (q * q + p * p) > p * p * r * r
}

fn shallower(u: i128, r: i128, p: i128, q: i128) -> bool {
    u * u * (q * q + p * p) < p * p * r * r
}

fn cap_levels(r: i128, b: i128) -> u32 {
    let mut count = 0u32;
    let mut scale = 1i128;
    while scale * b * b < 2 * r {
        count += 1;
        scale *= 3;
    }
    count
}

fn floor_levels(r: i128, a: i128, b: i128) -> u32 {
    if a < 1 || a * b + 1 > b * b {
        return 0;
    }
    let mut count = 0u32;
    let mut scale = 1i128;
    while 8 * scale * scale * b * b * b * b < r * r {
        count += 1;
        scale *= 3;
    }
    count
}

fn track(radius: u64, a: i128, b: i128) -> (u32, u64, [u64; 3]) {
    let r = radius as i128;
    let low = a * b - 1;
    let high = a * b + 1;
    let den = b * b;
    let deep = levels(radius);
    let mut run = 0u32;
    let mut boxes = 0u64;
    let mut seats = [0u64; 3];
    for level in 0..deep {
        let scale = 3i128.pow(level);
        let top = r / scale;
        let mut here = 0u64;
        let mut mine = [0u64; 3];
        for i in 0..=top {
            let left = scale * i;
            let right = scale * (i + 1);
            if 2 * right * right > r * r {
                break;
            }
            if steeper(left, r, low, den) && shallower(right, r, high, den) {
                here += 1;
                mine[(i % 3) as usize] += 1;
            }
        }
        if here > 0 {
            assert_eq!(run, level);
            run = level + 1;
            boxes = here;
            seats = mine;
        }
    }
    (run, boxes, seats)
}

fn resonance(label: &str, radius: u64) -> u64 {
    let r = radius as i128;
    let deep = levels(radius);
    let mut rows = 0u64;
    let mut live = 0u64;
    let mut at_cap = 0u64;
    let mut at_floor = 0u64;
    let mut budget = vec![0i128; deep as usize];
    for b in 1..=30i128 {
        for a in 0..=b {
            if !coprime(a, b) {
                continue;
            }
            let cap = cap_levels(r, b);
            let low = floor_levels(r, a, b);
            let (run, boxes, seats) = track(radius, a, b);
            assert!(run <= cap);
            assert!(run >= low);
            rows += 1;
            if cap < deep {
                live += 1;
            }
            if run == cap {
                at_cap += 1;
            }
            if low > 0 && run == low {
                at_floor += 1;
            }
            for seen in budget.iter_mut().take(run as usize) {
                *seen = (*seen).max(b);
            }
            println!(
                "track {label} r={radius} slope={a}/{b} eps=1/{den} run={run} cap={cap} floor={low} boxes={boxes} seats_mod3={s0},{s1},{s2}",
                den = b * b,
                s0 = seats[0],
                s1 = seats[1],
                s2 = seats[2]
            );
        }
    }
    println!(
        "track {label} r={radius} levels={deep} slopes={rows} live_caps={live} run_equals_cap={at_cap} run_equals_floor={at_floor}"
    );
    for (level, &seen) in budget.iter().enumerate() {
        let scale = 3i128.pow(level as u32);
        assert!(seen * seen * scale < 2 * r);
        println!(
            "budget {label} r={radius} level={level} rank={rank} biggest_denominator={seen} denominator_square_times_scale={weight} twice_r={twice}",
            rank = deep as usize - 1 - level,
            weight = seen * seen * scale,
            twice = 2 * r
        );
    }
    rows + deep as u64 + 1
}

fn secant(label: &str, radius: u64) -> u64 {
    let r = radius as i128;
    let square = r * r;
    let deep = levels(radius);
    let mut best = 0u32;
    let mut witness = 0i128;
    for level in 1..deep {
        let width = 3i128.pow(level);
        let step = (width - 1) / 2;
        let top = r / width;
        for i in 0..=top {
            let u = width * i;
            if u + 2 * step >= r {
                break;
            }
            let v0 = root_floor((square - u * u) as u64) as i128;
            let v1 = root_floor((square - (u + step) * (u + step)) as u64) as i128;
            let v2 = root_floor((square - (u + 2 * step) * (u + 2 * step)) as u64) as i128;
            if (v0 - 2 * v1 + v2).abs() <= 1 {
                best = level;
                witness = u;
                break;
            }
        }
    }
    let step = (3i128.pow(best) - 1) / 2;
    let next = (3i128.pow(best + 1) - 1) / 2;
    assert!(step * step <= 3 * r);
    assert!(next * next > 3 * r);
    println!(
        "secant {label} r={radius} levels={deep} line_depth={best} column={witness} spacing={step} spacing_square={sq} three_r={three} next_spacing={next} next_square={nsq}",
        sq = step * step,
        three = 3 * r,
        nsq = next * next
    );
    1
}

fn blind(sigma_num: i128, sigma_den: i128, back: u32) -> u64 {
    let scale = 3i128.pow(back);
    assert_eq!(scale % sigma_den, 0);
    let lift = scale / sigma_den;
    let slow = sigma_num * lift;
    let fast = slow + 1;
    let mut differ = 0i128;
    let samples = 3 * scale;
    for t in 0..samples {
        let mut a = [0i64; 4];
        let mut b = [0i64; 4];
        for k in 0..4i128 {
            let one = 2 * t + 1 - 2 * k * slow;
            let two = 2 * t + 1 - 2 * k * fast;
            a[k as usize] = one.div_euclid(2 * scale) as i64;
            b[k as usize] = two.div_euclid(2 * scale) as i64;
        }
        if mask_of(a, 3) != mask_of(b, 3) {
            differ += 1;
        }
    }
    assert!(differ * scale <= 6 * samples);
    println!(
        "blind slope={sigma_num}/{sigma_den} eps=1/{scale} samples={samples} differ={differ} fraction={frac:.9} bound={bound:.9}",
        frac = differ as f64 / samples as f64,
        bound = 6.0 / scale as f64
    );
    1
}

// BOX AGAINST BLOCK

fn heights(radius: u64) -> Vec<i128> {
    let r = radius as i128;
    let square = r * r;
    let mut out = Vec::with_capacity(radius as usize + 2);
    for x in 0..=r + 1 {
        if x > r {
            out.push(0);
        } else {
            out.push(root_floor((square - x * x) as u64) as i128);
        }
    }
    out
}

fn line_bracket(v: &[i128]) -> Option<((i128, i128), (i128, i128))> {
    let span = v.len() as i128;
    let mut low: Option<(i128, i128)> = None;
    let mut high: Option<(i128, i128)> = None;
    for a in 0..span {
        for b in 0..span {
            if a == b {
                continue;
            }
            let gap = b - a;
            let rise = v[b as usize] - v[a as usize] + 1;
            if gap > 0 {
                high = Some(match high {
                    None => (rise, gap),
                    Some(seen) => {
                        if rise * seen.1 < seen.0 * gap {
                            (rise, gap)
                        } else {
                            seen
                        }
                    }
                });
            } else {
                low = Some(match low {
                    None => (-rise, -gap),
                    Some(seen) => {
                        if -rise * seen.1 > seen.0 * -gap {
                            (-rise, -gap)
                        } else {
                            seen
                        }
                    }
                });
            }
        }
    }
    let (a, b) = (low?, high?);
    if a.0 * b.1 < b.0 * a.1 {
        Some((a, b))
    } else {
        None
    }
}

fn finest(radius: i128, left: i128, right: i128) -> Option<(i128, i128)> {
    let mut best: Option<(i128, i128)> = None;
    for b in 1..=30i128 {
        for a in 0..=b {
            if !coprime(a, b) {
                continue;
            }
            if steeper(left, radius, a * b - 1, b * b)
                && shallower(right, radius, a * b + 1, b * b)
                && best.is_none_or(|seen| b > seen.1)
            {
                best = Some((a, b));
            }
        }
    }
    best
}

fn blocks(label: &str, radius: u64, level: u32) -> u64 {
    let r = radius as i128;
    let scale = 3i128.pow(level);
    let v = heights(radius);
    let side = (r / scale + 2) as usize;
    let mut extent = vec![(i128::MAX, i128::MIN, i128::MAX, i128::MIN); side * side];
    for x in 0..=r {
        let top = v[x as usize];
        let bottom = v[(x + 1) as usize];
        let column = x / scale;
        let mut row = bottom;
        while row <= top {
            let band = row / scale;
            let slot = &mut extent[column as usize * side + band as usize];
            slot.0 = slot.0.min(x);
            slot.1 = slot.1.max(x);
            slot.2 = slot.2.min(row);
            slot.3 = slot.3.max(row);
            row = (band + 1) * scale;
        }
    }
    let mut boxes = 0u64;
    let mut across = 0u64;
    let mut short = 0u64;
    let mut lines = 0u64;
    let mut rows = 0u64;
    for column in 0..side as i128 {
        for band in 0..side as i128 {
            let slot = extent[column as usize * side + band as usize];
            if slot.1 < slot.0 {
                continue;
            }
            boxes += 1;
            let wide = slot.1 - slot.0 + 1;
            let tall = slot.3 - slot.2 + 1;
            let step = (wide - 1) / 2;
            let lift = (tall - 1) / 2;
            if wide == scale {
                across += 1;
            }
            if step * step <= 3 * r && lift * lift <= 3 * r {
                short += 1;
            }
            let stop = (slot.1 + 1).min(r);
            let cut: Vec<i128> = (slot.0..=stop).map(|x| v[x as usize]).collect();
            let Some((low, high)) = line_bracket(&cut) else {
                continue;
            };
            lines += 1;
            assert!(step * step <= 3 * r);
            let window = finest(r, slot.0, slot.1);
            let (a, b) = window.unwrap();
            let reach = t_of(r, slot.1) - t_of(r, slot.0);
            println!(
                "boxline {label} r={radius} level={level} box={column},{band} columns={from}..{to} wide={wide} tall={tall} t_in=({sl:.7},{sh:.7}) span={reach:.9} block_floor={floor:.9} share={share:.6} finest_window={a}/{b}",
                from = slot.0,
                to = slot.1,
                sl = -(high.0 as f64) / high.1 as f64,
                sh = -(low.0 as f64) / low.1 as f64,
                floor = scale as f64 / r as f64,
                share = reach / (scale as f64 / r as f64)
            );
            rows += 1;
        }
    }
    println!(
        "boxes {label} r={radius} level={level} boxes={boxes} side_to_side={across} both_extents_short={short} line_boxes={lines} block_cap_level={cap}",
        cap = 3i128.pow(level)
    );
    rows + 1
}

fn t_of(r: i128, u: i128) -> f64 {
    u as f64 / ((r * r - u * u) as f64).sqrt()
}

// PASS

pub fn derive() {
    let alphabet = line_alphabet(9);
    let masks: Vec<u16> = alphabet.iter().map(|entry| entry.0).collect();
    for span in [12i64, 15] {
        let other: Vec<u16> = line_alphabet(span).iter().map(|entry| entry.0).collect();
        assert_eq!(other, masks);
    }
    println!(
        "line_alphabet states={} list=[{}]",
        masks.len(),
        show(&masks)
    );
    for radius in [6560u64, 19682, 12345] {
        compare("carpet", radius);
    }
    for radius in [6560u64, 19682] {
        frozen("carpet", radius);
        centre("carpet", radius);
    }
    let unit: i128 = 300000000000;
    let slopes: [(&str, i128, i128); 10] = [
        ("half", unit / 2, 0),
        ("half_shifted", unit / 2, 123456789012),
        ("third", unit / 3, 0),
        ("third_shifted", unit / 3, 123456789012),
        ("third_again", unit / 3, 271828182845),
        ("seventh", unit / 7, 0),
        ("seventh_shifted", unit / 7, 123456789012),
        ("root2m1", 124264068712, 123456789012),
        ("golden", 185410196625, 123456789012),
        ("quarterpi", 235619449020, 123456789012),
    ];
    for (name, num, off) in slopes {
        for depth in 4..=12u32 {
            ladder(name, num, unit, off, depth);
        }
    }
    for step in 2..=9u32 {
        let drift = unit / 3i128.pow(step);
        for depth in 4..=12u32 {
            ladder(
                &format!("near_third_3^-{step}"),
                unit / 3 + drift,
                unit,
                123456789012,
                depth,
            );
        }
    }
    let mut rows = 0u64;
    for radius in [728u64, 2186, 6560, 19682] {
        rows += resonance("carpet", radius);
        rows += secant("carpet", radius);
    }
    for back in 2..=7u32 {
        rows += blind(1, 3, back);
    }
    rows += blocks("carpet", 19682, 6);
    println!("track carpet rows={rows}");
}
