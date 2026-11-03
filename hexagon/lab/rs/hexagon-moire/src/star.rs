use crate::lattice::{Family, Rule, FAMILIES};
use crate::sums::{mean, odds};
use mrly_math::six::{geometry::cut, FILL, GRID};
use mrly_num::series::{beta, dirichlet, CATALAN, EULER};
use num_rational::Ratio;

pub struct Slice {
    pub rows: usize,
    pub cols: usize,
    pub types: Vec<u8>,
}

pub fn slice(family: Family, number: usize) -> Slice {
    let hex = cut(&family.cube(number)).expect("a cut cube");
    let (rows, cols) = (hex.height(), hex.width());
    let types = (0..rows * cols)
        .map(|at| hex.cell.types().at(at) as u8)
        .collect();
    Slice { rows, cols, types }
}

impl Slice {
    pub fn ink(&self) -> Ratio<i64> {
        let filled = self.types.iter().filter(|cell| **cell == FILL).count();
        let inside = self.types.iter().filter(|cell| **cell != GRID).count();
        Ratio::new(filled as i64, inside as i64)
    }

    fn coordinates(&self, number: usize) -> Vec<[f64; 3]> {
        let (n, cols) = (number as i64, self.cols as i64);
        let mut out = vec![[0.0; 3]; self.rows * self.cols];
        for row in 0..self.rows as i64 {
            let z = 2 * row;
            let target = 6 * n - 2 - z;
            let low = 0.max(target - (4 * n - 1));
            let reach = (4 * n - 1).min(target) - low + 1;
            let offset = (cols - reach) / 2;
            for step in 0..reach {
                let x = low + step;
                let scale = 4.0 * n as f64;
                out[(row * cols + offset + step) as usize] = [
                    x as f64 / scale,
                    (target - x) as f64 / scale,
                    z as f64 / scale,
                ];
            }
        }
        out
    }
}

pub fn law(family: Family, number: usize) -> Ratio<i64> {
    let n = number as i64;
    let chi = if (3 * n - 1) / 2 % 2 == 0 { 1 } else { -1 };
    let r = |p: i64, q: i64| Ratio::new(p, q);
    let carpet = r(1, 2) + r(chi, 8) + r(1, 2 * n) - r(chi, 8 * n * n);
    match family {
        Family::Carpet => carpet,
        Family::Net => r(1, 1) - carpet,
        Family::Tree => r(1, 4) + (r(1, 3) - r(chi, 12)) / r(n, 1) + r(1 - chi, 6 * n * n),
        Family::Void => r(1, 4) - r(chi, 4 * n) + r(1, 2 * n * n),
    }
}

pub fn ink_laws(limit: usize) {
    println!("cut ink of every odd n <= {limit} against the closed forms, exact rationals");
    for family in FAMILIES {
        let mut matched = 0;
        let mut centre = 0;
        let mut layers = 0;
        for number in odds(limit) {
            let hex = slice(family, number);
            matched += usize::from(hex.ink() == law(family, number));
            centre += usize::from(hex.types[number * hex.cols + 2 * number - 1] == FILL);
            layers += 1;
        }
        println!(
            "  {}: {matched}/{layers} layers match, centre cell ink in {centre}/{layers}",
            family.name()
        );
    }
}

const HEIGHT: usize = 1200;
const WIDTH: usize = 2399;
const DEEPEST: usize = 111;

fn sample(count: usize, extent: usize) -> Vec<usize> {
    (0..count)
        .map(|slot| ((slot as f64 + 0.5) / count as f64 * extent as f64).floor() as usize)
        .collect()
}

fn raster<T: Copy>(hex: &Slice, cells: &[T]) -> Vec<T> {
    let lines = sample(HEIGHT, hex.rows);
    let slots = sample(WIDTH, hex.cols);
    let mut out = Vec::with_capacity(HEIGHT * WIDTH);
    for line in &lines {
        for slot in &slots {
            out.push(cells[line * hex.cols + slot]);
        }
    }
    out
}

fn masked(values: &[f64], keep: &[bool]) -> f64 {
    let kept: Vec<f64> = values
        .iter()
        .zip(keep)
        .filter(|(_, inside)| **inside)
        .map(|(value, _)| *value)
        .collect();
    mean(&kept)
}

pub fn ghost() {
    println!("carpet cut layers rendered on a {HEIGHT} by {WIDTH} raster, odd n <= {DEEPEST}");
    let layers: Vec<(Vec<bool>, Vec<bool>)> = odds(DEEPEST)
        .map(|number| {
            let hex = slice(Family::Carpet, number);
            let fill = raster(
                &hex,
                &hex.types.iter().map(|c| *c == FILL).collect::<Vec<bool>>(),
            );
            let inside = raster(
                &hex,
                &hex.types.iter().map(|c| *c != GRID).collect::<Vec<bool>>(),
            );
            (fill, inside)
        })
        .collect();
    let mut hexagon = vec![true; HEIGHT * WIDTH];
    for (_, inside) in &layers {
        for (keep, ok) in hexagon.iter_mut().zip(inside) {
            *keep &= ok;
        }
    }
    let deepest = slice(Family::Carpet, DEEPEST);
    let gaps = raster(&deepest, &deepest.coordinates(DEEPEST));
    let near = |cell: usize, width: f64| {
        let [x, y, z] = gaps[cell];
        (x - y).abs() < width || (y - z).abs() < width || (z - x).abs() < width
    };
    let star: Vec<bool> = (0..HEIGHT * WIDTH)
        .map(|cell| hexagon[cell] && near(cell, 0.004))
        .collect();
    let background: Vec<bool> = (0..HEIGHT * WIDTH)
        .map(|cell| hexagon[cell] && !near(cell, 0.02))
        .collect();
    println!(
        "  pixels: hexagon {}  star {}  background {}",
        hexagon.iter().filter(|k| **k).count(),
        star.iter().filter(|k| **k).count(),
        background.iter().filter(|k| **k).count()
    );
    let mut total = vec![0.0f64; HEIGHT * WIDTH];
    for (count, (fill, _)) in layers.iter().enumerate() {
        for (slot, value) in total.iter_mut().zip(fill) {
            *slot += f64::from(*value);
        }
        let stacked = count + 1;
        if [5, 28, 56].contains(&stacked) {
            let field: Vec<f64> = total.iter().map(|value| value / stacked as f64).collect();
            let (star_ink, background_ink) = (masked(&field, &star), masked(&field, &background));
            println!(
                "  layers {stacked:2}: star {star_ink:.5}  background {background_ink:.5}  star minus background {:+.5}",
                star_ink - background_ink
            );
        }
    }
}

pub fn fading_lattice(rule: &Rule) {
    println!("carpet lattice frame, band |x-y| <= 0.01 of the cube side, per-layer excess over the exact ink law");
    let counts = [14usize, 28, 56, 100, 200, 400];
    let mut excesses = Vec::new();
    for number in odds(2 * counts[counts.len() - 1] - 1) {
        let n = number as i64;
        let size = 4 * n;
        let step = 1.0 / size as f64;
        let reach = 2.max((0.01 * size as f64).ceil() as i64 + 2);
        let (mut total, mut inked) = (0.0f64, 0.0f64);
        for index in 0..2 * n {
            let z = 2 * index;
            let target = 6 * n - 2 - z;
            for offset in -reach..=reach {
                let x = target.div_euclid(2) + offset;
                let y = target - x;
                if !(0..size).contains(&x) || !(0..size).contains(&y) {
                    continue;
                }
                let d = (x - y) as f64 * step;
                let weight = ((d + step).min(0.01) - (d - step).max(-0.01)).max(0.0);
                total += weight;
                inked += weight * f64::from(rule.filled(x, y, z));
            }
        }
        let background = law(Family::Carpet, number);
        let background = *background.numer() as f64 / *background.denom() as f64;
        excesses.push(if total > 0.0 {
            inked / total - background
        } else {
            0.0
        });
    }
    let mut scaled = Vec::new();
    for count in counts {
        let excess = mean(&excesses[..count]);
        println!(
            "  L = {count:3}: excess {excess:+.6}  excess * L {:+.4}",
            excess * count as f64
        );
        scaled.push(excess * count as f64);
    }
    println!(
        "  slope of excess * L against ln L: {:+.4} (100 to 200)  {:+.4} (200 to 400)",
        (scaled[4] - scaled[3]) / 2f64.ln(),
        (scaled[5] - scaled[4]) / 2f64.ln()
    );
}

fn law_value(number: usize) -> f64 {
    let value = law(Family::Carpet, number);
    *value.numer() as f64 / *value.denom() as f64
}

fn arm_ink(rule: &Rule, n: i64, half: i64) -> f64 {
    let size = 4 * n;
    let (mut total, mut inked) = (0i64, 0i64);
    for index in 0..2 * n {
        let z = 2 * index;
        let target = 6 * n - 2 - z;
        let base = target / 2;
        for x in base - half - 1..=base + half + 1 {
            let y = target - x;
            if !(0..size).contains(&x) || !(0..size).contains(&y) || (x - y).abs() > half {
                continue;
            }
            total += 1;
            inked += i64::from(rule.filled(x, y, z));
        }
    }
    inked as f64 / total as f64
}

fn arm_count(rule: &Rule, n: i64, axis: usize) -> Ratio<i64> {
    let size = 4 * n;
    let (mut total, mut inked) = (0i64, 0i64);
    for step in 0..size {
        let (x, y, z) = match axis {
            0 => (step, step, 6 * n - 2 - 2 * step),
            1 => (6 * n - 2 - 2 * step, step, step),
            _ => (step, 6 * n - 2 - 2 * step, step),
        };
        if !(0..size).contains(&x)
            || !(0..size).contains(&y)
            || !(0..size).contains(&z)
            || z % 2 != 0
        {
            continue;
        }
        total += 1;
        inked += i64::from(rule.filled(x, y, z));
    }
    Ratio::new(inked, total)
}

fn arm_law(n: i64) -> Ratio<i64> {
    let rhythm = [0, 1, 0, -1, 0, -1, 0, 1][(n % 8) as usize];
    Ratio::new(1, 2) + Ratio::new(rhythm, 2 * n)
}

fn nearest(value: f64) -> (i64, i64) {
    let mut best = (0i64, 1i64, f64::INFINITY);
    for denominator in 1..=400i64 {
        let numerator = (value * denominator as f64).round() as i64;
        let gap = (value - numerator as f64 / denominator as f64).abs();
        if gap < best.2 - 1e-12 {
            best = (numerator, denominator, gap);
        }
    }
    (best.0, best.1)
}

const ARM_LIMIT: usize = 6400;
const LAW_LIMIT: usize = 2001;
const ARMS_LIMIT: usize = 221;
const WIDTH_LIMIT: usize = 800;
const HELD_LIMIT: usize = 1600;
const CHECK_LIMIT: usize = 3200;

pub fn cell_frame(rule: &Rule) {
    println!(
        "carpet cell frame: star the exact arm x = y of the cut, background the exact ink law"
    );
    let (mut matched, mut layers) = (0, 0);
    for number in odds(LAW_LIMIT) {
        matched += usize::from(arm_count(rule, number as i64, 0) == arm_law(number as i64));
        layers += 1;
    }
    let (mut agreed, mut triples) = (0, 0);
    for number in odds(ARMS_LIMIT) {
        let n = number as i64;
        let inks = [
            arm_count(rule, n, 0),
            arm_count(rule, n, 1),
            arm_count(rule, n, 2),
        ];
        agreed += usize::from(inks[0] == inks[1] && inks[1] == inks[2]);
        triples += 1;
    }
    println!(
        "  arm ink 1/2 + chi_8(n)/(2n), chi_8 the mod-8 character of Q(sqrt 2), exact rationals: {matched}/{layers} layers match at odd n <= {LAW_LIMIT}, three arms agree in {agreed}/{triples} at odd n <= {ARMS_LIMIT}"
    );
    let excesses: Vec<f64> = odds(2 * ARM_LIMIT - 1)
        .map(|number| arm_ink(rule, number as i64, 0) - law_value(number))
        .collect();
    let constant = (1.0 + 2f64.sqrt()).ln() / (2.0 * 2f64.sqrt())
        - CATALAN / 8.0
        - EULER / 4.0
        - 2f64.ln() / 2.0;
    let mut scaled = Vec::new();
    let mut residuals = Vec::new();
    println!(
        "  L = 0 mod 4: the -chi/8 term sums to 0, residual is excess * L + (ln L)/4 - C, against -23/192 = {:+.8}",
        -23.0 / 192.0
    );
    for count in [28usize, 56, 100, 200, 400, 800, 1600, 3200, 6400] {
        let excess = mean(&excesses[..count]);
        let scale = excess * count as f64;
        let logged = scale + (count as f64).ln() / 4.0;
        println!(
            "  L = {count:4}: excess {excess:+.8}  excess * L {scale:+.6}  excess * L + (ln L)/4 {logged:+.10}  residual * L^2 {:+.8}",
            (logged - constant) * (count as f64) * (count as f64)
        );
        scaled.push(scale);
        residuals.push(logged - constant);
    }
    print!("  residual ratio per doubling from L = 100:");
    for step in 2..residuals.len() - 1 {
        print!(" {:.2}", residuals[step] / residuals[step + 1]);
    }
    println!();
    println!(
        "  L = 2 mod 4: the -chi/8 term still sums to 0 and the character tail changes sign, against +25/192 = {:+.8}",
        25.0 / 192.0
    );
    for count in [102usize, 202, 402, 802, 1602, 3202] {
        let excess = mean(&excesses[..count]);
        let logged = excess * count as f64 + (count as f64).ln() / 4.0;
        println!(
            "  L = {count:4}: excess * L + (ln L)/4 {logged:+.10}  residual * L^2 {:+.8}",
            (logged - constant) * (count as f64) * (count as f64)
        );
    }
    println!(
        "  L odd: the -chi/8 term sums to +1/8, so the constant is C + 1/8 and the error is O(1/L)"
    );
    for count in [27usize, 99, 401, 1601, 6399] {
        let excess = mean(&excesses[..count]);
        let logged = excess * count as f64 + (count as f64).ln() / 4.0;
        println!(
            "  L = {count:4}: excess * L + (ln L)/4 - C {:+.6}  less 1/8 {:+.8}  times L {:+.6}",
            logged - constant,
            logged - constant - 0.125,
            (logged - constant - 0.125) * count as f64
        );
    }
    println!(
        "  slope of excess * L against ln L: {:+.6} (1600 to 3200)  {:+.6} (3200 to 6400)",
        (scaled[7] - scaled[6]) / 2f64.ln(),
        (scaled[8] - scaled[7]) / 2f64.ln()
    );
    println!(
        "  constant ln(1+sqrt2)/(2 sqrt2) - G/8 - gamma/4 - (ln2)/2 = {constant:.10}  G from its own series {:.10}  L(1, chi_8) series {:.8} against ln(1+sqrt2)/sqrt2 {:.8}",
        beta(2.0, 2_000_000),
        dirichlet(1.0, &[0, 1, 0, -1, 0, -1, 0, 1], 20_000_000),
        (1.0 + 2f64.sqrt()).ln() / 2f64.sqrt()
    );
}

struct Width {
    k: i64,
    b: i64,
    kappa: Ratio<i64>,
    m: Ratio<i64>,
    q: Ratio<i64>,
}

fn edge_run(k: i64) -> i64 {
    let hits = (-k..=k)
        .filter(|j| matches!(j.rem_euclid(8), 3..=5))
        .count() as i64;
    hits + (k + 2).div_euclid(4) - k
}

fn width_terms(half: i64) -> Width {
    let k = half.div_euclid(2);
    let b = i64::from(k.div_euclid(2) % 2 == 0);
    let run = edge_run(k);
    Width {
        k,
        b,
        kappa: Ratio::new(if k % 2 == 0 { -1 } else { 1 }, 8 * (2 * k + 1)),
        m: Ratio::new(-(k + b), 2 * (2 * k + 1)),
        q: Ratio::new(1 - 2 * run, 2 * (2 * k + 1)),
    }
}

fn band_count(rule: &Rule, n: i64, half: i64) -> Ratio<i64> {
    let size = 4 * n;
    let (mut total, mut inked) = (0i64, 0i64);
    for index in 0..2 * n {
        let z = 2 * index;
        let target = 6 * n - 2 - z;
        let base = target / 2;
        for x in base - half - 1..=base + half + 1 {
            let y = target - x;
            if !(0..size).contains(&x) || !(0..size).contains(&y) || (x - y).abs() > half {
                continue;
            }
            total += 1;
            inked += i64::from(rule.filled(x, y, z));
        }
    }
    Ratio::new(inked, total)
}

fn width_excess_law(half: i64, n: i64) -> Ratio<i64> {
    let w = width_terms(half);
    let chi = if (3 * n - 1) / 2 % 2 == 0 { 1 } else { -1 };
    let eight = [0i64, 1, 0, -1, 0, -1, 0, 1][(n % 8) as usize];
    w.kappa * chi + (w.m + w.q * eight) / n + Ratio::new(chi, 8 * n * n)
}

fn flat(value: Ratio<i64>) -> f64 {
    *value.numer() as f64 / *value.denom() as f64
}

const IDENTITY_LIMIT: usize = 201;
const LADDER_LIMIT: usize = 1602;

pub fn cell_width_law(rule: &Rule) {
    println!("carpet cell frame, one closed form for every band half-width: excess_W(n) = kappa chi + (m + q chi_8(n))/n + chi/(8 n^2), exact rationals against the counted band");
    println!("  K = W/2 since x - y is even on the cut, b = 1 when floor(K/2) is even, and E(K) = #(|j| <= K, j = 3,4,5 mod 8) + floor((K+2)/4) - K is the tail's chi_8 weight, derived and not fitted");
    let periodic = (0..=200i64)
        .filter(|k| edge_run(*k) == [0i64, -1, -1, 0, 1, 2, 2, 1][(k % 8) as usize])
        .count();
    println!("  E(K) against the 8-periodic run 0 -1 -1 0 1 2 2 1, which a block of eight leaves alone since 6 + 2 - 8 = 0: {periodic}/201 at K = 0..200");
    println!("  the match is broken out by n mod 8 so no residue class hides; odd W repeats its even neighbour, so only the first W at each K is counted in the total");
    let (mut seen, mut kept, mut total) = (Vec::new(), 0usize, 0usize);
    for half in [
        0i64, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 24, 32, 36, 64,
    ] {
        let w = width_terms(half);
        let fresh = !seen.contains(&w.k);
        if fresh {
            seen.push(w.k);
        }
        let (mut hit, mut all) = ([0usize; 4], [0usize; 4]);
        for number in odds(IDENTITY_LIMIT) {
            let n = number as i64;
            if n < w.k {
                continue;
            }
            let slot = ((n % 8) / 2) as usize;
            all[slot] += 1;
            hit[slot] += usize::from(
                band_count(rule, n, half) - law(Family::Carpet, number)
                    == width_excess_law(half, n),
            );
        }
        if fresh {
            kept += hit.iter().sum::<usize>();
            total += all.iter().sum::<usize>();
        }
        println!(
            "  W = {half:2}: K {:2} b {} kappa {} m {} q {} slope m/2 {} = {:+.8}  n = 1,3,5,7 mod 8 match {}/{} {}/{} {}/{} {}/{}",
            w.k,
            w.b,
            w.kappa,
            w.m,
            w.q,
            w.m / 2,
            flat(w.m) / 2.0,
            hit[0], all[0], hit[1], all[1], hit[2], all[2], hit[3], all[3]
        );
    }
    println!("  {} distinct half-widths, {kept}/{total} distinct checks; the nine odd rows above repeat their even neighbours and are not counted twice", seen.len());
}

pub fn cell_width_ladders(rule: &Rule) {
    println!("carpet cell frame, the width family summed: L * excess_L = (m/2) ln L + C_W + O(1/L^2) at even L, C_W = m (ln 2 + gamma/2) + q L(1, chi_8) - G/8 + Delta_W");
    println!("  Delta_W is the exact contribution of the clipped layers n < K, the 1/L^2 coefficient is -q/4 + 1/64 + m/48 at L = 0 mod 4 and +q/4 + 1/64 + m/48 at L = 2 mod 4, and the odd-L constant is C_W - kappa");
    let lvalue = (1.0 + 2f64.sqrt()).ln() / 2f64.sqrt();
    for half in [0i64, 2, 4, 6, 8, 12, 16] {
        let w = width_terms(half);
        let excesses = width_excesses(rule, half, LADDER_LIMIT);
        let mut drift = 0.0f64;
        for number in odds(w.k as usize) {
            if (number as i64) < w.k {
                drift += excesses[(number - 1) / 2] - flat(width_excess_law(half, number as i64));
            }
        }
        let (m, q, kappa) = (flat(w.m), flat(w.q), flat(w.kappa));
        let slope = m / 2.0;
        let constant = m * (2f64.ln() + EULER / 2.0) + q * lvalue - CATALAN / 8.0 + drift;
        let miss = |count: usize| {
            mean(&excesses[..count]) * count as f64 - slope * (count as f64).ln() - constant
        };
        let window = |count: usize| {
            (mean(&excesses[..count]) * count as f64
                - mean(&excesses[..count / 2]) * (count / 2) as f64)
                / 2f64.ln()
        };
        println!(
            "  W = {half:2}: slope {} = {slope:+.8}  C_W {constant:+.10}  Delta_W {drift:+.8}",
            w.m / 2
        );
        println!(
            "    L = 0 mod 4: residual * L^2 {:+.7} {:+.7} {:+.7} at L = 400, 800, 1600 against -q/4 + 1/64 + m/48 = {:+.7}",
            miss(400) * 160000.0,
            miss(800) * 640000.0,
            miss(1600) * 2560000.0,
            -q / 4.0 + 1.0 / 64.0 + m / 48.0
        );
        println!(
            "    L = 2 mod 4: residual * L^2 {:+.7} {:+.7} {:+.7} at L = 402, 802, 1602 against +q/4 + 1/64 + m/48 = {:+.7}",
            miss(402) * 402.0 * 402.0,
            miss(802) * 802.0 * 802.0,
            miss(1602) * 1602.0 * 1602.0,
            q / 4.0 + 1.0 / 64.0 + m / 48.0
        );
        println!(
            "    L odd: (residual + kappa) * L {:+.6} {:+.6} at L = 401, 1601, so the odd-L constant is C_W - kappa with an O(1/L) error",
            (miss(401) + kappa) * 401.0,
            (miss(1601) + kappa) * 1601.0
        );
        println!(
            "    sliding window L/2 to L: {:+.8} at L = 1600 against m/2 {:+.8}, {:+.8} at L = 1602 against m/2 + kappa/ln 2 {:+.8}",
            window(1600),
            slope,
            window(1602),
            slope + kappa / 2f64.ln()
        );
    }
}

fn width_excesses(rule: &Rule, half: i64, limit: usize) -> Vec<f64> {
    odds(2 * limit - 1)
        .map(|number| arm_ink(rule, number as i64, half) - law_value(number))
        .collect()
}

fn width_slope(excesses: &[f64], limit: usize) -> f64 {
    let near = mean(&excesses[..limit / 2]) * (limit / 2) as f64;
    let far = mean(&excesses[..limit]) * limit as f64;
    (far - near) / 2f64.ln()
}

fn edge_law(half: i64) -> f64 {
    let w = width_terms(half);
    -(w.k as f64 + w.b as f64) / (4.0 * (2 * w.k + 1) as f64)
}

pub fn cell_widths(rule: &Rule) {
    println!("carpet cell frame, band half-width W cells about the arm, slope of excess * L against ln L, K = W/2 and b = 1 when floor(K/2) is even");
    println!("  the twelve swept widths against the derived rational, with a best-rational search over denominators to 400 beside them as a blind reading");
    for half in [0i64, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32, 64] {
        let slope = width_slope(&width_excesses(rule, half, WIDTH_LIMIT), WIDTH_LIMIT);
        let (numerator, denominator) = nearest(slope);
        println!(
            "  W = {half:2} cells at L = {WIDTH_LIMIT}: slope {slope:+.6}  best rational {numerator}/{denominator} = {:+.6}  search residual {:.1e}  against -(K + b)/(4(2K + 1)) = {:+.6}",
            numerator as f64 / denominator as f64,
            (slope - numerator as f64 / denominator as f64).abs(),
            edge_law(half)
        );
    }
    println!("  odd W is not a new band: x - y is even on the arm, so W and W - 1 read one point set and one rational");
    for half in [0i64, 1, 2, 3] {
        let slope = width_slope(&width_excesses(rule, half, WIDTH_LIMIT), WIDTH_LIMIT);
        println!(
            "  W = {half:2} cells at L = {WIDTH_LIMIT}: slope {slope:+.6}  against -(K + b)/(4(2K + 1)) = {:+.6}",
            edge_law(half)
        );
    }
    println!("  the seven smallest widths recomputed to L = {CHECK_LIMIT}, four times the sweep, against the exact rational");
    for half in [0i64, 2, 4, 6, 8, 10, 12] {
        let excesses = width_excesses(rule, half, CHECK_LIMIT);
        let target = edge_law(half);
        let gaps: Vec<f64> = [WIDTH_LIMIT, HELD_LIMIT, CHECK_LIMIT]
            .iter()
            .map(|limit| (width_slope(&excesses, *limit) - target).abs())
            .collect();
        println!(
            "  W = {half:2} cells: slope {:+.8} at L = {CHECK_LIMIT} against {target:+.8}  gap {:.1e} (L = {WIDTH_LIMIT})  {:.1e} (L = {HELD_LIMIT})  {:.1e} (L = {CHECK_LIMIT})",
            width_slope(&excesses, CHECK_LIMIT),
            gaps[0],
            gaps[1],
            gaps[2]
        );
    }
    println!(
        "  four widths held out of the sweep, predicted before measuring, at L = {HELD_LIMIT}"
    );
    for half in [14i64, 18, 28, 36] {
        let slope = width_slope(&width_excesses(rule, half, HELD_LIMIT), HELD_LIMIT);
        let target = edge_law(half);
        println!(
            "  W = {half:2} cells: predicted {target:+.8}  measured {slope:+.8}  gap {:.1e}",
            (slope - target).abs()
        );
    }
    println!(
        "  the width law's own limit as W grows is -1/8 = {:+.6}, the value the fixed-fraction lattice band measures",
        -0.125
    );
}
