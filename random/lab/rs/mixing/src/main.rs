mod catalog;
mod cone;
mod kernel;
mod triple;

use catalog::{cells, cycles, four, three, word, Round, Spec, NEIGHBOURHOODS, SOURCES};
use kernel::{radius, rate, sigma2, tv, Kernel};
use std::env;
use std::time::Instant;

// CHECK

fn check() {
    let spec = Spec {
        n: 4,
        fill: word(7, 7),
        hood: 0,
        source: 0,
        spread: 3,
    };
    let round = spec.round();
    let r: Vec<u32> = round.table().iter().map(|y| y ^ spec.fill).collect();
    let lengths = cycles(&r);
    println!(
        "  {}: cycles {}, longest {}, fixed {}",
        spec.name(),
        lengths.len(),
        lengths.iter().max().unwrap(),
        lengths.iter().filter(|&&l| l == 1).count()
    );
    for (n, masks, want) in [(3, three(), 4320), (4, four(), 2856)] {
        let all = cells(n, &masks);
        let nonlinear: Vec<&Spec> = all.iter().filter(|s| !s.linear()).collect();
        let mean = nonlinear
            .iter()
            .map(|s| {
                let round = s.round();
                let r: Vec<u32> = round.table().iter().map(|y| y ^ s.fill).collect();
                cycles(&r).len() as f64
            })
            .sum::<f64>()
            / nonlinear.len() as f64;
        println!(
            "  n {n}: {} cells (want {want}), {} nonlinear, mean cycles {mean:.2}",
            all.len(),
            nonlinear.len()
        );
    }
}

// THREE

struct Curve {
    d: Vec<f64>,
    worst: Vec<usize>,
}

fn curve(k: &Kernel, starts: &[usize], eps: f64, most: usize) -> Curve {
    let m = k.m;
    let mut rows: Vec<Vec<f64>> = starts
        .iter()
        .map(|&s| {
            let mut v = vec![0.0; m];
            v[s] = 1.0;
            v
        })
        .collect();
    let mut buf = vec![0.0; m];
    let mut d = vec![1.0 - 1.0 / (m - 1) as f64];
    let mut worst = vec![starts[0]];
    while d.len() < most && *d.last().unwrap() > eps {
        let (mut top, mut at) = (0.0f64, 0);
        for (i, row) in rows.iter_mut().enumerate() {
            k.forward(row, &mut buf);
            std::mem::swap(row, &mut buf);
            let t = tv(row);
            if t > top {
                top = t;
                at = starts[i];
            }
        }
        d.push(top);
        worst.push(at);
    }
    Curve { d, worst }
}

fn sigma_of_power(k: &Kernel, j: usize) -> kernel::Spectrum {
    let pp = |v: &[f64], out: &mut [f64]| {
        let mut a = v.to_vec();
        let mut b = vec![0.0; v.len()];
        for _ in 0..j {
            k.forward(&a, &mut b);
            std::mem::swap(&mut a, &mut b);
        }
        for _ in 0..j {
            k.backward(&a, &mut b);
            std::mem::swap(&mut a, &mut b);
        }
        out.copy_from_slice(&a);
    };
    sigma2(&pp, k.m, 0x6d69_7869_6e67, 150)
}

fn lambda2(apply: &dyn Fn(&[f64], &mut [f64]), m: usize, seed: u64) -> kernel::Eigen {
    radius(apply, m, seed, 80, 1e-10, 40)
}

fn first_below(d: &[f64], eps: f64) -> Option<usize> {
    d.iter().position(|&x| x < eps).map(|i| i + 1)
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn census_three() {
    let all = cells(3, &three());
    let mut faults = 0;
    let mut irreducible = [0usize; 2];
    let mut floors: Vec<f64> = Vec::new();
    let mut sig = Vec::new();
    let mut lam = Vec::new();
    let mut quarter = Vec::new();
    let mut tight = Vec::new();
    let mut kstars = Vec::new();
    let mut drift: f64 = 0.0;
    let mut eff = Vec::new();
    for spec in &all {
        let round = spec.round();
        let h = round.table();
        let k = Kernel::new(&round);
        for a in 1..k.m {
            if k.direct(&h, a) != k.sparse(a) {
                faults += 1;
            }
        }
        let sizes = k.classes();
        let lin = spec.linear() as usize;
        if sizes.len() == 1 {
            irreducible[lin] += 1;
        } else if lin == 0 {
            floors.push(1.0 - *sizes.iter().min().unwrap() as f64 / (k.m - 1) as f64);
        }
        if sizes.len() > 1 || lin == 1 {
            if sizes.len() == 1 {
                let all: Vec<usize> = (1..k.m).collect();
                let c = curve(&k, &all, 1e-9, 12);
                println!(
                    "  affine irreducible: {}: d {:?}",
                    spec.name(),
                    &c.d[..4.min(c.d.len())]
                );
            }
            continue;
        }
        let mut powers: Vec<f64> = Vec::new();
        while powers.len() < 24 && powers.last().map_or(true, |&x| x > 1.0 - 1e-9) {
            powers.push(sigma_of_power(&k, powers.len() + 1).sigma);
        }
        let kstar = powers
            .iter()
            .position(|&x| x < 1.0 - 1e-9)
            .map_or(0, |i| i + 1);
        let s = sigma_of_power(&k, kstar.max(1));
        let fw = |v: &[f64], out: &mut [f64]| k.forward(v, out);
        let l = rate(&fw, k.m, 0x6d69_7869_6e67, 400, 200);
        let short = rate(&fw, k.m, 0x6d69_7869_6e67, 40, 40);
        drift = drift.max((short - l).abs());
        let all: Vec<usize> = (1..k.m).collect();
        let c = curve(&k, &all, 1e-6, 400);
        let q = first_below(&c.d, 0.25).unwrap_or(0);
        let t = first_below(&c.d, 0.01).unwrap_or(0);
        let light = c.worst[1..]
            .iter()
            .filter(|&&a| (a as u32 & spec.fill) == a as u32 || (a as u32).count_ones() == 1)
            .count();
        println!(
            "  {}: sigma2(P^j) {:.6?}, k* {kstar} (res {:.1e}, {} steps), |lambda2| {:.6}, rounds to 1/4 {q}, to 0.01 {t}, to 1e-6 {:?}, worst start fill-only or one cell {light} of {}, d {:.4?}",
            spec.name(),
            powers,
            s.residual,
            s.steps,
            l,
            first_below(&c.d, 1e-6),
            c.d.len() - 1,
            &c.d[..8.min(c.d.len())]
        );
        kstars.push(kstar as f64);
        eff.push(if kstar > 0 {
            s.sigma.powf(1.0 / kstar as f64)
        } else {
            1.0
        });
        sig.push(s.sigma);
        lam.push(l);
        quarter.push(q as f64);
        tight.push(t as f64);
    }
    println!("  |lambda2| from 40 + 40 power steps against 400 + 200: largest gap {drift:.2e}");
    println!(
        "  n 3: {} cells, {faults} kernel rows that differ from the direct table",
        all.len()
    );
    println!(
        "  irreducible: nonlinear {}, affine {}",
        irreducible[0], irreducible[1]
    );
    let lo = floors.iter().cloned().fold(1.0, f64::min);
    println!(
        "  reducible nonlinear {}: floor 1 - min class/511, smallest floor {lo:.6}",
        floors.len()
    );
    if !sig.is_empty() {
        let mn = |v: &[f64]| v.iter().cloned().fold(f64::MAX, f64::min);
        let mx = |v: &[f64]| v.iter().cloned().fold(f64::MIN, f64::max);
        println!(
            "  sigma2 min {:.6} median {:.6} max {:.6}; |lambda2| min {:.6} median {:.6} max {:.6}",
            mn(&sig),
            median(&mut sig.clone()),
            mx(&sig),
            mn(&lam),
            median(&mut lam.clone()),
            mx(&lam)
        );
        println!(
            "  k*: min {} median {} max {}; sigma2(P^k*)^(1/k*) min {:.6} median {:.6} max {:.6}",
            mn(&kstars),
            median(&mut kstars.clone()),
            mx(&kstars),
            mn(&eff),
            median(&mut eff.clone()),
            mx(&eff)
        );
        println!(
            "  rounds to 1/4: min {} median {} max {}; to 0.01: min {} median {} max {}",
            mn(&quarter),
            median(&mut quarter.clone()),
            mx(&quarter),
            mn(&tight),
            median(&mut tight.clone()),
            mx(&tight)
        );
    }
}

// TRIPLE

fn census_triple(every: usize) {
    let all = cells(3, &three());
    let mut picked = Vec::new();
    for spec in all.iter().filter(|s| !s.linear()) {
        let k = Kernel::new(&spec.round());
        if k.classes().len() == 1 {
            picked.push(spec);
        }
    }
    for spec in picked.iter().step_by(every) {
        let round = spec.round();
        let k = Kernel::new(&round);
        let t = triple::Triple::new(&round);
        let f3 = |v: &[f64], out: &mut [f64]| t.forward(v, out);
        let mut powers: Vec<f64> = Vec::new();
        while powers.len() < 8 && powers.last().map_or(true, |&x| x > 1.0 - 1e-9) {
            let j = powers.len() + 1;
            let pp = |v: &[f64], out: &mut [f64]| {
                let mut a = v.to_vec();
                let mut b = vec![0.0; v.len()];
                for _ in 0..j {
                    t.forward(&a, &mut b);
                    std::mem::swap(&mut a, &mut b);
                }
                for _ in 0..j {
                    t.backward(&a, &mut b);
                    std::mem::swap(&mut a, &mut b);
                }
                out.copy_from_slice(&a);
            };
            powers.push(sigma2(&pp, t.size, 0x7472_6970_6c65, 60).sigma);
        }
        let l3 = rate(&f3, t.size, 0x7472_6970_6c65, 150, 100);
        let fw = |v: &[f64], out: &mut [f64]| k.forward(v, out);
        let l2 = rate(&fw, k.m, 0x6d69_7869_6e67, 400, 200);
        let p0 = spec.fill.trailing_zeros() as usize;
        let starts: Vec<usize> = (0..9).filter(|&q| q != p0).collect();
        let mut rows: Vec<Vec<f64>> = starts
            .iter()
            .map(|&q| {
                let mut v = vec![0.0; t.size];
                v[t.index(1 << p0, 1 << q)] = 1.0;
                v
            })
            .collect();
        let two = curve(&k, &[1 << p0], 0.0, 40);
        let mut buf = vec![0.0; t.size];
        let mut d = vec![1.0 - 1.0 / (t.size - 1) as f64];
        let mut lump: f64 = 0.0;
        while d.len() < 120 && *d.last().unwrap() > 0.01 {
            let mut worst: f64 = 0.0;
            for row in rows.iter_mut() {
                t.forward(row, &mut buf);
                std::mem::swap(row, &mut buf);
                worst = worst.max(tv(row));
                if d.len() < two.d.len() {
                    let mg = t.marginal(row);
                    lump = lump.max((tv(&mg) - two.d[d.len()]).abs());
                }
            }
            d.push(worst);
        }
        println!(
            "  {}: 3-wise sigma2(P3^j) {:.6?}, |lambda2(P3)| {:.6}, |lambda2(P)| {:.6}, from 8 starts rounds to 1/4 {:?}, to 0.01 {:?}, marginal gap {lump:.1e}, d3 {:.4?}",
            spec.name(),
            powers,
            l3,
            l2,
            first_below(&d, 0.25),
            first_below(&d, 0.01),
            &d[..8.min(d.len())]
        );
    }
}

// STAMPS

fn law_gap(
    rows: &std::collections::BTreeMap<(u32, u32), u64>,
    against: &dyn Fn(u32) -> Vec<(u32, f64)>,
) -> f64 {
    let mut by: std::collections::BTreeMap<u32, Vec<(u32, u64)>> =
        std::collections::BTreeMap::new();
    for (&(a, b), &c) in rows {
        by.entry(a).or_default().push((b, c));
    }
    let mut worst: f64 = 0.0;
    for (a, row) in by {
        let total: u64 = row.iter().map(|r| r.1).sum();
        let mut diff: std::collections::BTreeMap<u32, f64> = std::collections::BTreeMap::new();
        for (b, c) in row {
            *diff.entry(b).or_default() += c as f64 / total as f64;
        }
        for (b, q) in against(a) {
            *diff.entry(b).or_default() -= q;
        }
        worst = worst.max(0.5 * diff.values().map(|x| x.abs()).sum::<f64>());
    }
    worst
}

fn stamps() {
    let spec = cells(3, &three())
        .into_iter()
        .find(|s| !s.linear() && Kernel::new(&s.round()).classes().len() == 1)
        .unwrap();
    let round = spec.round();
    let h = round.table();
    let k = Kernel::new(&round);
    let ddt = |a: u32| -> Vec<(u32, f64)> {
        k.sparse(a as usize)
            .into_iter()
            .map(|(b, c)| (b, c as f64 / k.m as f64))
            .collect()
    };
    let designs: Vec<u32> = NEIGHBOURHOODS.iter().map(|d| d.1).collect();
    let uniform: Vec<u32> = (0..512).collect();
    let void = 1u32 << round.void().trailing_zeros();
    let every: Vec<u32> = (1..512).collect();
    for (name, set, starts) in [
        ("uniform", &uniform, vec![void]),
        ("designs", &designs, vec![void]),
        ("designs", &designs, every),
    ] {
        let mut one: std::collections::BTreeMap<(u32, u32), u64> =
            std::collections::BTreeMap::new();
        let mut two: std::collections::BTreeMap<(u32, u32), u64> =
            std::collections::BTreeMap::new();
        let mut three_: std::collections::BTreeMap<(u32, u32, u32), u64> =
            std::collections::BTreeMap::new();
        let d0 = if starts.len() == 1 { starts[0] } else { 0 };
        for x0 in 0..512u32 {
            for &s0 in &starts {
                let d1 = h[x0 as usize] ^ h[(x0 ^ s0) as usize];
                for &c1 in set.iter() {
                    let x1 = h[x0 as usize] ^ c1;
                    let d2 = h[x1 as usize] ^ h[(x1 ^ d1) as usize];
                    *one.entry((d1, d2)).or_default() += 1;
                    if name == "uniform" {
                        continue;
                    }
                    for &c2 in set.iter() {
                        let x2 = h[x1 as usize] ^ c2;
                        let d3 = h[x2 as usize] ^ h[(x2 ^ d2) as usize];
                        *two.entry((d2, d3)).or_default() += 1;
                        *three_.entry((d1, d2, d3)).or_default() += 1;
                    }
                }
            }
        }
        let g1 = law_gap(&one, &ddt);
        if name == "uniform" {
            println!("  {}: stamps {name}: start {d0:#05x}, largest gap between the law of d2 given d1 and the kernel {g1:.3e}", spec.name());
            continue;
        }
        let mut pairs: std::collections::BTreeMap<u32, Vec<(u32, f64)>> =
            std::collections::BTreeMap::new();
        let mut totals: std::collections::BTreeMap<u32, u64> = std::collections::BTreeMap::new();
        for (&(b, _), &c) in &two {
            *totals.entry(b).or_default() += c;
        }
        for (&(b, d), &c) in &two {
            pairs
                .entry(b)
                .or_default()
                .push((d, c as f64 / totals[&b] as f64));
        }
        let mut hist: std::collections::BTreeMap<(u32, u32), u64> =
            std::collections::BTreeMap::new();
        let mut keyed: std::collections::BTreeMap<(u32, u32), Vec<(u32, u64)>> =
            std::collections::BTreeMap::new();
        for (&(a, b, d), &c) in &three_ {
            keyed.entry((a, b)).or_default().push((d, c));
            *hist.entry((a, b)).or_default() += c;
        }
        let mut worst: (f64, u32, u32) = (0.0, 0, 0);
        for ((a, b), row) in keyed {
            let total = hist[&(a, b)] as f64;
            let mut diff: std::collections::BTreeMap<u32, f64> = std::collections::BTreeMap::new();
            for (d, c) in row {
                *diff.entry(d).or_default() += c as f64 / total;
            }
            for &(d, q) in &pairs[&b] {
                *diff.entry(d).or_default() -= q;
            }
            let g = 0.5 * diff.values().map(|x| x.abs()).sum::<f64>();
            if g > worst.0 {
                worst = (g, a, b);
            }
        }
        println!(
            "  {}: stamps {name} {:?}: start {d0:#05x}, gap of d2 given d1 from the kernel {g1:.4}, gap of d3 given (d1, d2) from d3 given d2 {:.4} at d1 {:#05x}, d2 {:#05x}",
            spec.name(),
            set,
            worst.0,
            worst.1,
            worst.2
        );
    }
}

// CONE

fn sched_curve(ks: &[&Kernel], starts: &[usize], eps: f64, most: usize) -> Vec<f64> {
    let m = ks[0].m;
    let mut rows: Vec<Vec<f64>> = starts
        .iter()
        .map(|&s| {
            let mut v = vec![0.0; m];
            v[s] = 1.0;
            v
        })
        .collect();
    let mut buf = vec![0.0; m];
    let mut d = vec![1.0 - 1.0 / (m - 1) as f64];
    while d.len() < most && *d.last().unwrap() > eps {
        let k = ks[d.len() % ks.len()];
        let mut worst: f64 = 0.0;
        for row in rows.iter_mut() {
            k.forward(row, &mut buf);
            std::mem::swap(row, &mut buf);
            worst = worst.max(tv(row));
        }
        d.push(worst);
    }
    d
}

const FAMILY: [(&str, &[&str]); 2] = [
    (
        "carpet",
        &["slide", "scales", "row", "shear", "rotate", "cat"],
    ),
    (
        "board",
        &["slide", "shift", "shear", "rotate", "cat", "hyper"],
    ),
];

fn light() {
    for (design, spreads) in FAMILY {
        for level in 2..=6 {
            let g = cone::Grid::split(level, design, 0);
            let mut row = Vec::new();
            for &sch in spreads {
                let t = g
                    .cover(sch, 4 * g.n)
                    .iter()
                    .map(|x| x.unwrap_or(usize::MAX))
                    .max()
                    .unwrap();
                let hidden = if sch == "scales" { 0 } else { g.hidden(sch) };
                let stuck = if t < usize::MAX { 0 } else { g.stuck(sch) };
                let shown = if t < usize::MAX {
                    t.to_string()
                } else if stuck > 0 {
                    format!("never, {stuck} starts in a closed set")
                } else {
                    format!("over {}", 4 * g.n)
                };
                row.push(format!("{sch} {shown} ({hidden})"));
            }
            println!(
                "  cone {design} L {level} n {}, fill {} cells, every one of {} starts: rounds to cover the torus (spread orbits inside the fill) {}",
                g.n,
                g.reads.len(),
                g.cells,
                row.join(", ")
            );
        }
    }
    assert_eq!(cone::Grid::new(2).mask(), word(7, 7));
    assert_eq!(cone::Grid::new(2).rule, catalog::rule_of(SOURCES[0].1));
    let mut runs = vec![("carpet", "slide", 0), ("carpet", "scales", 0)];
    for source in 0..3 {
        for sch in ["shift", "shear", "hyper"] {
            runs.push(("board", sch, source));
        }
    }
    for (design, sch, source) in runs {
        let rule = catalog::rule_of(SOURCES[source].1);
        let g = cone::Grid::split(2, design, rule);
        let mask = g.mask();
        let make = |t: usize| {
            Kernel::new(&Round::with_perm(
                4,
                mask,
                NEIGHBOURHOODS[0].1,
                rule,
                cone::perm(2, sch, t),
            ))
        };
        let ks: Vec<Kernel> = if sch == "scales" {
            vec![make(0), make(1)]
        } else {
            vec![make(0)]
        };
        let refs: Vec<&Kernel> = ks.iter().collect();
        let starts: Vec<usize> = (0..16)
            .map(|p| 1usize << p)
            .chain([mask as usize, (!mask & 0xffff) as usize])
            .collect();
        let irreducible: Vec<bool> = ks.iter().map(|k| k.irreducible()).collect();
        let d = sched_curve(&refs, &starts, 0.01, 200);
        let e = if ks.len() == 1 {
            let fw = |v: &[f64], out: &mut [f64]| ks[0].forward(v, out);
            lambda2(&fw, ks[0].m, 0x6d69_7869_6e67)
        } else {
            let fw = |v: &[f64], out: &mut [f64]| {
                let mut t = vec![0.0; v.len()];
                ks[1].forward(v, &mut t);
                ks[0].forward(&t, out);
            };
            lambda2(&fw, ks[0].m, 0x6d69_7869_6e67)
        };
        let l = e.modulus.powf(1.0 / ks.len() as f64);
        println!(
            "  exact n 4 {design} {} {sch}: fill {mask:#06x}, irreducible {irreducible:?}, |lambda2| per round {l:.6} (residual {:.1e}), from 18 starts rounds to 1/2 {:?}, to 1/4 {:?}, to 0.05 {:?}, to 0.01 {:?}, d {:.4?}",
            SOURCES[source].0,
            e.residual,
            first_below(&d, 0.5),
            first_below(&d, 0.25),
            first_below(&d, 0.05),
            first_below(&d, 0.01),
            &d[..12.min(d.len())]
        );
    }
}

fn multi(runs: usize) {
    for level in 2..=4 {
        for (design, sch, source) in [
            ("carpet", "slide", 0),
            ("board", "shift", 1),
            ("board", "hyper", 1),
        ] {
            let g = cone::Grid::split(level, design, catalog::rule_of(SOURCES[source].1));
            let times: Vec<usize> = g
                .cover(sch, 4 * g.n)
                .iter()
                .map(|x| x.unwrap_or(usize::MAX))
                .collect();
            let worst = (0..g.cells)
                .max_by_key(|&p| (times[p], usize::MAX - p))
                .unwrap();
            let s = cone::sample(
                &g,
                sch,
                worst,
                runs,
                8 * g.n,
                0x6d75_6c74_6900 + level as u64,
            );
            let first = |tau: f64| s.bias.iter().position(|&b| b < tau).map(|i| i + 1);
            let pairs = (runs * (runs - 1) / 2) as f64;
            let expect = pairs / ((1u128 << g.cells.min(127)) as f64 - 1.0);
            let quiet = s
                .collisions
                .iter()
                .position(|&c| (c as f64) < 2.0 * expect + 10.0)
                .map(|i| i + 1);
            println!(
                "  sampled {design} {sch} L {level} n {}: start {worst}, {runs} runs, rounds to max bit bias below 0.2 {:?}, 0.1 {:?}, 0.05 {:?}; rounds to collisions below twice uniform plus 10 {:?}; bias {:.3?}",
                g.n,
                first(0.2),
                first(0.1),
                first(0.05),
                quiet,
                &s.bias[..14.min(s.bias.len())]
            );
        }
    }
}

// FOUR

fn census_four() {
    let all = cells(4, &four());
    let (mut closed, mut reducible, mut irreducible) = (0, 0, Vec::new());
    for spec in all.iter().filter(|s| !s.linear()) {
        let round = spec.round();
        if round.closed() {
            closed += 1;
            continue;
        }
        let k = Kernel::new(&round);
        if k.irreducible() {
            println!("  irreducible {}: {} edges", spec.name(), k.edges());
            irreducible.push(spec);
        } else {
            reducible += 1;
        }
    }
    println!(
        "  n 4 nonlinear: read-closed {closed}, other reducible {reducible}, irreducible {}",
        irreducible.len()
    );
    let mut rows = Vec::new();
    let (mut worst, mut gap): (f64, f64) = (0.0, 0.0);
    for (i, spec) in irreducible.iter().enumerate() {
        let figure =
            spec.fill == word(7, 7) && spec.hood == 0 && spec.source == 0 && spec.spread == 3;
        if i % 9 != 0 && !figure {
            continue;
        }
        let round = spec.round();
        let k = Kernel::new(&round);
        let h = round.table();
        let mut rng = 0x6d69_7869_6e67u64;
        for _ in 0..32 {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            let a = 1 + (rng % (k.m as u64 - 1)) as usize;
            assert_eq!(
                k.direct(&h, a),
                k.sparse(a),
                "a 4 x 4 kernel row that differs from the direct table"
            );
        }
        let mut powers: Vec<f64> = Vec::new();
        while powers.len() < 8 && powers.last().map_or(true, |&x| x > 1.0 - 1e-9) {
            powers.push(sigma_of_power(&k, powers.len() + 1).sigma);
        }
        let kstar = powers
            .iter()
            .position(|&x| x < 1.0 - 1e-9)
            .map_or(0, |i| i + 1);
        let fw = |v: &[f64], out: &mut [f64]| k.forward(v, out);
        let e = lambda2(&fw, k.m, 0x6d69_7869_6e67);
        let again = lambda2(&fw, k.m, 0x7365_636f_6e64);
        gap = gap.max((e.modulus - again.modulus).abs());
        let starts: Vec<usize> = (0..16)
            .map(|p| 1usize << p)
            .chain([spec.fill as usize, (!spec.fill & 0xffff) as usize])
            .collect();
        let c = curve(&k, &starts, 0.01, 120);
        let q = first_below(&c.d, 0.25).unwrap_or(0);
        let t = first_below(&c.d, 0.01).unwrap_or(0);
        println!(
            "  {}: |A| {}, sigma2(P^j) {:.6?}, k* {kstar}, |lambda2| {:.6} (lambda2 {:.6} {:+.6}i, residual {:.1e}, {} matvecs, {} restarts), from 18 starts rounds to 1/4 {q}, to 0.01 {t}, d {:.4?}",
            spec.name(),
            spec.fill.count_ones(),
            powers,
            e.modulus,
            e.value.0,
            e.value.1,
            e.residual,
            e.matvecs,
            e.restarts,
            &c.d[..10.min(c.d.len())]
        );
        worst = worst.max(e.residual).max(again.residual);
        rows.push((kstar, *powers.last().unwrap(), e.modulus, q, t));
    }
    let col = |f: &dyn Fn(&(usize, f64, f64, usize, usize)) -> f64| {
        let mut v: Vec<f64> = rows.iter().map(f).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        (v[0], v[v.len() / 2], v[v.len() - 1])
    };
    println!(
        "  sampled {} of {} irreducible cells (every 9th and the figure cell)",
        rows.len(),
        irreducible.len()
    );
    println!("  k* {:?}", col(&|r| r.0 as f64));
    println!("  sigma2(P^k*) {:.6?}", col(&|r| r.1));
    println!(
        "  |lambda2| {:.6?}, largest Arnoldi residual {worst:.1e}, a second seed within {gap:.1e}",
        col(&|r| r.2)
    );
    println!(
        "  rounds to 1/4 {:?}, to 0.01 {:?}",
        col(&|r| r.3 as f64),
        col(&|r| r.4 as f64)
    );
}

fn main() {
    let verb = env::args().nth(1).unwrap_or_else(|| "check".to_string());
    let start = Instant::now();
    match verb.as_str() {
        "check" => check(),
        "three" => census_three(),
        "four" => census_four(),
        "cone" => light(),
        "multi" => multi(env::args().nth(2).map_or(16384, |a| a.parse().unwrap())),
        "stamps" => stamps(),
        "triple" => census_triple(env::args().nth(2).map_or(25, |a| a.parse().unwrap())),
        _ => println!("verbs: check three"),
    }
    println!("# {:.1} s", start.elapsed().as_secs_f64());
}
