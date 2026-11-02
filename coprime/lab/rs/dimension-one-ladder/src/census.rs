use mrly_num::factor::gcd;

pub struct Ray {
    pub a: u64,
    pub b: u64,
    pub weight: u64,
}

fn octave(a: u64, b: u64) -> usize {
    let top = a.max(b);
    if top == 1 {
        return 0;
    }
    let mut j = 1;
    while 3u64.pow(j as u32) <= top {
        j += 1;
    }
    j
}

pub fn rays(level: u32) -> (u64, Vec<Ray>) {
    let size = 1usize << level;
    let value: Vec<u64> = (0..size)
        .map(|mask| {
            (0..level)
                .filter(|j| mask >> j & 1 == 1)
                .map(|j| 3u64.pow(j))
                .sum()
        })
        .collect();
    let mut keys: Vec<u64> = Vec::with_capacity(3usize.pow(level));
    let mut coprime = 0u64;
    for x in 0..size {
        let free = !x & (size - 1);
        let mut y = free;
        loop {
            let (a, b) = (value[x], value[y]);
            let g = gcd(a as u128, b as u128) as u64;
            if g == 1 {
                coprime += 1;
            }
            if let (Some(num), Some(den)) = (a.checked_div(g), b.checked_div(g)) {
                keys.push(num << 32 | den);
            }
            if y == 0 {
                break;
            }
            y = (y - 1) & free;
        }
    }
    keys.sort_unstable();
    let mut out = Vec::new();
    let mut index = 0;
    while index < keys.len() {
        let key = keys[index];
        let mut end = index;
        while end < keys.len() && keys[end] == key {
            end += 1;
        }
        out.push(Ray {
            a: key >> 32,
            b: key & 0xffff_ffff,
            weight: (end - index) as u64,
        });
        index = end;
    }
    (coprime, out)
}

pub fn report(level: u32, detail: bool) {
    let (coprime, rays) = rays(level);
    let points = 3f64.powi(level as i32);
    let mut bins = vec![0u64; level as usize + 1];
    let mut z = 0u128;
    let mut mass = 0u64;
    let mut peak = 0u64;
    let mut singles = 0u64;
    let mut light = 0u64;
    let mut light_z = 0u128;
    for ray in rays.iter() {
        bins[octave(ray.a, ray.b)] += 1;
        if ray.a == 0 || ray.b == 0 {
            continue;
        }
        z += (ray.weight as u128).pow(2);
        mass += ray.weight;
        peak = peak.max(ray.weight);
        if ray.weight == 1 {
            singles += 1;
        }
        if (2..=5).contains(&ray.weight) {
            light += 1;
            light_z += (ray.weight as u128).pow(2);
        }
    }
    println!(
        "census n {} gasket points {} primitive {}",
        level,
        3u64.pow(level),
        coprime
    );
    println!(
        "occupied rays {} with fibres {} without",
        rays.len(),
        rays.len() - 2
    );
    println!(
        "Z {} Z/3^n {:.6} sum M {} 3^n-2^(n+1)+1 {} max M {}",
        z,
        z as f64 / points,
        mass,
        3u64.pow(level) - 2u64.pow(level + 1) + 1,
        peak
    );
    let octaves: Vec<String> = bins
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0)
        .map(|(j, c)| format!("{}:{}", j, c))
        .collect();
    println!("octaves {}", octaves.join(" "));
    if let Some(occ) = bins.get(8) {
        println!("occ(8,n)/3^8 {:.3}", *occ as f64 / 6561.0);
    }
    let cut = |j: usize| (bins[..=j].iter().sum::<u64>() as f64).ln() / 3f64.ln() / level as f64;
    if level.is_multiple_of(2) {
        println!("band exponent j <= n/2 {:.4}", cut(level as usize / 2));
    } else {
        println!(
            "band exponent j <= floor(n/2) {:.4} or ceil {:.4}",
            cut(level as usize / 2),
            cut(level as usize / 2 + 1)
        );
    }
    if !detail {
        return;
    }
    let share = |part: u128| 100.0 * part as f64 / z as f64;
    println!(
        "M = 1 rays {} ({:.0}% of Z), M in [2,5] rays {} ({:.0}% of Z)",
        singles,
        share(singles as u128),
        light,
        share(light_z)
    );
    let mut heavy: Vec<&Ray> = rays.iter().filter(|r| r.a > 0 && r.b > 0).collect();
    heavy.sort_by(|p, q| q.weight.cmp(&p.weight).then(p.a.cmp(&q.a)));
    let top: Vec<String> = heavy[..10]
        .iter()
        .map(|r| format!("({},{}):{}", r.a, r.b, r.weight))
        .collect();
    let top_z: u128 = heavy[..10].iter().map(|r| (r.weight as u128).pow(2)).sum();
    println!(
        "ten heaviest {} carry {:.0}% of Z",
        top.join(" "),
        share(top_z)
    );
}
