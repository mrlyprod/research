use std::env;
use std::thread;
use std::time::Instant;

const P: u32 = 7;
const DIM: usize = 6;
const N: usize = 117649;
const THREADS: usize = 8;

// FIELDS

fn norm2(a: u32, b: u32) -> u32 {
    (a * a + 4 * b * b) % P
}

fn mul3(x: [u32; 3], y: [u32; 3]) -> [u32; 3] {
    let (a, b, c) = (x[0], x[1], x[2]);
    let (d, e, f) = (y[0], y[1], y[2]);
    [
        (a * d + 2 * (b * f + c * e)) % P,
        (a * e + b * d + 2 * c * f) % P,
        (a * f + b * e + c * d) % P,
    ]
}

fn det3(m: [[u32; 3]; 3]) -> u32 {
    let t = |i: usize, j: usize| m[i][j] as i64;
    let d = t(0, 0) * (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1))
        - t(0, 1) * (t(1, 0) * t(2, 2) - t(1, 2) * t(2, 0))
        + t(0, 2) * (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0));
    d.rem_euclid(P as i64) as u32
}

fn norm3_det(x: [u32; 3]) -> u32 {
    let cols = [mul3(x, [1, 0, 0]), mul3(x, [0, 1, 0]), mul3(x, [0, 0, 1])];
    let mut m = [[0u32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = cols[j][i];
        }
    }
    det3(m)
}

fn norm3(x: [u32; 3]) -> u32 {
    let (a, b, c) = (x[0], x[1], x[2]);
    (a * a * a + 2 * b * b * b + 4 * c * c * c + a * b * c) % P
}

fn field_checks() {
    for a in 0..P {
        for b in 0..P {
            assert!(norm2(a, b) != 0 || (a == 0 && b == 0));
        }
    }
    let mut hit2 = [false; 7];
    for a in 0..P {
        for b in 0..P {
            hit2[norm2(a, b) as usize] = true;
        }
    }
    assert!(hit2.iter().all(|&h| h));
    let mut hit3 = [false; 7];
    for i in 0..343 {
        let x = [i % 7, (i / 7) % 7, i / 49];
        assert_eq!(norm3(x), norm3_det(x));
        assert!(norm3(x) != 0 || i == 0);
        hit3[norm3(x) as usize] = true;
        for j in 0..343 {
            let y = [j % 7, (j / 7) % 7, j / 49];
            assert_eq!(norm3(mul3(x, y)), norm3(x) * norm3(y) % P);
        }
    }
    assert!(hit3.iter().all(|&h| h));
    for a in 0..P {
        assert!((a * a + 4) % P != 0);
        assert!((a * a * a + 5) % P != 0);
    }
}

// POINTS

fn coords(i: usize) -> [u32; DIM] {
    let mut c = [0u32; DIM];
    let mut x = i;
    for k in 0..DIM {
        c[k] = (x % P as usize) as u32;
        x /= P as usize;
    }
    c
}

fn pack(c: &[u32; DIM]) -> u32 {
    let mut v = 0u32;
    for k in 0..DIM {
        v |= c[k] << (4 * k);
    }
    v
}

fn unpack(v: u32) -> [u32; DIM] {
    let mut c = [0u32; DIM];
    for k in 0..DIM {
        c[k] = (v >> (4 * k)) & 15;
    }
    c
}

fn padd(a: u32, b: u32) -> u32 {
    let s = a + b;
    let m = (s + 0x111111) & 0x888888;
    s - (m >> 3) * 7
}

fn pneg(a: u32) -> u32 {
    let t = 0x777777 - a;
    let m = (t + 0x111111) & 0x888888;
    t - (m >> 3) * 7
}

fn poly(c: &[u32; DIM]) -> u32 {
    (c[0] + norm2(c[1], c[2]) + norm3([c[3], c[4], c[5]])) % P
}

struct Bits(Vec<u64>);

impl Bits {
    fn new() -> Bits {
        Bits(vec![0u64; 1 << 18])
    }
    fn set(&mut self, v: u32) {
        self.0[(v >> 6) as usize] |= 1u64 << (v & 63);
    }
    #[inline(always)]
    fn has(&self, v: u32) -> bool {
        self.0[(v >> 6) as usize] >> (v & 63) & 1 == 1
    }
}

fn level_set(u: u32) -> (Bits, Vec<u32>, Vec<u32>) {
    let mut bits = Bits::new();
    let mut inside = Vec::new();
    let mut outside = Vec::new();
    for i in 0..N {
        let c = coords(i);
        let v = pack(&c);
        if poly(&c) == u {
            bits.set(v);
            inside.push(v);
        } else {
            outside.push(v);
        }
    }
    (bits, inside, outside)
}

fn all_points() -> Vec<u32> {
    (0..N).map(|i| pack(&coords(i))).collect()
}

fn chunks<T: Clone>(v: &[T]) -> Vec<Vec<T>> {
    let size = (v.len() + THREADS - 1) / THREADS;
    v.chunks(size.max(1)).map(|c| c.to_vec()).collect()
}

// LINES

fn line_histogram(bits: &Bits) -> [u64; 8] {
    let points = all_points();
    let mut dirs = Vec::new();
    for &d in points.iter().skip(1) {
        let c = unpack(d);
        let k = (0..DIM).find(|&k| c[k] != 0).unwrap();
        if c[k] == 1 {
            dirs.push((d, k));
        }
    }
    assert_eq!(dirs.len(), (N - 1) / 6);
    let parts = chunks(&dirs);
    let hist = thread::scope(|s| {
        let handles: Vec<_> = parts
            .iter()
            .map(|part| {
                let points = &points;
                s.spawn(move || {
                    let mut h = [0u64; 8];
                    for &(d, k) in part {
                        for &x in points {
                            if (x >> (4 * k)) & 15 != 0 {
                                continue;
                            }
                            let mut n = 0;
                            let mut y = x;
                            for _ in 0..P {
                                n += bits.has(y) as usize;
                                y = padd(y, d);
                            }
                            h[n] += 1;
                        }
                    }
                    h
                })
            })
            .collect();
        let mut h = [0u64; 8];
        for handle in handles {
            let part = handle.join().unwrap();
            for i in 0..8 {
                h[i] += part[i];
            }
        }
        h
    });
    hist
}

// SATURATION

fn saturation(bits: &Bits, outside: &[u32]) -> Vec<(u32, u32, u32)> {
    let points = all_points();
    let dirs: Vec<u32> = points[1..].to_vec();
    let parts = chunks(outside);
    thread::scope(|s| {
        let handles: Vec<_> = parts
            .iter()
            .map(|part| {
                let dirs = &dirs;
                s.spawn(move || {
                    let mut out = Vec::with_capacity(part.len());
                    for &b in part {
                        let mut ends = 0u32;
                        let mut inner = 0u32;
                        for &d in dirs {
                            let b1 = padd(b, d);
                            if !bits.has(b1) {
                                continue;
                            }
                            let b2 = padd(b1, d);
                            if !bits.has(b2) {
                                continue;
                            }
                            if bits.has(padd(b2, d)) {
                                ends += 1;
                            }
                            if bits.has(padd(b, pneg(d))) {
                                inner += 1;
                            }
                        }
                        out.push((b, ends, inner));
                    }
                    out
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    })
}

fn summarize(name: &str, values: &[u32]) {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    let mut runs: Vec<(u32, usize)> = Vec::new();
    for &v in &sorted {
        match runs.last_mut() {
            Some((w, c)) if *w == v => *c += 1,
            _ => runs.push((v, 1)),
        }
    }
    println!(
        "{name}: min {} median {} max {} distinct {}",
        sorted[0],
        sorted[n / 2],
        sorted[n - 1],
        runs.len()
    );
    let shown: Vec<String> = runs.iter().map(|(v, c)| format!("{v}:{c}")).collect();
    println!("{name} histogram value:count {}", shown.join(" "));
}

// INDEPENDENT CHECK

fn digits_add(a: &[u8; DIM], b: &[u8; DIM]) -> [u8; DIM] {
    let mut c = [0u8; DIM];
    for k in 0..DIM {
        let s = a[k] + b[k];
        c[k] = if s >= 7 { s - 7 } else { s };
    }
    c
}

fn digits_index(c: &[u8; DIM]) -> usize {
    let mut i = 0usize;
    for k in (0..DIM).rev() {
        i = i * 7 + c[k] as usize;
    }
    i
}

fn digits_of(i: usize) -> [u8; DIM] {
    let mut c = [0u8; DIM];
    let mut x = i;
    for k in 0..DIM {
        c[k] = (x % 7) as u8;
        x /= 7;
    }
    c
}

fn progressions_by_digits(member: &[bool]) -> u64 {
    let set: Vec<usize> = (0..N).filter(|&i| member[i]).collect();
    let parts = chunks(&set);
    thread::scope(|s| {
        let handles: Vec<_> = parts
            .iter()
            .map(|part| {
                s.spawn(move || {
                    let mut found = 0u64;
                    for &x in part {
                        let xd = digits_of(x);
                        for d in 1..N {
                            let dd = digits_of(d);
                            let y = digits_add(&xd, &dd);
                            if !member[digits_index(&y)] {
                                continue;
                            }
                            let z = digits_add(&y, &dd);
                            if !member[digits_index(&z)] {
                                continue;
                            }
                            let w = digits_add(&z, &dd);
                            if member[digits_index(&w)] {
                                found += 1;
                            }
                        }
                    }
                    found
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    })
}

fn membership(inside: &[u32]) -> Vec<bool> {
    let mut member = vec![false; N];
    for &v in inside {
        let c = unpack(v);
        let d: [u8; DIM] = [
            c[0] as u8, c[1] as u8, c[2] as u8, c[3] as u8, c[4] as u8, c[5] as u8,
        ];
        member[digits_index(&d)] = true;
    }
    member
}

// MAIN

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let verb = args.first().map(|s| s.as_str()).unwrap_or("all");
    let u: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
    let t0 = Instant::now();
    field_checks();
    println!("fields: F_49 = F_7(sqrt 3), F_343 = F_7(cbrt 2), norms multiplicative, anisotropic, onto F_7");
    let (bits, inside, outside) = level_set(u);
    println!(
        "level set P = {u}: {} points of {N}, outside {}, {:.1}s",
        inside.len(),
        outside.len(),
        t0.elapsed().as_secs_f64()
    );
    assert_eq!(inside.len(), 16807);
    if verb == "lines" || verb == "all" {
        let h = line_histogram(&bits);
        let total: u64 = h.iter().sum();
        assert_eq!(total as usize, N * (N - 1) / 42);
        println!(
            "lines {total}: meeting the set in 0..7 points {:?}, largest {}, {:.1}s",
            h,
            (0..8).rev().find(|&i| h[i] > 0).unwrap(),
            t0.elapsed().as_secs_f64()
        );
    }
    if verb == "saturate" || verb == "all" {
        let sat = saturation(&bits, &outside);
        let ends: Vec<u32> = sat.iter().map(|s| s.1).collect();
        let inner: Vec<u32> = sat.iter().map(|s| s.2).collect();
        let total: Vec<u32> = sat.iter().map(|s| s.1 + s.2).collect();
        summarize("end", &ends);
        summarize("inner", &inner);
        summarize("all", &total);
        let addable = sat.iter().filter(|s| s.1 + s.2 == 0).count();
        println!(
            "outside points {}, addable (on no progression with three points of the set) {addable}, {:.1}s",
            sat.len(),
            t0.elapsed().as_secs_f64()
        );
        let sample = sat.iter().min_by_key(|s| s.1 + s.2).unwrap();
        println!(
            "least covered point {:?}: {} progressions as an end, {} as an inner term",
            unpack(sample.0),
            sample.1,
            sample.2
        );
    }
    if verb == "check" || verb == "all" {
        let member = membership(&inside);
        let found = progressions_by_digits(&member);
        println!(
            "independent check by digit arithmetic: {} four-term progressions with distinct terms in the set, {:.1}s",
            found,
            t0.elapsed().as_secs_f64()
        );
        assert_eq!(found, 0);
        let mut spoiled = member.clone();
        let b = outside[outside.len() / 3];
        let c = unpack(b);
        let d: [u8; DIM] = [
            c[0] as u8, c[1] as u8, c[2] as u8, c[3] as u8, c[4] as u8, c[5] as u8,
        ];
        spoiled[digits_index(&d)] = true;
        let found = progressions_by_digits(&spoiled);
        println!(
            "negative control, the set plus the outside point {:?}: {} progressions, {:.1}s",
            c,
            found,
            t0.elapsed().as_secs_f64()
        );
        assert!(found > 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(state: &mut u64) -> u32 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*state >> 33) as u32
    }

    #[test]
    fn packed_arithmetic_matches_digits() {
        let mut s = 7u64;
        for _ in 0..20000 {
            let i = lcg(&mut s) as usize % N;
            let j = lcg(&mut s) as usize % N;
            let a = coords(i);
            let b = coords(j);
            let mut c = [0u32; DIM];
            let mut neg = [0u32; DIM];
            for k in 0..DIM {
                c[k] = (a[k] + b[k]) % P;
                neg[k] = (P - a[k]) % P;
            }
            assert_eq!(padd(pack(&a), pack(&b)), pack(&c));
            assert_eq!(pneg(pack(&a)), pack(&neg));
            assert_eq!(unpack(pack(&a)), a);
        }
    }

    #[test]
    fn fields_are_fields() {
        field_checks();
    }

    #[test]
    fn largest_progression_free_subset_of_the_line_has_four_points() {
        let mut best = 0;
        for mask in 0u32..128 {
            let mut free = true;
            for s in 0..7 {
                for e in 1..7 {
                    if (0..4).all(|i| mask >> ((s + i * e) % 7) & 1 == 1) {
                        free = false;
                    }
                }
            }
            if free {
                best = best.max(mask.count_ones());
            }
        }
        assert_eq!(best, 4);
    }

    #[test]
    fn a_line_of_the_level_set_holds_at_most_three_points() {
        let (bits, _, _) = level_set(0);
        let mut s = 11u64;
        for _ in 0..200000 {
            let x = pack(&coords(lcg(&mut s) as usize % N));
            let d = pack(&coords(1 + lcg(&mut s) as usize % (N - 1)));
            let mut n = 0;
            let mut y = x;
            for _ in 0..P {
                n += bits.has(y) as u32;
                y = padd(y, d);
            }
            assert!(n <= 3);
        }
    }
}
