use crate::census::describe;
use crate::count::{arcs, blocks, mirrors, series};
use crate::design::Design;
use std::time::Instant;

fn law(n: u32) -> i128 {
    (8i128.pow(n) - 1) / 7 - 3i128.pow(n) + n as i128 + 1
}

fn lemma(n: u32) -> Vec<u32> {
    let side = 3usize.pow(n);
    let all: Vec<usize> = (0..=n).map(|j| (3usize.pow(j) - 1) / 2).collect();
    let mut turns: Vec<(usize, usize)> = Vec::new();
    for level in 0..n as usize {
        let s = 3usize.pow(level as u32);
        let mut next = Vec::new();
        for j in 0..3 {
            next.extend(turns.iter().map(|&(a, b)| (j * s + a, j * s + b)));
        }
        for j in 1..3 {
            next.extend(all[..level].iter().map(|&a| (j * s - 1 - a, j * s + a)));
        }
        turns = next;
    }
    let (b, r, t, l) = (0, side, 2 * side, 3 * side);
    let mut p = vec![u32::MAX; 4 * side];
    let mut tie = |x: usize, y: usize| {
        p[x] = y as u32;
        p[y] = x as u32;
    };
    for (k, &a) in all.iter().enumerate() {
        tie(b + a, l + a);
        tie(r + side - 1 - a, t + side - 1 - a);
        if k < n as usize {
            tie(b + side - 1 - a, r + a);
            tie(l + side - 1 - a, t + a);
        }
    }
    for &(x, y) in &turns {
        for o in [b, r, t, l] {
            tie(o + x, o + y);
        }
    }
    p
}

pub fn carpet() {
    let d = Design::new(7, 3, 2);
    println!("CARPET {}, kept {}", d.name(), d.kept());
    for level in 0..=7 {
        let clock = Instant::now();
        let (side, on) = d.cells(level);
        let (loops, strands) = arcs(side, &on);
        let euler = mirrors(side, &on);
        println!("  level {level}: side {side}, arcs loops {loops}, strands {strands}, mirror cycle rank {euler}, law {}, {:.2} s", law(level as u32), clock.elapsed().as_secs_f64());
    }
    let s = series(3, &d.tile, 15);
    let ok = s
        .iter()
        .enumerate()
        .filter(|(n, &x)| x as i128 == law(*n as u32))
        .count();
    let shown: Vec<String> = s.iter().map(|x| x.to_string()).collect();
    println!("  blocks levels 0..15: {}", shown.join(" "));
    println!(
        "  law (8^n - 1)/7 - 3^n + n + 1 holds at {ok} of {} levels; {}",
        s.len(),
        describe(&s)
    );
    let built = blocks(3, &d.tile, 10);
    let same = (0..=10)
        .filter(|&n| built[n].partner == lemma(n as u32))
        .count();
    println!("  matching lemma equals the glued matching at {same} of 11 levels 0..10");
    let splits = (0..15)
        .filter(|&n| {
            let u = (3i128.pow(n as u32) - 2 * n as i128 - 1) / 2;
            s[n + 1] as i128 - 8 * s[n] as i128 == 10 * u + 3 * n as i128
        })
        .count();
    println!("  new loops per gluing = 10 U_n + 3n with U_n = (3^n - 2n - 1)/2 at {splits} of 15 gluings");
}
