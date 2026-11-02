use crate::count::blocks;
use crate::design::Design;

fn lemma(code: u128, n: u32) -> Vec<u32> {
    let side = 1usize << n;
    let (b, r, t, l) = (0, side, 2 * side, 3 * side);
    let mut p = vec![u32::MAX; 4 * side];
    let mut tie = |x: usize, y: usize| {
        p[x] = y as u32;
        p[y] = x as u32;
    };
    match code {
        7 => {
            tie(b, l);
            tie(b + 1, l + 1);
            for k in 1..side / 2 {
                tie(b + 2 * k, b + 2 * k + 1);
                tie(l + 2 * k, l + 2 * k + 1);
            }
            for k in 0..side / 2 {
                tie(r + 2 * k, r + 2 * k + 1);
                tie(t + 2 * k, t + 2 * k + 1);
            }
        }
        11 => {
            tie(b, l);
            tie(b + side - 1, r);
            tie(r + side - 1, t + side - 1);
            for y in 1..side {
                tie(l + y, t + side - 1 - y);
            }
            for k in 1..side / 2 {
                tie(b + 2 * k - 1, b + 2 * k);
                tie(r + 2 * k - 1, r + 2 * k);
            }
        }
        _ => {
            tie(b, l);
            tie(r + side - 1, t + side - 1);
            for x in 1..side {
                tie(b + x, r + side - 1 - x);
                tie(l + x, t + side - 1 - x);
            }
        }
    }
    p
}

pub fn two() {
    println!("BASE 2 LEMMAS: glued matching against the stated matching, levels 1..14");
    for (code, gain) in [(7u128, "2^n - 2"), (11, "2^(n-1)"), (9, "1")] {
        let d = Design::full(code, 2);
        let built = blocks(2, &d.tile, 14);
        let same = (1..=14)
            .filter(|&n| built[n].partner == lemma(code, n as u32))
            .count();
        let law = |n: u32| -> u128 {
            match code {
                7 => (1 << n) - 2,
                11 => 1 << (n - 1),
                _ => 1,
            }
        };
        let gains = (1..14)
            .filter(|&n| built[n + 1].loops - d.kept() as u128 * built[n].loops == law(n as u32))
            .count();
        println!("  {}: lemma holds at {same} of 14 levels, new loops per gluing = {gain} at {gains} of 13 gluings", d.name());
    }
}
