use crate::design::BASE;
use mrly_math::bang::factory::residue_corners;

pub fn digits() -> Vec<[i64; 3]> {
    residue_corners(3, BASE)
        .into_iter()
        .filter(|corner| corner.iter().filter(|digit| **digit == 1).count() <= 1)
        .map(|corner| [corner[0] as i64, corner[1] as i64, corner[2] as i64])
        .collect()
}

pub fn cube_digits() -> Vec<[i64; 3]> {
    residue_corners(3, BASE)
        .into_iter()
        .map(|corner| [corner[0] as i64, corner[1] as i64, corner[2] as i64])
        .collect()
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

pub fn primitive(v: [i64; 3]) -> bool {
    gcd(gcd(v[0], v[1]), v[2]) == 1
}

pub fn shadow(level: usize, view: [i64; 3], digits: &[[i64; 3]]) -> usize {
    let count = digits.len();
    let total = count.pow(level as u32);
    let side = (BASE as i64).pow(level as u32);
    let span = 4 * side * 3 + 4;
    let mut keys: Vec<i64> = Vec::with_capacity(total);
    for index in 0..total {
        let mut left = index;
        let mut point = [0i64; 3];
        for _ in 0..level {
            let digit = digits[left % count];
            left /= count;
            for axis in 0..3 {
                point[axis] = point[axis] * BASE as i64 + digit[axis];
            }
        }
        let cross = [
            point[1] * view[2] - point[2] * view[1],
            point[2] * view[0] - point[0] * view[2],
            point[0] * view[1] - point[1] * view[0],
        ];
        let key = ((cross[0] + span) * 2 * span + cross[1] + span) * 2 * span + cross[2] + span;
        keys.push(key);
    }
    keys.sort_unstable();
    keys.dedup();
    keys.len()
}

pub fn views(bound: i64) -> Vec<[i64; 3]> {
    let mut out = Vec::new();
    for a in 0..=bound {
        for b in a..=bound {
            for c in b..=bound {
                if c == 0 || !primitive([a, b, c]) {
                    continue;
                }
                out.push([a, b, c]);
            }
        }
    }
    out
}
