use crate::rules::RULES;
use mrly_math::bang::symmetries;
use mrly_math::bang::universe::{apply, corner_index, corners, permutations};
use std::collections::BTreeSet;

pub type Elem = (Vec<usize>, Vec<u8>, bool);

pub fn act(code: usize, element: &Elem) -> usize {
    let cells = corners(3);
    let mut image = 0usize;
    for (i, cell) in cells.iter().enumerate() {
        if (code >> i) & 1 == 1 {
            image |= 1 << corner_index(&apply(&(element.0.clone(), element.1.clone()), cell));
        }
    }
    if element.2 {
        image ^= 0xff;
    }
    image
}

fn flip(bits: usize) -> Vec<u8> {
    (0..3).map(|j| ((bits >> (2 - j)) & 1) as u8).collect()
}

pub fn group(name: &str) -> Vec<Elem> {
    match name {
        "R" => vec![
            (vec![0, 1, 2], vec![0, 0, 0], false),
            (vec![2, 1, 0], vec![0, 0, 0], false),
        ],
        "H" => vec![
            (vec![0, 1, 2], vec![0, 0, 0], false),
            (vec![2, 1, 0], vec![0, 0, 0], false),
            (vec![0, 1, 2], vec![1, 1, 1], true),
            (vec![2, 1, 0], vec![1, 1, 1], true),
        ],
        "flips" => (0..8).map(|f| (vec![0, 1, 2], flip(f), false)).collect(),
        "perms" => permutations(3)
            .into_iter()
            .map(|p| (p, vec![0, 0, 0], false))
            .collect(),
        "B3" => symmetries(3)
            .into_iter()
            .map(|(p, f)| (p, f, false))
            .collect(),
        "B3xZ2" => symmetries(3)
            .into_iter()
            .flat_map(|(p, f)| [(p.clone(), f.clone(), false), (p.clone(), f.clone(), true)])
            .collect(),
        _ => unreachable!(),
    }
}

pub fn orbit(code: usize, elements: &[Elem]) -> BTreeSet<usize> {
    elements.iter().map(|e| act(code, e)).collect()
}

pub fn representatives(elements: &[Elem]) -> Vec<usize> {
    (0..RULES)
        .map(|code| {
            *orbit(code, elements)
                .iter()
                .next()
                .expect("an orbit is nonempty")
        })
        .collect()
}

fn walk_count(elements: &[Elem]) -> usize {
    representatives(elements)
        .into_iter()
        .collect::<BTreeSet<_>>()
        .len()
}

fn burnside(elements: &[Elem]) -> usize {
    let fixed: usize = elements
        .iter()
        .map(|e| (0..RULES).filter(|&c| act(c, e) == c).count())
        .sum();
    assert_eq!(fixed % elements.len(), 0, "Burnside sum is not divisible");
    fixed / elements.len()
}

fn action(elements: &[Elem]) -> BTreeSet<Vec<usize>> {
    elements
        .iter()
        .map(|e| (0..RULES).map(|c| act(c, e)).collect())
        .collect()
}

pub const NAMES: [&str; 6] = ["R", "H", "flips", "perms", "B3", "B3xZ2"];

pub fn report() {
    println!("GROUP LATTICE");
    println!("group order classes burnside");
    let mut counts = Vec::new();
    for name in NAMES {
        let elements = group(name);
        let walk = walk_count(&elements);
        let burn = burnside(&elements);
        assert_eq!(walk, burn, "{name} disagrees with its Burnside average");
        println!("{name} {} {walk} {burn}", elements.len());
        counts.push((name, walk));
    }
    for (name, want) in [("H", 88), ("B3", 22), ("B3xZ2", 14)] {
        let got = counts
            .iter()
            .find(|(n, _)| *n == name)
            .expect("the group is listed")
            .1;
        assert_eq!(got, want, "{name} class count is not {want}");
    }
    let actions: Vec<(&str, BTreeSet<Vec<usize>>)> =
        NAMES.iter().map(|n| (*n, action(&group(n)))).collect();
    println!("containments");
    for (a, sa) in &actions {
        let inside: Vec<&str> = actions
            .iter()
            .filter(|(b, sb)| b != a && sa.is_subset(sb))
            .map(|(b, _)| *b)
            .collect();
        let inside = if inside.is_empty() {
            "nothing".to_string()
        } else {
            inside.join(" ")
        };
        println!("{a} sits inside {inside}");
    }
    let h = &actions.iter().find(|(n, _)| *n == "H").expect("H").1;
    let b3 = &actions.iter().find(|(n, _)| *n == "B3").expect("B3").1;
    let r = &actions.iter().find(|(n, _)| *n == "R").expect("R").1;
    let meet: BTreeSet<Vec<usize>> = h.intersection(b3).cloned().collect();
    assert_eq!(&meet, r, "H meets B3 outside the reflection group");
    assert!(!h.is_subset(b3), "H is a subgroup of B3");
    assert!(!b3.is_subset(h), "B3 is a subgroup of H");
    println!("H and B3 are incomparable and meet in R of order 2");
    println!("the walk is 256 -> 88 under H and 256 -> 22 under B3, two branches, joined at 14 under B3xZ2");
}
