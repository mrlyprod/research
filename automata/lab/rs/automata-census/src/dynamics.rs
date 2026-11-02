use crate::groups::{group, orbit};
use crate::rules::{output, RULES};
use std::collections::BTreeSet;

fn arrows(rule: usize) -> Vec<(usize, usize, u8)> {
    let mut out = Vec::new();
    for a in 0..2u8 {
        for b in 0..2u8 {
            for c in 0..2u8 {
                out.push((
                    2 * a as usize + b as usize,
                    2 * b as usize + c as usize,
                    output(rule, a, b, c),
                ));
            }
        }
    }
    out
}

pub fn surjective(rule: usize) -> bool {
    let edges = arrows(rule);
    let start = 0b1111usize;
    let mut seen = BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(set) = stack.pop() {
        for label in 0..2u8 {
            let mut next = 0usize;
            for &(u, v, l) in &edges {
                if l == label && (set >> u) & 1 == 1 {
                    next |= 1 << v;
                }
            }
            if next == 0 {
                return false;
            }
            if seen.insert(next) {
                stack.push(next);
            }
        }
    }
    true
}

pub fn balanced_to(rule: usize, length: usize) -> bool {
    let mut counts = vec![0u32; 1 << length];
    for source in 0..1usize << (length + 2) {
        let mut word = 0usize;
        for i in 0..length {
            let a = ((source >> i) & 1) as u8;
            let b = ((source >> (i + 1)) & 1) as u8;
            let c = ((source >> (i + 2)) & 1) as u8;
            word |= (output(rule, a, b, c) as usize) << i;
        }
        counts[word] += 1;
    }
    counts.iter().all(|&c| c == 4)
}

pub fn injective(rule: usize) -> bool {
    let edges = arrows(rule);
    let mut adjacency = vec![BTreeSet::new(); 16];
    for &(u1, v1, l1) in &edges {
        for &(u2, v2, l2) in &edges {
            if l1 == l2 {
                adjacency[4 * u1 + u2].insert(4 * v1 + v2);
            }
        }
    }
    let mut core: BTreeSet<usize> = (0..16).collect();
    loop {
        let outs: BTreeSet<usize> = core
            .iter()
            .filter(|p| adjacency[**p].iter().any(|q| core.contains(q)))
            .copied()
            .collect();
        let ins: BTreeSet<usize> = outs
            .iter()
            .flat_map(|p| adjacency[*p].iter().filter(|q| outs.contains(q)).copied())
            .collect();
        let next: BTreeSet<usize> = outs.intersection(&ins).copied().collect();
        if next == core {
            break;
        }
        core = next;
    }
    !core.iter().any(|p| p / 4 != p % 4)
}

pub fn surjective_set() -> Vec<usize> {
    (0..RULES).filter(|&r| surjective(r)).collect()
}

pub fn injective_set() -> Vec<usize> {
    (0..RULES).filter(|&r| injective(r)).collect()
}

pub fn report() {
    println!("GEOMETRY IS NOT DYNAMICS");
    let surj = surjective_set();
    let inj = injective_set();
    println!(
        "de Bruijn subset walk: {} surjective rules {surj:?}",
        surj.len()
    );
    let mut previous: Vec<usize> = (0..RULES).collect();
    let mut settles = None;
    for length in 1..=12 {
        let balanced: Vec<usize> = (0..RULES).filter(|&r| balanced_to(r, length)).collect();
        assert!(
            balanced.iter().all(|r| previous.contains(r)),
            "the balanced set grows from length {} to {length}",
            length - 1
        );
        assert!(
            surj.iter().all(|r| balanced.contains(r)),
            "a surjective rule is not balanced at length {length}"
        );
        println!(
            "balanced on words of length {length}: {} rules",
            balanced.len()
        );
        if balanced == surj && settles.is_none() {
            settles = Some(length);
        }
        previous = balanced;
    }
    let settles = settles.expect("the balance test never reaches the de Bruijn set");
    println!("the exact balance test on all words of length n shrinks monotonically and equals the de Bruijn set from n = {settles} to n = 12");
    if surj.len() == 30 {
        println!("the count is 30, the published figure");
    } else {
        println!(
            "UNCLEAR the count is {} against the published 30",
            surj.len()
        );
    }
    println!("pair-graph core: {} injective rules {inj:?}", inj.len());
    assert_eq!(
        inj,
        vec![15, 51, 85, 170, 204, 240],
        "the reversible rules are not the expected six"
    );
    assert!(
        inj.iter().all(|r| surj.contains(r)),
        "an injective rule is not surjective"
    );
    let b3 = group("B3");
    let identity_class: Vec<usize> = orbit(204, &b3).into_iter().collect();
    assert_eq!(
        inj, identity_class,
        "the reversible rules are not the B3 orbit of 204"
    );
    println!(
        "the reversible six are exactly the B3 orbit of 204, the single-axis degree-1 designs"
    );
    let surjective_flags: Vec<bool> = (0..RULES).map(|r| surj.contains(&r)).collect();
    let mut mixed = Vec::new();
    let mut constant_yes = 0usize;
    let mut constant_no = 0usize;
    let mut reps: Vec<usize> = (0..RULES)
        .map(|c| *orbit(c, &b3).iter().next().expect("nonempty"))
        .collect();
    reps.sort();
    reps.dedup();
    for rep in &reps {
        let cls: Vec<usize> = orbit(*rep, &b3).into_iter().collect();
        let yes: Vec<usize> = cls
            .iter()
            .copied()
            .filter(|c| surjective_flags[*c])
            .collect();
        let no: Vec<usize> = cls
            .iter()
            .copied()
            .filter(|c| !surjective_flags[*c])
            .collect();
        if yes.is_empty() {
            constant_no += 1;
        } else if no.is_empty() {
            constant_yes += 1;
        } else {
            mixed.push((*rep, cls.len(), yes, no));
        }
    }
    println!(
        "B3 classes: {constant_yes} all surjective, {constant_no} none surjective, {} mixed",
        mixed.len()
    );
    for (rep, size, yes, no) in &mixed {
        println!("mixed class rep {rep} size {size} surjective {yes:?} not {no:?}");
    }
    assert!(
        !mixed.is_empty(),
        "surjectivity is constant on every B3 class"
    );
    let thirty = orbit(30, &b3);
    assert!(thirty.contains(&54), "30 and 54 do not share a B3 orbit");
    assert!(
        surjective_flags[30] && !surjective_flags[54],
        "30 and 54 do not split on surjectivity"
    );
    println!("witness: 30 and 54 share the B3 orbit of rep 30, 30 is surjective and 54 is not");
    println!("surjectivity is not a B3 invariant, witnessed by 30 and 54 in one orbit");
    let mut reps: Vec<usize> = (0..RULES)
        .map(|c| *orbit(c, &b3).iter().next().expect("nonempty"))
        .collect();
    reps.sort();
    reps.dedup();
    let mut mixed_rev = 0usize;
    for rep in &reps {
        let cls: Vec<usize> = orbit(*rep, &b3).into_iter().collect();
        if cls.iter().any(|c| injective(*c)) && cls.iter().any(|c| !injective(*c)) {
            mixed_rev += 1;
        }
    }
    assert_eq!(
        mixed_rev, 0,
        "reversibility is not constant on every B3 orbit"
    );
    println!("reversibility IS constant on every B3 orbit, 0 mixed of {}, the reversible six being exactly the orbit of 204", reps.len());
    let mut mixed_h = 0usize;
    let h = group("H");
    let mut hreps: Vec<usize> = (0..RULES)
        .map(|c| *orbit(c, &h).iter().next().expect("nonempty"))
        .collect();
    hreps.sort();
    hreps.dedup();
    for rep in &hreps {
        let cls: Vec<usize> = orbit(*rep, &h).into_iter().collect();
        if cls.iter().any(|c| surjective_flags[*c]) && cls.iter().any(|c| !surjective_flags[*c]) {
            mixed_h += 1;
        }
    }
    println!(
        "H classes with mixed surjectivity: {mixed_h} of {}",
        hreps.len()
    );
}
