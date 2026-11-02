use crate::rules::{evolve, mirrored, RULES};
use std::collections::{BTreeMap, BTreeSet};

pub const STEPS: usize = 256;

pub struct Census {
    pub class: Vec<usize>,
    pub occurring: Vec<u8>,
}

struct Pass {
    class: Vec<usize>,
    occurring: Vec<u8>,
    distinct: usize,
    reflected: usize,
    histogram: BTreeMap<usize, usize>,
}

fn classify(steps: usize, pad: usize) -> Pass {
    let mut seen: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut folded: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut class = vec![0usize; RULES];
    let mut occurring = vec![0u8; RULES];
    for rule in 0..RULES {
        let (diagram, met) = evolve(rule, steps, pad);
        occurring[rule] = met;
        let mirror = mirrored(&diagram);
        let next = seen.len();
        class[rule] = *seen.entry(diagram.cells.clone()).or_insert(next);
        folded.insert(diagram.cells.min(mirror));
    }
    let mut sizes: BTreeMap<usize, usize> = BTreeMap::new();
    for &c in &class {
        *sizes.entry(c).or_insert(0) += 1;
    }
    let mut histogram: BTreeMap<usize, usize> = BTreeMap::new();
    for size in sizes.values() {
        *histogram.entry(*size).or_insert(0) += 1;
    }
    Pass {
        distinct: seen.len(),
        reflected: folded.len(),
        class,
        occurring,
        histogram,
    }
}

pub fn report() -> Census {
    println!("SINGLE-SEED CENSUS");
    println!("one live cell on a line padded by pad cells beyond the 2T+1 window that is cropped and compared");
    let mut last: Option<Pass> = None;
    for steps in [64usize, 128, STEPS] {
        let narrow = classify(steps, steps);
        let wide = classify(steps, 2 * steps);
        assert_eq!(
            narrow.class, wide.class,
            "the diagrams move between pad = T and pad = 2T at T = {steps}"
        );
        assert_eq!(
            narrow.occurring, wide.occurring,
            "the occurring neighbourhoods move between pad = T and pad = 2T at T = {steps}"
        );
        println!(
            "T = {steps}: {} distinct diagrams, {} up to left-right reflection, class-size histogram {:?}, identical at pad = T and pad = 2T",
            narrow.distinct, narrow.reflected, narrow.histogram
        );
        if let Some(before) = &last {
            assert_eq!(
                before.occurring, narrow.occurring,
                "the occurring neighbourhoods move with T at {steps}"
            );
        }
        last = Some(narrow);
    }
    let pass = last.expect("three depths ran");
    let (class, occurring, distinct) = (pass.class, pass.occurring, pass.distinct);
    println!("the occurring neighbourhood set is the same at T = 64, 128 and {STEPS}");
    let mut failures = 0usize;
    for a in 0..RULES {
        for b in 0..RULES {
            let mask = occurring[a] as usize;
            if occurring[a] == occurring[b] && a & mask == b & mask && class[a] != class[b] {
                failures += 1;
            }
        }
    }
    assert_eq!(
        failures, 0,
        "rules sharing the occurring key do not share the diagram"
    );
    println!("law: equal key implies equal diagram, 0 failures over all {} ordered pairs, the key being the occurring set together with the rule restricted to it", RULES * RULES);
    let keys: BTreeSet<(u8, usize)> = (0..RULES)
        .map(|r| (occurring[r], r & occurring[r] as usize))
        .collect();
    let mut members: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (rule, &c) in class.iter().enumerate() {
        members.entry(c).or_default().push(rule);
    }
    let split: Vec<usize> = members
        .values()
        .filter(|rules| {
            rules
                .iter()
                .map(|r| (occurring[*r], r & occurring[*r] as usize))
                .collect::<BTreeSet<_>>()
                .len()
                > 1
        })
        .flatten()
        .copied()
        .collect();
    let split_classes: BTreeSet<usize> = split.iter().map(|r| class[*r]).collect();
    let split_keys: BTreeSet<(u8, usize)> = split
        .iter()
        .map(|r| (occurring[*r], r & occurring[*r] as usize))
        .collect();
    println!(
        "the key takes {} values against {distinct} distinct diagrams, so it separates strictly more than the diagram does and is not a complete invariant",
        keys.len()
    );
    println!("over-separated rules, {} of them spread over {} diagram classes carrying {} keys in all: {split:?}", split.len(), split_classes.len(), split_keys.len());
    assert_eq!(split_keys.len(), 11);
    assert!(
        keys.len() > distinct,
        "the key no longer over-separates the diagram"
    );
    assert_eq!(
        split,
        vec![23, 31, 55, 63, 87, 95, 119, 127, 151, 159, 183, 191, 215, 223, 247, 255],
        "the over-separated rules are not the expected 16"
    );
    assert_eq!(
        split_classes.len(),
        2,
        "the over-separated rules do not fall into two diagram classes"
    );
    let mut occ_sizes: BTreeMap<u32, usize> = BTreeMap::new();
    for o in &occurring {
        *occ_sizes.entry(o.count_ones()).or_insert(0) += 1;
    }
    println!("occurring-set sizes {occ_sizes:?}");
    Census { class, occurring }
}
