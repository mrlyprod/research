use crate::series::Rep;
use crate::word::{render, CODES};

pub fn family(prefix: u8, last: u8, length: usize) -> Vec<u8> {
    let mut word = vec![prefix; length - 1];
    word.push(last);
    word
}

pub fn components(word: &[u8]) -> u64 {
    render(word).components()
}

pub fn merging(word: &[u8]) -> u64 {
    render(word).merging()
}

pub fn max_components(length: usize) -> (u64, Vec<u8>) {
    let mut best = 0u64;
    let mut witness: Vec<u8> = Vec::new();
    let mut word = vec![0u8; length];
    walk(&mut word, 0, &mut |candidate: &[u8]| {
        let count = render(candidate).components();
        if count > best {
            best = count;
            witness = candidate.to_vec();
        }
    });
    (best, witness)
}

pub fn max_merging(length: usize) -> (u64, Vec<u8>) {
    let mut best = 0u64;
    let mut witness: Vec<u8> = Vec::new();
    let mut word = vec![0u8; length];
    walk(&mut word, 0, &mut |candidate: &[u8]| {
        let count = render(candidate).merging();
        if count > best {
            best = count;
            witness = candidate.to_vec();
        }
    });
    (best, witness)
}

fn walk(word: &mut Vec<u8>, at: usize, visit: &mut dyn FnMut(&[u8])) {
    if at == word.len() {
        visit(word);
        return;
    }
    for &code in CODES.iter() {
        word[at] = code;
        walk(word, at + 1, visit);
    }
}

pub fn max_predicted(rep: &Rep, length: usize) -> i64 {
    let state: Vec<i128> = rep.lambda.iter().map(|value| *value as i128).collect();
    let mut best = i64::MIN;
    descend(rep, &state, length, &mut best);
    best
}

fn descend(rep: &Rep, state: &[i128], left: usize, best: &mut i64) {
    if left == 0 {
        let value: i128 = state
            .iter()
            .zip(rep.gamma.iter())
            .map(|(a, b)| a * *b as i128)
            .sum();
        if value as i64 > *best {
            *best = value as i64;
        }
        return;
    }
    for &code in CODES.iter() {
        let matrix = &rep.matrices[&code];
        let mut next = vec![0i128; state.len()];
        for (i, weight) in state.iter().enumerate() {
            if *weight == 0 {
                continue;
            }
            for (j, slot) in next.iter_mut().enumerate() {
                *slot += weight * matrix[i][j] as i128;
            }
        }
        descend(rep, &next, left - 1, best);
    }
}

pub fn independent_bound(length: usize) -> u64 {
    let side = 1u64 << length;
    side * side / 2
}
