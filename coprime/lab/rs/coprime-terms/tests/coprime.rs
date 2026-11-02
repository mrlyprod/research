use coprime_terms::brute;
use coprime_terms::design::{CARPET, MENGER, VICSEK};
use coprime_terms::engine::{methods, terms, terms_with, Mode};

const MENGER_TERMS: [i128; 11] = [
    12,
    270,
    5916,
    123504,
    2538447,
    51497040,
    1038074187,
    20860210527,
    418429711224,
    8382927031902,
    167827226563374,
];

const CARPET_TERMS: [i128; 12] = [
    4,
    32,
    274,
    2320,
    19178,
    155392,
    1248416,
    10013432,
    80226680,
    642182854,
    5138968090,
    41117712068,
];

const VICSEK_TERMS: [i128; 14] = [
    5, 16, 90, 418, 2178, 10560, 54120, 266478, 1338422, 6658480, 33439772, 166998096, 835514606,
    4172530386,
];

#[test]
fn known_terms_reproduce() {
    let menger: Vec<i128> = terms(&MENGER, 11, 4).iter().map(|l| l.value).collect();
    assert_eq!(menger, MENGER_TERMS.to_vec());
    let carpet: Vec<i128> = terms(&CARPET, 12, 4).iter().map(|l| l.value).collect();
    assert_eq!(carpet, CARPET_TERMS.to_vec());
    let vicsek: Vec<i128> = terms(&VICSEK, 14, 4).iter().map(|l| l.value).collect();
    assert_eq!(vicsek, VICSEK_TERMS.to_vec());
}

#[test]
fn enumeration_agrees() {
    for design in [MENGER, CARPET, VICSEK] {
        let top = if design.dimension == 3 { 5 } else { 7 };
        for level in terms(&design, top, 4) {
            let found = brute::count(&design, level.level, 4);
            assert_eq!(
                level.value, found as i128,
                "{} level {}",
                design.name, level.level
            );
        }
    }
}

#[test]
fn inner_counters_agree() {
    for design in [MENGER, CARPET, VICSEK] {
        for level in 3..=7u32 {
            let span = 3u64.pow(level);
            for modulus in 1..span {
                if modulus % 3 == 0 {
                    continue;
                }
                let all = methods(&design, level, modulus);
                assert!(
                    all.iter().all(|v| *v == all[0]),
                    "{} level {} modulus {} {:?}",
                    design.name,
                    level,
                    modulus,
                    all
                );
            }
        }
    }
}

#[test]
fn pinned_counters_agree() {
    for design in [MENGER, CARPET, VICSEK] {
        let reference: Vec<i128> = terms(&design, 9, 4).iter().map(|l| l.value).collect();
        for (mode, top) in [
            (Mode::Direct, 9),
            (Mode::Zeta, 9),
            (Mode::Convolve, 9),
            (Mode::Bitset, 8),
            (Mode::Rows, 9),
            (Mode::Cube, 9),
        ] {
            let run: Vec<i128> = terms_with(&design, top, 4, mode)
                .iter()
                .map(|l| l.value)
                .collect();
            assert_eq!(run, reference[..top as usize], "{} {:?}", design.name, mode);
        }
    }
}
