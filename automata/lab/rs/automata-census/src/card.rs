use crate::dynamics::{injective, surjective};
use crate::groups::{group, orbit};
use crate::seed::Census;
use crate::table::{genus, levels};
use mrly_math::bang::universe::degree;
use mrly_math::bang::Code;
use mrly_num::boolean::walsh_spectrum;
use std::collections::BTreeSet;

fn line(rule: usize, census: &Census, b3: &[crate::groups::Elem]) -> String {
    let sigma = levels(rule);
    format!(
        "{rule} surj {} rev {} diagram {} deg {} pop {} genus {} occ {:08b} S {} {} {} {}",
        u8::from(surjective(rule)),
        u8::from(injective(rule)),
        census.class[rule],
        degree(Code::from(rule as u128), 3),
        (rule as u32).count_ones(),
        genus(rule, b3),
        census.occurring[rule],
        sigma[0],
        sigma[1],
        sigma[2],
        sigma[3]
    )
}

pub fn report(census: &Census) {
    println!("THE 110 CARD");
    let b3 = group("B3");
    let h = group("H");
    let big = orbit(110, &b3);
    let small = orbit(110, &h);
    let expected: BTreeSet<usize> = [
        61, 62, 91, 94, 103, 110, 118, 122, 124, 155, 157, 167, 173, 181, 185, 188, 199, 203, 211,
        217, 218, 227, 229, 230,
    ]
    .into_iter()
    .collect();
    assert_eq!(big, expected, "the B3 orbit of 110 is not the expected 24");
    assert!(
        !big.contains(&137) && !big.contains(&193),
        "137 or 193 sits in the B3 orbit of 110"
    );
    assert_eq!(
        small,
        BTreeSet::from([110, 124, 137, 193]),
        "the H class of 110 is not the expected four"
    );
    println!("B3 orbit of 110, {} members", big.len());
    for rule in &big {
        println!("{}", line(*rule, census, &b3));
    }
    println!("H class of 110, 4 members");
    for rule in &small {
        println!("{}", line(*rule, census, &b3));
    }
    let surj: BTreeSet<u8> = big.iter().map(|r| u8::from(surjective(*r))).collect();
    let revs: BTreeSet<u8> = big.iter().map(|r| u8::from(injective(*r))).collect();
    let degs: BTreeSet<i32> = big
        .iter()
        .map(|r| degree(Code::from(*r as u128), 3))
        .collect();
    let pops: BTreeSet<u32> = big.iter().map(|r| (*r as u32).count_ones()).collect();
    let diagrams: BTreeSet<usize> = big.iter().map(|r| census.class[*r]).collect();
    let occs: BTreeSet<u8> = big.iter().map(|r| census.occurring[*r]).collect();
    println!(
        "on the B3 orbit: surjective values {surj:?}, reversible {revs:?}, degree {degs:?}, popcount {pops:?}, diagram classes {}, occurring sets {}",
        diagrams.len(),
        occs.len()
    );
    let hdiagrams: BTreeSet<usize> = small.iter().map(|r| census.class[*r]).collect();
    let hpops: BTreeSet<u32> = small.iter().map(|r| (*r as u32).count_ones()).collect();
    println!(
        "on the H class: diagram classes {}, popcount {hpops:?}",
        hdiagrams.len()
    );
    let sigmas: BTreeSet<[i64; 4]> = big.iter().map(|r| levels(*r)).collect();
    println!(
        "signed Walsh level sums on the B3 orbit: {} distinct vectors",
        sigmas.len()
    );
    for rule in 0..256usize {
        let profile = amplitudes(rule);
        for mate in orbit(rule, &b3) {
            assert_eq!(
                profile,
                amplitudes(mate),
                "the Walsh amplitude profile of {rule} moves at {mate}"
            );
        }
    }
    println!("law: the multiset of |W| at each character weight is a B3 invariant on all 256 rules, while the signed level sums are not");
    println!("constant on the B3 orbit of 110: surjectivity, reversibility, degree, popcount, genus, the Walsh amplitude profile, so none of them separates 110 from a class-mate");
    for rule in [122usize, 218] {
        assert_eq!(
            census.occurring[rule], 0b0011_0111,
            "rule {rule} does not meet exactly the occurring set 00110111"
        );
        assert_eq!(
            census.occurring[rule].count_ones(),
            5,
            "rule {rule} does not meet exactly 5 neighbourhoods"
        );
    }
    for rule in big.iter().filter(|r| **r != 122 && **r != 218) {
        assert_eq!(
            census.occurring[*rule], 0xff,
            "rule {rule} on the B3 orbit of 110 does not meet all 8 neighbourhoods"
        );
    }
    println!("separating on that orbit: the single-seed diagram, distinct on all 24, and the occurring neighbourhood set, which splits 122 and 218, meeting the 5 neighbourhoods 00110111, off from the other 22, which meet all 8");
    println!("popcount is not an H invariant: the H class of 110 carries popcounts 5, 5, 3, 3, because conjugation complements the output");
}

fn amplitudes(rule: usize) -> Vec<(u32, i64)> {
    let mut out: Vec<(u32, i64)> = walsh_spectrum(rule as u128, 3)
        .into_iter()
        .enumerate()
        .map(|(s, w)| (s.count_ones(), w.abs()))
        .collect();
    out.sort();
    out
}
