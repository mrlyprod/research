mod alphabets;
mod cocycle;
mod growth;
mod order;
mod series;
mod unsolved;
mod word;

use mrly_core::rng::Rng;
use series::{build, product, Rep, Table};
use word::{render, spell, CODES, LIBRARY, SERIES};

const SEED: u64 = 1618033988;

fn line(cells: &[String], widths: &[usize]) -> String {
    cells
        .iter()
        .zip(widths.iter())
        .map(|(cell, width)| {
            let width = *width;
            format!("{cell:>width$}")
        })
        .collect::<Vec<_>>()
        .join("  ")
}

fn matrix_text(matrix: &[Vec<i64>]) -> String {
    let rows: Vec<String> = matrix
        .iter()
        .map(|row| {
            let cells: Vec<String> = row.iter().map(|value| format!("{value}")).collect();
            format!("[{}]", cells.join(","))
        })
        .collect();
    format!("[{}]", rows.join(","))
}

fn codes_text(codes: &[u8]) -> String {
    let cells: Vec<String> = codes.iter().map(|code| format!("{code}")).collect();
    format!("{{{}}}", cells.join(","))
}

fn alphabet() {
    println!("ALPHABET");
    println!("base 2, D = 2, bit i of a code is residue corner i in row-major order, code 3 the top row and code 5 the left column");
    let two = order::sensitivity(&CODES, 2);
    println!(
        "length 2: {} words over the 15 non-empty codes, {} multisets with two or more orderings",
        15 * 15,
        two.total
    );
    let three = order::sensitivity(&LIBRARY, 3);
    println!(
        "length 3: {} words over the library {}, {} multisets with two or more orderings",
        10 * 10 * 10,
        codes_text(&LIBRARY),
        three.total
    );
    println!("the library is every code of fill 2 or 3, the 15 less the four one-cell codes and the full tile");
    let (scanned, hits) = order::library_search();
    println!(
        "the length-3 row 36 188 188 100 is reproduced by {} of the {scanned} ten-code subsets of the 15, namely {}",
        hits.len(),
        hits.iter().map(|hit| codes_text(hit)).collect::<Vec<_>>().join(" ")
    );
    println!();
    println!("ORDER SENSITIVITY");
    let widths = [24usize, 17, 17];
    println!(
        "{}",
        line(
            &[
                "observable".into(),
                "length 2 of 105".into(),
                "length 3 of 210".into()
            ],
            &widths
        )
    );
    let rows: [(&str, usize, usize); 7] = [
        ("fill side density", two.fill, three.fill),
        ("main-diagonal count", two.diagonal, three.diagonal),
        ("boundary", two.boundary, three.boundary),
        ("components", two.components, three.components),
        ("Euler characteristic", two.euler, three.euler),
        ("holes", two.holes, three.holes),
        ("anti-diagonal profile", two.profile, three.profile),
    ];
    for (name, left, right) in rows.iter() {
        println!(
            "{}",
            line(
                &[(*name).into(), format!("{left}"), format!("{right}")],
                &widths
            )
        );
    }
    println!(
        "profile peak sensitive on {} of {} at length 2, profile support on {}",
        two.peak, two.total, two.support
    );
    println!(
        "boundary here is the count of cells with a void or exterior neighbour; the exposed-face count 4N - 2E is a second reading, sensitive on {} of 105 and {} of 210",
        two.perimeter, three.perimeter
    );
    assert_eq!(two.fill, 0, "fill is order-blind");
    assert_eq!(two.diagonal, 0, "the diagonal count is order-blind");
    assert_eq!(three.diagonal, 0, "the diagonal count is order-blind");
    assert_eq!(two.boundary, 0, "boundary is order-blind at length 2");
}

fn witnesses() {
    println!();
    println!("WITNESSES");
    let left = render(&[3, 6]).components();
    let right = render(&[6, 3]).components();
    println!("comp(A_3 (x) A_6) = {left}, comp(A_6 (x) A_3) = {right}");
    assert_eq!((left, right), (4, 2), "the minimal witness stands");
    let pairs = order::pair_component_table(&[3, 5, 10, 12, 6, 9]);
    let adjacent = [3u8, 5, 10, 12];
    let mut same = 0usize;
    let mut mixed = 0usize;
    for (a, b, one, two) in pairs.iter() {
        let split = adjacent.contains(a) != adjacent.contains(b);
        if split {
            mixed += 1;
            assert_eq!(
                (*one, *two),
                (4, 2),
                "adjacent against diagonal is 4 then 2"
            );
        } else {
            same += 1;
            assert_eq!(one, two, "a pair inside one class commutes");
        }
    }
    let twins = render(&[6, 9]).components();
    println!(
        "k = 2 codes: all {same} pairs inside one class commute of {} in all, diagonal against diagonal at {twins}, and all {mixed} adjacent-against-diagonal pairs give 4 against 2",
        pairs.len()
    );
    let big = order::pair_component_table(&[7, 11, 13, 14, 15]);
    let ones = big.iter().filter(|(_, _, a, b)| *a == 1 && *b == 1).count();
    println!(
        "k >= 3 codes: {} of {} pairs give one component in either order",
        ones,
        big.len()
    );
    assert_eq!(ones, big.len(), "every heavy pair is one component");
}

fn blind_laws() {
    println!();
    println!("ORDER-BLIND LAWS");
    let diagonal = order::diagonal_factors();
    println!(
        "the main-diagonal count factors on all {} words of length 3, mismatches {diagonal}",
        15 * 15 * 15
    );
    let contacts = order::contacts_multiply();
    println!(
        "row and column contacts multiply on all {} words of length 3, mismatches {contacts}",
        15 * 15 * 15
    );
    let (pairs, bad) = order::boundary_pairs();
    println!("boundary and interior agree under the factor swap on all {pairs} ordered code pairs, the 15 non-empty codes with the empty one, mismatches {bad}");
    assert_eq!(diagonal, 0, "the diagonal count factors");
    assert_eq!(contacts, 0, "contacts multiply");
    assert_eq!(bad, 0, "boundary is order-blind at length 2");
}

fn representation(table: &mut Table) -> Rep {
    println!();
    println!("RATIONAL SERIES");
    let widths = [24usize, 6, 20];
    println!(
        "{}",
        line(
            &[
                "observable".into(),
                "rank".into(),
                "distinct matrices".into()
            ],
            &widths
        )
    );
    let mut built: Vec<Rep> = Vec::new();
    for (which, &name) in SERIES.iter().enumerate() {
        let rep = build(which, table);
        let classes = rep.classes();
        println!(
            "{}",
            line(
                &[
                    name.into(),
                    format!("{}", rep.basis.len()),
                    format!("{}", classes.len())
                ],
                &widths
            )
        );
        built.push(rep);
    }
    for which in 0..4 {
        let names: Vec<String> = built[which].basis.iter().map(|word| spell(word)).collect();
        println!("{} basis words: {}", SERIES[which], names.join(", "));
    }
    let components = &built[0];
    println!(
        "components lambda = ({}), gamma = ({})^T",
        components
            .lambda
            .iter()
            .map(|v| format!("{v}"))
            .collect::<Vec<_>>()
            .join(","),
        components
            .gamma
            .iter()
            .map(|v| format!("{v}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let classes = components.classes();
    for (members, matrix) in classes.iter() {
        println!(
            "class {} matrix {}",
            codes_text(members),
            matrix_text(matrix)
        );
    }
    let mut commuting: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut pairs = 0usize;
    for i in 0..classes.len() {
        for j in i + 1..classes.len() {
            pairs += 1;
            if product(&classes[i].1, &classes[j].1) == product(&classes[j].1, &classes[i].1) {
                commuting.push((classes[i].0.clone(), classes[j].0.clone()));
            }
        }
    }
    println!(
        "{} of {} class pairs fail to commute, the commuting pair {}",
        pairs - commuting.len(),
        pairs,
        commuting
            .iter()
            .map(|(a, b)| format!("{} and {}", codes_text(a), codes_text(b)))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let light = components.matrices[&1].clone();
    let diagonal = components.matrices[&6].clone();
    let doubled: Vec<Vec<i64>> = light
        .iter()
        .map(|row| row.iter().map(|v| 2 * v).collect())
        .collect();
    println!("M(6,9) = 2 M(1,2,4,8): {}", doubled == diagonal);
    assert_eq!(
        doubled, diagonal,
        "the two zero-contact classes are proportional"
    );
    let mut checked = 0usize;
    let mut bad = 0usize;
    for length in 1..5 {
        for word in order::words(&CODES, length) {
            let obs = word::observe(&word);
            checked += 1;
            for (which, rep) in built.iter().enumerate() {
                if rep.predict(&word) != word::series_value(&obs, which) {
                    bad += 1;
                }
            }
        }
    }
    println!("all {checked} words of length 1 to 4 over the 15 codes, four observables, mismatches {bad}");
    assert_eq!(bad, 0, "the representation is exact to length 4");
    let mut rng = Rng::new(SEED);
    let mut long = 0usize;
    let mut wrong = 0usize;
    for length in 5..8 {
        for _ in 0..40 {
            let word: Vec<u8> = (0..length).map(|_| CODES[rng.below(15)]).collect();
            let obs = word::observe(&word);
            long += 1;
            for (which, rep) in built.iter().enumerate() {
                if rep.predict(&word) != word::series_value(&obs, which) {
                    wrong += 1;
                }
            }
        }
    }
    println!(
        "{long} words of length 5 to 7 drawn at seed {SEED}, four observables, mismatches {wrong}"
    );
    assert_eq!(wrong, 0, "the representation is exact on the drawn words");
    built.remove(0)
}

fn adjudication(rep: &Rep) {
    println!();
    println!("GROWTH");
    let widths = [3usize, 14, 11, 14, 9, 13, 11, 8];
    println!(
        "{}",
        line(
            &[
                "L".into(),
                "(15^(L-1),6)".into(),
                "2*4^(L-1)".into(),
                "(15^(L-1),3)".into(),
                "2^(L-1)".into(),
                "(7^(L-1),6)".into(),
                "2*3^(L-1)".into(),
                "kappa(3)".into(),
            ],
            &widths
        )
    );
    for length in 2..11usize {
        let power = (length - 1) as u32;
        let checker = growth::components(&growth::family(15, 6, length));
        let stripes = growth::components(&growth::family(15, 3, length));
        let gasket = growth::components(&growth::family(7, 6, length));
        let kappa = growth::merging(&growth::family(15, 3, length));
        let expected = [2 * 4u64.pow(power), 2u64.pow(power), 2 * 3u64.pow(power)];
        assert_eq!(checker, expected[0], "the checkerboard family is 2*4^(L-1)");
        assert_eq!(stripes, expected[1], "the striped family is 2^(L-1)");
        assert_eq!(gasket, expected[2], "the gasket family is 2*3^(L-1)");
        assert_eq!(
            kappa, expected[1],
            "the striped family merges on 2^(L-1) components"
        );
        assert_eq!(
            rep.predict(&growth::family(15, 6, length)) as u64,
            checker,
            "the series predicts the checkerboard family"
        );
        assert_eq!(
            rep.predict(&growth::family(15, 3, length)) as u64,
            stripes,
            "the series predicts the striped family"
        );
        assert_eq!(
            rep.predict(&growth::family(7, 6, length)) as u64,
            gasket,
            "the series predicts the gasket family"
        );
        println!(
            "{}",
            line(
                &[
                    format!("{length}"),
                    format!("{checker}"),
                    format!("{}", expected[0]),
                    format!("{stripes}"),
                    format!("{}", expected[1]),
                    format!("{gasket}"),
                    format!("{}", expected[2]),
                    format!("{kappa}"),
                ],
                &widths
            )
        );
    }
    for length in 2..5usize {
        let (best, witness) = growth::max_components(length);
        let bound = growth::independent_bound(length);
        println!(
            "length {length}: the largest component count over all {} words is {best}, attained by {}, against the independent-set bound {bound}",
            15usize.pow(length as u32),
            spell(&witness)
        );
        assert_eq!(best, bound, "the checkerboard bound is attained");
    }
    let mergers: Vec<String> = (1..5usize)
        .map(|length| {
            let (best, witness) = growth::max_merging(length);
            assert_eq!(
                best,
                1 << (length - 1),
                "the striped family maximises kappa"
            );
            format!("{best} at {}", spell(&witness))
        })
        .collect();
    println!(
        "the largest kappa over all words of length 1 to 4 is {}",
        mergers.join(", ")
    );
    let five = growth::max_predicted(rep, 5);
    println!(
        "length 5: the largest predicted component count over all {} words is {five}, against the independent-set bound {}",
        15usize.pow(5),
        growth::independent_bound(5)
    );
    assert_eq!(
        five as u64,
        growth::independent_bound(5),
        "the bound is attained at length 5"
    );
}

fn controls() {
    println!();
    println!("CONTROLS");
    for length in 2..4usize {
        let (checked, bad) = order::crate_agreement(length);
        println!("mrly_math::bang::magic agrees cell for cell on all {checked} words of length {length}, mismatches {bad}");
        assert_eq!(bad, 0, "the factory and the study draw one word");
    }
    let (cases, bad) = order::block_reduction();
    println!("a periodic word equals the self-Kronecker power of its one-period composite on {cases} cases at periods 2 and 3 and lengths to 6, mismatches {bad}");
    assert_eq!(bad, 0, "block reduction holds cell for cell");
}

fn main() {
    println!("MAGIC WORDS");
    println!();
    alphabet();
    witnesses();
    blind_laws();
    let mut table = Table::new();
    let rep = representation(&mut table);
    adjudication(&rep);
    cocycle::study(&rep);
    unsolved::study(&rep);
    alphabets::study();
    controls();
}
