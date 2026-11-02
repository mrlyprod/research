mod census;
mod collide;
mod counts;
mod tile;
mod witness;

use census::{
    cross_anatomy, rectangle_set, side_eight, side_four, side_nine, side_six, two_radix_lines,
    word_twelve, CrossAnatomy, SideSix,
};
use collide::{collisions, partners, residue_tile, self_power};
use counts::{line_brute, reducible_at, tile_total};
use mrly_math::bang::Code;
use num_bigint::BigUint;
use std::collections::{BTreeMap, BTreeSet};
use tile::{
    chain, cuts, factorisations, irreducible, kron, mask_tile, separable, split, totally_ordered,
    Tile,
};
use witness::{
    antidiagonal, from_cells, identity, line_closure, line_commuting, profiles_text, twelve_sweep,
};

fn tally(map: &BTreeMap<usize, usize>) -> String {
    map.iter()
        .map(|(key, value)| format!("{key}:{value}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn pair_tally(map: &BTreeMap<(usize, usize), usize>) -> String {
    map.iter()
        .map(|((a, b), value)| format!("({a},{b}):{value}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn set_text(list: &[usize]) -> String {
    format!(
        "{{{}}}",
        list.iter()
            .map(|value| format!("{value}"))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn alphabet() {
    println!("CODE FACTORISATION");
    println!("base-2 plane codes 1..15 render 0/1 2 x 2 tiles and base-3 codes 1..511 render 0/1 3 x 3 tiles, bit i of a code the cell (i / base, i mod base)");
    println!("a word folds by the Kronecker product, first letter outermost, so its composite is one 0/1 tile at side the product of the letter sides");
    let mut letters = 0usize;
    for (base, total) in [(2usize, 15u128), (3usize, 511u128)] {
        for code in 1..=total {
            let truth =
                mrly_math::bang::factory::create(Code::from(code), base, 2, base, 1).unwrap();
            let mine = mask_tile(code as u64, base);
            for r in 0..base {
                for c in 0..base {
                    assert_eq!(
                        truth.get(&[r, c]).is_ok_and(|v| v == 1),
                        mine.at(r, c),
                        "letter render"
                    );
                }
            }
            letters += 1;
        }
    }
    let mut products = 0usize;
    for a in 1u128..16 {
        let left = mrly_math::bang::factory::create(Code::from(a), 2, 2, 2, 1).unwrap();
        let mine_left = mask_tile(a as u64, 2);
        for b in 1u128..512 {
            let right = mrly_math::bang::factory::create(Code::from(b), 3, 2, 3, 1).unwrap();
            let truth = left.kron(&right);
            let mine = kron(&mine_left, &mask_tile(b as u64, 3));
            for r in 0..6 {
                for c in 0..6 {
                    assert_eq!(
                        truth.get(&[r, c]).is_ok_and(|v| v == 1),
                        mine.at(r, c),
                        "product render"
                    );
                }
            }
            products += 1;
        }
    }
    println!("{letters} letters and all {products} shape-(2,3) products agree cell for cell with mrly_math::bang::factory::create and mrly_core::Tensor::kron");
    println!();
}

fn block_test() {
    println!("THE BLOCK TEST");
    println!("a tile C of side N cuts at d | N when every nonzero d-block of C is one tile B, and then C = A (x) B with A the 0/1 indicator of the nonzero blocks");
    println!("the test is exact integer comparison, it names both factors outright, and it never touches a rearrangement singular value");
    let four = side_four();
    println!(
        "side 4: {} ordered pairs of the 15 base-2 codes, {} distinct products, largest preimage {}",
        four.products, four.distinct, four.max_preimage
    );
    println!(
        "side 4 by the block test over all 65535 nonempty masks: {} reducible, and the two sets are equal: {}",
        four.by_block.len(),
        four.by_block == four.by_product
    );
    assert_eq!(four.by_block, four.by_product);
    let total = 65535usize;
    println!(
        "side 4: reducible {}, irreducible {}, reducible share {:.6}%",
        four.distinct,
        total - four.distinct,
        100.0 * four.distinct as f64 / total as f64
    );
    println!();
}

fn side_six_report(six: &SideSix, anatomy: &CrossAnatomy) {
    println!("SIDE 6");
    println!(
        "shape (2,3): {} products, {} distinct, {} internal collisions",
        15 * 511,
        six.image23.len(),
        six.collisions23
    );
    println!(
        "shape (3,2): {} products, {} distinct, {} internal collisions",
        511 * 15,
        six.image32.len(),
        six.collisions32
    );
    assert_eq!(six.collisions23, 0);
    assert_eq!(six.collisions32, 0);
    println!("cross-shape, a tile in both images: {}", six.cross.len());
    let universe = tile_total(6);
    let reducible = BigUint::from(six.reducible.len());
    println!(
        "reducible {} + {} - {} = {} of {}, de-duplicated once before printing; irreducible {}",
        six.image23.len(),
        six.image32.len(),
        six.cross.len(),
        six.reducible.len(),
        universe,
        &universe - &reducible
    );
    assert_eq!(
        six.reducible.len(),
        six.image23.len() + six.image32.len() - six.cross.len()
    );
    println!(
        "reducible share {:.7}%",
        100.0 * six.reducible.len() as f64 / 68719476735.0
    );
    println!(
        "the {}: axis-separable {}, not separable {}",
        six.cross.len(),
        anatomy.separable,
        anatomy.non_separable
    );
    println!(
        "the {}: commutations {}, rewritings {}",
        six.cross.len(),
        anatomy.commuting,
        anatomy.rewriting
    );
    println!(
        "not separable and not commuting {}, not separable and commuting {}, the second exactly the diagonal and the antidiagonal",
        anatomy.non_separable_non_commuting,
        anatomy.non_separable_commuting.len()
    );
    println!("fill over the 171: {}", tally(&anatomy.fills));
    println!(
        "fill over the 50 not separable: {}",
        tally(&anatomy.fills_non_separable)
    );
    println!(
        "outer-fill signature of the 50, sorted pair: {}",
        pair_tally(&anatomy.outer_fills)
    );
    println!(
        "a one-cell letter anywhere in a reading, over the 171: {}; over the 50: {}",
        tally(&anatomy.one_cell_any),
        tally(&anatomy.one_cell_any_non_separable)
    );
    println!(
        "a one-cell OUTER factor, over the 171: {}; over the 50: {}",
        tally(&anatomy.one_cell_outer),
        tally(&anatomy.one_cell_outer_non_separable)
    );
    println!(
        "three readings of 48 over the 50, not one statistic: not commuting {}, a one-cell letter in at least one reading {}, a one-cell outer factor in at least one reading {}, and the three sets are equal: {}",
        anatomy.set_non_commuting.len(),
        anatomy.set_one_any.len(),
        anatomy.set_one_outer.len(),
        anatomy.set_non_commuting == anatomy.set_one_any
            && anatomy.set_one_any == anatomy.set_one_outer
    );
    assert_eq!(anatomy.set_non_commuting, anatomy.set_one_any);
    assert_eq!(anatomy.set_one_any, anatomy.set_one_outer);
    println!();
}

fn two_radix(anatomy: &CrossAnatomy) {
    println!("THE TWO-RADIX LINES");
    let lines = two_radix_lines();
    let text: Vec<String> = lines.iter().map(|line| set_text(line)).collect();
    println!(
        "a line of side 6 is a subset of {{0..5}}; {} of the 63 nonempty subsets factor in both radix orders: {}",
        lines.len(),
        text.join(" ")
    );
    let rectangles = rectangle_set();
    println!(
        "the {} axis-separable tiles are exactly the products R x C of two such lines, checked as set equality: {}",
        anatomy.separable,
        rectangles == anatomy.separable_set
    );
    assert_eq!(rectangles, anatomy.separable_set);
    println!(
        "so {} = {} x {} is arithmetic with a checked bijection behind it, never a numeric coincidence",
        anatomy.separable,
        lines.len(),
        lines.len()
    );
    println!();
}

fn commutation(anatomy: &CrossAnatomy) {
    println!("COMMUTATION");
    let pairs: Vec<String> = anatomy
        .commuting_pairs
        .iter()
        .map(|(a, b)| format!("({a},{b})"))
        .collect();
    println!(
        "the {} commuting (base-2, base-3) code pairs: {}",
        anatomy.commuting,
        pairs.join(" ")
    );
    let odd: Vec<String> = anatomy
        .non_separable_commuting
        .iter()
        .map(|(a, b)| format!("c{a} (x) c{b}.q3"))
        .collect();
    println!(
        "9 of them are a commuting row line against a commuting column line, 3 x 3, and the other 2 are {}, the diagonal and the antidiagonal",
        odd.join(" and ")
    );
    let partnered: BTreeSet<u32> = anatomy.commuting_pairs.iter().map(|(a, _)| *a).collect();
    let orphans: Vec<String> = (1u32..16)
        .filter(|code| !partnered.contains(code))
        .map(|code| format!("{code}"))
        .collect();
    println!(
        "base-2 codes {} have no commuting base-3 partner, so the carpet code itself does not commute",
        orphans.join(", ")
    );
    println!("one-cell letters commute exactly when a(n-1) = b(m-1), which gives gcd(m-1,n-1)+1 singleton pairs and, where no common power exists, gcd(m-1,n-1)+2 commuting pairs in one dimension");
    println!("m  n  singletons  gcd+1  commuting  gcd+2");
    for (m, n) in [
        (2usize, 3usize),
        (3, 5),
        (4, 7),
        (5, 9),
        (3, 7),
        (5, 7),
        (4, 5),
        (6, 11),
        (5, 13),
        (3, 9),
    ] {
        let found = line_commuting(m, n);
        let singles = found
            .iter()
            .filter(|(a, b)| a.count_ones() == 1 && b.count_ones() == 1)
            .count();
        let g = tile::gcd(m - 1, n - 1);
        println!("{m}  {n}  {singles}  {}  {}  {}", g + 1, found.len(), g + 2);
    }
    let a = from_cells(3, &[(1, 1)]);
    let b = from_cells(5, &[(2, 2)]);
    let left = kron(&a, &b);
    let right = kron(&b, &a);
    println!(
        "the cells [3]{{(1,1)}} and [5]{{(2,2)}} commute at side 15, both readings giving {}: {}",
        left.text(),
        left == right
    );
    assert_eq!(left, right);
    println!("neither is a corner cell, a full row, a full column, the full tile, a diagonal or an antidiagonal, so the side-6 picture of four corner cells is an artifact of gcd(1,2) = 1");
    println!();
}

fn higher_sides() {
    println!("SIDE 8 AND SIDE 9");
    let eight = side_eight();
    println!(
        "side 8, shape (2,4): {} tiles; shape (4,2): {} tiles; intersection {}",
        eight.image24, eight.image42, eight.intersection
    );
    println!(
        "the intersection is exactly the set of triple products X (x) Y (x) Z of nonempty base-2 codes: {}, so 3375 = 15^3 is associativity and nothing more",
        eight.triples_match
    );
    assert!(eight.triples_match);
    let universe8 = tile_total(8);
    let reducible8 = BigUint::from(eight.image24 + eight.image42 - eight.intersection);
    println!(
        "side 8: reducible {} + {} - {} = {} of {}, irreducible {}",
        eight.image24,
        eight.image42,
        eight.intersection,
        reducible8,
        universe8,
        &universe8 - &reducible8
    );
    let nine = side_nine();
    let universe9 = tile_total(9);
    let reducible9 = BigUint::from(nine.distinct);
    println!(
        "side 9: {} ordered base-3 pairs, {} distinct, {} collisions; reducible {} of {}, irreducible {}",
        nine.products,
        nine.distinct,
        nine.collisions,
        nine.distinct,
        universe9,
        &universe9 - &reducible9
    );
    assert_eq!(nine.collisions, 0);
    println!("3375 never stands beside 171: 8 is a prime power, where no order of the sides can change, and the two counts measure different things");
    println!();
}

fn prime_powers() {
    println!("COUNTING AT PRIME-POWER SIDE");
    println!("at a prime-power side the divisors form a chain, factorisation is unique, and the irreducible series is I = T/(1+T) over the grading");
    for (prime, power) in [(2usize, 2usize), (2, 3), (3, 2), (2, 4), (5, 2)] {
        let side = prime.pow(power as u32);
        println!(
            "side {side} = {prime}^{power}: reducible {}",
            reducible_at(prime, power, true)
        );
    }
    let square = BigUint::from(33554431u64) * BigUint::from(33554431u64);
    println!(
        "side 25 reads (2^25 - 1)^2 = 33554431^2 = {square}, which is not 65535^2 = 4294836225"
    );
    assert_eq!(square, reducible_at(5, 2, true));
    println!("one dimension, the same series against exhaustive brute force:");
    for (prime, power) in [(2usize, 2usize), (2, 3), (2, 4), (3, 2)] {
        let side = prime.pow(power as u32);
        let formula = reducible_at(prime, power, false);
        let brute = BigUint::from(line_brute(side));
        println!(
            "N = {side}: formula {formula}, brute force {brute}, agree {}",
            formula == brute
        );
        assert_eq!(formula, brute);
    }
    println!();
}

fn report_chains(name: &str, tile: &Tile) {
    let words = factorisations(tile);
    let chains: Vec<String> = words.iter().map(|word| set_text(&chain(word))).collect();
    let mut union: BTreeSet<usize> = BTreeSet::new();
    for word in &words {
        union.extend(chain(word));
    }
    let ordered = totally_ordered(&union);
    println!(
        "{name}: cut chains {}, union {} totally ordered by divisibility {}, so the readings share no common refinement",
        chains.join(" and "),
        set_text(&union.iter().copied().collect::<Vec<usize>>()),
        ordered
    );
    assert!(!ordered);
}

fn witnesses(six: &SideSix) {
    println!("WITNESSES");
    let w1 = from_cells(6, &[(0, 0), (2, 2)]);
    let (a2, b3) = split(&w1, 2).unwrap();
    let (x3, y2) = split(&w1, 3).unwrap();
    println!(
        "W1 {}: cut at 2 gives ({} = c{}, {} = c{}.q3); cut at 3 gives ({} = c{}.q3, {} = c{})",
        w1.text(),
        a2.text(),
        a2.pack(),
        b3.text(),
        b3.pack(),
        x3.text(),
        x3.pack(),
        y2.text(),
        y2.pack()
    );
    println!(
        "W1: the four factors are irreducible {} {} {} {}; axis-separable {}; in the cross-shape 171 {}; cut set {}",
        irreducible(&a2),
        irreducible(&b3),
        irreducible(&x3),
        irreducible(&y2),
        separable(&w1),
        six.cross.contains(&w1.pack()),
        set_text(&cuts(&w1))
    );
    assert!(irreducible(&a2) && irreducible(&b3) && irreducible(&x3) && irreducible(&y2));
    assert!(!separable(&w1));
    assert!(six.cross.contains(&w1.pack()));
    report_chains("W1", &w1);
    let w2 = from_cells(12, &[(0, 0), (3, 3)]);
    println!("W2 {}: profiles {}", w2.text(), profiles_text(&w2));
    let inner4 = from_cells(4, &[(0, 0), (3, 3)]);
    let words = factorisations(&w2);
    let lengths: BTreeSet<usize> = words.iter().map(|word| word.len()).collect();
    println!(
        "W2: the side-4 factor {} is irreducible {}; cut set {}, which holds 2 and 3 and not 6, so cut sets are not closed under lcm",
        inner4.text(),
        irreducible(&inner4),
        set_text(&cuts(&w2))
    );
    println!(
        "W2: lengths {:?} and side multisets {{2,2,3}} and {{3,4}}, so neither the length nor the side multiset is an invariant of the tile; axis-separable {}",
        lengths,
        separable(&w2)
    );
    assert!(irreducible(&inner4));
    assert_eq!(lengths.len(), 2);
    assert!(!separable(&w2));
    report_chains("W2", &w2);
    println!("W2 is minimal: every side below 12 is a prime power, where factorisation is unique, or a product of two distinct primes, where every factorisation has prime-side factors only");
    let mut family = 0usize;
    for (m, n) in [
        (2usize, 3usize),
        (3, 2),
        (2, 5),
        (3, 5),
        (4, 3),
        (5, 7),
        (2, 9),
    ] {
        let (im, in_) = (identity(m), identity(n));
        let (em, en) = (antidiagonal(m), antidiagonal(n));
        assert_eq!(kron(&im, &in_), identity(m * n));
        assert_eq!(kron(&in_, &im), identity(m * n));
        assert_eq!(kron(&em, &en), antidiagonal(m * n));
        assert_eq!(kron(&en, &em), antidiagonal(m * n));
        family += 1;
    }
    println!("the infinite family: I_m (x) I_n = I_mn = I_n (x) I_m and E_m (x) E_n = E_mn = E_n (x) E_m, checked at {family} side pairs, so the failure is not a side-6 accident but lives at every side with two distinct prime factors");
    let d2 = from_cells(2, &[(0, 0)]);
    let d3 = from_cells(3, &[(1, 1)]);
    let e2 = from_cells(2, &[(1, 1)]);
    let e3 = from_cells(3, &[(0, 0)]);
    let one = kron(&d2, &d3);
    let two = kron(&e3, &e2);
    println!(
        "not a trace monoid: {} (x) {} = {} = {} (x) {}, four pairwise distinct irreducible letters, so the relation is no commutation of a letter pair: {}",
        d2.text(),
        d3.text(),
        one.text(),
        e3.text(),
        e2.text(),
        one == two
    );
    assert_eq!(one, two);
    println!();
}

fn cut_sets(six: &SideSix) {
    println!("CUT SETS");
    println!("L(C) is the set of d | N at which C cuts; two factorisations share a common refinement exactly when the union of their cut chains is a divisor chain");
    let line = line_closure(20);
    let (side, mask, list) = line.first_lcm.clone().unwrap();
    println!(
        "one dimension, every nonempty subset of {{0..N-1}} at N = 1..20: {} failures of gcd closure, {} of lcm closure",
        line.gcd_violations, line.lcm_violations
    );
    assert_eq!(line.gcd_violations, 0);
    let diagonal = 0b1001u128;
    println!(
        "the first lcm failure is at N = {side} on {} with L = {}, and the line {} that lifts to W2 carries the same cut set {}",
        tile::line_text(mask, side),
        set_text(&list),
        tile::line_text(diagonal, 12),
        set_text(&tile::line_cuts(diagonal, 12))
    );
    assert_eq!(tile::line_cuts(diagonal, 12), list);
    let sweep = twelve_sweep(&six.reducible);
    println!(
        "two dimensions, the {} side-12 composites of three-letter plane-code words counted above: {} failures of gcd closure, {} of lcm closure",
        sweep.tiles, sweep.gcd_violations, sweep.lcm_violations
    );
    assert_eq!(sweep.gcd_violations, 0);
    println!(
        "on the same family the criterion holds with {} mismatches: a tile carries two or more irreducible factorisations exactly when L(C) holds two incomparable divisors",
        sweep.mismatches
    );
    assert_eq!(sweep.mismatches, 0);
    println!(
        "{} of them carry two or more factorisations, {} carry factorisations of unequal length, and the largest number of irreducible factorisations is {}",
        sweep.multiple, sweep.unequal_length, sweep.max_factorisations
    );
    println!("gcd closure of L(C) is Conjecture and is the one missing structural fact; lcm closure is Refuted outright by the N = 12 witness above");
    println!();
}

fn word_census() {
    println!("THE WORD CENSUS AT SIDE 12");
    println!("words over the plane-code alphabet only, a strictly smaller universe than the reducible side-12 tiles, because every plane code has prime side");
    let w2 = from_cells(12, &[(0, 0), (3, 3)]);
    let census = word_twelve(&w2);
    println!(
        "each of the three shapes (2,2,3), (2,3,2), (3,2,2) holds {} words with {} {} {} distinct composites",
        census.per_shape, census.distinct[0], census.distinct[1], census.distinct[2]
    );
    println!(
        "pairwise {} {} {}, triple {}, union {}, and {} = 15 x 171",
        census.pairs[0],
        census.pairs[1],
        census.pairs[2],
        census.triple,
        census.union,
        census.pairs[0]
    );
    assert_eq!(census.pairs[0], 15 * 171);
    println!(
        "W2 lies in the (2,2,3) image {} and in the other two {} {}, because its length-2 reading needs the irreducible side-4 letter [4]{{(0,0),(3,3)}}, which is no plane code",
        census.witness_shapes[0], census.witness_shapes[1], census.witness_shapes[2]
    );
    assert!(census.witness_shapes[0] && !census.witness_shapes[1] && !census.witness_shapes[2]);
    println!("so length and the side multiset are invariants inside the magic-word submonoid, where every letter has prime side, and both fail in the full tile monoid");
    println!();
}

fn render_collisions() {
    println!("RENDER COLLISIONS");
    println!("a base-2 code and a base-3 code render one tile at side n when the two residue rules agree cell for cell; the universe is 15 x 511 = 7665 code pairs per side");
    let mut rows: Vec<String> = Vec::new();
    for side in [2usize, 3, 4, 5, 6, 7, 8, 9, 12, 18] {
        let found = collisions(side);
        rows.push(format!("{side}:{}", found.len()));
        if side >= 4 {
            assert_eq!(found, vec![(15u32, 511u32)]);
        }
    }
    println!("collisions by side: {}", rows.join(" "));
    println!("480 = 15 x 2^5 at side 2, one per nonempty base-2 code at side 3, and only the full tile (15, 511) at every side from 4 up");
    let match_list = partners(7, 3);
    let text: Vec<String> = match_list.iter().map(|code| format!("{code}")).collect();
    println!(
        "the carpet code 7 has match list [{}] at side 3, a unique partner of fill {}",
        text.join(","),
        residue_tile(495, 3, 3).fill()
    );
    assert_eq!(match_list, vec![495u32]);
    let carpet9 = residue_tile(7, 2, 9);
    let residue9 = residue_tile(495, 3, 9);
    let power9 = self_power(495, 3, 2);
    println!(
        "at side 9 the three readings separate: carpet residue fill {}, c495 residue fill {}, c495 self-power fill {}, and the three tiles are pairwise distinct: {}",
        carpet9.fill(),
        residue9.fill(),
        power9.fill(),
        carpet9 != residue9 && residue9 != power9 && carpet9 != power9
    );
    assert!(carpet9 != residue9 && residue9 != power9 && carpet9 != power9);
    println!("the two self-power ladders 2, 4, 8, 16 and 3, 9, 27 share no side, so the pure fractal reading carries no collision at all");
    println!();
}

fn main() {
    alphabet();
    block_test();
    let six = side_six();
    let anatomy = cross_anatomy(&six);
    side_six_report(&six, &anatomy);
    two_radix(&anatomy);
    commutation(&anatomy);
    higher_sides();
    prime_powers();
    witnesses(&six);
    word_census();
    cut_sets(&six);
    render_collisions();
}
