use super::{Axis, Key, Measure, Tag};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Record {
    pub id: &'static str,
    pub name: &'static str,
    pub offset: i32,
    pub shift: i32,
    pub terms: &'static str,
    pub status: Tag,
    pub key: Option<Key>,
    pub formula: &'static str,
    pub witness: &'static str,
}

const fn record(
    id: &'static str,
    name: &'static str,
    offset: i32,
    terms: &'static str,
    status: Tag,
    formula: &'static str,
    witness: &'static str,
) -> Record {
    Record {
        id,
        name,
        offset,
        shift: 0,
        terms,
        status,
        key: None,
        formula,
        witness,
    }
}

impl Record {
    const fn keyed(
        self,
        code: u128,
        dimension: usize,
        measure: Measure,
        axis: Axis,
        shift: i32,
    ) -> Record {
        Record {
            shift,
            key: Some(Key::new(code, dimension, 2, measure, axis)),
            ..self
        }
    }
}

const FILLS: &str = "mrly_math::counts::fill";
const SIDES: &str = "ledger::terms, the odd-side law";
const EXPOSURE: &str = "mrly_math::counts::exposure";

pub static RECORDS: &[Record] = &[
    record("A000029", "Number of necklaces with n beads of 2 colors, allowing turning over (these are also called bracelets)", 0, "1, 2, 3, 4, 6, 8, 13, 18, 30, 46, 78, 126", Tag::Verified, "none", "mrly_math::bang::baseq::bracelets, the base line at dim 1"),
    record("A000070", "a(n) = Sum_{k=0..n} p(k) where p(k) = number of partitions of k", 0, "1, 2, 4, 7, 12, 19, 30, 45, 67, 97, 139, 195", Tag::Verified, "none", "REFS.md"),
    record("A000244", "Powers of 3: a(n) = 3^n", 0, "1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147", Tag::Proved, "3^level", "mrly_math::three::diagonal, every admissible cut of bang dim 3, code 126 at side 2"),
    record("A000290", "The squares: a(n) = n^2", 0, "0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121", Tag::Proved, "k^2", SIDES).keyed(1, 2, Measure::Fills, Axis::Side, 0),
    record("A000351", "Powers of 5: a(n) = 5^n", 0, "1, 5, 25, 125, 625, 3125, 15625, 78125, 390625, 1953125, 9765625, 48828125", Tag::Proved, "5^level", FILLS).keyed(9, 2, Measure::Fills, Axis::Level, 0),
    record("A000370", "Number of NPN-equivalence classes of Boolean functions of n or fewer variables", 0, "1, 2, 4, 14, 222, 616126, 200253952527184", Tag::Verified, "none", "bijection.md, the NPN sibling at dim 1..4"),
    record("A000384", "Hexagonal numbers: a(n) = n*(2*n-1)", 0, "0, 1, 6, 15, 28, 45, 66, 91, 120, 153, 190, 231", Tag::Proved, "k(2k - 1)", SIDES).keyed(3, 2, Measure::Fills, Axis::Side, 0),
    record("A000420", "Powers of 7: a(n) = 7^n", 0, "1, 7, 49, 343, 2401, 16807, 117649, 823543, 5764801, 40353607, 282475249, 1977326743", Tag::Proved, "7^level", FILLS).keyed(11, 2, Measure::Fills, Axis::Level, 0),
    record("A000567", "Octagonal numbers: n*(3*n-2)", 0, "0, 1, 8, 21, 40, 65, 96, 133, 176, 225, 280, 341", Tag::Proved, "k(3k - 2)", SIDES).keyed(7, 2, Measure::Fills, Axis::Side, 0),
    record("A000578", "The cubes: a(n) = n^3", 0, "0, 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 1331", Tag::Proved, "k^3", SIDES).keyed(1, 3, Measure::Fills, Axis::Side, 0),
    record("A000616", "a(-1)=1 by convention; for n >= 0, a(n) = number of irreducible Boolean functions of n variables", -1, "1, 2, 3, 6, 22, 402, 1228158, 400507806843728", Tag::Proved, "none", "mrly_math::bang::baseq::sequence, bijection.md"),
    record("A001018", "Powers of 8: a(n) = 8^n", 0, "1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728, 1073741824, 8589934592", Tag::Proved, "8^level", FILLS).keyed(7, 2, Measure::Fills, Axis::Level, 0),
    record("A001024", "Powers of 15: a(n) = 15^n", 0, "1, 15, 225, 3375, 50625, 759375, 11390625, 170859375, 2562890625, 38443359375", Tag::Verified, "15^n", "REFS.md"),
    record("A001316", "Gould's sequence: number of odd entries in row n of Pascal's triangle", 0, "1, 2, 2, 4, 2, 4, 4, 8, 2, 4, 4, 8", Tag::Verified, "2^wt(n)", "REFS.md"),
    record("A001481", "Numbers that are the sum of 2 squares", 1, "0, 1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18", Tag::Verified, "none", "spin.md, the ring radii of a spun square lattice"),
    record("A001844", "Centered square numbers: a(n) = 2*n*(n+1)+1", 0, "1, 5, 13, 25, 41, 61, 85, 113, 145, 181, 221, 265", Tag::Proved, "2k^2 - 2k + 1", SIDES).keyed(9, 2, Measure::Fills, Axis::Side, -1),
    record("A002407", "Cuban primes: primes which are the difference of two consecutive cubes", 1, "7, 19, 37, 61, 127, 271, 331, 397, 547, 631, 919, 1657", Tag::Verified, "none", "mrly_math::counts::six theorems, slices.md"),
    record("A003136", "Loeschian numbers: numbers of the form x^2 + xy + y^2", 1, "0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, 25", Tag::Verified, "none", "spin.md, the ring radii of a spun hexagonal lattice"),
    record("A003180", "Number of equivalence classes of Boolean functions of n variables under action of symmetric group", 0, "2, 4, 12, 80, 3984, 37333248, 25626412338274304", Tag::Verified, "none", "bijection.md, the axis-permutation orbits at a one-term shift"),
    record("A003215", "Hex (or centered hexagonal) numbers: 3*n*(n+1)+1", 0, "1, 7, 19, 37, 61, 91, 127, 169, 217, 271, 331, 397", Tag::Proved, "3k^2 - 3k + 1", SIDES).keyed(11, 2, Measure::Fills, Axis::Side, -1),
    record("A003463", "a(n) = (5^n - 1)/4", 0, "0, 1, 6, 31, 156, 781, 3906, 19531, 97656, 488281, 2441406, 12207031", Tag::Verified, "(5^n - 1)/4", "REFS.md"),
    record("A004016", "Theta series of planar hexagonal lattice A_2", 0, "1, 6, 0, 6, 6, 0, 0, 12, 0, 6, 0, 0", Tag::Verified, "none", "mrly_num::gauss::shells, spin.md"),
    record("A004018", "Theta series of square lattice: number of ways of writing n as a sum of 2 squares", 0, "1, 4, 4, 0, 4, 8, 0, 0, 4, 4, 8, 0", Tag::Verified, "4(d_1(n) - d_3(n))", "mrly_num::gauss::shells, spin.md"),
    record("A004662", "Powers of 3 written in base 8", 0, "1, 3, 11, 33, 121, 363, 1331, 4213, 14641, 46343, 163251, 531773", Tag::Verified, "none", "this page, a near miss of A396934 - 1"),
    record("A005418", "Number of (n-1)-bead black-white reversible strings", 1, "1, 2, 3, 6, 10, 20, 36, 72, 136, 272, 528, 1056", Tag::Verified, "none", "REFS.md"),
    record("A005728", "Number of fractions in Farey series of order n", 0, "1, 2, 3, 5, 7, 11, 13, 19, 23, 29, 33, 43", Tag::Verified, "1 + sum of phi(k) for k <= n", "mrly_num::lattice::farey, the lit nodes of the stack in farey.md"),
    record("A005898", "Centered cube numbers: n^3 + (n+1)^3", 0, "1, 9, 35, 91, 189, 341, 559, 855, 1241, 1729, 2331, 3059", Tag::Proved, "k^3 + (k - 1)^3", SIDES).keyed(129, 3, Measure::Fills, Axis::Side, -1),
    record("A009964", "Powers of 20", 0, "1, 20, 400, 8000, 160000, 3200000, 64000000, 1280000000, 25600000000, 512000000000", Tag::Proved, "20^level", FILLS).keyed(23, 3, Measure::Fills, Axis::Level, 0),
    record("A009971", "Powers of 27", 0, "1, 27, 729, 19683, 531441, 14348907, 387420489, 10460353203, 282429536481, 7625597484987", Tag::Proved, "27^level", FILLS).keyed(255, 3, Measure::Fills, Axis::Level, 0),
    record("A011934", "a(n) = abs(1^3 - 2^3 + 3^3 - 4^3 + ... + (-1)^(n+1)*n^3)", 0, "0, 1, 7, 20, 44, 81, 135, 208, 304, 425, 575, 756", Tag::Verified, "none", "this page, the parent of the bisections A103532 and A395241"),
    record("A014105", "Second hexagonal numbers: a(n) = n*(2*n + 1)", 0, "0, 3, 10, 21, 36, 55, 78, 105, 136, 171, 210, 253", Tag::Proved, "2k^2 - 3k + 1", SIDES).keyed(3, 2, Measure::Voids, Axis::Side, -1),
    record("A016185", "a(n) = 9^n - 8^n", 0, "0, 1, 17, 217, 2465, 26281, 269297, 2685817, 26269505, 253202761, 2413042577, 22791125017", Tag::Proved, "9^level - 8^level", "mrly_math::counts::void").keyed(7, 2, Measure::Voids, Axis::Level, 0),
    record("A016754", "Odd squares: a(n) = (2n+1)^2, also centered octagonal numbers", 0, "1, 9, 25, 49, 81, 121, 169, 225, 289, 361, 441, 529", Tag::Proved, "(2k - 1)^2", SIDES).keyed(15, 2, Measure::Fills, Axis::Side, -1),
    record("A016755", "Odd cubes: a(n) = (2*n + 1)^3", 0, "1, 27, 125, 343, 729, 1331, 2197, 3375, 4913, 6859, 9261, 12167", Tag::Proved, "(2k - 1)^3", SIDES).keyed(255, 3, Measure::Fills, Axis::Side, -1),
    record("A018413", "Divisors of 363", 1, "1, 3, 11, 33, 121, 363", Tag::Verified, "none", "this page, a near miss of A396934 - 1"),
    record("A034474", "a(n) = 5^n + 1", 0, "2, 6, 26, 126, 626, 3126, 15626, 78126, 390626, 1953126, 9765626, 48828126", Tag::Verified, "5^n + 1", "REFS.md"),
    record("A045944", "Rhombic matchstick numbers: a(n) = n*(3*n+2)", 0, "0, 5, 16, 33, 56, 85, 120, 161, 208, 261, 320, 385", Tag::Proved, "3k^2 - 4k + 1", SIDES).keyed(1, 2, Measure::Voids, Axis::Side, -1),
    record("A046092", "4 times triangular numbers: a(n) = 2*n*(n+1)", 0, "0, 4, 12, 24, 40, 60, 84, 112, 144, 180, 220, 264", Tag::Proved, "2k^2 - 2k", SIDES).keyed(6, 2, Measure::Fills, Axis::Side, -1),
    record("A047999", "Sierpinski's triangle (or gasket): Pascal's triangle read by rows mod 2", 0, "1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0", Tag::Verified, "none", "coprime.md, the points i AND j = 0"),
    record("A048883", "a(n) = 3^wt(n), where wt(n) = A000120(n)", 0, "1, 3, 3, 9, 3, 9, 9, 27, 3, 9, 9, 27", Tag::Verified, "3^wt(n)", "REFS.md"),
    record("A054247", "Number of n X n binary matrices under action of dihedral group of the square D_4", 0, "1, 2, 6, 102, 8548, 4211744, 8590557312, 70368882591744, 2305843028004192256", Tag::Verified, "none", "README.md, the rigid hypercube census at dim 2"),
    record("A065473", "Decimal expansion of the strongly carefree constant: Product_{p prime} (1 - (3*p-2)/(p^3))", 0, "2, 8, 6, 7, 4, 7, 4, 2, 8, 4, 3, 4", Tag::Verified, "none", "coprime.md"),
    record("A069403", "a(n) = 2*Fibonacci(2*n+1) - 1", 0, "1, 3, 9, 25, 67, 177, 465, 1219, 3193, 8361, 21891, 57313", Tag::Verified, "2 F(2n + 1) - 1", "REFS.md"),
    record("A103532", "Number of divisors of 240^n", 0, "1, 20, 81, 208, 425, 756, 1225, 1856, 2673, 3700, 4961, 6480", Tag::Proved, "(4k - 3)k^2", SIDES).keyed(23, 3, Measure::Fills, Axis::Side, -1),
    record("A125833", "Numbers whose base-5 representation is 333333.......3", 0, "0, 3, 18, 93, 468, 2343, 11718, 58593, 292968, 1464843, 7324218, 36621093", Tag::Verified, "3(5^n - 1)/4", "REFS.md"),
    record("A128625", "Expansion of (1+3*x)/(1-5*x)", 0, "1, 8, 40, 200, 1000, 5000, 25000, 125000, 625000, 3125000, 15625000, 78125000", Tag::Verified, "8 5^(n-1) for n >= 1", "REFS.md"),
    record("A129824", "a(n) = Product_{k=0..n} (1 + binomial(n,k))", 0, "2, 4, 12, 64, 700, 17424, 1053696, 160579584, 62856336636, 63812936890000, 168895157342195152, 1169048914836855865344", Tag::Proved, "Prod_{k=0}^{n} (1 + C(n,k))", "lab/rs/design-census, the fill classes of base 2"),
    record("A141148", "Number of aperiodic ternary necklaces with n beads of each color and no adjacent beads of the same color", 1, "2, 3, 14, 65, 346, 1929, 11442, 70310, 445928, 2896239, 19186738, 129184583", Tag::Verified, "none", "REFS.md"),
    record("A154105", "a(n) = 12*n^2 + 18*n + 7", 0, "7, 37, 91, 169, 271, 397, 547, 721, 919, 1141, 1387, 1657", Tag::Verified, "12n^2 + 18n + 7", "slices.md, the centered hexagonal vertices at n = k - 1"),
    record("A192908", "Constant term in the reduction by (x^2 -> x + 1) of a polynomial family; a(n) = 2*Fibonacci(2n-2) + 1", 0, "1, 1, 3, 7, 17, 43, 111, 289, 755, 1975, 5169, 13531", Tag::Verified, "2 F(2n - 2) + 1", "REFS.md"),
    record("A229896", "Sizes of logical groups of the same integer in A229895", 1, "1, 1, 4, 1, 5, 27, 1, 7, 37, 256, 1, 9, 61, 369, 3125, 1, 11, 91, 671, 4651, 46656, 1, 13, 127, 1105, 9031, 70993, 823543, 1, 15, 169, 1695, 15961, 144495, 1273609, 16777216, 1, 17, 217, 2465, 26281, 269297, 2685817, 26269505, 387420489", Tag::Verified, "none", "this page, an interior window holding the carpet voids"),
    record("A255016", "Number of toroidal n X n binary arrays, allowing rotation and/or reflection of rows and/or columns as well as matrix transposition", 0, "1, 2, 6, 26, 805, 172112, 239123150, 1436120190288, 36028817512382026", Tag::Verified, "none", "mrly_math::bang::baseq::sequence at dim 2, bijection.md"),
    record("A268240", "Pascal's tetrahedron of trinomial coefficients read mod 2", 0, "1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1", Tag::Verified, "none", "REFS.md"),
    record("A299916", "a(n) = A299914(2n+1); the six-pointed-star holes of the Menger slice, by a comment", 0, "1, 6, 42, 306, 2250, 16578, 122202, 900882, 6641514, 48963042, 360969210, 2661166386", Tag::Verified, "a(n) = 9a(n-1) - 12a(n-2)", "mrly_math::counts::cut_fills, slices.md").keyed(23, 3, Measure::Triangles, Axis::Level, 1),
    record("A332705", "Number of unit square faces (or surface area) of a stage-n Menger sponge", 0, "6, 72, 1056, 18048, 336384, 6531072, 129048576, 2568388608, 51267108864, 1024536870912", Tag::Proved, "2*20^level + 4*8^level", EXPOSURE).keyed(23, 3, Measure::Surface, Axis::Level, 0),
    record("A347825", "Number of ways to cut a 2 X n rectangle into rectangles with integer sides up to symmetries of the rectangle", 0, "1, 2, 6, 17, 61, 220, 883, 3597, 15232, 65130, 282294, 1229729", Tag::Verified, "none", "this page, a near miss of A396934/2"),
    record("A361870", "Array read by downward antidiagonals: nonequivalent 2-colorings of the cells of an n-dimensional hypercube with edges k cells long", 0, "2, 2, 1, 2, 2, 1, 2, 3, 2, 1, 2, 6", Tag::Verified, "none", "README.md, the rigid hypercube census"),
    record("A381517", "Perimeter of the Sierpinski carpet at iteration n", 0, "4, 16, 80, 496, 3536, 26992, 212048, 1684720, 13442768, 107437168, 859182416, 6872514544", Tag::Proved, "(4*8^level + 16*3^level)/5", EXPOSURE).keyed(7, 2, Measure::Surface, Axis::Level, 0),
    record("A395134", "Decimal expansion of the probability that the line that passes through two points selected independently and uniformly at random in a half-disk intersects the arc at two points.", 0, "4, 5, 9, 6, 2, 0, 3, 5, 3, 9, 0, 7", Tag::Verified, "1 - 16/(3 Pi^2)", "coprime.md, the complement of the A396934 density"),
    record("A395241", "a(n) = n^2*(4*n + 3)", 0, "0, 7, 44, 135, 304, 575, 972, 1519, 2240, 3159, 4300, 5687", Tag::Verified, "a(n) = n^2*(4*n + 3)", "lab/rs/oeis-terms, the b-file to n = 10000").keyed(23, 3, Measure::Voids, Axis::Side, -1),
    record("A396922", "E.g.f. A(x) satisfies A( x / A(log(A(log(A(log(A(x))))))) ) = exp(x)", 0, "1, 1, 3, 40, 1421, 87896, 7921207, 951512332, 144407735033, 26715045346048", Tag::Verified, "none", "REFS.md"),
    record("A396934", "Number of pairs (i,j) with 0 <= i,j < 2^n, i AND j = 0, and gcd(i,j) = 1", 0, "0, 2, 4, 12, 34, 122, 362, 1130, 3406, 10506, 31550, 95260", Tag::Verified, "none", "lab/rs/oeis-terms, the b-file to n = 20"),
    record("A398348", "Number of toroidal n X n X n binary arrays, allowing rotation and/or reflection of the layers along each axis as well as all permutations of the axes", 1, "2, 22, 111618, 6005363762644688, 7089215977519836239803174210135872, 10157435539019790383692007859901914095646506996125324171134976", Tag::Verified, "none", "lab/rs/oeis-terms, the b-file to n = 14"),
    record("A399018", "Number of solid unit triangles in the hexagonal cross-section of the first iteration of the generalized Menger sponge with subdivision 2n+1, cut through the center perpendicular to a space diagonal", 0, "6, 42, 72, 204, 210, 486, 420, 888, 702, 1410, 1056, 2052", Tag::Verified, "3(n+1)(3n+2) for n even, 3(n+1)(5n+2) for n odd", "mrly_math::counts::cut_fills").keyed(23, 3, Measure::Triangles, Axis::Side, -1),
    record("A399019", "Number of empty unit triangles in the hexagonal cross-section of the first iteration of the generalized Menger sponge with subdivision 2n+1, cut through the center perpendicular to a space diagonal", 0, "0, 12, 78, 90, 276, 240, 594, 462, 1032, 756, 1590, 1122", Tag::Verified, "3n(5n+3) for n even, 3n(3n+1) for n odd", "mrly_math::counts::cut_voids").keyed(23, 3, Measure::Gaps, Axis::Side, -1),
];

pub fn record_by_id(id: &str) -> Option<&'static Record> {
    RECORDS.iter().find(|record| record.id == id)
}

fn parsed() -> &'static [Vec<i128>] {
    static PARSED: OnceLock<Vec<Vec<i128>>> = OnceLock::new();
    PARSED.get_or_init(|| {
        RECORDS
            .iter()
            .map(|record| {
                record
                    .terms
                    .split(',')
                    .map_while(|term| term.trim().parse().ok())
                    .collect()
            })
            .collect()
    })
}

impl Record {
    pub fn parsed(&self) -> &'static [i128] {
        let at = RECORDS
            .iter()
            .position(|record| record.id == self.id)
            .expect("every record is listed");
        &parsed()[at]
    }
}

pub fn identify(terms: &[i128]) -> Vec<(&'static Record, i32)> {
    if terms.is_empty() {
        return Vec::new();
    }
    let mut found: Vec<(&'static Record, i32)> = RECORDS
        .iter()
        .zip(parsed())
        .filter_map(|(record, known)| {
            known
                .windows(terms.len())
                .position(|window| window == terms)
                .map(|at| (record, record.offset + at as i32))
        })
        .collect();
    found.sort_by_key(|(record, shift)| (record.key.is_none(), shift.abs(), record.id));
    found
}

#[cfg(test)]
mod tests {
    use super::super::{terms, BUDGET, TERMS};
    use super::*;

    #[test]
    fn the_records_are_sorted_distinct_and_parse() {
        for pair in RECORDS.windows(2) {
            assert!(
                pair[0].id < pair[1].id,
                "{} before {}",
                pair[0].id,
                pair[1].id
            );
        }
        for record in RECORDS {
            assert!(record.parsed().len() >= 5, "{}", record.id);
            assert!(!record.name.contains('|'), "{}", record.id);
        }
        assert_eq!(RECORDS.len(), 65);
        assert_eq!(record_by_id("A398348").unwrap().parsed().len(), 5);
    }

    #[test]
    fn every_keyed_record_reads_its_generator_at_its_shift() {
        let mut keyed = 0;
        for record in RECORDS.iter().filter(|record| record.key.is_some()) {
            let key = record.key.unwrap();
            let (got, _) = terms(&key, TERMS, BUDGET).unwrap();
            let known = record.parsed();
            let mut compared = 0;
            for (index, &term) in got.iter().enumerate() {
                let at = key.axis.start() + index as i32 + record.shift - record.offset;
                if let Some(&expected) = usize::try_from(at).ok().and_then(|at| known.get(at)) {
                    assert_eq!(term, expected, "{} at {index}", record.id);
                    compared += 1;
                }
            }
            assert!(compared >= 3, "{} compared {compared}", record.id);
            keyed += 1;
        }
        assert_eq!(keyed, 25);
    }

    #[test]
    fn every_design_family_reads_its_record_by_value() {
        use mrly_life::source::{sequence, Source};
        for family in Source::designs() {
            let id = family.oeis().expect("every design family carries an id");
            let known = record_by_id(&id)
                .expect("the family's id is a record")
                .parsed();
            let limit = usize::try_from(*known.last().unwrap()).unwrap();
            let values: Vec<i128> = sequence(family, limit)
                .unwrap()
                .into_iter()
                .map(|value| value as i128)
                .collect();
            assert!(values.len() >= 5, "{family:?}");
            assert_eq!(
                identify(&values).iter().filter(|(r, _)| r.id == id).count(),
                1,
                "{family:?}"
            );
        }
    }

    #[test]
    fn identify_returns_the_shift() {
        let ids = |terms: &[i128]| -> Vec<(&str, i32)> {
            identify(terms).iter().map(|(r, s)| (r.id, *s)).collect()
        };
        assert_eq!(ids(&[8, 21, 40, 65]), [("A000567", 2)]);
        assert_eq!(ids(&[6, 15, 28]), [("A000384", 2)]);
        assert_eq!(ids(&[72, 1056, 18048]), [("A332705", 1)]);
        assert_eq!(ids(&[1, 17, 217]), [("A016185", 1), ("A229896", 37)]);
        assert_eq!(ids(&[3, 6, 22, 402]), [("A000616", 1)]);
        assert_eq!(ids(&[4, 12, 80, 3984]), [("A003180", 1)]);
        assert_eq!(ids(&[5, 7, 11, 13]), [("A005728", 3)]);
        assert!(ids(&[5, 7, 11, 14]).is_empty());
        assert!(ids(&[]).is_empty());
    }
}
