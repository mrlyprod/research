use crate::order;
use crate::series::Rep;
use crate::word::{render, CODES};

pub const UNIT: [u8; 4] = [1, 2, 4, 8];
pub const ROW: [u8; 2] = [3, 12];
pub const COLUMN: [u8; 2] = [5, 10];
pub const DIAGONAL: [u8; 2] = [6, 9];
pub const GASKET: [u8; 4] = [7, 11, 13, 14];
pub const FULL: u8 = 15;
pub const REACH: usize = 14;
pub const DRAWN: usize = 7;
pub const HORIZON: usize = 1 << 20;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Form {
    UnitDomino,
    UnitDiagonal,
    RowColumn,
    DominoDiagonal,
    DominoFull,
    Empty,
    All,
}

pub struct Family {
    pub name: &'static str,
    pub rule: &'static str,
    pub form: Form,
    pub pairs: Vec<(u8, u8)>,
}

fn gcd(a: u64, b: u64) -> u64 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 {
        1
    } else {
        a
    }
}

fn ratio(num: u64, den: u64) -> String {
    let g = gcd(num, den);
    if den / g == 1 {
        return format!("{}", num / g);
    }
    format!("{}/{}", num / g, den / g)
}

pub fn value(rep: &Rep, word: &[u8]) -> i128 {
    let mut state: Vec<i128> = rep.lambda.iter().map(|entry| *entry as i128).collect();
    for code in word {
        let matrix = &rep.matrices[code];
        let mut next = vec![0i128; state.len()];
        for (i, weight) in state.iter().enumerate() {
            if *weight == 0 {
                continue;
            }
            for (j, slot) in next.iter_mut().enumerate() {
                *slot += weight * matrix[i][j] as i128;
            }
        }
        state = next;
    }
    state
        .iter()
        .zip(rep.gamma.iter())
        .map(|(a, b)| a * *b as i128)
        .sum()
}

pub fn exponent(form: Form, word: &[u8], marked: u8) -> u32 {
    match form {
        Form::UnitDomino => {
            let mut last = 0usize;
            let mut count = 0usize;
            for (i, code) in word.iter().enumerate() {
                if *code == marked {
                    last = i + 1;
                    count += 1;
                }
            }
            (last - count) as u32
        }
        Form::UnitDiagonal => word.iter().filter(|code| **code == marked).count() as u32,
        Form::RowColumn => {
            let last = *word.last().expect("a word has a last letter");
            let run = word.iter().rev().take_while(|code| **code == last).count();
            (word.len() - run) as u32
        }
        Form::DominoDiagonal => {
            let mut last = 0usize;
            for (i, code) in word.iter().enumerate() {
                if *code == marked {
                    last = i + 1;
                }
            }
            last as u32
        }
        Form::DominoFull => {
            let count = word.iter().filter(|code| **code == marked).count();
            let run = word
                .iter()
                .rev()
                .take_while(|code| **code == marked)
                .count();
            (count - run) as u32
        }
        Form::Empty => 0,
        Form::All => word.len() as u32,
    }
}

fn pair_word(a: u8, b: u8, mask: usize, length: usize) -> Vec<u8> {
    (0..length)
        .map(|i| if (mask >> i) & 1 == 1 { b } else { a })
        .collect()
}

fn sweep(rep: &Rep, form: Form, a: u8, b: u8, reach: usize) -> (usize, usize) {
    let mut checked = 0usize;
    let mut bad = 0usize;
    for length in 1..=reach {
        for mask in 0..(1usize << length) {
            let word = pair_word(a, b, mask, length);
            let want = 1i128 << exponent(form, &word, b);
            if value(rep, &word) != want {
                bad += 1;
            }
            checked += 1;
        }
    }
    (checked, bad)
}

fn drawn(form: Form, a: u8, b: u8, reach: usize) -> (usize, usize) {
    let mut checked = 0usize;
    let mut bad = 0usize;
    for length in 1..=reach {
        for mask in 0..(1usize << length) {
            let word = pair_word(a, b, mask, length);
            let want = 1u64 << exponent(form, &word, b);
            if render(&word).components() != want {
                bad += 1;
            }
            checked += 1;
        }
    }
    (checked, bad)
}

pub fn families() -> Vec<Family> {
    let dominoes: Vec<u8> = ROW.iter().chain(COLUMN.iter()).copied().collect();
    let mut out: Vec<Family> = Vec::new();
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    for &domino in dominoes.iter() {
        for &unit in UNIT.iter() {
            pairs.push((domino, unit));
        }
    }
    out.push(Family {
        name: "unit and domino",
        rule: "2^(k - m), k the last unit place, m the unit count",
        form: Form::UnitDomino,
        pairs,
    });
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    for &unit in UNIT.iter() {
        for &diagonal in DIAGONAL.iter() {
            pairs.push((unit, diagonal));
        }
    }
    out.push(Family {
        name: "unit and diagonal",
        rule: "2^(diagonal count)",
        form: Form::UnitDiagonal,
        pairs,
    });
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    for &row in ROW.iter() {
        for &column in COLUMN.iter() {
            pairs.push((row, column));
        }
    }
    out.push(Family {
        name: "crossed dominoes",
        rule: "2^(L - r), r the terminal run",
        form: Form::RowColumn,
        pairs,
    });
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    for &domino in dominoes.iter() {
        for &diagonal in DIAGONAL.iter() {
            pairs.push((domino, diagonal));
        }
    }
    out.push(Family {
        name: "domino and diagonal",
        rule: "2^k, k the last diagonal place",
        form: Form::DominoDiagonal,
        pairs,
    });
    let pairs: Vec<(u8, u8)> = dominoes.iter().map(|domino| (*domino, FULL)).collect();
    out.push(Family {
        name: "domino and full",
        rule: "2^(n - j), n the full count, j the terminal full run",
        form: Form::DominoFull,
        pairs,
    });
    let pairs: Vec<(u8, u8)> = GASKET.iter().map(|code| (*code, FULL)).collect();
    out.push(Family {
        name: "gasket and full",
        rule: "1",
        form: Form::Empty,
        pairs,
    });
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    for class in [
        UNIT.to_vec(),
        ROW.to_vec(),
        COLUMN.to_vec(),
        GASKET.to_vec(),
    ]
    .iter()
    {
        for i in 0..class.len() {
            for j in i + 1..class.len() {
                pairs.push((class[i], class[j]));
            }
        }
    }
    out.push(Family {
        name: "inside one flat class",
        rule: "1",
        form: Form::Empty,
        pairs,
    });
    out.push(Family {
        name: "inside the diagonal class",
        rule: "2^L",
        form: Form::All,
        pairs: vec![(DIAGONAL[0], DIAGONAL[1])],
    });
    out
}

fn constants(rep: &Rep) {
    println!("CONSTANT WORDS");
    let mut doubling: Vec<u8> = Vec::new();
    let mut bad = 0usize;
    for &code in CODES.iter() {
        let one = value(rep, &[code]);
        for length in 1..8usize {
            let word = vec![code; length];
            if value(rep, &word) != one.pow(length as u32) {
                bad += 1;
            }
        }
        if one == 2 {
            doubling.push(code);
        }
    }
    println!("comp(c^L) = comp(c)^L for all 15 codes at L = 1..7, mismatches {bad}");
    println!(
        "comp(c) = 2 exactly for the codes {doubling:?} and 1 for the other {}",
        15 - doubling.len()
    );
    println!("so the constant-word rate is log 2 on the diagonal class and 0 elsewhere, and the frequency functional is Phi(f) = (f_6 + f_9) log 2");
    assert_eq!(bad, 0, "a constant word is a power");
    assert_eq!(doubling, vec![6, 9], "only the diagonal class doubles");
}

fn forms(rep: &Rep) {
    println!();
    println!("CLOSED FORMS");
    for family in families().iter() {
        let mut checked = 0usize;
        let mut bad = 0usize;
        let mut seen = 0usize;
        let mut wrong = 0usize;
        for (a, b) in family.pairs.iter() {
            let (one, two) = sweep(rep, family.form, *a, *b, REACH);
            checked += one;
            bad += two;
            let (three, four) = drawn(family.form, *a, *b, DRAWN);
            seen += three;
            wrong += four;
        }
        println!(
            "{}: {} pairs, comp = {}, {checked} words to L = {REACH}, mismatches {bad}, {seen} words drawn to L = {DRAWN}, mismatches {wrong}",
            family.name,
            family.pairs.len(),
            family.rule
        );
        assert_eq!(
            bad, 0,
            "the closed form is exact against the representation"
        );
        assert_eq!(wrong, 0, "the closed form is exact against the drawn word");
    }
    let named: usize = families().iter().map(|family| family.pairs.len()).sum();
    let all = CODES.len() * (CODES.len() - 1) / 2;
    println!("{named} of the {all} letter pairs carry a closed form, 9 of the 15 pairs of distinct classes and all 6 pairs inside one class; the {} pairs left open are the gasket class against the unit, domino and diagonal classes and the full tile against the unit and diagonal classes", all - named);
}

fn cut_law(rep: &Rep) {
    println!();
    println!("THE ZERO-CONTACT CUT");
    let mut across = [0u64; 16];
    let mut down = [0u64; 16];
    for &code in CODES.iter() {
        let (rows, columns) = render(&[code]).contacts();
        across[code as usize] = rows;
        down[code as usize] = columns;
    }
    let text = |table: &[u64; 16]| {
        CODES
            .iter()
            .map(|code| format!("{code}:{}", table[*code as usize]))
            .collect::<Vec<_>>()
            .join(" ")
    };
    println!("h {}", text(&across));
    println!("v {}", text(&down));
    let mut total = 0usize;
    let mut applies = 0usize;
    let mut bad = 0usize;
    for length in 1..5usize {
        for word in order::words(&CODES, length) {
            total += 1;
            let mut cut: Option<usize> = None;
            for start in 0..length {
                let rows: u64 = word[start..]
                    .iter()
                    .map(|code| across[*code as usize])
                    .product();
                let columns: u64 = word[start..]
                    .iter()
                    .map(|code| down[*code as usize])
                    .product();
                if rows == 0 && columns == 0 {
                    cut = Some(start);
                }
            }
            if let Some(start) = cut {
                applies += 1;
                let ahead: i128 = word[..start]
                    .iter()
                    .map(|code| code.count_ones() as i128)
                    .product();
                if value(rep, &word) != ahead * value(rep, &word[start..]) {
                    bad += 1;
                }
            }
        }
    }
    println!("comp = fill(prefix) * comp(suffix) at the last zero-contact suffix: {applies} of {total} words of length 1 to 4 admit the cut, mismatches {bad}");
    assert_eq!(bad, 0, "the cut law is exact where it applies");
}

pub fn morse(length: usize) -> Vec<u8> {
    (0..length).map(|n| (n.count_ones() & 1) as u8).collect()
}

fn doubling(length: usize) -> Vec<u8> {
    let mut word: Vec<u8> = vec![1];
    while word.len() < length {
        let mut next: Vec<u8> = Vec::with_capacity(2 * word.len());
        for letter in word.iter() {
            next.push(1);
            next.push(1 - letter);
        }
        word = next;
    }
    word.truncate(length);
    word
}

fn word_structure() {
    println!();
    println!("THE THUE-MORSE WORD");
    let bits = morse(HORIZON);
    let mut ones = 0usize;
    let mut uneven = 0usize;
    for (i, bit) in bits.iter().enumerate() {
        ones += *bit as usize;
        let length = i + 1;
        if length % 2 == 0 && ones * 2 != length {
            uneven += 1;
        }
    }
    println!("the prefix of every even length L <= {HORIZON} carries exactly L/2 of each letter, exceptions {uneven}");
    let mut longest = 0usize;
    let mut run = 0usize;
    let mut zeros = [0usize; 8];
    let mut zero_run = 0usize;
    for (i, bit) in bits.iter().enumerate() {
        if i > 0 && bits[i - 1] == *bit {
            run += 1;
        } else {
            run = 1;
        }
        zero_run = if *bit == 0 { zero_run + 1 } else { 0 };
        longest = longest.max(run);
        zeros[zero_run.min(7)] += 1;
    }
    println!("the longest terminal run over all L <= {HORIZON} is {longest}, either letter");
    println!(
        "the terminal run of the first letter takes the value 0 on {} prefixes, 1 on {}, 2 on {}, and never more",
        zeros[0], zeros[1], zeros[2]
    );
    let mut complete = [0usize; 8];
    let mut run = 1usize;
    for i in 1..bits.len() {
        if bits[i] == bits[i - 1] {
            run += 1;
        } else {
            complete[run.min(7)] += 1;
            run = 1;
        }
    }
    println!(
        "the first {HORIZON} letters hold {} complete runs of length 1 and {} of length 2, none longer, plus one unfinished run of length {run} at the cut",
        complete[1], complete[2]
    );
    let pd = doubling(bits.len());
    let mut bad = 0usize;
    for i in 0..bits.len() - 1 {
        if (bits[i] ^ bits[i + 1]) != pd[i] {
            bad += 1;
        }
    }
    println!("the run-boundary word t_n xor t_(n+1) is the period-doubling word on the first {} terms, mismatches {bad}", bits.len() - 1);
    assert_eq!(
        uneven, 0,
        "the letter counts are exactly equal at even length"
    );
    assert_eq!(longest, 2, "no three equal letters run together");
    assert_eq!(bad, 0, "the boundary word is the period-doubling word");
}

struct Reading {
    pair: (u8, u8),
    num: u64,
    den: u64,
    rate: &'static str,
    prediction: &'static str,
    fill: &'static str,
}

fn morse_rates(rep: &Rep) {
    println!();
    println!("THE THUE-MORSE EXPONENT");
    let bits = morse(HORIZON);
    let readings = [
        Reading {
            pair: (3, 1),
            num: 1,
            den: 2,
            rate: "(1/2) log 2",
            prediction: "0",
            fill: "(1/2) log 2",
        },
        Reading {
            pair: (1, 6),
            num: 1,
            den: 2,
            rate: "(1/2) log 2",
            prediction: "(1/2) log 2",
            fill: "(1/2) log 2",
        },
        Reading {
            pair: (3, 5),
            num: 1,
            den: 1,
            rate: "log 2",
            prediction: "0",
            fill: "log 2",
        },
        Reading {
            pair: (3, 6),
            num: 1,
            den: 1,
            rate: "log 2",
            prediction: "(1/2) log 2",
            fill: "log 2",
        },
        Reading {
            pair: (3, 15),
            num: 1,
            den: 2,
            rate: "(1/2) log 2",
            prediction: "0",
            fill: "(3/2) log 2",
        },
        Reading {
            pair: (7, 15),
            num: 0,
            den: 1,
            rate: "0",
            prediction: "0",
            fill: "(log 3 + log 4)/2",
        },
    ];
    for (family, reading) in families().iter().take(6).zip(readings.iter()) {
        let (a, b) = reading.pair;
        assert!(
            family.pairs.contains(&(a, b)),
            "the pick sits in its family"
        );
        let mut rates: Vec<String> = Vec::new();
        for swap in 0..2usize {
            let word: Vec<u8> = bits
                .iter()
                .map(|bit| if (*bit as usize) == swap { b } else { a })
                .collect();
            let mut short = 0usize;
            for length in 1..=REACH {
                let want = 1i128 << exponent(family.form, &word[..length], b);
                if value(rep, &word[..length]) != want {
                    short += 1;
                }
            }
            assert_eq!(short, 0, "the closed form reads the prefix");
            let power = exponent(family.form, &word, b) as u64;
            let gap = (power * reading.den) as i64 - (reading.num * HORIZON as u64) as i64;
            assert!(
                gap.unsigned_abs() <= 2 * reading.den,
                "the prefix rate sits within 2/L of the limit"
            );
            rates.push(ratio(power, HORIZON as u64));
        }
        println!(
            "{}: pair ({a},{b}), rate {}, prediction {}, fill rate {}, prefix rate at L = {HORIZON} is {} and {} under the two letter readings",
            family.name, reading.rate, reading.prediction, reading.fill, rates[0], rates[1]
        );
    }
    println!("the rate is the fill rate on every family above except domino and full, where it lies strictly between the per-letter value 0 and the fill rate, and gasket and full, where the word is connected");
}

fn boundary(rep: &Rep) {
    println!();
    println!("THE SIMPLEX BOUNDARY");
    let flat = 3u8;
    let mark = 6u8;
    let squares = |i: usize| {
        let root = (i as f64).sqrt() as usize;
        (root.saturating_sub(1)..root + 2).any(|r| r * r == i)
    };
    let mut bad = 0usize;
    for length in 1..=REACH {
        let powers: Vec<u8> = (1..=length)
            .map(|i| if i.is_power_of_two() { mark } else { flat })
            .collect();
        let boxes: Vec<u8> = (1..=length)
            .map(|i| if squares(i) { mark } else { flat })
            .collect();
        let flats = vec![flat; length];
        for word in [powers, boxes, flats] {
            let want = 1i128 << exponent(Form::DominoDiagonal, &word, mark);
            if value(rep, &word) != want {
                bad += 1;
            }
        }
    }
    println!("three words over the pair (3,6) whose marked letter has frequency 0, closed form checked to L = {REACH}, mismatches {bad}");
    let mut low: Vec<String> = Vec::new();
    for k in 2..14usize {
        let length = (1usize << (k + 1)) - 1;
        low.push(ratio(1u64 << k, length as u64));
    }
    println!("the marked letter at the powers of 2: the prefix rate is 1 at every L = 2^k and {} at L = 2^(k+1) - 1, so the rate has no limit, upper 1 and lower 1/2", low.join(", "));
    let mut squares_low: Vec<String> = Vec::new();
    for n in [2usize, 4, 8, 16, 32, 64] {
        squares_low.push(ratio((n * n) as u64, (n * n + 2 * n) as u64));
    }
    println!("the marked letter at the squares: the prefix rate is 1 at every L = n^2 and {} at L = (n+1)^2 - 1, so the rate is 1", squares_low.join(", "));
    println!("the constant word 3^L has rate 0, is periodic and hence minimal, and carries the same letter frequencies as both");
    println!("so at a frequency vector on the boundary of the simplex the exponent is 0, is 1, and fails to exist, and the closed forms give a rate only where both letters have positive frequency");
    assert_eq!(bad, 0, "the closed form covers the three boundary words");
}

pub fn study(rep: &Rep) {
    println!();
    println!("THE COMPONENT COCYCLE");
    constants(rep);
    forms(rep);
    cut_law(rep);
    word_structure();
    morse_rates(rep);
    boundary(rep);
}
