use crate::groups::{group, orbit, Elem};
use crate::rules::RULES;
use mrly_math::bang::universe::degree;
use mrly_math::bang::Code;
use mrly_num::boolean::walsh_spectrum;
use std::collections::BTreeSet;

pub fn levelset(code: usize) -> bool {
    let mut value = [-1i8; 4];
    for i in 0..8usize {
        let bit = ((code >> i) & 1) as i8;
        let w = i.count_ones() as usize;
        if value[w] != -1 && value[w] != bit {
            return false;
        }
        value[w] = bit;
    }
    true
}

pub fn pin(code: usize) -> bool {
    let cells: Vec<usize> = (0..8).filter(|i| (code >> i) & 1 == 1).collect();
    if cells.is_empty() {
        return false;
    }
    let fixed = (0..3)
        .filter(|axis| {
            cells
                .iter()
                .map(|c| (c >> axis) & 1)
                .collect::<BTreeSet<_>>()
                .len()
                == 1
        })
        .count();
    cells.len() == 1 << (3 - fixed)
}

pub fn genus(code: usize, b3: &[Elem]) -> &'static str {
    let cls = orbit(code, b3);
    if cls.iter().any(|&c| levelset(c)) {
        "iso"
    } else if cls.iter().any(|&c| pin(c)) {
        "axis"
    } else {
        "comp"
    }
}

pub fn levels(code: usize) -> [i64; 4] {
    let spectrum = walsh_spectrum(code as u128, 3);
    let mut out = [0i64; 4];
    for (s, w) in spectrum.iter().enumerate() {
        out[s.count_ones() as usize] += w;
    }
    out
}

pub struct Row {
    pub code: usize,
    pub rep_b3: usize,
    pub rep_h: usize,
    pub rep_both: usize,
    pub size_b3: usize,
    pub size_h: usize,
    pub size_both: usize,
    pub pop: u32,
    pub degree: i32,
    pub genus: &'static str,
    pub levels: [i64; 4],
    pub affine: bool,
}

pub fn rows() -> Vec<Row> {
    let b3 = group("B3");
    let h = group("H");
    let both = group("B3xZ2");
    (0..RULES)
        .map(|code| {
            let ob = orbit(code, &b3);
            let oh = orbit(code, &h);
            let oz = orbit(code, &both);
            let d = degree(Code::from(code as u128), 3);
            Row {
                code,
                rep_b3: *ob.iter().next().expect("nonempty"),
                rep_h: *oh.iter().next().expect("nonempty"),
                rep_both: *oz.iter().next().expect("nonempty"),
                size_b3: ob.len(),
                size_h: oh.len(),
                size_both: oz.len(),
                pop: (code as u32).count_ones(),
                degree: d,
                genus: genus(code, &b3),
                levels: levels(code),
                affine: d <= 1,
            }
        })
        .collect()
}

pub fn report() {
    println!("PER-RULE TABLE");
    println!("code repB3 repH repB3xZ2 sizeB3 sizeH sizeB3xZ2 pop deg genus S0 S1 S2 S3 affine");
    let rows = rows();
    for row in &rows {
        println!(
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            row.code,
            row.rep_b3,
            row.rep_h,
            row.rep_both,
            row.size_b3,
            row.size_h,
            row.size_both,
            row.pop,
            row.degree,
            row.genus,
            row.levels[0],
            row.levels[1],
            row.levels[2],
            row.levels[3],
            if row.affine { 1 } else { 0 }
        );
    }
    for row in &rows {
        assert_eq!(
            row.levels[0],
            8 - 2 * row.pop as i64,
            "the weight-zero Walsh sum of {} is not 8 - 2 pop",
            row.code
        );
    }
    let affine = rows.iter().filter(|r| r.affine).count();
    let iso = rows.iter().filter(|r| r.genus == "iso").count();
    let axis = rows.iter().filter(|r| r.genus == "axis").count();
    let comp = rows.iter().filter(|r| r.genus == "comp").count();
    println!("affine rules {affine}; genus counts iso {iso} axis {axis} comp {comp}");
    let mut classes: Vec<&Row> = rows.iter().filter(|r| r.code == r.rep_b3).collect();
    classes.sort_by_key(|r| r.code);
    let ci = classes.iter().filter(|r| r.genus == "iso").count();
    let ca = classes.iter().filter(|r| r.genus == "axis").count();
    let cc = classes.iter().filter(|r| r.genus == "comp").count();
    println!(
        "B3 classes {} with genus iso {ci} axis {ca} comp {cc}",
        classes.len()
    );
    let sizes: BTreeSet<usize> = classes.iter().map(|r| r.size_b3).collect();
    println!("B3 orbit sizes {sizes:?}");
}
