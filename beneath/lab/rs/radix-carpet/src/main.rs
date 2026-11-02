use mrly_math::radix::tile;

const BOX_CODE: u128 = 0b1_1110_1111;

fn main() {
    let carpet = tile(3, BOX_CODE).unwrap();
    let residues = carpet.base().residues().unwrap();
    let class_code = carpet.code().unwrap();
    let canonical = carpet.canonical().unwrap();

    println!("CARPET  base 3 on Z[i], the box digits minus the centre");
    println!("  box code {BOX_CODE}, bit r m + c: the cell at row r and column c");
    println!(
        "  size {}  fill(2) {}  distinct(2) {}  dimension {:.6}",
        carpet.size(),
        carpet.fill(2),
        carpet.distinct(2),
        carpet.dimension()
    );
    println!("  canonical digits {canonical}");
    println!("  class code {class_code}, bit i: the class of canonical residue i");
    println!("  canonical residues {residues:?}");
    println!("  codes differ {}", class_code != BOX_CODE);

    assert_eq!(carpet.size(), 8);
    assert_eq!(carpet.fill(2), 64);
    assert_eq!(carpet.distinct(2), 64);
    assert_eq!(format!("{:.6}", carpet.dimension()), "1.892789");
    assert!(!canonical);
    assert_eq!(class_code, 479);
    assert_ne!(class_code, BOX_CODE);
}
