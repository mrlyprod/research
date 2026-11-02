use crate::rules::{output, RULES};
use mrly_math::bang::code_to_corners;
use mrly_math::bang::Code;

pub fn report() {
    println!("IDENTITY");
    let mut checked = 0usize;
    for rule in 0..RULES {
        let corners =
            code_to_corners(Code::from(rule as u128), 3, 2).expect("the code fits three axes");
        let mut from_design = [false; 8];
        for corner in &corners {
            from_design[4 * corner[0] as usize + 2 * corner[1] as usize + corner[2] as usize] =
                true;
        }
        for l in 0..2u8 {
            for c in 0..2u8 {
                for r in 0..2u8 {
                    let index = 4 * l as usize + 2 * c as usize + r as usize;
                    assert_eq!(
                        from_design[index],
                        output(rule, l, c, r) == 1,
                        "rule {rule} disagrees with bang dim 3, code {rule} at ({l},{c},{r})"
                    );
                    checked += 1;
                }
            }
        }
    }
    println!("rule N and bang dim 3, code N agree on all {checked} rule cells over {RULES} codes");
    println!("the dictionary is (x0,x1,x2) = (l,c,r) and corner index i = 4 x0 + 2 x1 + x2");
}
