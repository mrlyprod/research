use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "paper-slice-recurrence-order";

fn main() -> Result<()> {
    let dims: Vec<usize> = (2..=14).collect();
    let free: Vec<usize> = dims.iter().map(|d| 2 * d + 1).collect();
    let proved: Vec<usize> = dims.iter().map(|d| d.div_ceil(2)).collect();
    assert_eq!(free.first(), Some(&5));
    assert_eq!(free.last(), Some(&29));
    assert_eq!(proved, vec![1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]);
    save(NAME, &json!({"free": free, "proved": proved}))
}
