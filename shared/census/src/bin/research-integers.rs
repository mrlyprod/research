use census::{integers, save};
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "research-integers";
const CEILING: i128 = 100_000;
const WINDOW: usize = 10_000;

fn main() -> Result<()> {
    let counts = integers::census(CEILING, WINDOW, &[usize::MAX]).remove(0);
    let never = counts[1..].iter().filter(|&&c| c == 0).count();
    let once = counts[1..].iter().filter(|&&c| c == 1).count();
    let many = counts[1..].iter().filter(|&&c| c > 1).count();
    assert_eq!((never, once, many), (3589, 765, 5646));
    let peak = *counts.iter().max().unwrap();
    assert_eq!(peak, counts[16]);
    assert_eq!(counts.len(), WINDOW + 1);
    save(NAME, &json!(counts))
}
