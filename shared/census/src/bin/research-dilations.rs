use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "research-dilations";
const LEVEL: usize = 4;
const LAGS: usize = 31;
const PEAK: u64 = 5;

fn stern_field(level: usize) -> Vec<u64> {
    let mut r = vec![1u64];
    for j in 0..level {
        let s = 1usize << j;
        let mut out = vec![0u64; r.len() + 2 * s];
        for (n, &v) in r.iter().enumerate() {
            for k in 0..3 {
                out[n + k * s] += v;
            }
        }
        r = out;
    }
    r
}

fn main() -> Result<()> {
    let field = stern_field(LEVEL);
    assert_eq!(field.len(), LAGS);
    assert_eq!(field.iter().sum::<u64>(), 81);
    assert_eq!(field.iter().max(), Some(&PEAK));
    let peaks: Vec<usize> = (0..LAGS).filter(|&n| field[n] == PEAK).collect();
    assert_eq!(peaks, vec![10, 12, 18, 20]);
    assert_eq!(field.iter().filter(|&&r| r == 1).count(), 2 * LEVEL + 1);
    assert_eq!(peaks.len() * peaks.len(), 16);
    save(NAME, &json!({"field": field}))
}
