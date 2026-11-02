use crate::checks;

#[derive(Clone, Copy)]
pub enum Cost {
    Cheap,
    Dear,
}

#[derive(Clone, Copy)]
pub enum Verdict {
    Checked,
    Failed,
    Unchecked,
    Skipped,
    Orphan,
}

pub type Check = fn() -> Result<(), String>;

pub type Entry = (&'static str, Cost, Check);

pub fn entries() -> &'static [Entry] {
    checks::table()
}
