//! The census: every number writer of the figures, one bin each, its json beside it.
#![deny(missing_docs)]

/// The ledger sweep the two integers writers share.
pub mod integers;

use mrly_core::error::{value_error, Result};
use mrly_core::Json;
use mrly_core::Error;
use std::path::PathBuf;

fn path(name: &str) -> Result<PathBuf> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let found: Vec<PathBuf> = std::fs::read_dir(&root)
        .map_err(|e| Error::Value(format!("cannot read {root:?}: {e}")))?
        .flatten()
        .map(|entry| entry.path().join("figures").join(format!("{name}.ts")))
        .filter(|figure| figure.is_file())
        .collect();
    match found.as_slice() {
        [figure] => Ok(figure.with_extension("json")),
        _ => value_error(format!(
            "census {name} needs exactly one figures/{name}.ts, found {}.",
            found.len()
        )),
    }
}

fn has_null(value: &Json) -> bool {
    match value {
        Json::Null => true,
        Json::Array(items) => items.iter().any(has_null),
        Json::Object(map) => map.values().any(has_null),
        _ => false,
    }
}

/// Saves the value as `<name>.json` beside the one `research/*/figures/<name>.ts`, whatever the cwd: compact, keys in written order, floats in their shortest round-trip form, one trailing newline.
///
/// # Errors
///
/// Errs on a name outside `a-z`, `0-9` and `-`, on a null anywhere in the value (a NaN or an infinity became one), when no single figure of that name exists, or when the file cannot be written.
pub fn save(name: &str, value: &Json) -> Result<()> {
    if name.is_empty()
        || !name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
    {
        return value_error(format!("census name {name:?} is not a-z, 0-9 and -."));
    }
    if has_null(value) {
        return value_error(format!("census {name} holds a null, a NaN or an infinity."));
    }
    let file = path(name)?;
    let text = format!("{value}\n");
    std::fs::write(&file, &text)
        .map_err(|e| Error::Value(format!("cannot write {file:?}: {e}")))?;
    println!("census {name} {} bytes", text.len());
    Ok(())
}
