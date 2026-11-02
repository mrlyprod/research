use std::fs;
use std::path::Path;

pub fn write_csv(path: &Path, header: &[String], rows: &[Vec<String>]) {
    let mut out = header.join(",");
    out.push('\n');
    for row in rows {
        let fields: Vec<String> = row.iter().map(|field| quote(field)).collect();
        out.push_str(&fields.join(","));
        out.push('\n');
    }
    fs::write(path, out).expect("the table is writable");
}

fn quote(field: &str) -> String {
    if field.contains([',', '"', '\n']) {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}
