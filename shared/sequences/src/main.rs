use mrly_core::error::Result;
use sequences::merge::{self, LONG};
use sequences::oeis::{self, Search};
use sequences::{io, page_rows, root, walk, write};
use std::time::Instant;

const BEAUTY: [(&str, &str); 8] = [
    (
        "ladder_base=3_dim=4",
        "the dimension 4 diagonal slice ladder of the base 3 parity solid, whose order 2 recurrence the cuts topic proves, sits in no record under the searched transforms.",
    ),
    (
        "ladder_base=5_dim=3",
        "the base 5 sibling of A299916: the same carry automaton one odd base up.",
    ),
    (
        "ladder_base=3_dim=5",
        "the dimension 5 ladder, order 3 by the proved order law ceil(dim/2), with no record.",
    ),
    (
        "sequence_dim=3_code=23_measure=holes_axis=side",
        "the sponge's hex slice holes are centered hexagonal at sides 4m+3 and zero at sides 4m+1, a parity switch no transform undoes.",
    ),
    (
        "sequence_dim=2_base=3_code=495_measure=surface_axis=side",
        "the base 3 carpet's perimeter at odd sides has second differences of period 3 in its 16 terms, a quasi-polynomial the dump lacks.",
    ),
    (
        "sequence_dim=3_code=126_measure=surface_axis=level",
        "code 126's exposure obeys a(level) = 23 a(level-1) - 90 a(level-2), roots 18 and 5, with no record.",
    ),
    (
        "press_dim=1_base=6_code=34",
        "the integers whose base 6 digits are all units mod 6, the digits 1 and 5.",
    ),
    (
        "sequence_dim=3_code=126_measure=holes_axis=side",
        "the holes of code 126's hex slice along the odd side are two quadratics interleaved by the parity of k in its 16 terms.",
    ),
];

fn main() -> Result<()> {
    let start = Instant::now();
    let mut rows = walk()?;
    let generated = rows.len();
    let walked = start.elapsed().as_secs_f64();
    rows.extend(page_rows()?);
    let classes = merge::classes(&rows);
    let relations = merge::relations(&rows, &classes);
    let search = match (std::env::var("OEIS_STRIPPED"), std::env::var("OEIS_NAMES")) {
        (Ok(stripped), Ok(names)) => oeis::search(
            &rows,
            &classes,
            &io(std::fs::read_to_string(stripped))?,
            &io(std::fs::read_to_string(names))?,
        )?,
        _ => {
            println!("OEIS_STRIPPED and OEIS_NAMES unset: no dump search");
            Search {
                hits: vec![Vec::new(); classes.len()],
                totals: vec![0; classes.len()],
                searched: vec![false; classes.len()],
                ..Search::default()
            }
        }
    };
    let out = root().join("data/research/shared/sequences");
    io(std::fs::create_dir_all(&out))?;
    io(std::fs::write(
        out.join("atlas.json"),
        write::json(&rows, &classes, &relations, &search),
    ))?;
    io(std::fs::write(
        out.join("atlas.md"),
        write::markdown(&rows, &classes, &search),
    ))?;
    let owner = write::class_of(&classes, rows.len());
    let collisions = classes.iter().filter(|c| c.members.len() > 1).count();
    let count = |kind: &str| classes.iter().filter(|c| c.kind == kind).count();
    let absent: Vec<usize> = (0..classes.len()).filter(|&c| search.absent(c)).collect();
    let long_absent: Vec<usize> = absent
        .iter()
        .copied()
        .filter(|&c| {
            let terms = &rows[classes[c].rep].terms;
            terms.len() >= LONG && !merge::simple(terms)
        })
        .collect();
    println!(
        "rows {} ({} generated, {} from the ledger page)",
        rows.len(),
        generated,
        rows.len() - generated
    );
    println!(
        "rows under {LONG} terms {} (the cell budget cuts the rendered level rows)",
        rows.iter().filter(|r| r.terms.len() < LONG).count()
    );
    println!("distinct sequences {}", classes.len());
    println!(
        "collisions {collisions}: explained {}, trivial {}, one design two measures {}, across topics {}, across codes and measures {}",
        count("explained"),
        count("trivial"),
        count("design"),
        count("topic"),
        count("codes")
    );
    println!(
        "searched {}, absent {} ({} with {LONG}+ terms, not constant, linear or geometric)",
        search.searched.iter().filter(|&&s| s).count(),
        absent.len(),
        long_absent.len()
    );
    println!(
        "stated ids failing their value check {}",
        search.failed.len()
    );
    for (i, id) in &search.failed {
        println!("  failed {} {}", rows[*i].name, id);
    }
    let names = |c: usize, n: usize| -> String {
        let names: Vec<&str> = classes[c]
            .members
            .iter()
            .filter(|m| rows[m.0].measure != "ledger")
            .take(n)
            .map(|m| rows[m.0].name.as_str())
            .collect();
        names.join(" = ")
    };
    for kind in ["design", "topic"] {
        for c in (0..classes.len())
            .filter(|&c| classes[c].kind == kind)
            .take(25)
        {
            let terms: Vec<String> = rows[classes[c].rep]
                .terms
                .iter()
                .take(6)
                .map(i128::to_string)
                .collect();
            println!("  {kind} [{}]: {}", terms.join(","), names(c, 3));
        }
    }
    let mut disagree = 0;
    for (i, row) in rows.iter().enumerate().skip(generated) {
        if row.generator != "the ledger page" {
            continue;
        }
        if let Some(j) = rows[..generated].iter().position(|r| r.name == row.name) {
            if owner[i].0 != owner[j].0 {
                disagree += 1;
                println!("  ledger disagrees {}", row.name);
            }
        }
    }
    println!("ledger rows disagreeing with the walk {disagree}");
    let mut laws = 0;
    for (i, row) in rows.iter().enumerate() {
        if let Some(rest) = row.name.strip_prefix("arcs_law_") {
            let drawn = format!("arcs_{rest}");
            if let Some(j) = rows.iter().position(|r| r.name == drawn) {
                if owner[i].0 != owner[j].0 {
                    laws += 1;
                    println!("  law disagrees {}", row.name);
                }
            }
        }
    }
    println!("loop laws disagreeing with the drawn counts {laws}");
    println!(
        "transform relations across measures {} ({} within one design or across topics)",
        relations.len(),
        relations.iter().filter(|r| r.near).count()
    );
    for r in relations.iter().filter(|r| r.near).take(30) {
        println!(
            "  {}({}) at offset {} = {}",
            r.transform,
            names(r.from, 1),
            r.offset,
            names(r.to, 1)
        );
    }
    let mut candidates = String::new();
    for &c in &long_absent {
        let rep = &rows[classes[c].rep];
        let terms: Vec<String> = rep.terms.iter().take(10).map(i128::to_string).collect();
        candidates += &format!(
            "{} | {} | {} | {} | {}\n",
            c,
            rep.name,
            rep.closed.as_deref().unwrap_or("none"),
            terms.join(","),
            classes[c].members.len()
        );
    }
    io(std::fs::write(out.join("absences.txt"), candidates))?;
    println!("the eight most beautiful absences");
    for (name, sentence) in BEAUTY {
        let state = rows
            .iter()
            .position(|r| r.name == name)
            .map(|i| owner[i].0)
            .map(|c| if search.absent(c) { "absent" } else { "found" })
            .unwrap_or("missing");
        println!("  {name} ({state}): {sentence}");
    }
    println!(
        "seconds {:.1} ({walked:.1} for the walk)",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
