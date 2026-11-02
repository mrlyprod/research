# Sequences

- The sequence atlas: one bin walks every integer-sequence generator in the tree and writes the census the sequence work reads from.
- Generators: the fill law, the exposure recurrence, the renders and the hex slice over every design of the nine spaces on both axes, orbit mates included where their sequences differ; the carry ladder, its fill and sign; the solid slice closed forms; the diagonal profile; the design counts, bracelets and fill classes; the press; the arc loops and their laws; the 23 life sources; the rows of `../sequences.md` as rendered.
- Budget: 16 terms, 600000 cells a rendered term, 2 seconds a sequence.
- Dedup: equal terms under a shift of up to 3 register once with every name; a row under 8 terms joins a class only when one class contains it.
- Kinds of collision: explained (one measure, a source read twice, a law and its drawing, plane surface and perimeter), trivial (constant, linear, geometric), one design two measures, across topics, across codes and measures.
- Relations: one class equal to `-a`, the first differences, the partial sums, `a-1`, `a+1`, `a/2` or a bisection of another.
- OEIS: four 5-term windows at offsets 0 to 3 for `a`, `-a`, the differences, the partial sums, `a-1`, `a+1`, `a/2` and `a/gcd`; a hit agrees with the record on every term both carry; an absence is Conjecture.
- Run: `OEIS_STRIPPED=<stripped dump> OEIS_NAMES=<names file> cargo run --release --manifest-path research/shared/sequences/Cargo.toml` from the repo root.
- Output: `data/research/shared/sequences/atlas.json` (rows, classes, relations, failed ids), `atlas.md` (one table sorted by generator), `absences.txt`.
- Test: `cargo test --manifest-path research/shared/sequences/Cargo.toml` pins the walk at 99152 rows and one orbit-mate dedup case.
