# Tree Walks

- The crate's own render and graph of `bang dim 2, base 3, code 79` (levels 1 to 6) and `bang dim 2, code 7` (levels 1 to 9): nodes, edges, and FNV-1a digests of the sorted filled cell list and of the sorted edge list over flat cell indices.
- It exists to pin `../../py/tree-walks/trees.py`, which rebuilds both designs from the parity rule with no shared code: the two print identical lines, so the Python graphs are the graphs the walk-dimension census measured.

## RUN

- `cargo run --release --manifest-path research/walks/lab/rs/tree-walks/Cargo.toml`
- About 2 s warm; prints only.

## WITNESSES

- trees.md, section "Two designs are trees": every line, `edges = nodes - 1` at all 15 levels.
