# Shared

The machinery every topic leans on: the sequence ledger, the sequence atlas, the claims checker, the figure census and the Lean package.

- [ledger/](ledger/): walks `research/*/notes` and `research/*/claims` and writes `sequences.md`.
- [claims/](claims/): re-runs claim witnesses.
- [sequences/](sequences/): the sequence atlas, one walk of every integer-sequence generator, de-duplicated, related and searched against the OEIS.
- [notes/atlas.md](notes/atlas.md): the atlas chapter, its census, the laws behind its collisions and the ranked absences.
- [claims/atlas.md](claims/atlas.md): the atlas claims.
- [figures/research-sequences-atlas.ts](figures/research-sequences-atlas.ts): the atlas figure, every distinct sequence one cell, ordered by growth.
- [census/](census/): one bin a figure, writing the JSON beside its figure.
- [figures/research-sequences.ts](figures/research-sequences.ts): the sequences figure, six codes by four odd sides, each checked against its fill count.
- [refs.md](refs.md): the OEIS rows the ledger cites.
- [sequences.md](sequences.md): the generated OEIS ledger every topic links.
- [lean/](lean/): the one Lake package with Mathlib; each topic's proofs live in its own `lean/` folder.
- [constants.md](constants.md): the constants ledger, each row the digits its generator prints and the OEIS entry if any.
