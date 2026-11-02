---
title: Beneath a design
lead: The object beneath a design as three slots, accept, place and glue: the memory dial priced by a transfer matrix with `kappa` zero on exactly the non-empty product rules over `19563` live classes, the Collatz carry landed on that dial as a four-state transducer of unbounded radius whose zero-carry sets are the golden and supergolden rules at `3n + 1` and `7n + 1`, and the radix dial that keeps every count, carries the gasket, the dragons and the Koch curve as codes with a twist vector, and loses the fill law the moment a digit turns; the Mobius meter of every width-three rule to `2^30` against its own mass, the memory zeta as a matrix ladder with one pole comb per eigenvalue, and the design census up to similarity and up to affine conjugacy.
figure: research-beneath
slug: beneath
---

A design as [core](../../core/notes/core.md) defines it freezes three independent choices at once, and only one of them is visible in the definition. This page names all three and turns the two that carry a parameter.

A design is three slots.

- **Accept.** Which digit words survive. Today: the full shift on the digit set `F`, every word over `F` allowed and no word tested against its neighbours.
- **Place.** Where each word lands. Today: place value at base `base` on `Z^dim`, the word `d_1 d_2 d_3 ...` landing on the point `sum_j d_j base^(-j)`.
- **Glue.** Which words name one point. Never chosen: two words name one point exactly when they land on the same point of the ambient.

Accept and place do not interfere, and the split is exact: the number of accepted words of length `level` is a function of the accept slot alone, and the point set drawn at level `level` is the image of those words under the place slot. (**Proved**, from the definitions: neither the set of accepted words nor its cardinality mentions the place map, and the drawn set is by definition the image.) So the counting side of this tree is carried by accept, with place entering a count only through the contraction ratio and through which words collide; the drawing side is carried entirely by place.

One sentence explains why the counts close. Today's accept slot tests each digit alone, so the accepted words of length `level` are the words over `F` and there are `card F^level` of them, a product of `level` equal factors. (**Proved**, from the definition of the full shift.) That product is the property the two dials below are measured against. The memory dial breaks it in a controlled way and prices the break with a number; the radix dial leaves it untouched and moves only the picture.

Throughout, a digit vector `d` in `{0,1}^dim` is read as the corner integer `c(d) = sum_i d_i 2^i`, axis `i` on bit `i`, so at `dim 2` bit `0` is the column `x` and bit `1` the row `y`. That is the crate's reading too: `mrly_math::bang::universe::corners(dim)` emits the corner vector row first and `corner_index` folds it most significant first, so a crate design's corner integer at `dim 2` is `c = x + 2y`. (**Verified**, `lab/py/memory-census`, pinned cell for cell by the memory demo's host test on codes `11` and `13`, which the axis swap exchanges and which are drawn differently.)

## Chapters

- [The memory dial](memory-dial.md): The accept slot dialled by a width-`k` window rule: its transfer matrix, the coupling `kappa` that is zero on exactly the non-empty product rules, the census of live classes, the famous constants one notch in, and the Collatz carry landed on the dial as a four-state transducer.
- [The memory meter and zeta](memory-meter.md): The Mobius meter of every width-three rule to `2^30` read against the accepted set's own mass, and the memory zeta, the Dirichlet series of a rule's accepted integers, as a matrix ladder with one pole comb per eigenvalue.
- [The question mark](question-mark.md): The Stern-Brocot tree as a base-2 design with a rank-2 carry, and Minkowski's `?` as the swap of the place slot that turns restricted partial quotients into run-length rules of the memory dial, with their two dimensions.
- [Rules as continued fraction sets](graph-directed.md): The dimension of the accepted set of any rule with a cycle, and every rule read as a graph-directed continued fraction set whose dimension is the zero of its pressure.
- [The radix dial](radix-dial.md): The place slot dialled by moving addresses to the Gaussian or Eisenstein integers: the fill law and where it stops, the twist law, the gasket, dragons and Koch curve as codes, the digit codes a base carries up to its unit group, the glue slot, and where every number is generated.
- [The chaos game](chaos-game.md): The game read as a shift register on digits: exact forgetting in `k` steps, de Bruijn painting in `card F^k + k - 1` draws, the random painting time as an alternating sum over restricted games, restricted games as the memory dial with the census of the square, and the families of the tree as games.
