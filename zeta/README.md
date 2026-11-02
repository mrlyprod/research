---
title: The zeta function of a digit design
lead: The design's own zeta function: an explicit zero-free half plane and the census it closes, a residue comb whose teeth are Rouche-certified one by one, a second family that is no critical line and obeys no counting, symmetry, contraction or gain law, the ordinate shadow as a constant-free Newton step from each zeta zero, the multiplicativity wall with the three products that stand where the Euler product does not, and the identity zeta_F N_F = 1 whose Mertens function runs the wrong way.
figure: research-zeta
slug: zeta
date: 2026-09-07
---

# The zeta function of a digit design

The design's own zeta function: an explicit zero-free half plane and the census it closes, a residue comb whose teeth are Rouche-certified one by one, a second family that is no critical line and obeys no counting, symmetry, contraction or gain law, the ordinate shadow as a constant-free Newton step from each zeta zero, the multiplicativity wall with the three products that stand where the Euler product does not, and the identity zeta_F N_F = 1 whose Mertens function runs the wrong way.

- [notes/zeta.md](notes/zeta.md): The zeta function of a digit design
- [notes/design-zeta-zeros.md](notes/design-zeta-zeros.md): The zeros of the design zeta
- [notes/wall-and-identity.md](notes/wall-and-identity.md): The wall, the products and the identity
- [notes/meter-and-closure.md](notes/meter-and-closure.md): The Mobius meter and the closure problem
- [claims/digit-designs-and-the-euler-product.md](claims/digit-designs-and-the-euler-product.md): Digit designs and the Euler product
- [claims/zeros-of-the-design-zeta.md](claims/zeros-of-the-design-zeta.md): The zeros of the design zeta
- [papers/design-dirichlet-inverse.md](papers/design-dirichlet-inverse.md): The Dirichlet Inverse of a Digit Design and the Transport of Zeros
- [refs.md](refs.md): the references this topic cites
- [lab/py/design-zeta/](lab/py/design-zeta/): The Dirichlet series of a digit design, `zeta_F(s) = sum_(n in S_F) n^(-s)` over the base-`base` integers whose digits all lie in a digit set `F`, continued to the whole plane and read for its zeros.
- [lab/py/mrly-euler/](lab/py/mrly-euler/): The multiplicative structure of a digit design `S_F` over the digit set `F`
- [lab/py/nyman-beurling/](lab/py/nyman-beurling/): The Nyman-Beurling distance restricted to a digit design
- [lab/py/transport-census/](lab/py/transport-census/): The transport census
- [lab/py/zeta-family/](lab/py/zeta-family/): The SECOND FAMILY of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)`
- [lab/py/zeta-locus/](lab/py/zeta-locus/): The locus of the zeros of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)`, with `F` the digits of the code and `S_F` the integers they write
- [lab/py/zeta-shadow/](lab/py/zeta-shadow/): The MECHANISM of the ordinate shadow
- [figures/paper-design-dirichlet-inverse.ts](figures/paper-design-dirichlet-inverse.ts): paper-design-dirichlet-inverse
- [figures/research-zeta.ts](figures/research-zeta.ts): research-zeta
- [figures/research-zeta.json](figures/research-zeta.json): research-zeta census

## Lost

- A functional equation for `zeta_F`, or Weil's theorem carried to it: the Hurwitz reflection moves the kernel and never relates `zeta_F(s)` to `zeta_F(1 - s)`, and the evaluation `P -> P(base)` is no bridge (THE WALL AND THE PRODUCTS; [claims/digit-designs-and-the-euler-product.md](claims/digit-designs-and-the-euler-product.md)).
- An Euler product on a proper design, or the conjecture on the design's own Mobius `nu_F`: the indicator of `S_F` is multiplicative only at the full digit set, and a boxed zero right of the abscissa, resolved not certified, refutes every square-root bound for `nu_F`, so the conjecture lives on `mu` restricted to `S_F` ([claims/digit-designs-and-the-euler-product.md](claims/digit-designs-and-the-euler-product.md)).
- Gluing `zeta_F M_F = 1 + D_F`: the abscissa of `D_F` is exactly `alpha`, so the glue gains nothing (THE IDENTITY).
- The per-denominator split of the position pairing: its `l^1` mass sits at the top level, the large sieve does not rescue it on the major arcs, and the exact second moment prices the whole grid at `(alpha + 1)/2` (THE WALL AND THE PRODUCTS).
- Laws for the second family of zeros: no curve, no symmetry axis, no counting law in `alpha` or the fill, and the rightmost zero's gain is no function of `alpha` and `fill/base` (THE ZEROS; the Refuted rows in [claims/zeros-of-the-design-zeta.md](claims/zeros-of-the-design-zeta.md)).
- The design's pole lattice and the quadratic `L`-zeros as frequencies of the Mobius meter: the meter oscillates at the zeta ordinates and carries neither (THE METER).
- Nyman-Beurling on a proper design: the closure reaches `chi` only when the set holds every squarefree integer, so a design's distance plateaus (at base 3 `{0,1}` its limit lies in `[0.101926, 0.104579]`) and the criterion is a floor, not a route (THE CLOSURE PROBLEM).
