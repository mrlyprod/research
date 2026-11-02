# carpet-stack-address

- Brightness of the odd parity-carpet stack at a rational point, in closed form, with no loop over layers anywhere.
- What layer `n` does at `x = (a_1/q, a_2/q)` depends only on `r = n mod 2q`, so `B_N(x) = ceil(N/2) - Sum_{r in S(x), r <= N} (floor((N-r)/2q) + 1)` with `S(x)` the odd residues whose two coordinates both land dark.
- Per-point cost is one pass over the `2q` residues and never depends on `N`; the form is set against literal stacking on 48 points at `q = 9973`, and against a second route that counts the bright residues directly instead of subtracting the dark ones.
- Farey line-stack at `N = 55`: the node value `floor(N/b)` against literal membership in each of the 55 scales, hashed to the same digest.
- The `512 x 512` render at `N = 55` by three routes sharing no inner loop: column bit-masks, layer accumulation, per-pixel triple loop.

## RUN

- `uv run python research/farey/lab/py/carpet-stack-address/stack.py`
- From the repo root. One core, about a second.
- Domain is the full source domain: 48 points at `N = 55`, `5555` and `10^18`, the sweep `N in {1, 2, 55, 5555, 19945, 19946, 19947, 40001}`, the whole `512 x 512` render at `N = 55`.
- Nothing is written to disk; the render is hashed in memory.

## WITNESSES

- `farey.md:214` - the closed form itself, the only expression the point evaluator computes.
- `farey.md:220-221` - 940 Farey nodes at `N = 55`, brightness sum `1540 = N(N+1)/2` landed by count, both routes on digest `f42186ea5b3670ac677a5a4d42dc1acaed16cf34d10884de000277165c810834`.
- `farey.md:221-222` - the `512 x 512` render sha256-identical by three routes, `6cddd85b6298cffdeb374df76365be324ab96f37c212110f683e07c7d4e60239`.
- `farey.md:222` - 48 of 48 points equal, closed form against literal stacking, at `N = 55` and `N = 5555`.
- `farey.md:223` - 384 comparisons, 0 mismatches, over the 8 sweep values, both sides of the `2q = 19946` period boundary.
- `farey.md:224-225` - `N = 10^18`, a stack of `5*10^17` layers, both routes agreeing on all 48 values in well under a tenth of a second; first four `375564022861726680`, `374862127744911263`, `374862127744911259`, `375012533841371703`.
- the claims line - the same numbers on the ledger line for this section.
