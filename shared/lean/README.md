# Lean

One Lake package for the whole research tree: it requires Mathlib, and every topic's proofs live in that topic's own `lean/` folder.

- `lakefile.toml`: the package `mrly`, Mathlib pinned at tag `v4.35.0`, one `lean_lib` per topic with `srcDir = "../../<topic>/lean"`.
- `lean-toolchain`: `leanprover/lean4:v4.35.0`, the toolchain the pinned Mathlib tag uses.
- Dependencies and build outputs go to `data/shared/lean/` at the repo root, never here; `.lake/` holds only the lakefile cache and is ignored.

## BUILD

- First time: `lake -d research/shared/lean update`, which clones the dependencies and pulls the Mathlib oleans from its cache; then delete `research/shared/lean/.lake/packages`, a stray copy the cache unpacks at the default path.
- Every time, from the repo root: `lake -d research/shared/lean build`.
- One topic: `lake -d research/shared/lean build <Topic>`.

## A TOPIC ADDS A FILE

- Create `research/<topic>/lean/<Topic>.lean`, starting with `module`; import only what the proof needs (core Lean often suffices, else `public import Mathlib`).
- Add to `lakefile.toml` a `[[lean_lib]]` with `name = "<Topic>"` and `srcDir = "../../<topic>/lean"`, and the name to `defaultTargets`.
- A second file in the same topic is a module `<Topic>.<Name>` at `research/<topic>/lean/<Topic>/<Name>.lean`, and `<Topic>.lean` must import it, or `lake build` never compiles it.
- `core/lean/Core.lean` is the template: the generator's functions under their own names in an `@[expose] public section`, no import when core Lean suffices, `decide` on the small cases, `native_decide` in a `_native` theorem past the kernel's reach.

## RULES

- A Lean file proves one finite headline, and its proof is a computation: `decide` preferred; `native_decide` allowed and marked in the theorem name.
- Never analytic number theory, never an infinite statement dressed as finite, never a `sorry`-free trick (an axiom, an `unsafe` definition, a `native_decide` on a hand-fed instance).
- A built file earns its claims row the tag Proved with witness `../lean/<File>.lean`.
- The statement names the generator it mirrors, so the row's number and the theorem's number are read off one object.
