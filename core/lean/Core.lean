module

@[expose] public section

def axisMaps (base : Nat) : List (List Nat) :=
  (List.range base).foldl (fun out b =>
    let rot := (List.range base).map fun r => (r + b) % base
    let ref := (List.range base).map fun r => (base + b - r) % base
    let out := if out.contains rot then out else out ++ [rot]
    if out.contains ref then out else out ++ [ref]) []

def residueCorners (dimension base : Nat) : List (List Nat) :=
  (List.range (base ^ dimension)).map fun i =>
    (List.range dimension).map fun j => i / base ^ (dimension - 1 - j) % base

def arrangements : Nat → List Nat → List (List Nat)
  | 0, _ => [[]]
  | k + 1, items => items.flatMap fun x => (arrangements k (items.erase x)).map (x :: ·)

def permutations (n : Nat) : List (List Nat) :=
  arrangements n (List.range n)

def choices (axis : List (List Nat)) (dimension : Nat) : List (List Nat) :=
  (List.range dimension).foldl (fun out _ =>
    out.flatMap fun pre => (List.range axis.length).map fun i => pre ++ [i]) [[]]

def group (base dimension : Nat) : List (List Nat) :=
  let cells := residueCorners dimension base
  let axis := axisMaps base
  (permutations dimension).flatMap fun perm =>
    (choices axis dimension).map fun choice =>
      cells.map fun cell =>
        (List.range dimension).foldl (fun acc i =>
          acc * base + ((axis.getD (choice.getD i 0) []).getD (cell.getD (perm.getD i 0) 0) 0)) 0

def carry (element : List Nat) (code : Nat) : Nat :=
  element.zipIdx.foldl (fun acc (image, index) =>
    if code.testBit index then acc + 2 ^ image else acc) 0

def insertSorted (x : Nat) : List Nat → List Nat
  | [] => [x]
  | y :: ys => if x < y then x :: y :: ys else if x = y then y :: ys else y :: insertSorted x ys

def orbit (group : List (List Nat)) (code : Nat) : List Nat :=
  group.foldl (fun set element => insertSorted (carry element code) set) []

def representatives (base dimension : Nat) : List (Nat × Nat) :=
  let group := group base dimension
  let walk := (List.range (2 ^ base ^ dimension)).foldl (fun (seen, out) bits =>
    if seen.testBit bits then (seen, out)
    else
      let members := orbit group bits
      (members.foldl (fun s m => s ||| 1 <<< m) seen, (bits, members.length) :: out)) (0, [])
  walk.2.reverse

set_option maxRecDepth 100000 in
theorem representatives_two_counts :
    (representatives 2 1).length = 3 ∧ (representatives 2 2).length = 6 ∧
      (representatives 2 3).length = 22 := by
  decide

theorem representatives_two_four_native : (representatives 2 4).length = 402 := by
  native_decide
