module

@[expose] public section

def VOID : Nat := 0

def FILL : Nat := 1

def GRID : Nat := 2

abbrev Point := Nat × Nat

abbrev Edge := Point × Point

def ones (number : Nat) : Nat × (Nat → Nat → Nat → Nat) :=
  (number, fun _ _ _ => 1)

def cut (cell : Nat × (Nat → Nat → Nat → Nat)) : List (List Nat) :=
  let scale := 4
  let grid := fun x y z => cell.2 (x / scale) (y / scale) (z / scale)
  let size := cell.1 * scale
  let k := 3 * (size - 1) / 2
  let rows := ((List.range size).filter (· % 2 = 0)).filterMap fun z =>
    let target := k - z
    let minX := target - (size - 1)
    let maxX := min (size - 1) target
    if minX > maxX then none
    else some ((List.range' minX (maxX + 1 - minX)).map fun x => grid x (target - x) z)
  let width := rows.foldl (fun m (row : List Nat) => max m row.length) 0
  rows.map fun (row : List Nat) =>
    let offset := (width - row.length) / 2
    (List.range width).map fun c =>
      if offset ≤ c ∧ c < offset + row.length then
        if row.getD (c - offset) 0 = 1 then FILL else VOID
      else GRID

def north (x y : Nat) : List Point :=
  [(x, 2 * y + 2), (x + 1, 2 * y), (x + 2, 2 * y + 2)]

def south (x y : Nat) : List Point :=
  [(x, 2 * y), (x + 1, 2 * y + 2), (x + 2, 2 * y)]

def corners (x y start : Nat) : List Point :=
  if (x + y + start) % 2 = 0 then north x y else south x y

def sorted (a b : Point) : Edge :=
  if a.1 < b.1 ∨ (a.1 = b.1 ∧ a.2 ≤ b.2) then (a, b) else (b, a)

def edgesOf (c : List Point) : List Edge :=
  [sorted (c.getD 0 (0, 0)) (c.getD 1 (0, 0)), sorted (c.getD 1 (0, 0)) (c.getD 2 (0, 0)),
    sorted (c.getD 0 (0, 0)) (c.getD 2 (0, 0))]

def bump (e : Edge) : List (Edge × Nat) → List (Edge × Nat)
  | [] => [(e, 1)]
  | (f, n) :: rest => if e = f then (f, n + 1) :: rest else (f, n) :: bump e rest

structure Mesh where
  vertices : List Point
  edges : List (Edge × Nat)

structure Census where
  triangles : Nat
  vertices : Nat
  edges : Nat
  boundaryEdges : Nat
  interiorEdges : Nat
  euler : Int

def Mesh.add (mesh : Mesh) (x y start : Nat) : Mesh :=
  let c := corners x y start
  { vertices := c.foldl (fun vs p => if vs.contains p then vs else p :: vs) mesh.vertices,
    edges := (edgesOf c).foldl (fun es e => bump e es) mesh.edges }

def Mesh.tally (mesh : Mesh) (triangles : Nat) : Census :=
  let edges := mesh.edges.length
  let boundary := (mesh.edges.filter fun (_, n) => n = 1).length
  { triangles, vertices := mesh.vertices.length, edges, boundaryEdges := boundary,
    interiorEdges := edges - boundary,
    euler := (mesh.vertices.length : Int) - edges + triangles }

def census (types : List (List Nat)) (start : Nat) : Census :=
  let cells := types.zipIdx.flatMap fun (row, y) =>
    (row.zipIdx.filter fun (v, _) => v ≠ GRID).map fun (_, x) => (x, y)
  (cells.foldl (fun mesh (x, y) => mesh.add x y start) (⟨[], []⟩ : Mesh)).tally cells.length

def solidSlice (number : Nat) : Census :=
  census (cut (ones number)) 0

def readings (rec : Census) : List Nat :=
  [rec.triangles, rec.boundaryEdges, rec.edges, rec.interiorEdges, rec.vertices]

def closedForms (k : Nat) : List Nat :=
  let number := 2 * k - 1
  [6 * number * number, 6 * number, 6 * number * (3 * k - 1), 6 * number * (3 * k - 2),
    3 * number ^ 2 + 3 * number + 1]

theorem solid_slice_census :
    ∀ k ∈ List.range' 1 3,
      readings (solidSlice (2 * k - 1)) = closedForms k ∧ (solidSlice (2 * k - 1)).euler = 1 := by
  decide +kernel

theorem solid_slice_census_native :
    ∀ k ∈ List.range' 1 10,
      readings (solidSlice (2 * k - 1)) = closedForms k ∧ (solidSlice (2 * k - 1)).euler = 1 := by
  native_decide

end
