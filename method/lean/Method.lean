module

public section

def solidCells (mask s : Nat) : List (Nat × Nat × Nat) :=
  (List.range s).flatMap fun x => (List.range s).flatMap fun y =>
    ((List.range s).filter fun z =>
      (mask >>> ((x &&& 1) ||| ((y &&& 1) <<< 1) ||| ((z &&& 1) <<< 2))) &&& 1 = 1).map
      fun z => (x, y, z)

def fills (mask s : Nat) : Nat := (solidCells mask s).length

theorem menger_odd_side_fill :
    ∀ k ∈ List.range' 1 6, fills 23 (2 * k - 1) = 4 * k ^ 3 - 3 * k ^ 2 := by
  decide +kernel

end
