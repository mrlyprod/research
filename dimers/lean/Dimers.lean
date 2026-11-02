module

@[expose] public section

def digits (base code : Nat) : List (Nat × Nat) :=
  (List.range base).flatMap fun i =>
    (List.range base).filterMap fun j => if code.testBit (base * i + j) then some (i, j) else none

def level (base code : Nat) : Nat → List (Nat × Nat)
  | 0 => [(0, 0)]
  | n + 1 => (level base code n).flatMap fun (r, c) =>
      (digits base code).map fun (i, j) => (base * r + i, base * c + j)

def signed (base code : Nat) : Int :=
  ((digits base code).map fun (i, j) => (-1 : Int) ^ (i + j)).sum

def imbalanceLaw (base code n : Nat) : Int :=
  let s := signed base code
  let fill : Int := (digits base code).length
  if base % 2 = 1 then s ^ n else fill ^ (n - 1) * s

def direct (base code n : Nat) : Int :=
  ((level base code n).map fun (r, c) => if (r + c) % 2 = 0 then (1 : Int) else -1).sum

def codes (base : Nat) : List Nat :=
  List.range' 1 (2 ^ (base * base) - 1)

def imbalance (base top : Nat) : Nat :=
  (codes base).foldl (fun wrong code =>
    (List.range' 1 top).foldl (fun wrong n =>
      if direct base code n = imbalanceLaw base code n then wrong else wrong + 1) wrong) 0

def balanced (base : Nat) : Nat :=
  (codes base).countP fun code => signed base code = 0

def bruteCount (cells : List (Nat × Nat)) (side : Nat) : Nat :=
  let grid := cells.foldl (fun g (r, c) => g ||| 1 <<< (side * r + c)) 0
  let inCells := fun r c => r < side && c < side && grid.testBit (side * r + c)
  let empty : Array Nat := Array.replicate (2 ^ side) 0
  let states := (List.range (side * side)).foldl (fun states p =>
    let r := p / side
    let c := p % side
    (List.range states.size).foldl (fun nxt mask =>
      let ways := states[mask]!
      let covered := mask % 2 = 1
      let rest := mask / 2
      if ways = 0 then nxt
      else if !inCells r c || covered then
        if !inCells r c && covered then nxt else nxt.modify rest (· + ways)
      else
        let nxt := if c + 1 < side && inCells r (c + 1) && rest % 2 = 0
          then nxt.modify (rest ||| 1) (· + ways) else nxt
        if inCells (r + 1) c then nxt.modify (rest ||| 1 <<< (side - 1)) (· + ways) else nxt)
      empty) (empty.set! 0 1)
  states[0]!

def unbalancedTilings (base top : Nat) : Nat :=
  (codes base).foldl (fun total code =>
    if signed base code = 0 then total
    else (List.range' 1 top).foldl (fun total n =>
      total + bruteCount (level base code n) (base ^ n)) total) 0

theorem imbalance_two_small : imbalance 2 3 = 0 ∧ balanced 2 = 5 := by
  decide

theorem imbalance_native :
    imbalance 2 6 = 0 ∧ imbalance 3 4 = 0 ∧ balanced 2 = 5 ∧ balanced 3 = 125 ∧
      unbalancedTilings 2 4 = 0 ∧ unbalancedTilings 3 2 = 0 := by
  native_decide

theorem brute_count_native :
    bruteCount (level 3 495 1) 3 = 2 ∧ bruteCount (level 3 495 2) 9 = 6724 ∧
      bruteCount (level 3 255 1) 3 = 4 ∧ bruteCount (level 3 255 2) 9 = 1291616 ∧
      bruteCount (level 2 15 2) 4 = 36 := by
  native_decide
