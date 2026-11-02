import { type Color, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";

export const units = { num };

const P = 7;
const LEVEL = 6;
const SPAN = Math.floor((2 * P - 1) / 3) + 1;
const SIDE = SPAN ** (LEVEL / 2);
const LOW = Math.floor((P - 1) / 3);
const NAME = "paper-digit-automaton-primes";

type Cell = { x: number; y: number; kind: number };

let memo: Cell[] | undefined;

function board() {
  if (memo) return memo;
  const carry = num.kummer.carry(P, 3, 1);
  const { moves, closes } = carry;
  const fill = (P + 1) / 2;
  for (let s = 0; s < 2; s++) {
    const kept = moves[s].filter((to) => to !== undefined).length;
    if (kept !== fill) throw new Error(`${NAME}: state ${s} keeps ${kept} digits, want ${fill}`);
    for (let d = SPAN; d < P; d++)
      if (moves[s][d] !== undefined) throw new Error(`${NAME}: digit ${d} is allowed in state ${s}, the crop drops it`);
    for (let d = 0; d <= LOW; d++)
      if (s === 0 && moves[s][d] !== 0) throw new Error(`${NAME}: digit ${d} leaves state 0`);
  }
  const witness = num.kummer.witness(P, 3, 1);
  if (witness?.[0] !== 4n || witness[1] !== 14n) throw new Error(`${NAME}: witness ${witness}, want 4 and 14`);
  const cells: Cell[] = [];
  let inK = 0;
  let design = 0;
  let endOne = 0;
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      const digits: number[] = [];
      for (let i = 0; i < LEVEL / 2; i++)
        digits.push(Math.floor(x / SPAN ** i) % SPAN, Math.floor(y / SPAN ** i) % SPAN);
      let state: number | undefined = 0;
      for (const d of digits) {
        state = moves[state][d];
        if (state === undefined) break;
      }
      const k = digits.reduce((sum, d, i) => sum + d * P ** i, 0);
      const accepted = state !== undefined && closes[state];
      if (accepted !== num.kummer.member(k, P, 3, 1))
        throw new Error(`${NAME}: ${k} read ${accepted} against the crate`);
      if (!accepted) continue;
      inK++;
      const low = digits.every((d) => d <= LOW);
      if (low) design++;
      if (state === 1) endOne++;
      cells.push({ x, y, kind: low ? 0 : state === 1 ? 1 : 2 });
    }
  }
  const total = BigInt(num.kummer.count(P, 3, 1, LEVEL));
  if (BigInt(inK) !== total || inK !== fill ** LEVEL)
    throw new Error(`${NAME}: ${inK} in K, want ${fill ** LEVEL} and the crate's ${total}`);
  if (design !== (LOW + 1) ** LEVEL) throw new Error(`${NAME}: ${design} in the design, want ${(LOW + 1) ** LEVEL}`);
  if (endOne !== (fill ** LEVEL - 1) / 3)
    throw new Error(`${NAME}: ${endOne} end in state 1, want ${(fill ** LEVEL - 1) / 3}`);
  memo = cells;
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const cells = board();
  const grid = new Grid(pen.frame(0.08), SIDE, SIDE, 0.14);
  const colors: Color[] = [ink.blue, ink.orange, ink.dim];
  for (const { x, y, kind } of cells) grid.fill(pen, x, SIDE - 1 - y, colors[kind]);
}
