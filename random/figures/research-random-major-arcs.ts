import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 216;
const MOD = SIDE * SIDE;
const NAME = "research-random-major-arcs";

let memo: { plus: boolean[][]; minus: boolean[][] } | undefined;

function sign(m: number, five: boolean) {
  let s = 1;
  for (let i = 0; i < 6; i++) {
    const d = m % 6;
    if (five ? d === 5 : d !== 0) s = -s;
    m = Math.floor(m / 6);
  }
  return s;
}

function place(m: number) {
  const d: number[] = [];
  for (let i = 0; i < 6; i++) {
    d.push(m % 6);
    m = Math.floor(m / 6);
  }
  return { row: 36 * d[1] + 6 * d[3] + d[5], col: 36 * d[0] + 6 * d[2] + d[4] };
}

function masks() {
  if (memo) return memo;
  const roots = new Array<number>(MOD).fill(0);
  for (let n = 0; n < MOD; n++) roots[(n * n) % MOD]++;
  let squares = 0;
  let five = 0;
  let rest = 0;
  for (let m = 0; m < MOD; m++) {
    if (!roots[m]) continue;
    squares++;
    five += roots[m] * sign(m, true);
    rest += roots[m] * sign(m, false);
  }
  if (squares !== 3288) throw new Error(`${NAME}: ${squares} squares mod 6^6, want 3288`);
  if (five !== 14168 || rest !== 2192) throw new Error(`${NAME}: local sums ${five}, ${rest}, want 14168, 2192`);
  const grid = () => Array.from({ length: SIDE }, () => new Array<boolean>(SIDE).fill(false));
  const plus = grid();
  const minus = grid();
  for (let m = 0; m < MOD; m++) {
    if (!roots[m]) continue;
    const { row, col } = place(m);
    if (col === row) continue;
    const s = sign(m, col > row);
    (s > 0 ? plus : minus)[row][col] = true;
  }
  memo = { plus, minus };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { plus, minus } = masks();
  const area = pen.frame(0.07);
  const side = Math.min(area.w, area.h);
  const cell = Math.floor(side / SIDE);
  const span = cell >= 1 ? cell * SIDE : side;
  const box = frame(area.x + (area.w - span) / 2, area.y + (area.h - span) / 2, span, span);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, plus, ink.blue);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, minus, ink.orange);
}
