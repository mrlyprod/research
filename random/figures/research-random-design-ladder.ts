import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as life from "@mrlyprod/life";

export const units = { life };

const STEPS = 63;
const WIDTH = 2 * STEPS + 1;
const NAME = "research-random-design-ladder";

let memo: { cell: boolean[][]; thue: boolean[][] } | undefined;

function parity(n: number) {
  let x = n;
  let p = 0;
  while (x > 0) {
    p ^= x & 1;
    x >>= 1;
  }
  return p;
}

function mask(cells: life.Tensor, width: number) {
  if (cells.shape.length !== 2 || cells.shape[0] !== STEPS + 1 || cells.shape[1] !== width)
    throw new Error(`${NAME}: diagram ${cells.shape}, want ${STEPS + 1},${width}`);
  const out: boolean[][] = [];
  for (let t = 0; t <= STEPS; t++) {
    out.push([]);
    for (let col = 0; col < WIDTH; col++) out[t].push(cells.data[t * width + col] !== 0);
  }
  return out;
}

function seeded(rule: number) {
  const width = 3 * STEPS + 1;
  const row = Array.from({ length: width }, (_, i) => (i >= STEPS ? parity(i - STEPS) : 0));
  return mask(life.history(row, rule, STEPS, false), width);
}

function diagrams() {
  if (memo) return memo;
  const cell = mask(life.single_seed(18, STEPS), WIDTH);
  const gasket = mask(life.single_seed(90, STEPS), WIDTH);
  for (let t = 0; t <= STEPS; t++) {
    for (let col = 0; col < WIDTH; col++) {
      const k = Math.abs(col - STEPS);
      const pascal = k <= t && (t + k) % 2 === 0 && (((t + k) / 2) & ((t - k) / 2)) === 0;
      if (cell[t][col] !== pascal) throw new Error(`${NAME}: rule 18 leaves Pascal mod 2 at ${t},${col}`);
      if (gasket[t][col] !== pascal) throw new Error(`${NAME}: rule 90 leaves Pascal mod 2 at ${t},${col}`);
    }
  }
  const thue = seeded(18);
  const linear = seeded(90);
  let apart = 0;
  for (let t = 0; t <= STEPS; t++) for (let col = 0; col < WIDTH; col++) if (thue[t][col] !== linear[t][col]) apart++;
  if (apart === 0) throw new Error(`${NAME}: rule 18 never leaves rule 90 on the Thue-Morse seed`);
  memo = { cell, thue };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { cell, thue } = diagrams();
  const area = pen.frame(0.07);
  const rows = 2 * (STEPS + 1);
  const gutter = area.h * 0.04;
  const size = Math.floor(Math.min(area.w / WIDTH, (area.h - gutter) / rows));
  const w = size * WIDTH;
  const h = size * (STEPS + 1);
  const left = area.x + (area.w - w) / 2;
  const top = area.y + (area.h - 2 * h - gutter) / 2;
  new Grid(frame(left, top, w, h), WIDTH, STEPS + 1, 0).carpet(pen, cell, ink.blue);
  new Grid(frame(left, top + h + gutter, w, h), WIDTH, STEPS + 1, 0).carpet(pen, thue, ink.orange);
}
