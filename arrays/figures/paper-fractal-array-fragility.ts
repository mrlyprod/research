import type { Color, Ink, Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const GENERATOR = [0, 1, 2, 4];
const BASE = 9;
const ORDERS = 4;
const WANT = [4, 9, 27, 81];
const GAPS = [0.022, 0.009, 0.0025, 0];
const ROW = 0.19;
const CORNER = 0.16;

type Cell = { x: number; w: number; essential: boolean; digits: number[] };

function flags(order: number): Map<number, boolean> {
  const sensors = Array.from(math.arrays.fractal(GENERATOR, BASE, order), Number);
  const essential = math.arrays.essential(sensors);
  if (sensors.length !== GENERATOR.length ** order)
    throw new Error(`paper-fractal-array-fragility: ${sensors.length} sensors at order ${order}`);
  return new Map(sensors.map((s, i) => [s, essential[i]]));
}

function rows(width: number): Cell[][] {
  const paired = Array.from(math.arrays.paired(GENERATOR), Number);
  if (paired.join() !== "0,1,4") throw new Error(`paper-fractal-array-fragility: U(G) = ${paired}, want 0,1,4`);
  const out: Cell[][] = [];
  let parents: { x: number; w: number; digits: number[] }[] = [{ x: 0, w: width, digits: [] }];
  for (let order = 1; order <= ORDERS; order++) {
    const gap = GAPS[order - 1] * width;
    const flag = flags(order);
    const row: Cell[] = [];
    for (const parent of parents) {
      const w = (parent.w - gap * (GENERATOR.length - 1)) / GENERATOR.length;
      GENERATOR.forEach((digit, k) => {
        const digits = [...parent.digits, digit];
        const sensor = digits.reduce((sum, d, i) => sum + d * BASE ** i, 0);
        const essential = flag.get(sensor);
        if (essential === undefined)
          throw new Error(`paper-fractal-array-fragility: ${sensor} is no sensor of order ${order}`);
        const rule = order === 1 || digits.every((d) => paired.includes(d));
        if (essential !== rule)
          throw new Error(
            `paper-fractal-array-fragility: sensor ${sensor} essential ${essential}, Theorem A says ${rule}`,
          );
        row.push({ x: parent.x + k * (w + gap), w, essential, digits });
      });
    }
    const count = row.filter((c) => c.essential).length;
    if (count !== WANT[order - 1])
      throw new Error(`paper-fractal-array-fragility: ${count} essential at order ${order}, want ${WANT[order - 1]}`);
    out.push(row);
    parents = row;
  }
  return out;
}

function runs(row: Cell[]): Cell[] {
  const out: Cell[] = [];
  for (const cell of row) {
    const last = out.at(-1);
    if (last && last.essential === cell.essential && Math.abs(last.x + last.w - cell.x) < 1e-9) last.w += cell.w;
    else out.push({ ...cell });
  }
  return out;
}

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const all = rows(frame.w);
  const band = frame.h * ROW;
  const step = (frame.h - band) / (ORDERS - 1);
  all.forEach((row, i) => {
    const y = frame.y + i * step;
    for (const cell of runs(row)) {
      const color: Color = cell.essential ? ink.orange : ink.blue;
      const r = GAPS[i] > 0 ? Math.min(cell.w, band) * CORNER : 0;
      pen.round_rect(frame.x + cell.x, y, cell.w, band, r, color);
    }
  });
}
