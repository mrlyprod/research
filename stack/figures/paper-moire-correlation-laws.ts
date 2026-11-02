import type { Ink, Pen } from "@mrlyprod/draw";
import { field, frame } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const SCALES = 31;
const CELL = 28;

let memo: Float64Array | undefined;

function values() {
  if (memo) return memo;
  const odds: number[] = [];
  for (let m = 3; m <= 63; m += 2) odds.push(m);
  if (odds.length !== SCALES)
    throw new Error(`paper-moire-correlation-laws: ${odds.length} odd scales, want ${SCALES}`);
  const correlation = math.moire.pairs.correlation;
  if (correlation(3, 5) !== 0)
    throw new Error(`paper-moire-correlation-laws: correlation(3, 5) is ${correlation(3, 5)}, want 0`);
  if (!(Math.abs(correlation(9, 9) - 1) < 1e-12))
    throw new Error(`paper-moire-correlation-laws: correlation(9, 9) is ${correlation(9, 9)}, want 1`);
  const grid = new Float64Array(SCALES * SCALES);
  let at = 0;
  for (const m of odds) {
    for (const n of odds) {
      const r = correlation(m, n);
      grid[at++] = r <= 0 ? 0 : 0.5 + 0.5 * Math.sqrt(r);
    }
  }
  const clear = grid.reduce((count, v) => count + (v === 0 ? 1 : 0), 0);
  if (clear !== 762) throw new Error(`paper-moire-correlation-laws: ${clear} clear pairs, want 762`);
  memo = grid;
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const ramp = new ink.Ramp([ink.ground, ink.blue, ink.yellow]);
  const side = CELL * SCALES;
  const edge = (pen.width - side) / 2;
  field.draw_range(pen, frame(edge, edge, side, side), SCALES, SCALES, values(), [0, 1], ramp);
}
