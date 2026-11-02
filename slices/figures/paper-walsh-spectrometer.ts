import type { Ink, Pen } from "@mrlyprod/draw";
import { hex } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const SIDE = 7;

type Slice = ReturnType<typeof math.six.cut_design>;

let memo: Slice | undefined;

function slice() {
  if (memo) return memo;
  const cut = math.six.cut_design(105, SIDE, 1, 2);
  const [height, width] = cut.cell.shape;
  if (height !== 2 * SIDE || width < hex.row_len(SIDE, SIDE - 1))
    throw new Error(`paper-walsh-spectrometer: the cut is ${height} by ${width}, want ${2 * SIDE} rows`);
  if (hex.count(SIDE) !== 294)
    throw new Error(`paper-walsh-spectrometer: the hexagon holds ${hex.count(SIDE)} triangles, want 294`);
  if (math.six.fills(cut) !== 72)
    throw new Error(`paper-walsh-spectrometer: the cut holds ${math.six.fills(cut)} fills, want 72`);
  memo = cut;
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const cut = slice();
  const FILL = math.six.FILL();
  const frame = pen.frame(0.08);
  hex.hexagon(pen, frame, SIDE, 0, () => ink.line);
  hex.hexagon(pen, frame, SIDE, 2.6, (row, col) => (hex.at(cut, SIDE, row, col) === FILL ? null : ink.ground));
  hex.hexagon(pen, frame, SIDE, 2.6, (row, col) => (hex.at(cut, SIDE, row, col) === FILL ? ink.blue : null));
}
