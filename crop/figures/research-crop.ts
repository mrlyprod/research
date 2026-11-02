import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const CUT = 1;
const IN = 2;

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const cells = math.two.create(7, 3, 4, 0, 2);
  const side = cells.shape[1];
  if (side !== 81) throw new Error(`research-crop: side ${side}, want 81`);
  if (cells.types.length !== side * side)
    throw new Error(`research-crop: ${cells.types.length} cell bytes, want ${side * side}`);
  const types = { shape: cells.shape, data: cells.types };
  let sum = 0;
  for (const kind of cells.types) sum += kind;
  if (sum !== 4096) throw new Error(`research-crop: ${sum} filled, want 4096`);
  const ball = math.shape.named("ball", 2, new math.shape.Frac(1, 2));
  const tally = math.shape.census(ball, types);
  if (tally.filled[IN] !== 2908) throw new Error(`research-crop: ${tally.filled[IN]} filled inside, want 2908`);
  if (tally.filled[CUT] !== 204) throw new Error(`research-crop: ${tally.filled[CUT]} filled cut, want 204`);
  const map = math.shape.regions(ball, [side, side]);
  if (map.data.length !== side * side)
    throw new Error(`research-crop: ${map.data.length} regions, want ${side * side}`);
  const grid = new Grid(area, side, side, 0.1);
  for (let row = 0; row < side; row++) {
    for (let col = 0; col < side; col++) {
      if (cells.types[row * side + col] === 0) continue;
      const region = map.data[row * side + col];
      if (region === IN) grid.fill(pen, col, row, ink.yellow);
      else if (region === CUT) grid.fill(pen, col, row, ink.orange);
    }
  }
}
