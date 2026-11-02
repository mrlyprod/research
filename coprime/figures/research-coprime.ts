import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";

export const units = { math, num };

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const cells = math.two.create(495, 3, 4, 0, 3);
  const side = cells.shape[1];
  if (side !== 81 || cells.shape[0] !== side) throw new Error(`research-coprime: shape ${cells.shape}, want 81 square`);
  let sum = 0;
  for (const kind of cells.types) sum += kind;
  if (sum !== 4096) throw new Error(`research-coprime: ${sum} filled, want 4096`);
  const grid = new Grid(area, side, side, 0.1);
  let lit = 0;
  for (let row = 0; row < side; row++) {
    for (let col = 0; col < side; col++) {
      if (cells.types[row * side + col] === 0) continue;
      const coprime = num.factor.gcd(col + 1, side - row) === "1";
      if (coprime) lit++;
      grid.fill(pen, col, row, coprime ? ink.yellow : ink.dim);
    }
  }
  if (!(lit > 0 && lit < 4096)) throw new Error(`research-coprime: ${lit} coprime of 4096, want some but not all`);
}
