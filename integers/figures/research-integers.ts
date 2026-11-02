import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-integers.json" with { type: "json" };

const WINDOW = 10000;
const SIDE = 100;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const counts: number[] = census;
  if (counts.length !== WINDOW + 1) throw new Error(`research-integers: ${counts.length} counts, want ${WINDOW + 1}`);
  const peak = Math.max(...counts);
  if (!(peak > 0)) throw new Error("research-integers: no row writes any integer");
  const ramp = ink.Ramp.tone(ink.dim, ink.yellow);
  const grid = new Grid(pen.frame(0.08), SIDE, SIDE, 0.12);
  const scale = Math.log(1 + peak);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      const count = counts[row * SIDE + col + 1];
      if (count === 0) continue;
      grid.fill(pen, col, row, ramp.at(Math.log(1 + count) / scale));
    }
  }
}
