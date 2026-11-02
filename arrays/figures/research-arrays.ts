import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-arrays.json" with { type: "json" };

const PX = 20;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { side, cells } = census;
  for (const cell of cells) {
    const [x, y, essential] = cell;
    if (cell.length !== 3 || !(x >= 0 && x < side && y >= 0 && y < side && (essential === 0 || essential === 1))) {
      throw new Error(`research-arrays: cell ${cell} off the ${side} by ${side} lattice`);
    }
  }
  const edge = side * PX;
  const area = frame(Math.round((pen.width - edge) / 2), Math.round((pen.height - edge) / 2), edge, edge);
  const grid = new Grid(area, side, side, 0.16);
  for (const [x, y, essential] of cells) grid.fill(pen, x, side - 1 - y, essential ? ink.orange : ink.blue);
}
