import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-sequences-atlas.json" with { type: "json" };

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { side, cells } = census;
  if (!(side > 0 && cells.length <= side * side && /^[0-3]*$/.test(cells))) {
    throw new Error(`research-sequences-atlas: ${cells.length} cells do not fit a ${side} by ${side} square`);
  }
  const inks = [ink.blue, ink.dim, ink.orange];
  const edge = Math.min(pen.width, pen.height) * 0.86;
  const area = frame((pen.width - edge) / 2, (pen.height - edge) / 2, edge, edge);
  const grid = new Grid(area, side, side, 0.12);
  for (let i = 0; i < cells.length; i++) {
    const state = cells.charCodeAt(i) - 48;
    if (state < 3) grid.fill(pen, i % side, Math.floor(i / side), inks[state]);
  }
}
