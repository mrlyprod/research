import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const universe = math.bang.bang(2);
  if (universe.total !== 16) throw new Error(`research-core: ${universe.total} codes, want 16`);
  if (universe.distinct() !== 6) throw new Error(`research-core: ${universe.distinct()} orbits, want 6`);
  const canonical = universe.canonical().map((d) => String(d.i));
  if (canonical.join() !== "0,1,3,6,7,15") throw new Error(`research-core: canonical ${canonical}, want 0,1,3,6,7,15`);
  const gutter = area.w * 0.03;
  const side = (area.w - 3 * gutter) / 4;
  for (let code = 0; code < 16; code++) {
    const color = universe.design(code).canonical ? ink.yellow : ink.blue;
    const x = area.x + (code % 4) * (side + gutter);
    const y = area.y + Math.floor(code / 4) * (side + gutter);
    const cell = frame(x, y, side, side);
    pen.rect(cell.x, cell.y, cell.w, cell.h, ink.panel);
    const cells = math.two.create(code, 2, 4, 0, 2);
    if (cells.shape[1] !== 16) throw new Error(`research-core: design ${code} is ${cells.shape[1]} wide, want 16`);
    new Grid(cell, 16, 16, 0).paint(pen, cells, (kind) => (kind !== 0 ? color : null));
  }
}
