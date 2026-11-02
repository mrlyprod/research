import type { Ink, Pen, Point } from "@mrlyprod/draw";
import census from "./research-toolpaths.json" with { type: "json" };

const STEPS = 2401;
const SHARP = 686;
const MARGIN = 0.08;
const PATH = 2.2;
const DOT = 3.1;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { path, sharp } = census as { path: Point[]; sharp: number[] };
  if (path.length !== STEPS + 1) throw new Error(`research-toolpaths: ${path.length} points, want ${STEPS + 1}`);
  if (sharp.length !== SHARP) throw new Error(`research-toolpaths: ${sharp.length} sharp vertices, want ${SHARP}`);
  let x0 = Infinity;
  let y0 = Infinity;
  let x1 = -Infinity;
  let y1 = -Infinity;
  for (const [x, y] of path) {
    x0 = Math.min(x0, x);
    y0 = Math.min(y0, y);
    x1 = Math.max(x1, x);
    y1 = Math.max(y1, y);
  }
  const area = pen.frame(MARGIN);
  const [cx, cy] = area.center();
  const scale = (2 * area.radius()) / Math.max(x1 - x0, y1 - y0);
  const [mx, my] = [(x0 + x1) / 2, (y0 + y1) / 2];
  const at = ([x, y]: Point): Point => [cx + (x - mx) * scale, cy - (y - my) * scale];
  pen.polyline(path.map(at), PATH, ink.dim);
  for (const i of sharp) {
    const [x, y] = at(path[i]);
    pen.disc(x, y, DOT, ink.orange);
  }
}
