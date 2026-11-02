import type { Ink, Pen, Point } from "@mrlyprod/draw";
import census from "./paper-slice-sign-even-half.json" with { type: "json" };

const DIMS = 20;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { three, others } = census as { three: number[]; others: number[][] };
  if (three.length !== DIMS) throw new Error(`paper-slice-sign-even-half: ${three.length} depths, want ${DIMS}`);
  if (others.length !== 4 || others.some((row) => row.length !== DIMS))
    throw new Error("paper-slice-sign-even-half: others are not four rows of 20");
  const frame = pen.frame(0.08);
  const peak = three[DIMS - 1];
  const at = (dim: number, depth: number): Point => [
    frame.x + (frame.w * (dim - 2)) / 38,
    frame.y + frame.h * (1 - depth / peak),
  ];
  pen.rect(frame.x, frame.y + frame.h, frame.w, 2, ink.line);

  others.forEach((row, order) => {
    const steps: Point[] = [];
    row.forEach((value, index) => {
      const dim = 2 + 2 * index;
      steps.push(at(dim, value));
      steps.push(at(Math.min(dim + 2, 40), value));
    });
    pen.polyline(steps, 3, ink.fade(ink.dim, 0.95 - 0.1 * order));
  });

  const curve = three.map((k, index) => at(2 + 2 * index, k));
  pen.polyline(curve, 3, ink.orange);
  for (const [x, y] of curve) pen.disc(x, y, 8, ink.orange);
}
