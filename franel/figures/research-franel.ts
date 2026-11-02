import type { Color, Ink, Pen } from "@mrlyprod/draw";
import { Grid } from "@mrlyprod/draw";
import census from "./research-franel.json" with { type: "json" };

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { order, low, high, kernel, terms } = census as {
    order: number;
    low: number;
    high: number;
    kernel: number[];
    terms: [number, number, number][];
  };
  if (kernel.length !== (order * (order + 1)) / 2)
    throw new Error(
      `research-franel: ${kernel.length} kernel values at order ${order}, want ${(order * (order + 1)) / 2}`,
    );
  if (terms.some(([d, e]) => d < 1 || e < 1 || d > order || e > order))
    throw new Error("research-franel: a term off the lattice");
  const weight = (d: number, e: number) =>
    d >= e ? kernel[((d - 1) * d) / 2 + e - 1] : kernel[((e - 1) * e) / 2 + d - 1];
  const wash = (k: number): Color => ink.mix(ink.ground, ink.dim, 0.015 + 0.7 * Math.sqrt(k));
  const lit = (v: number): Color => {
    const s = (Math.log(Math.abs(v)) - Math.log(low)) / (Math.log(high) - Math.log(low));
    return ink.mix(ink.ground, v > 0 ? ink.blue : ink.orange, 0.55 + 0.45 * s);
  };
  const frame = pen.frame(0.08);
  const fine = new Grid(frame, order, order, 0.34);
  const bold = new Grid(frame, order, order, 0.08);
  for (let d = 1; d <= order; d++) {
    for (let e = 1; e <= order; e++) fine.fill(pen, e - 1, d - 1, wash(weight(d, e)));
  }
  for (const [d, e, v] of terms) bold.fill(pen, e - 1, d - 1, lit(v));
}
