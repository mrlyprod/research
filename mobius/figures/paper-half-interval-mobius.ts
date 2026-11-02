import type { Color, Ink, Pen, Point } from "@mrlyprod/draw";
import census from "./paper-half-interval-mobius.json" with { type: "json" };

export const units = {};

const LOW = 2;
const HIGH = 7;
const FIFTH = 94939;
const QUARTER = 3789;
const BASES = 1551;

export default function draw(pen: Pen, ink: Ink) {
  const { bases, q25, q20 } = census as { bases: number[]; q25: number[]; q20: number[] };
  if (bases.length !== BASES || q25.length !== BASES || q20.length !== BASES)
    throw new Error(
      `paper-half-interval-mobius: ${bases.length} bases, ${q25.length} and ${q20.length} margins, want ${BASES}`,
    );
  if (bases.some((b, i) => i > 0 && b <= bases[i - 1]))
    throw new Error("paper-half-interval-mobius: bases do not ascend");
  if (!bases.includes(QUARTER) || !bases.includes(FIFTH))
    throw new Error("paper-half-interval-mobius: a wall is not among the bases");
  const frame = pen.frame(0.08);
  const top = q25[BASES - 1];
  const bottom = q20[0];
  const pad = 0.04 * (top - bottom);
  const hi = top + pad;
  const lo = bottom - pad;
  const at = (b: number, m: number): Point => {
    const u = (Math.log10(b) - LOW) / (HIGH - LOW);
    return [frame.x + frame.w * u, frame.y + (frame.h * (hi - m)) / (hi - lo)];
  };

  for (let decade = 3; decade < 7; decade++) {
    const x = frame.x + (frame.w * (decade - LOW)) / (HIGH - LOW);
    pen.segment([x, frame.y], [x, frame.y + frame.h], 1.4, ink.line);
  }
  const zero = at(bases[0], 0)[1];
  pen.segment([frame.x, zero], [frame.x + frame.w, zero], 2.2, ink.dim);

  const bars: [number[], number, Color][] = [
    [q25, QUARTER, ink.orange],
    [q20, FIFTH, ink.blue],
  ];
  for (const [margins, wall, color] of bars) {
    pen.polyline(
      bases.map((b, i) => at(b, margins[i])),
      5,
      color,
    );
    const [x, y] = at(wall, 0);
    pen.disc(x, y, 15, ink.ground);
    pen.disc(x, y, 9.5, color);
  }
}
