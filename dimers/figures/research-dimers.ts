import type { Color, Ink, Pen } from "@mrlyprod/draw";
import census from "./research-dimers.json" with { type: "json" };

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { side, dominoes } = census as { side: number; dominoes: [number, number, number][] };
  for (const [r, c, flat] of dominoes) {
    if (r < 0 || c < 0 || r + 1 - flat >= side || c + flat >= side)
      throw new Error(`research-dimers: domino ${r}, ${c} off the board`);
  }
  const frame = pen.frame(0.08);
  const unit = frame.w / side;
  const pad = unit * 0.11;
  const thick = unit - 2 * pad;
  for (const [r, c, flat] of dominoes) {
    const x = frame.x + c * unit;
    const y = frame.y + r * unit;
    const [w, h, color]: [number, number, Color] =
      flat === 1 ? [2 * unit, unit, ink.blue] : [unit, 2 * unit, ink.orange];
    pen.round_rect(x + pad, y + pad, w - 2 * pad, h - 2 * pad, thick * 0.3, color);
  }
}
