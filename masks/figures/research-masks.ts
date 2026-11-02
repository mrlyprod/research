import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./research-masks.json" with { type: "json" };

const MARGIN = 0.08;
const FILL = 0.82;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { side, masks } = census as { side: number; masks: { family: string; rows: number[][] }[] };
  const mask = masks.find((m) => m.family === "Digit") ?? masks[0];
  if (mask.rows.length !== side) throw new Error(`research-masks: ${mask.rows.length} rows, want ${side}`);
  const frame = pen.frame(MARGIN);
  const cell = frame.w / side;
  const size = cell * FILL;
  const inset = (cell - size) / 2;
  for (let i = 0; i < side; i++) {
    for (let j = 0; j < side; j++) {
      if (!mask.rows[i][j]) continue;
      pen.rect(frame.x + j * cell + inset, frame.y + i * cell + inset, size, size, j <= i ? ink.blue : ink.yellow);
    }
  }
}
