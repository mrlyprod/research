import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./research-dilations.json" with { type: "json" };

const LAGS = 31;
const MARGIN = 0.08;
const FILL = 0.46;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { field } = census as { field: number[] };
  if (field.length !== LAGS) throw new Error(`research-dilations: ${field.length} lags, want ${LAGS}`);
  const peak = Math.max(...field);
  const frame = pen.frame(MARGIN);
  const cell = frame.w / LAGS;
  const top = cell * FILL;
  for (let m = 0; m < LAGS; m++) {
    for (let n = 0; n < LAGS; n++) {
      const count = field[m] * field[n];
      const r = top * Math.sqrt(count / (peak * peak));
      pen.disc(
        frame.x + (n + 0.5) * cell,
        frame.y + (m + 0.5) * cell,
        r,
        count === peak * peak ? ink.yellow : ink.blue,
      );
    }
  }
}
