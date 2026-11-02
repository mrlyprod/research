import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./paper-first-base-below-a-quarter.json" with { type: "json" };

export const units = {};

const LOW = 10;
const HIGH = 40;
const QUARTER = 0.25;
const RUNGS = 395;

export default function draw(pen: Pen, ink: Ink) {
  const { rungs } = census as { rungs: [number, number][] };
  if (rungs.length !== RUNGS) throw new Error(`paper-first-base-below-a-quarter: ${rungs.length} rungs, want ${RUNGS}`);
  if (rungs.some(([q]) => !Number.isInteger(q) || q < LOW || q > HIGH))
    throw new Error(`paper-first-base-below-a-quarter: a rung base is outside ${LOW}..${HIGH}`);
  const lo = Math.min(...rungs.map((r) => r[1]));
  const hi = Math.max(...rungs.map((r) => r[1]));
  const pad = (hi - lo) * 0.06;
  const foot = lo - pad;
  const head = hi + pad;
  if (!(foot < QUARTER && QUARTER < head))
    throw new Error("paper-first-base-below-a-quarter: the quarter line is outside");

  const frame = pen.frame(0.08);
  const slot = frame.w / (HIGH - LOW + 1);
  const dash = slot * 0.72;
  const place = (e: number) => frame.y + (frame.h * (head - e)) / (head - foot);

  const line = place(QUARTER);
  pen.rect(frame.x, line - 1, frame.w, 2, ink.fade(ink.blue, 0.9));

  for (const [q, e] of rungs) {
    const x = frame.x + (q - LOW + 0.5) * slot - dash / 2;
    const y = place(e);
    const under = e < QUARTER;
    const thick = under ? 3.8 : 2.6;
    pen.rect(x, y - thick / 2, dash, thick, under ? ink.yellow : ink.fade(ink.dim, 0.7));
  }
}
