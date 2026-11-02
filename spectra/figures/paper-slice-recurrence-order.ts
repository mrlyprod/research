import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./paper-slice-recurrence-order.json" with { type: "json" };

export const units = {};

const DIMS = 13;

export default function draw(pen: Pen, ink: Ink) {
  const { free, proved } = census as { free: number[]; proved: number[] };
  if (free.length !== DIMS || proved.length !== DIMS)
    throw new Error(
      `paper-slice-recurrence-order: ${free.length} free and ${proved.length} proved orders, want ${DIMS}`,
    );
  const frame = pen.frame(0.08);
  const peak = free[DIMS - 1];
  const slot = frame.w / DIMS;
  const foot = frame.y + frame.h;
  const unit = frame.h / peak;
  free.forEach((order, index) => {
    const wide = slot * 0.62;
    const tall = order * unit;
    pen.rect(frame.x + index * slot + (slot - wide) / 2, foot - tall, wide, tall, ink.fade(ink.dim, 0.45));
  });
  proved.forEach((order, index) => {
    const wide = slot * 0.26;
    const tall = order * unit;
    pen.rect(frame.x + index * slot + (slot - wide) / 2, foot - tall, wide, tall, ink.blue);
  });
  pen.rect(frame.x, foot, frame.w, 2, ink.line);
}
