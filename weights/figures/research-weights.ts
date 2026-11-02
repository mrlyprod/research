import { type Color, type Frame, frame, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import census from "./research-weights.json" with { type: "json" };

const SAMPLES = 720;

export const units = {};

// PANEL

function center(ys: number[]) {
  let sum = 0;
  for (const y of ys) sum += y;
  const mean = sum / ys.length;
  return ys.map((y) => y - mean);
}

function reach(first: number[], second: number[]) {
  let top = 0;
  for (const y of [...first, ...second]) top = Math.max(top, Math.abs(y));
  return top * 1.14;
}

function trace(pen: Pen, area: Frame, ys: number[], span: number, thick: number, color: Color) {
  const last = ys.length - 1;
  const pts = ys.map((y, i): Point => [area.x + (area.w * i) / last, area.y + area.h * (0.5 - (0.5 * y) / span)]);
  pen.polyline(pts, thick, color);
}

function panel(pen: Pen, ink: Ink, area: Frame, first: number[], second: number[]) {
  plot.axis(pen, area, ink.line);
  const inner = area.inset(18);
  const mid = inner.y + inner.h / 2;
  pen.segment([inner.x, mid], [inner.x + inner.w, mid], 1, ink.line);
  const span = reach(first, second);
  trace(pen, inner, second, span, 3.4, ink.orange);
  trace(pen, inner, first, span, 2.1, ink.blue);
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { lattice, smooth, thin, fat } = census as Record<string, number[]>;
  for (const ys of [lattice, smooth, thin, fat]) {
    if (ys?.length !== SAMPLES) throw new Error(`research-weights: a curve of ${ys?.length} samples, want ${SAMPLES}`);
  }
  const box = pen.frame(0.08);
  const half = box.h / 2;
  panel(pen, ink, frame(box.x, box.y, box.w, half).inset(14), center(lattice), center(smooth));
  panel(pen, ink, frame(box.x, box.y + half, box.w, half).inset(14), center(thin), center(fat));
}
