import { type Color, type Frame, frame, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";
import census from "./paper-novelty-meter.json" with { type: "json" };

const LOW = 8;
const HIGH = 20.5;
const PER_OCTAVE = 16;
const ZEROS = 29;
const CURVE = 2400;

export const units = { num };

let memo: number[] | undefined;

function wave() {
  if (memo) return memo;
  const { smooth, rough } = census;
  const samples = Math.round((HIGH - LOW) * PER_OCTAVE) + 1;
  if (smooth.length !== samples || rough.length !== samples)
    throw new Error(`paper-novelty-meter: ${smooth.length} smooth and ${rough.length} rough samples, want ${samples}`);
  const line = new num.zeta.Line();
  const gammas = line.zeros(ZEROS);
  if (!(Math.abs(gammas[0] - 14.134725) < 1e-6))
    throw new Error(`paper-novelty-meter: first zero ${gammas[0]}, want 14.134725`);
  if (!(Math.abs(gammas[ZEROS - 1] - 98.831194) < 1e-5))
    throw new Error(`paper-novelty-meter: last zero ${gammas[ZEROS - 1]}, want 98.831194`);
  const coef = line.novelty_coefficients(gammas);
  if (!(Math.abs(num.zeta.Complex.abs(coef[0]) - 0.1879) < 5e-4))
    throw new Error(`paper-novelty-meter: first coefficient ${num.zeta.Complex.abs(coef[0])}, want 0.1879`);
  if (!(Math.abs(num.zeta.Complex.abs(coef[9]) - 4.286e-3) < 5e-6))
    throw new Error(`paper-novelty-meter: tenth coefficient ${num.zeta.Complex.abs(coef[9])}, want 0.004286`);
  const at = (j: number) => num.zeta.novelty_wave(gammas, coef, -j * Math.LN2);
  const peak = smooth.reduce((a, v) => Math.max(a, Math.abs(v)), 0);
  const miss = smooth.reduce((m, v, k) => Math.max(m, Math.abs(v - at(LOW + k / PER_OCTAVE))), 0);
  if (!(miss / peak < 1e-2)) throw new Error(`paper-novelty-meter: the zeros miss ${miss / peak}, want under 0.01`);
  memo = Array.from({ length: CURVE + 1 }, (_, k) => at(LOW + (HIGH - LOW) * (k / CURVE)));
  return memo;
}

function beads(pen: Pen, area: Frame, values: number[], span: number, r: number, color: Color) {
  const last = values.length - 1;
  plot.dots(
    pen,
    values.map((v, i): Point => [area.x + (area.w * i) / last, area.y + area.h * (0.5 - (0.5 * v) / span)]),
    r,
    color,
  );
}

function reach(values: number[]) {
  return values.reduce((a, v) => Math.max(a, Math.abs(v)), 0) * 1.12;
}

export default function draw(pen: Pen, ink: Ink) {
  const { smooth, rough } = census;
  const values = wave();
  const area = pen.frame(0.08);
  const gap = 28;
  const topH = area.h * 0.62 - gap / 2;
  const top = frame(area.x, area.y, area.w, topH);
  const foot = frame(area.x, area.y + topH + gap, area.w, area.h - topH - gap);
  for (const box of [top, foot]) {
    plot.axis(pen, box, ink.line);
    const mid = box.y + box.h / 2;
    pen.segment([box.x, mid], [box.x + box.w, mid], 1, ink.line);
  }
  const innerTop = top.inset(16);
  const innerFoot = foot.inset(16);
  const spanTop = reach(smooth);
  const spanFoot = reach(rough);
  const curve = values.map((v, k): Point => {
    const t = k / CURVE;
    return [innerTop.x + innerTop.w * t, innerTop.y + innerTop.h * (0.5 - (0.5 * v) / spanTop)];
  });
  pen.polyline(curve, 2, ink.orange);
  beads(pen, innerTop, smooth, spanTop, 3.4, ink.blue);
  beads(pen, innerFoot, rough, spanFoot, 3.4, ink.blue);
}
