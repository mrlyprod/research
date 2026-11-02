import { type Color, type Frame, frame, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";
import census from "./paper-franel-converse.json" with { type: "json" };

const BASE = 3;
const DIGITS = [0, 1];
const ORDER = 81;

export const units = { num, math };

type Node = [number, number];

let memo: { thin: Node[]; full: Node[] } | undefined;

// FACTS

function nodes() {
  if (memo) return memo;
  const design = Array.from(math.design.elements(BASE, DIGITS, 5), Number).filter((n) => n <= ORDER);
  if (design.length !== 16) throw new Error(`paper-franel-converse: ${design.length} design elements, want 16`);
  if (design.slice(0, 5).join() !== "1,3,4,9,10")
    throw new Error(`paper-franel-converse: design starts ${design.slice(0, 5)}, want 1,3,4,9,10`);
  const held = new Set(design);
  const full: Node[] = num.lattice
    .farey(ORDER)
    .filter((n) => n.num > 0)
    .map((n) => [n.num, n.den]);
  const thin = full.filter(([, b]) => held.has(b));
  if (thin.length !== 241 || census.thin.length !== 241)
    throw new Error(`paper-franel-converse: ${thin.length} thin nodes and ${census.thin.length} thin deltas, want 241`);
  if (full.length !== 2020 || census.full.length !== 2020)
    throw new Error(
      `paper-franel-converse: ${full.length} full nodes and ${census.full.length} full deltas, want 2020`,
    );
  memo = { thin, full };
  return memo;
}

// MARKS

function ticks(pen: Pen, area: Frame, list: Node[], thick: number, color: Color) {
  const reach = Math.log(ORDER);
  const foot = area.y + area.h;
  for (const [a, b] of list) {
    const x = area.x + (area.w * a) / b;
    const h = area.h * (1 - (0.82 * Math.log(b)) / reach);
    pen.segment([x, foot], [x, foot - h], thick, color);
  }
}

function trace(area: Frame, list: Node[], delta: number[], span: number): Point[] {
  return list.map(([a, b], i) => [area.x + (area.w * a) / b, area.y + area.h * (0.5 - (0.5 * delta[i]) / span)]);
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { thin, full } = nodes();
  const area = pen.frame(0.08);
  const gap = 26;
  const row = (area.h - 2 * gap) * 0.21;
  const top = frame(area.x, area.y, area.w, row);
  const mid = frame(area.x, area.y + row + gap, area.w, row);
  const footY = area.y + 2 * (row + gap);
  const foot = frame(area.x, footY, area.w, area.y + area.h - footY);

  plot.baseline(pen, top, ink.line);
  plot.baseline(pen, mid, ink.line);
  plot.axis(pen, foot, ink.line);
  const zero = foot.y + foot.h / 2;
  pen.segment([foot.x, zero], [foot.x + foot.w, zero], 1, ink.line);

  ticks(pen, top.inset(2), thin, 2.2, ink.blue);
  ticks(pen, mid.inset(2), full, 1, ink.fade(ink.dim, 0.55));

  const inner = foot.inset(16);
  const span = census.thin.reduce((a, v) => Math.max(a, Math.abs(v)), 0) * 1.12;
  const beads = trace(inner, thin, census.thin, span);
  pen.polyline(trace(inner, full, census.full, span), 1.2, ink.fade(ink.dim, 0.7));
  pen.polyline(beads, 2.4, ink.orange);
  plot.dots(pen, beads, 2.6, ink.blue);
}
