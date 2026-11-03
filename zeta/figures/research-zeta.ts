import { type Frame, frame, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import census from "./research-zeta.json" with { type: "json" };

type Panel = { columns: number[]; ordinates: number[]; family: number[][]; teeth: number[][]; hollow: number[][] };

const HEIGHT = 40;
const WIDE = 1.2;
const REACH = 0.7;

export const units = {};

function check(p: Panel, sizes: [number, number, number]) {
  const counts = [p.family.length, p.teeth.length, p.hollow.length];
  if (counts.some((n, i) => n !== sizes[i])) throw new Error(`research-zeta: zeros ${counts}, want ${sizes}`);
  if (p.columns.length !== 2) throw new Error(`research-zeta: ${p.columns.length} pole columns, want 2`);
  const [alpha] = p.columns;
  for (const [re, im] of [...p.family, ...p.teeth, ...p.hollow]) {
    if (!(re > alpha - WIDE && re < alpha + REACH && im > 0 && im < HEIGHT))
      throw new Error(`research-zeta: zero ${re} ${im} off the strip`);
  }
}

function panel(pen: Pen, ink: Ink, box: Frame, p: Panel) {
  const [alpha, below] = p.columns;
  const lo = alpha - WIDE;
  const hi = alpha + REACH;
  const at = (re: number, im: number): Point => [
    box.x + (box.w * (re - lo)) / (hi - lo),
    box.y + box.h * (1 - im / HEIGHT),
  ];
  pen.rect(box.x, box.y, box.w, box.h, ink.panel);
  plot.axis(pen, box, ink.line);
  pen.segment(at(alpha, 0), at(alpha, HEIGHT), 1.6, ink.fade(ink.dim, 0.5));
  pen.segment(at(below, 0), at(below, HEIGHT), 1.6, ink.fade(ink.dim, 0.3));
  for (const t of p.ordinates) {
    if (t < 0.4) continue;
    for (const line of p.columns) {
      const [x, y] = at(line, t);
      pen.ring(x, y, 11, 1.6, ink.fade(ink.dim, 0.9));
    }
  }
  for (const [re, im] of p.family) pen.disc(...at(re, im), 7, ink.blue);
  for (const [re, im] of p.teeth) pen.disc(...at(re, im), 7, ink.yellow);
  for (const [re, im] of p.hollow) pen.disc(...at(re, im), 7, ink.fade(ink.fg, 0.9));
}

export default function draw(pen: Pen, ink: Ink) {
  const { design, full } = census as { design: Panel; full: Panel };
  check(design, [7, 6, 0]);
  check(full, [6, 0, 4]);
  const area = pen.area(0.08);
  const wide = area.w * 0.45;
  const tall = area.h * 0.8;
  panel(pen, ink, frame(area.x, area.y, wide, tall), design);
  panel(pen, ink, frame(area.x + area.w - wide, area.y + area.h - tall, wide, tall), full);
}
