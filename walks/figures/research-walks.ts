import { type Color, type Frame, frame, type Ink, type Pen, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./research-walks.json" with { type: "json" };

const LEVEL = 4;
const SIDE = 3 ** LEVEL;
const SITES = 7 ** LEVEL;

export const units = { math };

function panel(pen: Pen, ink: Ink, area: Frame, code: number, trail: number[], color: Color) {
  const cell = math.two.create(code, 3, LEVEL, 0, 3);
  if (cell.shape.length !== 2 || cell.shape[0] !== SIDE || cell.shape[1] !== SIDE)
    throw new Error(`research-walks: code ${code} is ${cell.shape}, want ${SIDE} by ${SIDE}`);
  const mask = Array.from(cell.types, (kind) => kind !== 0);
  const sites = mask.filter(Boolean).length;
  if (sites !== SITES) throw new Error(`research-walks: code ${code} holds ${sites} sites, want ${SITES}`);
  const seen = new Set<number>();
  for (const at of trail) {
    if (!Number.isInteger(at) || at < 0 || at >= mask.length || !mask[at])
      throw new Error(`research-walks: code ${code} visits ${at}, not a site`);
    if (seen.has(at)) throw new Error(`research-walks: code ${code} visits ${at} twice`);
    seen.add(at);
  }
  const step = area.w / SIDE;
  const dot = step * 0.7;
  const faint = ink.fade(ink.dim, 0.55);
  mask.forEach((on, at) => {
    if (!on) return;
    const row = Math.floor(at / SIDE);
    const col = at % SIDE;
    pen.rect(area.x + col * step + (step - dot) / 2, area.y + row * step + (step - dot) / 2, dot, dot, faint);
  });
  for (const at of trail) {
    const row = Math.floor(at / SIDE);
    const col = at % SIDE;
    pen.rect(area.x + col * step, area.y + row * step, step, step, color);
  }
}

export default function draw(pen: Pen, ink: Ink) {
  const walks = census as { "127": number[]; "239": number[] };
  const box = pen.frame(0.08);
  const half = box.w / 2;
  const first = frame(box.x, box.y, half, half).inset(8);
  const second = frame(box.x + half, box.y + half, half, half).inset(8);
  for (const area of [first, second]) plot.axis(pen, area, ink.line);
  panel(pen, ink, first.inset(10), 127, walks["127"], ink.blue);
  panel(pen, ink, second.inset(10), 239, walks["239"], ink.orange);
}
