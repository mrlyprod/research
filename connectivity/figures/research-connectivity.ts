import * as core from "@mrlyprod/core";
import { frame as box, type Color, type Frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { core, math };

const SIDE = 64;
const CELLS = 729;
const SEED = 20260902;

// MASKS

function scatter() {
  const mask = new Array<boolean>(SIDE * SIDE).fill(false);
  for (const slot of new core.Rng(SEED).sample_indices(SIDE * SIDE, CELLS)) mask[slot] = true;
  const on = mask.filter(Boolean).length;
  if (on !== CELLS) throw new Error(`research-connectivity: scatter holds ${on} cells, want ${CELLS}`);
  return mask;
}

function gasket() {
  const cells = math.two.create(7, 2, 6, 0, 2);
  if (cells.shape[1] !== SIDE) throw new Error(`research-connectivity: gasket ${cells.shape[1]} wide, want ${SIDE}`);
  const mask = Array.from(cells.types, (kind) => kind !== 0);
  const on = mask.filter(Boolean).length;
  if (on !== CELLS) throw new Error(`research-connectivity: gasket holds ${on} cells, want ${CELLS}`);
  return mask;
}

// DRAW

function plate(pen: Pen, frame: Frame, mask: boolean[], color: Color, rule: number, ink: Ink) {
  const edge = frame.inset(-rule * 2);
  pen.rect(edge.x, edge.y, edge.w, rule, ink.line);
  pen.rect(edge.x, edge.y + edge.h - rule, edge.w, rule, ink.line);
  pen.rect(edge.x, edge.y, rule, edge.h, ink.line);
  pen.rect(edge.x + edge.w - rule, edge.y, rule, edge.h, ink.line);
  const grid = new Grid(frame, SIDE, SIDE, 0);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      if (mask[row * SIDE + col]) grid.fill(pen, col, row, color);
    }
  }
}

export default function draw(pen: Pen, ink: Ink) {
  const rule = pen.width / 320;
  const lay = pen.frame(0.08).inset(rule * 2);
  const gap = rule * 4;
  const side = (lay.w - gap) / 2;
  plate(pen, box(lay.x, lay.y, side, side), gasket(), ink.green, rule, ink);
  plate(pen, box(lay.x + side + gap, lay.y + side + gap, side, side), scatter(), ink.pink, rule, ink);
}
