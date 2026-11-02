import { frame as box, type Color, type Frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const SIDE = 5;
const GAP = 0.14;
const CODES = [1, 3, 7, 9, 11, 15];
const COUNTS = [9, 15, 21, 13, 19, 25];

function outline(pen: Pen, frame: Frame, thick: number, ink: Ink) {
  pen.polyline(
    [
      [frame.x, frame.y],
      [frame.x + frame.w, frame.y],
      [frame.x + frame.w, frame.y + frame.h],
      [frame.x, frame.y + frame.h],
      [frame.x, frame.y],
    ],
    thick,
    ink.line,
  );
}

function panel(pen: Pen, frame: Frame, code: number, tone: Color, ink: Ink) {
  const design = math.two.create(code, SIDE, 1, 0, 2);
  if (design.shape.join() !== `${SIDE},${SIDE}`)
    throw new Error(`paper-sequence-census: design ${design.shape}, want ${SIDE} square`);
  const grid = new Grid(frame, SIDE, SIDE, GAP);
  let filled = 0;
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      if (design.types[row * SIDE + col] === 0) continue;
      filled++;
      const [x, y, w, h] = grid.cell(col, row);
      pen.rect(x, y, w, h, tone);
    }
  }
  const pitch = frame.w / SIDE;
  outline(pen, box(frame.x + pitch, frame.y + pitch, 3 * pitch, 3 * pitch), 2, ink);
  return filled;
}

export default function draw(pen: Pen, ink: Ink) {
  const tones = [ink.blue, ink.blue, ink.blue, ink.indigo, ink.indigo, ink.indigo];
  const frame = pen.frame(0.08);
  const height = frame.h * 0.8;
  const block = box(frame.x, frame.y + (frame.h - height) / 2, frame.w, height);
  const counts: number[] = [];
  block.rows(2).forEach((row, band) => {
    row.cols(3).forEach((slot, index) => {
      const place = band * 3 + index;
      counts.push(panel(pen, slot.inset(16), CODES[place], tones[place], ink));
    });
  });
  if (counts.join() !== COUNTS.join()) throw new Error(`paper-sequence-census: fills ${counts}, want ${COUNTS}`);
}
