import { frame as box, type Color, type Frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

function letter(code: number) {
  return math.bang.MagicLayer.new(new math.name.Bang(code, 2, 2), 2);
}

function lattice(pen: Pen, frame: Frame, side: number, thick: number, ink: Ink) {
  const step = frame.w / side;
  for (let k = 0; k <= side; k++) {
    const offset = k * step;
    pen.segment([frame.x + offset, frame.y], [frame.x + offset, frame.y + frame.h], thick, ink.line);
    pen.segment([frame.x, frame.y + offset], [frame.x + frame.w, frame.y + offset], thick, ink.line);
  }
}

function panel(pen: Pen, frame: Frame, order: [number, number], tone: Color, ink: Ink) {
  const picture = math.bang.magic([letter(order[0]), letter(order[1])]);
  if (picture.shape.join() !== "4,4")
    throw new Error(`paper-order-sensitivity-of-kronecker-words: picture ${picture.shape}, want 4 square`);
  lattice(pen, frame, 4, 2.5, ink);
  const grid = new Grid(frame, 4, 4, 0);
  let filled = 0;
  for (let row = 0; row < 4; row++) {
    for (let col = 0; col < 4; col++) {
      if (picture.data[row * 4 + col] === 0) continue;
      filled++;
      const [x, y, w, h] = grid.cell(col, row);
      pen.rect(x, y, w, h, tone);
    }
  }
  return filled;
}

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const half = frame.w / 2;
  const top = box(frame.x, frame.y, half, half);
  const low = box(frame.x + half, frame.y + half, half, half);
  const first = panel(pen, top, [3, 6], ink.blue, ink);
  const second = panel(pen, low, [6, 3], ink.orange, ink);
  if (first !== 4 || second !== 4)
    throw new Error(`paper-order-sensitivity-of-kronecker-words: panels fill ${first} and ${second}, want 4 and 4`);
  const forward = math.bang.word.components([letter(3), letter(6)]);
  const backward = math.bang.word.components([letter(6), letter(3)]);
  if (forward !== "4" || backward !== "2")
    throw new Error(`paper-order-sensitivity-of-kronecker-words: components ${forward} and ${backward}, want 4 and 2`);
}
