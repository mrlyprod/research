import { type Color, type Frame, frame, Grid, type Ink, type Pen, type Tensor } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const SIDE = 15;

function letters(first: [number, number], second: [number, number]) {
  return [first, second].map(([code, number]) => math.bang.MagicLayer.new(new math.name.Bang(code, 2, 2), number));
}

function panel(pen: Pen, box: Frame, cells: Tensor, color: Color) {
  new Grid(box, SIDE, SIDE, 0.09).paint(pen, cells, (kind) => (kind !== 0 ? color : null));
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const gutter = area.w * 0.055;
  const side = (area.w - gutter) / 2;
  const top = area.y + (area.h - side) / 2;

  const forward = letters([7, 3], [14, 5]);
  const reverse = letters([14, 5], [7, 3]);
  if (math.bang.word.side(forward) !== String(SIDE))
    throw new Error(`research-magic: forward side ${math.bang.word.side(forward)}, want ${SIDE}`);
  if (math.bang.word.side(reverse) !== String(SIDE))
    throw new Error(`research-magic: reverse side ${math.bang.word.side(reverse)}, want ${SIDE}`);
  if (math.bang.word.fill(forward) !== math.bang.word.fill(reverse))
    throw new Error("research-magic: the two orders fill differently");

  const left = math.bang.magic(forward);
  const right = math.bang.magic(reverse);
  for (const tensor of [left, right]) {
    if (
      tensor.shape.length !== 2 ||
      tensor.shape[0] !== SIDE ||
      tensor.shape[1] !== SIDE ||
      tensor.data.length !== SIDE * SIDE
    )
      throw new Error(`research-magic: a word is ${tensor.shape}, want ${SIDE} square`);
  }
  const sums = [0, 0];
  let differs = false;
  for (let i = 0; i < SIDE * SIDE; i++) {
    sums[0] += left.data[i];
    sums[1] += right.data[i];
    if (left.data[i] !== right.data[i]) differs = true;
  }
  if (sums[0] !== sums[1]) throw new Error(`research-magic: fills ${sums[0]} and ${sums[1]} differ`);
  if (!differs) throw new Error("research-magic: the two orders are the same word");

  panel(pen, frame(area.x, top, side, side), left, ink.blue);
  panel(pen, frame(area.x + side + gutter, top, side, side), right, ink.orange);
}
