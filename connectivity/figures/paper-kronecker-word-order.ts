import { frame as box, type Color, type Frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LENGTH = 6;
const SIDE = 64;
const CELLS = 216;
const GASKET = 7;
const DOMINO = 3;
const PIECES = [14, 45];

// WORDS

function thueMorse(index: number) {
  let bits = 0;
  for (let n = index; n > 0; n >>= 1) bits += n & 1;
  return bits & 1;
}

function letter(code: number) {
  return math.bang.MagicLayer.new(new math.name.Bang(code, 2, 2), 2);
}

function reading(gasketAt: number) {
  return Array.from({ length: LENGTH }, (_, place) => letter(thueMorse(place) === gasketAt ? GASKET : DOMINO));
}

// PIECES

function pieces(cells: ArrayLike<number>) {
  const label = new Int32Array(SIDE * SIDE).fill(-1);
  let count = 0;
  let filled = 0;
  for (let start = 0; start < SIDE * SIDE; start++) {
    if (!cells[start] || label[start] >= 0) continue;
    const stack = [start];
    label[start] = count;
    while (stack.length) {
      const at = stack.pop() as number;
      filled++;
      const row = Math.floor(at / SIDE);
      const col = at % SIDE;
      for (const [dr, dc] of [
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
      ]) {
        const r = row + dr;
        const c = col + dc;
        if (r < 0 || r >= SIDE || c < 0 || c >= SIDE) continue;
        const next = r * SIDE + c;
        if (cells[next] && label[next] < 0) {
          label[next] = count;
          stack.push(next);
        }
      }
    }
    count++;
  }
  return { label, count, filled };
}

// PANEL

function panel(pen: Pen, frame: Frame, letters: ReturnType<typeof letter>[], inks: Color[], rule: number, ink: Ink) {
  const picture = math.bang.magic(letters);
  if (picture.shape.join() !== `${SIDE},${SIDE}`)
    throw new Error(`paper-kronecker-word-order: picture ${picture.shape}, want ${SIDE} square`);
  const { label, count, filled } = pieces(picture.data);
  if (filled !== CELLS) throw new Error(`paper-kronecker-word-order: ${filled} cells, want ${CELLS}`);
  const claimed = String(math.bang.word.components(letters));
  if (claimed !== String(count))
    throw new Error(`paper-kronecker-word-order: crate counts ${claimed} pieces, the drawing ${count}`);
  const edge = frame.inset(-rule * 2);
  pen.rect(edge.x, edge.y, edge.w, rule, ink.line);
  pen.rect(edge.x, edge.y + edge.h - rule, edge.w, rule, ink.line);
  pen.rect(edge.x, edge.y, rule, edge.h, ink.line);
  pen.rect(edge.x + edge.w - rule, edge.y, rule, edge.h, ink.line);
  const grid = new Grid(frame, SIDE, SIDE, 0);
  for (let at = 0; at < SIDE * SIDE; at++) {
    if (label[at] >= 0) grid.fill(pen, at % SIDE, Math.floor(at / SIDE), inks[label[at] % inks.length]);
  }
  return count;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const rule = pen.width / 320;
  const lay = pen.frame(0.06).inset(rule * 2);
  const gap = rule * 4;
  const side = (lay.w - gap) / 2;
  const inks = [ink.blue, ink.orange];
  const first = panel(pen, box(lay.x, lay.y, side, side), reading(0), inks, rule, ink);
  const second = panel(pen, box(lay.x + side + gap, lay.y + side + gap, side, side), reading(1), inks, rule, ink);
  if (first !== PIECES[0] || second !== PIECES[1])
    throw new Error(`paper-kronecker-word-order: ${first} and ${second} pieces, want ${PIECES.join(" and ")}`);
}
