import { frame as boxed, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const CODES = [1, 3, 7, 9, 11, 15];
const SIDES = [3, 5, 7, 9];

// DESIGNS

type Cells = ReturnType<typeof math.two.create>;

let memo: Cells[][] | undefined;

function designs() {
  if (memo) return memo;
  memo = CODES.map((code) =>
    SIDES.map((number) => {
      const cells = math.two.create(code, number, 1, 0, 2);
      if (cells.shape[1] !== number)
        throw new Error(`research-sequences: code ${code} is ${cells.shape[1]} wide, want ${number}`);
      let sum = 0n;
      for (const kind of cells.types) sum += BigInt(kind);
      const fill = BigInt(math.counts.fill(code, number, 2, 1, 2));
      if (sum !== fill)
        throw new Error(`research-sequences: code ${code} at ${number} holds ${sum} cells, want ${fill}`);
      return cells;
    }),
  );
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const cells = designs();
  const area = pen.frame(0.08);
  const pitch = area.h / CODES.length;
  const block = boxed(area.x + (area.w - pitch * SIDES.length) / 2, area.y, pitch * SIDES.length, area.h);
  const rule = 1.5;
  for (let col = 0; col <= SIDES.length; col++) {
    const x = block.x + col * pitch;
    pen.rect(x - rule / 2, block.y, rule, block.h, ink.line);
  }
  for (let row = 0; row <= CODES.length; row++) {
    const y = block.y + row * pitch;
    pen.rect(block.x, y - rule / 2, block.w, rule, ink.line);
  }
  CODES.forEach((_, row) => {
    SIDES.forEach((number, col) => {
      const tile = boxed(block.x + col * pitch, block.y + row * pitch, pitch, pitch).inset(pitch * 0.1);
      new Grid(tile, number, number, 0.1).paint(pen, cells[row][col], (kind) => (kind !== 0 ? ink.yellow : null));
    });
  });
}
