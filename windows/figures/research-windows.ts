import type { Color, Ink, Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./research-windows.json" with { type: "json" };

export const units = { math };

const LEVEL = 3;
const SIDE = 27;
const RING = 11;
const LAST = RING - 1;
const SLOTS = 4 * RING - 4;
const CORNERS: [number, number][] = [
  [0, 0],
  [0, LAST],
  [LAST, LAST],
  [LAST, 0],
];
const CELL = 22;
const HAIR = 2;
const GUTTER = 15;

let memo: ArrayLike<number> | undefined;

function carpet() {
  if (memo) return memo;
  const design = math.bang.factory.create("495", 3, 2, 3, LEVEL);
  if (design.shape.length !== 2 || design.shape[0] !== SIDE || design.shape[1] !== SIDE)
    throw new Error(`research-windows: the carpet is ${design.shape}, want ${SIDE} by ${SIDE}`);
  let filled = 0;
  for (const kind of design.data) filled += kind;
  if (filled !== 8 ** LEVEL) throw new Error(`research-windows: the carpet holds ${filled} cells, want ${8 ** LEVEL}`);
  memo = design.data;
  return memo;
}

function seat(slot: number): [number, number] {
  if (slot < 4) return CORNERS[slot];
  const n = Math.floor((slot - 4) / 4);
  const sides: [number, number][] = [
    [0, n + 1],
    [n + 1, LAST],
    [LAST, LAST - 1 - n],
    [LAST - 1 - n, 0],
  ];
  return sides[(slot - 4) % 4];
}

function cells(pen: Pen, x: number, y: number, side: number, lit: ArrayLike<number>, on: Color, off: Color | null) {
  const size = CELL - HAIR;
  for (let i = 0; i < lit.length; i++) {
    const color = lit[i] === 1 ? on : off;
    if (color) pen.rect(x + (i % side) * CELL, y + Math.floor(i / side) * CELL, size, size, color);
  }
}

export default function draw(pen: Pen, ink: Ink) {
  const { windows } = census;
  if (windows.length !== SLOTS) throw new Error(`research-windows: ${windows.length} windows, want ${SLOTS}`);
  const tile = 3 * CELL - HAIR;
  const span = RING * tile + (RING - 1) * GUTTER;
  const origin = Math.floor((pen.width - span) / 2);
  const inner = Math.floor((pen.width - (SIDE * CELL - HAIR)) / 2);
  cells(pen, inner, inner, SIDE, carpet(), ink.fade(ink.dim, 0.75), null);
  windows.forEach((mask, slot) => {
    if (!(Number.isInteger(mask) && mask >= 0 && mask < 512))
      throw new Error(`research-windows: window ${mask} is not nine bits`);
    const [r, c] = seat(slot);
    const bits = Array.from({ length: 9 }, (_, i) => (mask >> i) & 1);
    cells(pen, origin + c * (tile + GUTTER), origin + r * (tile + GUTTER), 3, bits, ink.blue, ink.line);
  });
}
