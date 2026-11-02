import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./research-information.json" with { type: "json" };

export const units = { math };

const SPLIT = 9;
const SIDE = SPLIT * SPLIT;

function whole() {
  const cells = math.two.create(495, 3, 4, 0, 3);
  if (cells.shape[0] !== SIDE || cells.shape[1] !== SIDE)
    throw new Error(`research-information: level 4 is ${cells.shape}, want ${SIDE} square`);
  const mask: boolean[][] = [];
  let ones = 0;
  for (let row = 0; row < SIDE; row++) {
    mask.push([]);
    for (let col = 0; col < SIDE; col++) {
      const on = cells.types[row * SIDE + col] !== 0;
      mask[row].push(on);
      if (on) ones++;
    }
  }
  if (ones !== 4096) throw new Error(`research-information: ${ones} cells at level 4, want 4096`);
  return mask;
}

function rearranged() {
  const rows = census.rearranged;
  if (rows.length !== SIDE || rows.some((row) => row.length !== SIDE))
    throw new Error(`research-information: the rearrangement is not ${SIDE} square`);
  return rows.map((row) => [...row].map((c) => c === "1"));
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.07);
  const gutter = area.w * 0.06;
  const cell = Math.floor((area.w - gutter) / 2 / SIDE);
  const side = cell * SIDE;
  const left = area.x + (area.w - 2 * side - gutter) / 2;
  const top = area.y + (area.h - side) / 2;
  new Grid(frame(left, top, side, side), SIDE, SIDE, 0).carpet(pen, whole(), ink.blue);
  new Grid(frame(left + side + gutter, top, side, side), SIDE, SIDE, 0).carpet(pen, rearranged(), ink.orange);
}
