import { frame as boxed, type Color, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const M = 5;
const K = 7;
const PX = 6;
const SIDE = 2 ** K;

// GASKETS

type Cells = ReturnType<typeof math.two.create>;

let memo: { two: Cells; three: Cells } | undefined;

function gaskets() {
  if (memo) return memo;
  const two = math.two.create(7, 2, K, 0, 2);
  const three = math.two.create(11, 3, M, 0, 3);
  if (two.shape[0] !== SIDE || two.shape[1] !== SIDE)
    throw new Error(`research-cobham: the base 2 gasket is ${two.shape}, want ${SIDE},${SIDE}`);
  if (three.shape[0] < SIDE || three.shape[1] < SIDE)
    throw new Error(`research-cobham: the base 3 gasket is ${three.shape}, want at least ${SIDE},${SIDE}`);
  memo = { two, three };
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { two, three } = gaskets();
  if (!(SIDE > (3 ** M - 1) / 2)) throw new Error(`research-cobham: side ${SIDE} does not pass ${(3 ** M - 1) / 2}`);
  const edge = SIDE * PX;
  const area = boxed(Math.round((pen.width - edge) / 2), Math.round((pen.height - edge) / 2), edge, edge);
  const grid = new Grid(area, SIDE, SIDE, 0);
  const wide = three.shape[1];
  let b = 0;
  let both = 0;
  let axes = 0;
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      const y = SIDE - 1 - row;
      const a = two.types[y * SIDE + col] !== 0;
      const inB = three.types[y * wide + col] !== 0;
      if (inB) b++;
      if (a && inB) both++;
      if (a && inB && (col === 0 || y === 0)) axes++;
      let tone: Color;
      if (a && inB) tone = ink.yellow;
      else if (inB) tone = ink.blue;
      else if (a) tone = ink.fade(ink.dim, 0.45);
      else continue;
      grid.fill(pen, col, row, tone);
    }
  }
  if (b !== 3 ** M) throw new Error(`research-cobham: ${b} base 3 cells, want ${3 ** M}`);
  if (both !== 111) throw new Error(`research-cobham: ${both} shared cells, want 111`);
  if (axes !== 2 ** (M + 1) - 1) throw new Error(`research-cobham: ${axes} axis cells, want ${2 ** (M + 1) - 1}`);
}
