import { Grid, type Ink, type Pen, type Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LEVEL = 3;
const SIDE = 27;
const LOOPS = 50;
const STRANDS = 2 * SIDE;
const HALF = Math.PI / 2;

// ARCS

let memo: number[] | undefined;

function arcs() {
  if (memo) return memo;
  const design = math.bang.factory.create(495, 3, 2, 3, LEVEL);
  if (design.shape.length !== 2 || design.shape[0] !== SIDE || design.shape[1] !== SIDE)
    throw new Error(`research-arcs: the design is ${design.shape}, want ${SIDE},${SIDE}`);
  let sum = 0;
  for (const kind of design.data) sum += kind;
  if (sum !== 8 ** LEVEL) throw new Error(`research-arcs: the design fills ${sum} cells, want ${8 ** LEVEL}`);
  const parity = math.bang.factory.create(7, 3, 2, 2, LEVEL);
  if (parity.data.length !== design.data.length || parity.data.some((kind, i) => kind !== design.data[i]))
    throw new Error("research-arcs: code 7 at base 2 differs from the carpet");

  const drawn = math.arcs.draw(495, 3, LEVEL);
  if (drawn.side !== SIDE) throw new Error(`research-arcs: the arcs are ${drawn.side} wide, want ${SIDE}`);
  if (drawn.cells.length !== SIDE * SIDE)
    throw new Error(`research-arcs: ${drawn.cells.length} arc cells, want ${SIDE * SIDE}`);
  for (let f = 0; f < SIDE * SIDE; f++) {
    if ((drawn.cells[f] & 1) !== (design.data[f] === 1 ? 1 : 0))
      throw new Error(`research-arcs: arc cell ${f} disagrees with the design`);
  }
  if (drawn.strands !== STRANDS) throw new Error(`research-arcs: ${drawn.strands} strands, want ${STRANDS}`);
  if (drawn.loops !== LOOPS) throw new Error(`research-arcs: ${drawn.loops} loops, want ${LOOPS}`);
  const law = math.arcs.law(495, 3, LEVEL);
  if (law?.loops !== LOOPS) throw new Error(`research-arcs: the loop law gives ${law?.loops}, want ${LOOPS}`);
  memo = drawn.cells;
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const cells = arcs();
  const frame = pen.frame(0.08);
  const cell = frame.w / SIDE;
  const thick = cell * 0.15;
  const radius = cell / 2;
  const grid = new Grid(frame, SIDE, SIDE, 0);
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      if (cells[y * SIDE + x] & 1) grid.fill(pen, x, SIDE - 1 - y, ink.line);
    }
  }
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      const byte = cells[y * SIDE + x];
      const left = frame.x + x * cell;
      const top = frame.y + (SIDE - 1 - y) * cell;
      const sweeps: [Point, [number, number]][] =
        byte & 1
          ? [
              [
                [left, top + cell],
                [-HALF, 0],
              ],
              [
                [left + cell, top],
                [HALF, Math.PI],
              ],
            ]
          : [
              [
                [left + cell, top + cell],
                [Math.PI, Math.PI + HALF],
              ],
              [
                [left, top],
                [0, HALF],
              ],
            ];
      sweeps.forEach(([centre, angles], arc) => {
        pen.arc(centre, radius, angles, thick, byte & (2 << arc) ? ink.orange : ink.blue);
      });
    }
  }
}
