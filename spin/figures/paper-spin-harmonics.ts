import { Grid, type Ink, type Pen, type Point } from "@mrlyprod/draw";
import census from "./paper-spin-harmonics.json" with { type: "json" };

const LONG = 4;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { side, cells, chords } = census;
  const inside = (n: number) => Number.isInteger(n) && n >= 0 && n < side;
  if (cells.some((cell) => cell.length !== 2 || !cell.every(inside)))
    throw new Error(`paper-spin-harmonics: a cell off the ${side} by ${side} lattice`);
  if (
    chords.length !== (cells.length * (cells.length - 1)) / 2 ||
    chords.some((c) => c.length !== 5 || !c.slice(0, 4).every(inside))
  ) {
    throw new Error(`paper-spin-harmonics: ${chords.length} chords, want one for each pair of ${cells.length} cells`);
  }
  const frame = pen.frame(0.08);
  const step = frame.w / side;
  const lattice = new Grid(frame, side, side, 0.1);
  const centre = (row: number, col: number): Point => [frame.x + (col + 0.5) * step, frame.y + (row + 0.5) * step];

  for (let row = 0; row < side; row++) {
    for (let col = 0; col < side; col++) lattice.fill(pen, col, row, ink.fade(ink.dim, 0.16));
  }
  for (const [row, col] of cells) lattice.fill(pen, col, row, ink.blue);

  const thin = step * 0.017;
  const thick = step * 0.038;
  const weight = (d2: number) => (d2 === LONG ? thick : thin);
  const casing = (d2: number) => weight(d2) + step * (d2 === LONG ? 0.016 : 0.008);
  for (const [r1, c1, r2, c2, d2] of chords) pen.segment(centre(r1, c1), centre(r2, c2), casing(d2), ink.ground);
  for (const [r1, c1, r2, c2, d2] of chords) {
    if (d2 !== LONG) pen.segment(centre(r1, c1), centre(r2, c2), thin, ink.dim);
  }
  for (const [r1, c1, r2, c2, d2] of chords) {
    if (d2 === LONG) pen.segment(centre(r1, c1), centre(r2, c2), thick, ink.pink);
  }
}
