import { type Frame, field, frame, Grid, type Ink, type Pen, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LEVEL = 4;
const SIDE = 81;
const EPS = 1 / 36;
const PLATE = 600;
const BANDS = [
  [1 / 12, 2.122718, 2.122723],
  [1 / 8, 2.134668, 2.135742],
  [1 / 6, 2.135019, 2.136794],
];
const TUBES = [
  [1 / 12, 0.180947086, 0.180947093],
  [1 / 8, 0.234186414, 0.234701259],
];

const sponge = math.three.sponge;
const point = new Float64Array([0, 0, 0.5]);

function dist(x: number, y: number) {
  point[0] = x;
  point[1] = y;
  return sponge.distance(point);
}

let memo: { cells: ArrayLike<number>; covers: Float64Array } | undefined;

function facts() {
  if (memo) return memo;
  const cube = math.three.carpet(3, LEVEL);
  let filled = 0;
  for (const kind of cube.types) filled += kind;
  if (filled !== 160000) throw new Error(`paper-sponge-measurability: carpet sums to ${filled}, want 160000`);
  const dust = math.three.slice(cube, 2, (SIDE - 1) / 2);
  const cells = dust.types;
  if (cells.length !== SIDE * SIDE)
    throw new Error(`paper-sponge-measurability: slice of ${cells.length} cells, want ${SIDE * SIDE}`);
  let lit = 0;
  for (const kind of cells) if (kind !== 0) lit++;
  if (lit !== 256) throw new Error(`paper-sponge-measurability: slice lights ${lit} cells, want 256`);

  const ident = (Math.PI + 8) / 36 - Math.sqrt(2) / 27;
  if (!(Math.abs((sponge.tube(1 / 6) ?? 0) - ident) < 1e-9))
    throw new Error("paper-sponge-measurability: tube(1/6) is off the closed form");
  for (const [delta, lo, hi] of TUBES) {
    const t = sponge.tube(delta) ?? 0;
    if (!(lo <= t && t <= hi))
      throw new Error(`paper-sponge-measurability: tube(${delta}) is ${t}, want ${lo} to ${hi}`);
  }
  const values = [0, 0, 0];
  BANDS.forEach(([eps, lo, hi], k) => {
    values[k] = sponge.profile(eps) ?? 0;
    if (!(lo <= values[k] && values[k] <= hi))
      throw new Error(`paper-sponge-measurability: profile(${eps}) is ${values[k]}, want ${lo} to ${hi}`);
  });
  if (!(values[2] - values[0] >= 0.012296))
    throw new Error("paper-sponge-measurability: the profile gap is under 0.012296");
  if (!(Math.abs(dist(0.5, 0.5) - sponge.COVER()) < 1e-15))
    throw new Error("paper-sponge-measurability: the centre is off the covering radius");
  if (!(Math.abs(dist(1 / 6, 1 / 6) - Math.sqrt(2) / 18) < 1e-15))
    throw new Error("paper-sponge-measurability: the arm is off sqrt(2)/18");

  const covers = new Float64Array(PLATE * PLATE);
  for (let j = 0; j < PLATE; j++) {
    for (let i = 0; i < PLATE; i++) {
      const cover = 0.5 + (EPS - dist((i + 0.5) / PLATE, (j + 0.5) / PLATE)) * PLATE;
      covers[j * PLATE + i] = Math.min(Math.max(cover, 0), 1);
    }
  }
  memo = { cells, covers };
  return memo;
}

function bands(pen: Pen, ink: Ink, box: Frame) {
  plot.axis(pen, box, ink.line);
  const area = box.inset(24);
  const [lo, hi] = [2.12, 2.14];
  const yOf = (v: number) => area.y + (area.h * (hi - v)) / (hi - lo);
  const width = area.w * 0.3;
  const [a, b] = [BANDS[0], BANDS[2]];
  [a, b].forEach((band, k) => {
    const x = area.x + area.w * (0.25 + 0.5 * k) - width / 2;
    const [top, foot] = [yOf(band[2]), yOf(band[1])];
    pen.rect(x, top, width, Math.max(foot - top, 4), ink.orange);
  });
  const [cx, gapLo, gapHi] = [area.x + area.w / 2, yOf(a[2]) - 2, yOf(b[1]) + 2];
  pen.segment([cx, gapLo], [cx, gapHi], 2, ink.dim);
  pen.segment([cx - 10, gapLo], [cx + 10, gapLo], 2, ink.dim);
  pen.segment([cx - 10, gapHi], [cx + 10, gapHi], 2, ink.dim);
}

export default function draw(pen: Pen, ink: Ink) {
  const { cells, covers } = facts();
  const margin = Math.round(pen.width * 0.08);
  const sheet = frame(margin, margin, PLATE, PLATE);
  pen.rect(sheet.x, sheet.y, sheet.w, sheet.h, ink.panel);
  const tube = field.patch(sheet.x, sheet.y, PLATE, PLATE);
  for (let j = 0; j < PLATE; j++) {
    for (let i = 0; i < PLATE; i++) {
      const cover = covers[j * PLATE + i];
      if (cover > 0) tube.blend(sheet.x + i, sheet.y + j, ink.blue, cover);
    }
  }
  tube.paint(pen);
  const lattice = new Grid(sheet, SIDE, SIDE, 0);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) if (cells[row * SIDE + col] !== 0) lattice.fill(pen, col, row, ink.fg);
  }
  bands(pen, ink, frame(pen.width - margin - 240, pen.height - margin - 240, 240, 240));
}
