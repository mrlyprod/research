import type { Ink, Pen } from "@mrlyprod/draw";
import { iso } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const ROWS = 9;
const WIDEST = 6.5;
const COUNTS = [1, 1, 3, 3, 6, 3, 3, 1, 1];

type Quad = ReturnType<typeof math.three.quads>[number];

let memo: Quad[][][] | undefined;

function rows() {
  if (memo) return memo;
  const codes: string[][] = Array.from({ length: ROWS }, () => []);
  for (const design of math.bang.bang(3).canonical()) codes[design.rule().length].push(design.i);
  if (codes.some((row, r) => row.length !== COUNTS[r]))
    throw new Error(`research-method: rows of ${codes.map((row) => row.length)}, want ${COUNTS}`);
  memo = codes.map((row, r) =>
    row.map((code) => {
      const cube = math.three.create(code, 2, 1, 2);
      if (math.three.fills(cube) !== r)
        throw new Error(`research-method: design ${code} fills ${math.three.fills(cube)}, want ${r}`);
      return math.three.quads(cube);
    }),
  );
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const shade = [ink.blue, ink.mix(ink.blue, ink.ground, 0.4), ink.mix(ink.blue, ink.ground, 0.65)];
  const wire = ink.mix(ink.line, ink.dim, 0.3);
  const [mx, my] = frame.center();
  const pitch = frame.h / ROWS;
  const step = frame.w / WIDEST;
  const s = pitch * 0.225;
  rows().forEach((row, r) => {
    const cy = my + (r - (ROWS - 1) / 2) * pitch;
    row.forEach((quads, j) => {
      const cx = mx + (j - (row.length - 1) / 2) * step;
      iso.cage(pen, cx, cy, s, 1, s / 20, wire);
      iso.stamp(pen, quads, cx, cy, s, shade, ink.ground, s / 14);
    });
  });
}
