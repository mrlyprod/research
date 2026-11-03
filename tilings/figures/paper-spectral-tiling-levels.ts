import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./paper-spectral-tiling-levels.json" with { type: "json" };

type Column = [number[], number[]];
type Census = {
  base: number;
  digits: number[];
  depth: number;
  low: number;
  high: number;
  points: number;
  columns: number;
  gap: number;
  thin: number;
  heights: number[];
  stacks: Column[][];
};

export const units = {};

const BASE = 12;
const DIGITS = [0, 3, 6, 9];
const DEPTH = 5;
const LOW = -1.25;
const HIGH = 2;
const GAP = 4;

function check(c: Census, cols: number, side: number) {
  if (c.base !== BASE || c.digits.join() !== DIGITS.join() || c.depth !== DEPTH || c.low !== LOW || c.high !== HIGH)
    throw new Error("paper-spectral-tiling-levels: the census draws another measure");
  if (c.points !== 4 ** DEPTH)
    throw new Error(`paper-spectral-tiling-levels: ${c.points} lattice points, want ${4 ** DEPTH}`);
  if (
    c.columns !== cols ||
    c.heights.length !== cols ||
    c.stacks.length !== 2 ||
    c.stacks.some((s) => s.length !== cols)
  )
    throw new Error(`paper-spectral-tiling-levels: the census holds ${c.columns} columns, the frame ${cols}`);
  if (c.gap !== GAP || Math.abs(c.thin * side - GAP) > 1e-9)
    throw new Error(`paper-spectral-tiling-levels: bands thinner than ${c.thin * side} px merged, want ${GAP}`);
  const sum = (v: number[]) => v.reduce((a, b) => a + b, 0);
  for (let k = 0; k < cols; k++) {
    const q = c.heights[k];
    if (!(q <= 1 + 1e-9 && q > 1 - 1e-2)) throw new Error(`paper-spectral-tiling-levels: Psi at column ${k} = ${q}`);
    let total = 0;
    for (const [bands, runs] of c.stacks.map((s) => s[k])) {
      if (runs.length !== bands.length + 1)
        throw new Error(
          `paper-spectral-tiling-levels: column ${k} holds ${bands.length} bands and ${runs.length} runs`,
        );
      if (bands.some((b) => b < c.thin) || runs.some((r) => r < 0))
        throw new Error(`paper-spectral-tiling-levels: column ${k} holds a thin band or a negative run`);
      total += sum(bands) + sum(runs);
    }
    if (Math.abs(total - q) > 1e-6)
      throw new Error(`paper-spectral-tiling-levels: column ${k} stacks to ${total}, Psi = ${q}`);
  }
}

export default function draw(pen: Pen, ink: Ink) {
  const data = census as Census;
  const frame = pen.frame(0.08);
  const cols = Math.round(frame.w);
  check(data, cols, frame.h);
  const tones = [ink.blue, ink.orange];
  for (let c = 0; c < cols; c++) {
    let y = frame.y + frame.h;
    data.stacks.forEach((stack, g) => {
      const [bands, runs] = stack[c];
      const run = (q: number) => {
        const h = q * frame.h;
        if (h > 0) pen.rect(frame.x + c, y - h, 1, h, tones[g]);
        y -= h;
      };
      bands.forEach((q, i) => {
        run(runs[i]);
        const h = q * frame.h;
        pen.rect(frame.x + c, y - h, 1, h - 1.5, tones[g]);
        pen.rect(frame.x + c, y - 1.5, 1, 1.5, ink.ground);
        y -= h;
      });
      run(runs[bands.length]);
    });
  }
}
