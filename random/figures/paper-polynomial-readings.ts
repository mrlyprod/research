import type { Ink, Pen, Point } from "@mrlyprod/draw";

const NAME = "paper-polynomial-readings";
const DEPTH = 12;
const LAGS = 1 << DEPTH;
const TOTAL = 2 / 3;
const SLACK = 1e-9;

let memo: { walks: [number, number][][]; lo: number; hi: number } | undefined;

function correlations() {
  const c = new Float64Array(LAGS + 1);
  c[0] = 1;
  c[1] = -1 / 3;
  for (let h = 2; h <= LAGS; h++) c[h] = h & 1 ? -(c[h >> 1] + c[(h >> 1) + 1]) / 2 : c[h >> 1];
  return c;
}

function walks() {
  if (memo) return memo;
  const c = correlations();
  const sums = [0, 0, 0];
  const out: [number, number][][] = [[[0, 0]], [[0, 0]], [[0, 0]]];
  let lo = 0;
  let hi = 0;
  for (let h = 0; h < LAGS; h++) {
    sums[h % 3] += c[h];
    const H = h + 1;
    for (let r = 0; r < 3; r++) {
      out[r].push([H, sums[r]]);
      lo = Math.min(lo, sums[r]);
      hi = Math.max(hi, sums[r]);
    }
    if (H >= 2 && (H & h) === 0) {
      const k = Math.log2(H);
      const closed = (16 / 45) * 1.5 ** k + (k % 2 ? 7 / 15 : 8 / 15);
      if (Math.abs(sums[0] - closed) > SLACK)
        throw new Error(`${NAME}: U_0(2^${k}) = ${sums[0]}, the closed form gives ${closed}`);
      if (Math.abs(sums[0] + sums[1] + sums[2] - TOTAL) > SLACK)
        throw new Error(`${NAME}: the correlations below 2^${k} sum to ${sums[0] + sums[1] + sums[2]}, want 2/3`);
    }
  }
  memo = { walks: out, lo, hi };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { walks: lines, lo, hi } = walks();
  const area = pen.frame(0.08);
  const at = ([H, v]: [number, number]): Point => [
    area.x + (H / LAGS) * area.w,
    area.y + area.h - ((v - lo) / (hi - lo)) * area.h,
  ];
  const thick = Math.max(1.5, area.w * 0.005);
  const zero = at([0, 0])[1];
  pen.segment([area.x, zero], [area.x + area.w, zero], thick * 0.5, ink.line);
  pen.polyline(lines[2].map(at), thick, ink.fade(ink.orange, 0.5));
  pen.polyline(lines[1].map(at), thick, ink.orange);
  pen.polyline(lines[0].map(at), thick, ink.blue);
}
