import type { Ink, Pen, Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./paper-erdos-125-upper-density.json" with { type: "json" };

const LEVEL = 6;
const SAMPLES = 3072;

export const units = { math };

const { curves, orbit: chain } = census as { curves: number[][]; orbit: [number, number][] };

let memo: [number, number][] | undefined;

function orbit() {
  if (memo) return memo;
  memo = chain.map(([four, h], three) => {
    if (three === 0) return [1, h];
    const pair = new math.sumset.Pair(three, four);
    const tau = pair.scale();
    pair.free();
    if (!(tau >= 1 && tau < 4))
      throw new Error(`paper-erdos-125-upper-density: scale ${tau} of pair (${three}, ${four}) off [1, 4)`);
    return [tau, h];
  });
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  if (curves.length !== LEVEL + 1 || curves.some((hs) => hs.length !== SAMPLES + 1))
    throw new Error(
      `paper-erdos-125-upper-density: ${curves.length} curves of ${curves.map((hs) => hs.length)}, want ${LEVEL + 1} of ${SAMPLES + 1}`,
    );
  if (chain.length !== LEVEL + 1)
    throw new Error(`paper-erdos-125-upper-density: ${chain.length} orbit points, want ${LEVEL + 1}`);
  const frame = pen.frame(0.08);
  const lo = 0.5 * 0.92;
  const hi = curves[LEVEL].reduce((a, h) => Math.max(a, h), 0) * 1.03;
  const at = (u: number, h: number): Point => [frame.x + frame.w * u, frame.y + (frame.h * (hi - h)) / (hi - lo)];
  const marks = orbit().map(([tau, h]) => at(Math.log(tau) / Math.log(4), h));
  const foot = frame.y + frame.h;
  for (const [x, y] of marks) pen.segment([x, foot], [x, y], 1.4, ink.line);
  pen.segment([frame.x, foot], [frame.x + frame.w, foot], 1.4, ink.line);
  curves.forEach((hs, k) => {
    const t = k / LEVEL;
    pen.polyline(
      hs.map((h, i) => at(i / SAMPLES, h)),
      1.3 + 1.3 * t,
      ink.mix(ink.dim, ink.blue, t * t),
    );
  });
  for (const [x, y] of marks) {
    pen.disc(x, y, 9, ink.ground);
    pen.disc(x, y, 6.5, ink.orange);
  }
}
