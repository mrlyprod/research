import type { Ink, Pen, Point } from "@mrlyprod/draw";
import { plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const CODE = 126;
const LEVEL = 7;
const HEIGHTS = [190, 191];
const HALF = 6;

let memo: Point[][] | undefined;

// FACTS

function place(point: ArrayLike<number>): Point {
  const [u, v] = math.three.shadow(point);
  return [Number(u) / Math.SQRT2, Number(v) / Math.sqrt(6)];
}

function pieces(points: Uint32Array[]) {
  const tally = Array<number>(8).fill(0);
  for (const p of points) {
    const key = (p[0] >> HALF) * 4 + (p[1] >> HALF) * 2 + (p[2] >> HALF);
    if (key >= tally.length) throw new Error(`research-cuts: piece ${key} of point ${p}, want under ${tally.length}`);
    tally[key]++;
  }
  return tally.filter((n) => n > 0);
}

function cuts() {
  if (memo) return memo;
  const profile = math.three.profile(CODE, 2, LEVEL, 2);
  memo = HEIGHTS.map((height) => {
    if (profile[height] !== "2187")
      throw new Error(`research-cuts: profile ${profile[height]} at ${height}, want 2187`);
    const points = math.three.diagonal_slice(CODE, 2, LEVEL, 2, height);
    if (points.length !== 2187) throw new Error(`research-cuts: ${points.length} points at ${height}, want 2187`);
    const tally = pieces(points);
    if (tally.join() !== "729,729,729")
      throw new Error(`research-cuts: pieces ${tally} at ${height}, want 729,729,729`);
    return points.map(place);
  });
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const slices = cuts();
  const flat = slices.flat();
  if (flat.length !== 4374) throw new Error(`research-cuts: ${flat.length} points, want 4374`);
  let lo: Point = [Infinity, Infinity];
  let hi: Point = [-Infinity, -Infinity];
  for (const [u, v] of flat) {
    lo = [Math.min(lo[0], u), Math.min(lo[1], v)];
    hi = [Math.max(hi[0], u), Math.max(hi[1], v)];
  }
  const frame = pen.frame(0.08);
  const span = [hi[0] - lo[0], hi[1] - lo[1]];
  const scale = Math.min(frame.w / span[0], frame.h / span[1]);
  const [cx, cy] = frame.center();
  const screen = ([u, v]: Point): Point => [
    cx + (u - (lo[0] + hi[0]) / 2) * scale,
    cy - (v - (lo[1] + hi[1]) / 2) * scale,
  ];
  for (const [i, cut] of slices.entries()) plot.dots(pen, cut.map(screen), 2.2, [ink.blue, ink.orange][i]);
}
