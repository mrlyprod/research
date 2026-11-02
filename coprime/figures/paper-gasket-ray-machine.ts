import type { Ink, Pen, Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";

export const units = { math, num };

const LEVEL = 5;

type Gasket = { points: Point[]; rays: Map<string, { ray: Point; on: Point[] }> };

let memo: Gasket | undefined;

function gasket() {
  if (memo) return memo;
  const cell = math.two.from_corners(
    [
      [0, 0],
      [1, 0],
      [0, 1],
    ],
    3,
    LEVEL,
    0,
    3,
  );
  const side = 3 ** LEVEL;
  if (cell.shape.length !== 2 || cell.shape[0] !== side || cell.shape[1] !== side)
    throw new Error(`paper-gasket-ray-machine: gasket ${cell.shape}, want ${side} square`);
  const points: Point[] = [];
  for (let row = 0; row < side; row++) {
    for (let col = 0; col < side; col++) {
      if (cell.types[row * side + col] !== 0) points.push([row, col]);
    }
  }
  if (points.length !== 3 ** LEVEL)
    throw new Error(`paper-gasket-ray-machine: ${points.length} points, want ${3 ** LEVEL}`);
  const rays = new Map<string, { ray: Point; on: Point[] }>();
  for (const point of points) {
    if (point[0] === 0 && point[1] === 0) continue;
    const step = Number(num.factor.gcd(point[0], point[1]));
    const ray: Point = [point[0] / step, point[1] / step];
    const key = `${ray[0]},${ray[1]}`;
    const held = rays.get(key);
    if (held) held.on.push(point);
    else rays.set(key, { ray, on: [point] });
  }
  const ladder = num.series.fibonacci(1000);
  const shift = rays.get("1,3")?.on.length;
  if (shift !== ladder[LEVEL] - 1)
    throw new Error(`paper-gasket-ray-machine: ray 1,3 holds ${shift}, want ${ladder[LEVEL] - 1}`);
  if (rays.get("3,1")?.on.length !== shift)
    throw new Error(`paper-gasket-ray-machine: ray 3,1 holds ${rays.get("3,1")?.on.length}, want ${shift}`);
  memo = { points, rays };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { points, rays } = gasket();
  const frame = pen.frame(0.08);
  let reach = 0;
  for (const point of points) reach = Math.max(reach, point[0], point[1]);
  const step = frame.w / reach;
  const at = (point: Point): Point => [frame.x + point[1] * step, frame.y + frame.h - point[0] * step];
  const origin = at([0, 0]);
  const order = [...rays.values()].sort(
    (a, b) => a.on.length - b.on.length || a.ray[0] - b.ray[0] || a.ray[1] - b.ray[1],
  );
  for (const { ray, on } of order) {
    let far = on[0];
    for (const point of on) if (point[0] + point[1] >= far[0] + far[1]) far = point;
    const mass = on.length;
    const gold = (ray[0] === 1 && ray[1] === 3) || (ray[0] === 3 && ray[1] === 1);
    const color = gold ? ink.yellow : ink.fade(ink.blue, Math.min(0.16 + 0.13 * (mass - 1), 1));
    pen.segment(origin, at(far), gold ? 3 : 1.6, color);
  }
  for (const point of points) {
    const [x, y] = at(point);
    pen.disc(x, y, 2.8, ink.fg);
  }
}
