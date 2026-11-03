import type { Color, Ink, Pen, Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";

export const units = { num, math };

const BASE = 10;
const LEVEL = 3;
const Z = 8;
const COUNTS = [732, 202, 26, 40];
const NAMES = ["A", "B", "C1", "C2"];

let memo: number[] | undefined;

function regions() {
  if (memo) return memo;
  const y = BASE ** LEVEL;
  const found = math.dissection.regions(BASE, LEVEL, Z).map((name) => NAMES.indexOf(name));
  if (found.length !== y || found.includes(-1))
    throw new Error(`paper-unconditional-mertens-at-large-base: ${found.length} regions, want ${y} named ${NAMES}`);
  const cap = math.dissection.cap(y);
  for (let a = 0; a < y; a++) {
    const [l, d] = math.dissection.fraction(a, y, cap);
    const h = BigInt(a) * d - l * BigInt(y);
    if (d > cap || (h < 0n ? -h : h) * cap > BigInt(y) || num.factor.gcd(l, d) !== "1")
      throw new Error(
        `paper-unconditional-mertens-at-large-base: frequency ${a} has fraction ${l}/${d} off the Dirichlet bound`,
      );
  }
  const counts = [0, 0, 0, 0];
  for (const r of found) counts[r]++;
  if (counts.some((count, r) => count !== COUNTS[r]))
    throw new Error(`paper-unconditional-mertens-at-large-base: counts ${counts}, want ${COUNTS}`);
  memo = found;
  return memo;
}

function band(pen: Pen, center: Point, radii: Point, angles: Point, color: Color) {
  const at = (r: number, t: number): Point => [center[0] + r * Math.sin(t), center[1] - r * Math.cos(t)];
  const steps = Math.max(Math.ceil((angles[1] - angles[0]) / 0.004), 1);
  const pts: Point[] = [];
  for (let i = 0; i <= steps; i++) pts.push(at(radii[1], angles[0] + ((angles[1] - angles[0]) * i) / steps));
  for (let i = steps; i >= 0; i--) pts.push(at(radii[0], angles[0] + ((angles[1] - angles[0]) * i) / steps));
  pen.polygon(pts, color);
}

function runs(found: number[], pass: number) {
  const n = found.length;
  let start = found.findIndex((r) => r !== pass);
  if (start < 0) start = 0;
  const out: Point[] = [];
  let open: number | undefined;
  for (let step = 1; step <= n; step++) {
    const a = start + step;
    const inside = found[a % n] === pass;
    if (open === undefined && inside) open = a;
    else if (open !== undefined && !inside) {
      out.push([open, a]);
      open = undefined;
    }
  }
  if (open !== undefined) out.push([open, start + n + 1]);
  return out;
}

export default function draw(pen: Pen, ink: Ink) {
  const found = regions();
  const y = BASE ** LEVEL;
  const frame = pen.frame(0.08);
  const center = frame.center();
  const outer = frame.radius();
  const slot = (2 * Math.PI) / y;
  const width = 0.12 * outer;
  const tracks = [0.28, 0.48, 0.68, 0.88];
  const inks = [ink.blue, ink.dim, ink.orange, ink.yellow];
  for (let pass = 0; pass < 4; pass++) {
    const mid = tracks[pass] * outer + width / 2;
    pen.ring(center[0], center[1], mid, 1, ink.line);
    let held = 0;
    for (const [first, end] of runs(found, pass)) {
      held += end - first;
      band(pen, center, [mid - width / 2, mid + width / 2], [(first - 0.5) * slot, (end - 0.5) * slot], inks[pass]);
    }
    if (held !== COUNTS[pass])
      throw new Error(`paper-unconditional-mertens-at-large-base: pass ${pass} holds ${held}, want ${COUNTS[pass]}`);
  }
}
