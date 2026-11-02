import type { Frame, Ink, Pen, Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LEVEL = 6;
const POINTS = 4096;
const DIGITS = 4;
const NORM = 9n;
const PROBE = 2;
const STEP = 1 / 9;
const TOL = 1e-12;
const MARGIN = 0.08;
const THIN = 1.4;

// THE STEP

function steps(pts: Point[]) {
  return pts.slice(1).map(([x, y], i) => Math.sqrt((x - pts[i][0]) ** 2 + (y - pts[i][1]) ** 2));
}

// THE LAYOUT

function fit(pts: Point[], frame: Frame): Point[] {
  let low: Point = [Number.MAX_VALUE, Number.MAX_VALUE];
  let high: Point = [-Number.MAX_VALUE, -Number.MAX_VALUE];
  for (const [x, y] of pts) {
    low = [Math.min(low[0], x), Math.min(low[1], y)];
    high = [Math.max(high[0], x), Math.max(high[1], y)];
  }
  const span = [high[0] - low[0], high[1] - low[1]];
  const scale = Math.min(frame.w / span[0], frame.h / span[1]);
  const mid = [(low[0] + high[0]) / 2, (low[1] + high[1]) / 2];
  const [cx, cy] = frame.center();
  return pts.map(([x, y]) => [cx + (x - mid[0]) * scale, cy - (y - mid[1]) * scale]);
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const design = math.radix.koch();
  if (design.size() !== DIGITS) throw new Error(`research-beneath: ${design.size()} digits, want ${DIGITS}`);
  if (design.base().norm() !== NORM)
    throw new Error(`research-beneath: base norm ${design.base().norm()}, want ${NORM}`);
  if (design.fill(LEVEL) !== String(POINTS))
    throw new Error(`research-beneath: fill ${design.fill(LEVEL)}, want ${POINTS}`);

  const probe = design.plane(PROBE);
  if (probe.length !== DIGITS * DIGITS)
    throw new Error(`research-beneath: ${probe.length} probe points, want ${DIGITS * DIGITS}`);
  if (!steps(probe).every((s) => Math.abs(s - STEP) < TOL))
    throw new Error(`research-beneath: the probe steps are not ${STEP}`);

  const curve = design.plane(LEVEL);
  if (curve.length !== POINTS) throw new Error(`research-beneath: ${curve.length} points, want ${POINTS}`);

  const frame = pen.frame(MARGIN);
  pen.polyline(fit(curve, frame), THIN, ink.blue);
}
