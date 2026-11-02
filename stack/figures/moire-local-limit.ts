import { type Frame, field, frame, type Ink, type Pen, type Point } from "@mrlyprod/draw";
import census from "./moire-local-limit.json" with { type: "json" };

const SIDE = 410;
const HALF = SIDE / 2;
const EDGE = 82;
const STEPS = 400;

export const units = {};

// FIELDS

function check() {
  const { octant, axis, arc, scale } = census;
  const folded = (HALF * (HALF + 1)) / 2;
  if (octant.length !== folded) throw new Error(`moire-local-limit: ${octant.length} octant values, want ${folded}`);
  if (axis.length !== SIDE) throw new Error(`moire-local-limit: ${axis.length} axis values, want ${SIDE}`);
  if (arc.length !== STEPS + 1) throw new Error(`moire-local-limit: ${arc.length} arc points, want ${STEPS + 1}`);
  if (!(scale > 0) || octant.some((n) => n < 0 || n > scale))
    throw new Error(`moire-local-limit: octant values off 0..${scale}`);
}

function overlay(): Float64Array {
  const slot = new Uint32Array(HALF * HALF);
  let next = 0;
  for (let i = 0; i < HALF; i++) {
    for (let j = i; j < HALF; j++) {
      slot[i * HALF + j] = next;
      slot[j * HALF + i] = next;
      next++;
    }
  }
  const values = new Float64Array(SIDE * SIDE);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      const i = Math.min(row, SIDE - 1 - row);
      const j = Math.min(col, SIDE - 1 - col);
      values[row * SIDE + col] = census.octant[slot[i * HALF + j]] / census.scale;
    }
  }
  return values;
}

function smooth(): Float64Array {
  const values = new Float64Array(SIDE * SIDE);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) values[row * SIDE + col] = 0.5 * (1 - census.axis[col] * census.axis[row]);
  }
  return values;
}

// ARCS

function arcs(pen: Pen, ink: Ink, box: Frame, thick: number, chords: boolean) {
  const kappa = census.kappa;
  for (const [sx, sy] of [
    [1, 1],
    [1, -1],
    [-1, 1],
    [-1, -1],
  ]) {
    if (chords) {
      const a = box.at(0.5 * (1 - sx), 0.5 * (1 - sy * kappa));
      const b = box.at(0.5 * (1 - sx * kappa), 0.5 * (1 - sy));
      pen.segment(a, b, 2, ink.dim);
    }
    const line: Point[] = census.arc.map(([x, y]) => box.at(0.5 * (1 - sx * x), 0.5 * (1 - sy * y)));
    pen.polyline(line, thick, ink.orange);
  }
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  check();
  const ramp = ink.Ramp.tone(ink.panel, ink.blue);
  const first = frame(EDGE, EDGE, SIDE, SIDE);
  const far = pen.width - EDGE - SIDE;
  const second = frame(far, far, SIDE, SIDE);
  field.draw_range(pen, first, SIDE, SIDE, overlay(), [0, 1], ramp);
  field.draw_range(pen, second, SIDE, SIDE, smooth(), [0, 0.5], ramp);
  arcs(pen, ink, second, 4, true);
  arcs(pen, ink, first, 3, false);
}
