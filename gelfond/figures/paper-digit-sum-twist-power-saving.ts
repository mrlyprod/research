import { field, type Ink, type Pen, type Point } from "@mrlyprod/draw";

const BASE = 10;
const DIGITS = [0, 1, 2, 3, 4, 5, 6, 8, 9];
const FILL = DIGITS.length;
const EXCLUDED = BASE - FILL;
const TWIST = 1 / 2;
const PATH = 200000;
const GAMMA = 0.85;

// FACTS

function digit(v: number) {
  let re = 0;
  let im = 0;
  for (const d of DIGITS) {
    re += Math.cos(2 * Math.PI * d * v);
    im += Math.sin(2 * Math.PI * d * v);
  }
  return Math.hypot(re, im);
}

function norm(v: number) {
  return Math.abs(v - Math.round(v));
}

function kernel(x: number) {
  return Math.max(Math.sin(Math.PI * BASE * x) / Math.sin(Math.PI * x), 1 / Math.sin(Math.PI / BASE));
}

function pair(x: number) {
  return FILL - 4 * Math.sin((Math.PI * x) / 2) ** 2;
}

function constant(beta: number) {
  const mean = DIGITS.reduce((s, d) => s + d, 0) / FILL;
  const variance = DIGITS.reduce((s, d) => s + (d - mean) ** 2, 0) / FILL;
  const x0 = beta / (BASE + 1);
  const xf = 1 / (2 * (BASE - 1));
  const g3 = Math.max(FILL * Math.sqrt(1 - 16 * x0 * x0 * variance), Math.min(EXCLUDED + kernel(xf), pair(xf)));
  return 1 - Math.min(pair(x0), EXCLUDED + kernel(x0), g3) / FILL;
}

function facts() {
  const beta = norm((BASE - 1) * TWIST);
  const reach = beta / (BASE + 1);
  const c = constant(beta);
  if (Math.abs(reach - 1 / 22) > 1e-15)
    throw new Error(`paper-digit-sum-twist-power-saving: reach ${reach}, want 1/22`);
  if (!(c >= 0.15)) throw new Error(`paper-digit-sum-twist-power-saving: c = ${c}, want at least 0.15`);
  if (Math.abs(digit(0) * digit(0) - FILL * FILL) > 1e-9)
    throw new Error("paper-digit-sum-twist-power-saving: corner not full");
  let far = 1;
  let share = 0;
  for (let i = 0; i < PATH; i++) {
    const u = i / PATH;
    const v1 = u + TWIST;
    const v2 = BASE * u + TWIST;
    far = Math.min(far, Math.max(norm(v1), norm(v2)));
    share = Math.max(share, (digit(v1) * digit(v2)) / (FILL * FILL));
  }
  if (far < reach - 1e-12)
    throw new Error(`paper-digit-sum-twist-power-saving: path within ${far} of the corner, want ${reach}`);
  if (share > 1 - c + 1e-12)
    throw new Error(`paper-digit-sum-twist-power-saving: path keeps ${share}, above 1 - c = ${1 - c}`);
}

function segments(shift: number) {
  const cuts = [-0.5, 0.5];
  for (let k = -BASE - 1; k <= BASE + 1; k++) {
    const v = (k + shift + 0.5) / BASE;
    if (v > -0.5 + 1e-12 && v < 0.5 - 1e-12) cuts.push(v);
  }
  cuts.sort((a, b) => a - b);
  const out: [number, number][][] = [];
  for (let i = 0; i + 1 < cuts.length; i++) {
    const [p, q] = [cuts[i], cuts[i + 1]];
    const n = Math.round(BASE * ((p + q) / 2) - shift);
    out.push([
      [p, BASE * p - shift - n],
      [q, BASE * q - shift - n],
    ]);
  }
  return out;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  facts();
  const frame = pen.frame(0.08);
  const side = Math.round(frame.w);
  const table = new Float64Array(side);
  for (let i = 0; i < side; i++) table[i] = digit((i + 0.5) / side - 0.5) / FILL;
  const values = new Float64Array(side * side);
  for (let row = 0; row < side; row++) {
    const up = table[side - 1 - row];
    for (let col = 0; col < side; col++) values[row * side + col] = (up * table[col]) ** GAMMA;
  }
  field.draw_range(pen, frame, side, side, values, [0, 1], new ink.Ramp([ink.ground, ink.blue, ink.yellow]));
  const at = ([v1, v2]: [number, number]): Point => [frame.x + frame.w * (v1 + 0.5), frame.y + frame.h * (0.5 - v2)];
  const lattice = segments(0);
  const twisted = segments(((BASE - 1) * TWIST) % 1);
  if (lattice.length !== BASE + 1 || twisted.length !== BASE)
    throw new Error(
      `paper-digit-sum-twist-power-saving: ${lattice.length} and ${twisted.length} segments, want ${BASE + 1} and ${BASE}`,
    );
  const through = lattice.some(
    ([p, q]) => p[0] < 0 && q[0] > 0 && Math.abs(p[1] - (p[0] * (q[1] - p[1])) / (q[0] - p[0])) < 1e-12,
  );
  if (!through) throw new Error("paper-digit-sum-twist-power-saving: no lattice segment through the centre");
  for (const [p, q] of lattice) pen.segment(at(p), at(q), 2.6, ink.dim);
  for (const [p, q] of twisted) pen.segment(at(p), at(q), 4.4, ink.orange);
}
