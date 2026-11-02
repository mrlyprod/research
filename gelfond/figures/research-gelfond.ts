import type { Ink, Pen, Point } from "@mrlyprod/draw";

const BASE = 10;
const DIGITS = [0, 1, 2, 3, 4, 5, 6, 8, 9];
const FILL = DIGITS.length;
const EXCLUDED = BASE - FILL;
const SAMPLES = 1800;
const DEPTHS = [2, 4, 8, 16, 32];
const LATTICE = 10;

// FACTS

const mean = DIGITS.reduce((s, d) => s + d, 0) / FILL;
const variance = DIGITS.reduce((s, d) => s + (d - mean) ** 2, 0) / FILL;

function kernel(x: number) {
  return Math.max(Math.sin(Math.PI * BASE * x) / Math.sin(Math.PI * x), 1 / Math.sin(Math.PI / BASE));
}

function pair(x: number) {
  return FILL - 4 * Math.sin((Math.PI * x) / 2) ** 2;
}

function constant(alpha: number) {
  const y = (BASE - 1) * alpha;
  const beta = Math.abs(y - Math.round(y));
  if (beta < 1e-12) return 0;
  const x0 = beta / (BASE + 1);
  const xf = 1 / (2 * (BASE - 1));
  const g3 = Math.max(FILL * Math.sqrt(1 - 16 * x0 * x0 * variance), Math.min(EXCLUDED + kernel(xf), pair(xf)));
  return 1 - Math.min(pair(x0), EXCLUDED + kernel(x0), g3) / FILL;
}

function constants() {
  const cs = Array.from({ length: SAMPLES + 1 }, (_, i) => constant(i / SAMPLES));
  const full = cs.filter((c) => c === 0).length;
  if (full !== LATTICE) throw new Error(`research-gelfond: ${full} lattice samples, want ${LATTICE}`);
  if (!(cs[SAMPLES / 2] >= 0.15)) throw new Error(`research-gelfond: c(1/2) = ${cs[SAMPLES / 2]}, want at least 0.15`);
  if (cs.some((c) => c < 0 || c >= 1)) throw new Error("research-gelfond: a constant outside [0, 1)");
  return cs;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const cs = constants();
  const frame = pen.frame(0.08);
  const at = (alpha: number, v: number): Point => [frame.x + frame.w * alpha, frame.y + frame.h * (1 - v)];
  for (let r = 0; r < LATTICE; r++) pen.segment(at(r / 9, 1), at(r / 9, 0), 1.4, ink.line);
  pen.segment(at(0, 0), at(1, 0), 1.4, ink.line);
  DEPTHS.forEach((k, j) => {
    const tone = ink.mix(ink.blue, ink.ground, 0.62 * (1 - j / (DEPTHS.length - 1)));
    const pts = cs.map((c, i): Point => at(i / SAMPLES, (1 - c) ** (k / 2)));
    pen.polyline(pts, 3.2, tone);
  });
  for (let r = 0; r < LATTICE; r++) {
    const [x, y] = at(r / 9, 1);
    pen.disc(x, y, 11, ink.orange);
  }
}
