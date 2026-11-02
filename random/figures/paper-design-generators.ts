import type { Ink, Pen, Point } from "@mrlyprod/draw";

export const units = {};

const NAME = "paper-design-generators";
const SHRINK = 103 / 200;
const RUNS = 12;
const STEPS = 256;

type Mat = [number, number, number, number];

let memo: { outer: Point[]; inner: Point[]; images: Point[][] } | undefined;

function run(t: number): Mat {
  return [0, 1, t, (t - 1) / 3];
}

function form(x: number, y: number) {
  return 20 * x * x + 20 * x * y + 41 * y * y;
}

function ellipse(level: number): Point[] {
  const pts: Point[] = [];
  for (let i = 0; i <= STEPS; i++) {
    const a = (2 * Math.PI * i) / STEPS;
    const c = Math.cos(a);
    const s = Math.sin(a);
    const r = Math.sqrt(level / form(c, s));
    pts.push([r * c, r * s]);
  }
  return pts;
}

function apply(m: Mat, [x, y]: Point): Point {
  return [m[0] * x + m[1] * y, m[2] * x + m[3] * y];
}

function worst([p, q, r, s]: Mat) {
  const n11 = 20 * p * p + 20 * p * r + 41 * r * r;
  const n22 = 20 * q * q + 20 * q * s + 41 * s * s;
  const n12 = 20 * p * q + 10 * (p * s + q * r) + 41 * r * s;
  const lin = 41 * n11 + 20 * n22 - 20 * n12;
  const det = n11 * n22 - n12 * n12;
  return (lin + Math.sqrt(Math.max(0, lin * lin - 4 * 720 * det))) / (2 * 720);
}

function certificate() {
  const a = (t: number) => 103 / 10 - 41 * t * t;
  const d = (t: number) => 5807 / 1800 + (22 / 9) * t - (41 / 9) * t * t;
  const b = (t: number) => 103 / 10 + (22 / 3) * t - (82 / 3) * t * t;
  const delta = (t: number) => 30179 / 1125 - (2266 / 45) * t - (93367 / 450) * t * t;
  const want: [number, number][] = [
    [a(-1 / 2), 1 / 20],
    [a(1 / 4), 619 / 80],
    [d(-1 / 2), 173 / 200],
    [d(1 / 4), 1421 / 400],
    [delta(-1 / 2), 133 / 1000],
    [delta(1 / 4), 5077 / 4000],
  ];
  for (const [got, exact] of want)
    if (Math.abs(got - exact) > 1e-12) throw new Error(`${NAME}: certificate value ${got}, want ${exact}`);
  for (const t of [-0.5, -0.3, -0.1, 0, 0.1, 0.25]) {
    const gap = 4 * a(t) * d(t) - b(t) * b(t) - delta(t);
    if (Math.abs(gap) > 1e-9) throw new Error(`${NAME}: discriminant mismatch ${gap} at t = ${t}`);
  }
  const [x, y] = apply(run(-1 / 2), [1, -1 / 3]);
  if (Math.abs(x + 1 / 3) > 1e-15 || Math.abs(y + 1 / 3) > 1e-15) throw new Error(`${NAME}: A(-1/2) e_0 is ${x}, ${y}`);
}

function objects() {
  if (memo) return memo;
  certificate();
  const images: Point[][] = [];
  let tight = 0;
  let most = 0;
  for (let b = 1; b <= RUNS; b++) {
    const m = run((-1 / 2) ** b);
    const ratio = worst(m);
    if (ratio > SHRINK) throw new Error(`${NAME}: run of length ${b} stretches the form by ${ratio}`);
    if (ratio > most) [most, tight] = [ratio, b];
    images.push(ellipse(1).map((p) => apply(m, p)));
  }
  if (tight !== 2) throw new Error(`${NAME}: the tightest run has length ${tight}, not 2`);
  for (let i = 0; i <= 96; i++) {
    const ratio = worst(run(-1 / 2 + (3 / 4) * (i / 96)));
    if (ratio > SHRINK) throw new Error(`${NAME}: the form fails at sample ${i}, ratio ${ratio}`);
  }
  memo = { outer: ellipse(1), inner: ellipse(SHRINK), images };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { outer, inner, images } = objects();
  const box = pen.frame(0.08);
  let reach = 0;
  for (const [x, y] of outer) reach = Math.max(reach, Math.abs(x), Math.abs(y));
  const scale = box.w / (2 * reach);
  const [cx, cy] = box.center();
  const place = (pts: Point[]) => pts.map(([x, y]) => [cx + scale * x, cy - scale * y] as Point);
  const thin = Math.max(1, box.w / 400);
  for (const curve of images) pen.polyline(place(curve), thin, ink.orange);
  pen.polyline(place(inner), 2.5 * thin, ink.blue);
  pen.polyline(place(outer), 2.5 * thin, ink.fg);
}
