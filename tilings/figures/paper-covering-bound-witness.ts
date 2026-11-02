import type { Ink, Pen, Point } from "@mrlyprod/draw";

export const units = {};

const F = [0, 10, 11, 20, 23, 35, 45, 47, 50, 55, 59];
const RINGS: { m: number; fibres: number[][] }[] = [
  {
    m: 15,
    fibres: [
      [2, 5, 8, 11, 14],
      [0, 5, 10],
      [0, 5, 10],
    ],
  },
  {
    m: 18,
    fibres: [
      [5, 11, 17],
      [0, 9],
      [1, 10],
      [2, 11],
      [5, 14],
    ],
  },
  {
    m: 20,
    fibres: [
      [3, 7, 11, 15, 19],
      [0, 10],
      [0, 10],
      [5, 15],
    ],
  },
];
const CENTERS: Point[] = [
  [0.27, 0.3],
  [0.73, 0.3],
  [0.5, 0.71],
];

function vanishes(residues: number[], m: number): boolean {
  let re = 0;
  let im = 0;
  for (const r of residues) {
    re += Math.cos((2 * Math.PI * r) / m);
    im += Math.sin((2 * Math.PI * r) / m);
  }
  return Math.hypot(re, im) < 1e-9;
}

function multiset(values: number[]): Map<number, number> {
  const counts = new Map<number, number>();
  for (const v of values) counts.set(v, (counts.get(v) ?? 0) + 1);
  return counts;
}

function check({ m, fibres }: { m: number; fibres: number[][] }) {
  const residues = F.map((f) => f % m);
  if (!vanishes(residues, m)) throw new Error(`paper-covering-bound-witness: the residues mod ${m} do not vanish`);
  const want = multiset(residues);
  const got = multiset(fibres.flat());
  for (const [r, n] of want)
    if (got.get(r) !== n)
      throw new Error(`paper-covering-bound-witness: residue ${r} mod ${m} drawn ${got.get(r) ?? 0} times, want ${n}`);
  if (got.size !== want.size) throw new Error(`paper-covering-bound-witness: stray residue mod ${m}`);
  for (const fibre of fibres) {
    const p = fibre.length;
    const step = m / p;
    if (m % p !== 0 || fibre.some((r, i) => (r - fibre[0] - i * step) % m !== 0))
      throw new Error(`paper-covering-bound-witness: ${fibre} is no ${p}-fibre mod ${m}`);
  }
  return want;
}

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const radius = frame.w * 0.21;
  const tones: Record<number, [number, number, number, number]> = { 2: ink.fg, 3: ink.orange, 5: ink.blue };
  RINGS.forEach((ring, k) => {
    const counts = check(ring);
    const [cx, cy] = frame.at(CENTERS[k][0], CENTERS[k][1]);
    const at = (r: number, scale: number): Point => {
      const angle = -Math.PI / 2 + (2 * Math.PI * r) / ring.m;
      return [cx + scale * radius * Math.cos(angle), cy + scale * radius * Math.sin(angle)];
    };
    pen.ring(cx, cy, radius, 1.5, ink.fg);
    const seen = new Map<string, number>();
    for (const fibre of ring.fibres) {
      const key = fibre.join(",");
      const times = seen.get(key) ?? 0;
      seen.set(key, times + 1);
      const scale = 1 - 0.1 * times;
      const pts = fibre.map((r) => at(r, scale));
      if (pts.length === 2) pen.segment(pts[0], pts[1], 2.5 + 3 * times, tones[2]);
      else pen.polyline([...pts, pts[0]], 2.5, tones[fibre.length]);
    }
    for (const [r, n] of counts) {
      const [x, y] = at(r, 1);
      pen.disc(x, y, n > 1 ? 9 : 6, ink.fg);
    }
  });
}
