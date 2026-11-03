import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 10;
const CELLS = SIDE * SIDE;
const BINARY = [1, 2, 4, 8].reduce((r, v) => r | (1 << v), 0);
const DEPTHS = [1, 4, 8, 12, 16, 18, 16, 12, 8, 4, 1];

const at = (i: number, j: number) => (((i % SIDE) + SIDE) % SIDE) * SIDE + (((j % SIDE) + SIDE) % SIDE);
const bit = (t: number) => (BINARY >> t) & 1;
const ones = (v: number) => {
  let c = 0;
  for (; v; v &= v - 1) c++;
  return c;
};

function reduce(rows: number[], v: number) {
  for (const r of rows) if ((v ^ r) < v) v ^= r;
  return v;
}

function insert(rows: number[], v: number) {
  const r = reduce(rows, v);
  if (r) {
    rows.push(r);
    rows.sort((a, b) => b - a);
  }
}

function reads(fill: boolean[]) {
  const out: [number, number[]][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!fill[a]) continue;
    const seen: number[] = [];
    for (let di = -1; di <= 1; di++)
      for (let dj = -1; dj <= 1; dj++) {
        const q = at(Math.floor(a / SIDE) + di, (a % SIDE) + dj);
        if ((di || dj) && !fill[q]) seen.push(q);
      }
    const live = seen.some((_, t) => bit(t + 1) !== bit(0));
    if (live) out.push([a, seen]);
  }
  return out;
}

function graph(colour: boolean[]) {
  const out: number[][] = Array.from({ length: CELLS }, () => []);
  for (const fill of [colour, colour.map((c) => !c)]) {
    const written = reads(fill);
    for (let p = 0; p < CELLS; p++) {
      const readers = written.filter(([, q]) => q.includes(p));
      if (!readers.length) continue;
      const window = [...new Set(readers.flatMap(([, q]) => q))].filter((q) => q !== p);
      const local = readers.map(([, q]) => window.reduce((m, c, k) => (q.includes(c) ? m | (1 << k) : m), 0));
      const vector = (z: number) =>
        local.reduce((v, m, i) => {
          const t = ones(z & m);
          return v | ((bit(t) ^ bit(t + 1)) << i);
        }, 0);
      const base = vector(0);
      const rows: number[] = [];
      for (let z = 1; z < 1 << window.length; z++) insert(rows, vector(z) ^ base);
      readers.forEach(([a], i) => {
        if (reduce(rows, 1 << i) === 0) out[p].push(a);
      });
    }
  }
  return out;
}

function depths() {
  const colour = Array.from({ length: CELLS }, (_, p) => (Math.floor(p / SIDE) + (p % SIDE)) % 2 === 0);
  const edges = graph(colour);
  const depth = new Array<number>(CELLS).fill(-1);
  const queue = [0];
  depth[0] = 0;
  for (let k = 0; k < queue.length; k++)
    for (const q of edges[queue[k]])
      if (depth[q] < 0) {
        depth[q] = depth[queue[k]] + 1;
        queue.push(q);
      }
  const count = new Array<number>(DEPTHS.length).fill(0);
  for (const d of depth) {
    if (d < 0 || d >= DEPTHS.length) throw new Error(`research-random-schedules: a cell at depth ${d}`);
    count[d]++;
  }
  if (count.join() !== DEPTHS.join())
    throw new Error(`research-random-schedules: depths ${count.join()}, want ${DEPTHS.join()}`);
  return { colour, depth };
}

export default function draw(pen: Pen, ink: Ink) {
  const { colour, depth } = depths();
  const area = pen.frame(0.07);
  const cell = Math.floor(Math.min(area.w, area.h) / SIDE);
  const left = area.x + (area.w - cell * SIDE) / 2;
  const top = area.y + (area.h - cell * SIDE) / 2;
  const grid = new Grid(frame(left, top, cell * SIDE, cell * SIDE), SIDE, SIDE, 0.12);
  const far = DEPTHS.length - 1;
  const shift = (p: number) => at(Math.floor(p / SIDE) + SIDE / 2, (p % SIDE) + SIDE / 2);
  for (let p = 0; p < CELLS; p++) {
    const q = shift(p);
    const alpha = 1 - (0.82 * depth[p]) / far;
    grid.fill(pen, q % SIDE, Math.floor(q / SIDE), ink.fade(colour[p] ? ink.blue : ink.orange, alpha));
  }
}
