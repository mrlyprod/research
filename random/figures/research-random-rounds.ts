import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 3;
const CELLS = SIDE * SIDE;
const STATES = 1 << CELLS;
const FILL = 0x5f;
const RULE = [2, 3, 5, 7].reduce((r, v) => r | (1 << v), 0);
const NET: [number, number][] = [
  [-1, 0],
  [1, 0],
  [0, -1],
  [0, 1],
];
const COLS = 32;
const ROWS = 16;
const NAME = "research-random-rounds";

let memo: { order: number[]; cycle: number[] } | undefined;

function machine() {
  const reads: [number, number][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!((FILL >> a) & 1)) continue;
    let mask = 0;
    for (const [di, dj] of NET) {
      const i = (Math.floor(a / SIDE) + di + SIDE) % SIDE;
      const j = ((a % SIDE) + dj + SIDE) % SIDE;
      mask |= 1 << (SIDE * i + j);
    }
    reads.push([a, mask & ~FILL & (STATES - 1)]);
  }
  if (reads.length !== 6) throw new Error(`${NAME}: the split writes ${reads.length} cells, want 6`);
  const slide = Array.from(
    { length: CELLS },
    (_, p) => SIDE * ((Math.floor(p / SIDE) + 1) % SIDE) + (((p % SIDE) + 1) % SIDE),
  );
  const ones = (v: number) => {
    let c = 0;
    for (; v; v &= v - 1) c++;
    return c;
  };
  return (x: number) => {
    let y = x;
    for (const [a, mask] of reads) y ^= ((RULE >> ones(x & mask)) & 1) << a;
    let z = 0;
    for (let p = 0; p < CELLS; p++) z |= ((y >> p) & 1) << slide[p];
    return z ^ FILL;
  };
}

function census() {
  if (memo) return memo;
  const round = machine();
  const next = new Uint16Array(STATES);
  const seen = new Uint8Array(STATES);
  for (let x = 0; x < STATES; x++) {
    next[x] = round(x);
    if (seen[next[x]]++) throw new Error(`${NAME}: the round is not a bijection at ${x}`);
  }
  seen.fill(0);
  const cycles: number[][] = [];
  for (let s = 0; s < STATES; s++) {
    if (seen[s]) continue;
    const run: number[] = [];
    for (let x = s; !seen[x]; x = next[x]) {
      seen[x] = 1;
      run.push(x);
    }
    cycles.push(run);
  }
  cycles.sort((a, b) => b.length - a.length || a[0] - b[0]);
  const longest = cycles[0].length;
  const fixed = cycles.filter((c) => c.length === 1).length;
  if (cycles.length !== 62 || longest !== 36 || fixed !== 4)
    throw new Error(`${NAME}: ${cycles.length} cycles, longest ${longest}, ${fixed} fixed; want 62, 36, 4`);
  if (cycles.reduce((s, c) => s + c.length - 1, 0) % 2) throw new Error(`${NAME}: the round is an odd permutation`);
  const order: number[] = [];
  const cycle: number[] = [];
  cycles.forEach((run, k) => {
    for (const x of run) {
      order.push(x);
      cycle.push(k);
    }
  });
  if (order.length !== COLS * ROWS) throw new Error(`${NAME}: ${order.length} states laid out, want ${COLS * ROWS}`);
  memo = { order, cycle };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { order, cycle } = census();
  const area = pen.frame(0.07);
  const gap = area.w * 0.004;
  const tile = Math.floor((area.w - (COLS - 1) * gap) / COLS / SIDE) * SIDE;
  const wide = COLS * tile + (COLS - 1) * gap;
  const high = ROWS * tile + (ROWS - 1) * gap;
  const left = area.x + (area.w - wide) / 2;
  const top = area.y + (area.h - high) / 2;
  for (let t = 0; t < order.length; t++) {
    const x = order[t];
    const box = frame(left + (t % COLS) * (tile + gap), top + Math.floor(t / COLS) * (tile + gap), tile, tile);
    pen.rect(box.x, box.y, box.w, box.h, cycle[t] % 2 ? ink.orange : ink.blue);
    const live = Array.from({ length: SIDE }, (_, i) =>
      Array.from({ length: SIDE }, (_, j) => ((x >> (SIDE * i + j)) & 1) === 1),
    );
    new Grid(box, SIDE, SIDE, 0.45).carpet(pen, live, ink.fg);
  }
}
