import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 16;
const CELLS = SIDE * SIDE;
const LETTER = 14;
const RULE = [1, 3, 6, 8];
const SEED = 0x6d616368696e65n;
const LOAD = [11, 8];
const TILES = 4;
const OFFSETS = [
  [-1, -1],
  [-1, 0],
  [1, 1],
];
const SLIDE = [6, 12, 18, 24];
const BRANCH = [4, 8, 11, 15];
const LINE = [17, 34, 51, 68];
const RING = [44, 88, 133, 177];

function fill() {
  const mask = (1n << 64n) - 1n;
  let state = SEED;
  const cells = new Uint8Array(CELLS);
  for (let w = 0; w < CELLS / 64; w++) {
    state = (state + 0x9e3779b97f4a7c15n) & mask;
    let z = state;
    z = ((z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n) & mask;
    z = ((z ^ (z >> 27n)) * 0x94d049bb133111ebn) & mask;
    z ^= z >> 31n;
    for (let b = 0; b < 64; b++) cells[64 * w + b] = Number((z >> BigInt(b)) & 1n);
  }
  return cells;
}

const at = (i: number, j: number) => (((i % SIDE) + SIDE) % SIDE) * SIDE + (((j % SIDE) + SIDE) % SIDE);
const slide = (p: number) => at((p >> 4) + 1, (p & 15) + 1);

function feistel() {
  const side = new Uint8Array(CELLS);
  for (let p = 0; p < CELLS; p++) {
    let on = 1;
    for (let l = 0; l < 4; l++) on &= (LETTER >> (2 * (((p >> 4) >> l) & 1) + (((p & 15) >> l) & 1))) & 1;
    side[p] = on;
  }
  const reads: [number, number[]][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!side[a]) continue;
    const near: number[] = [];
    for (let di = -1; di <= 1; di++)
      for (let dj = -1; dj <= 1; dj++) {
        const q = at((a >> 4) + di, (a & 15) + dj);
        if ((di || dj) && !side[q]) near.push(q);
      }
    reads.push([a, near]);
  }
  const total = reads.reduce((s, [, near]) => s + near.length, 0);
  if (reads.length !== 81 || total !== 346)
    throw new Error(`research-random-machine: ${reads.length} written cells reading ${total}, want 81 and 346`);
  return (x: Uint8Array) => {
    for (const [a, near] of reads) {
      let count = 0;
      for (const q of near) count += x[q];
      if (RULE.includes(count)) x[a] ^= 1;
    }
  };
}

function machine(sum: boolean) {
  const mess = feistel();
  return (x: Uint8Array) => {
    mess(x);
    const y = new Uint8Array(CELLS);
    for (let p = 0; p < CELLS; p++) {
      if (!sum) y[slide(p)] = x[p];
      else for (const [a, b] of OFFSETS) y[p] ^= x[at((p >> 4) - a, (p & 15) - b)];
    }
    return y;
  };
}

function rule30(moved: boolean) {
  return (x: Uint8Array) => {
    const y = new Uint8Array(CELLS);
    for (let p = 0; p < CELLS; p++) {
      const v = x[(p + CELLS - 1) % CELLS] ^ (x[p] | x[(p + 1) % CELLS]);
      y[moved ? slide(p) : p] = v;
    }
    return y;
  };
}

function differences(step: (x: Uint8Array) => Uint8Array, flip: number, keep: number[]) {
  let a: Uint8Array = fill();
  let b: Uint8Array = fill();
  b[flip] ^= 1;
  const out: Uint8Array[] = [];
  for (let r = 1; r <= keep[keep.length - 1]; r++) {
    a = step(a);
    b = step(b);
    if (keep.includes(r)) out.push(a.map((v, p) => v ^ b[p]));
  }
  return out;
}

export default function draw(pen: Pen, ink: Ink) {
  const block = (p: number) => p >> 4 < 8 && (p & 15) < 8;
  const run = (p: number) => p < 64;
  const load = LOAD[0] * SIDE + LOAD[1];
  const tiles = [
    ...differences(machine(false), load, SLIDE).map((d) => ({ d, out: block })),
    ...differences(machine(true), load, BRANCH).map((d) => ({ d, out: block })),
    ...differences(rule30(true), load, LINE).map((d) => ({ d, out: block })),
    ...differences(rule30(false), CELLS / 2 + 16, RING).map((d) => ({ d, out: run })),
  ];
  const area = pen.frame(0.07);
  const gap = area.w * 0.04;
  const tile = Math.floor((area.w - (TILES - 1) * gap) / TILES / SIDE) * SIDE;
  const span = TILES * tile + (TILES - 1) * gap;
  const left = area.x + (area.w - span) / 2;
  const top = area.y + (area.h - span) / 2;
  tiles.forEach(({ d, out }, t) => {
    const mask = (want: (p: number) => boolean) =>
      Array.from({ length: SIDE }, (_, i) => Array.from({ length: SIDE }, (_, j) => want(SIDE * i + j)));
    const grid = new Grid(
      frame(left + (t % TILES) * (tile + gap), top + Math.floor(t / TILES) * (tile + gap), tile, tile),
      SIDE,
      SIDE,
      0,
    );
    grid.carpet(
      pen,
      mask((p) => out(p) && !d[p]),
      ink.line,
    );
    grid.carpet(
      pen,
      mask((p) => !out(p) && d[p] === 1),
      ink.blue,
    );
    grid.carpet(
      pen,
      mask((p) => out(p) && d[p] === 1),
      ink.orange,
    );
  });
}
