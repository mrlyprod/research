import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 4;
const CELLS = SIDE * SIDE;
const STATES = 1 << CELLS;
const CARPET = 7;
const RULE = [2, 3, 5, 7].reduce((r, v) => r | (1 << v), 0);
const TILES = 8;

function word(outer: number, inner: number) {
  let mask = 0;
  for (let i = 0; i < SIDE; i++)
    for (let j = 0; j < SIDE; j++) {
      const a = (outer >> (2 * (i >> 1) + (j >> 1))) & 1;
      const b = (inner >> (2 * (i & 1) + (j & 1))) & 1;
      mask |= (a & b) << (SIDE * i + j);
    }
  return mask;
}

function machine() {
  const fill = word(CARPET, CARPET);
  if (fill !== 0x135f) throw new Error(`research-random: the level-2 carpet is ${fill}, want 0x135f`);
  const reads: [number, number][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!((fill >> a) & 1)) continue;
    let mask = 0;
    for (let di = -1; di <= 1; di++)
      for (let dj = -1; dj <= 1; dj++) {
        if (!di && !dj) continue;
        const q = ((((a >> 2) + di + SIDE) % SIDE) << 2) | (((a & 3) + dj + SIDE) % SIDE);
        mask |= 1 << q;
      }
    reads.push([a, mask & ~fill & (STATES - 1)]);
  }
  const slide = Array.from({ length: CELLS }, (_, p) => ((((p >> 2) + 1) % SIDE) << 2) | (((p & 3) + 1) % SIDE));
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
    return z ^ fill;
  };
}

function census(round: (x: number) => number) {
  const next = new Uint16Array(STATES);
  const seen = new Uint8Array(STATES);
  for (let x = 0; x < STATES; x++) {
    next[x] = round(x);
    if (seen[next[x]]++) throw new Error(`research-random: the round is not a bijection at ${x}`);
  }
  seen.fill(0);
  let cycles = 0;
  let longest = 0;
  let fixed = 0;
  for (let s = 0; s < STATES; s++) {
    if (seen[s]) continue;
    let len = 0;
    for (let x = s; !seen[x]; x = next[x]) {
      seen[x] = 1;
      len++;
    }
    cycles++;
    longest = Math.max(longest, len);
    if (len === 1) fixed++;
  }
  if (cycles !== 258 || longest !== 24865 || fixed !== 3)
    throw new Error(`research-random: ${cycles} cycles, longest ${longest}, ${fixed} fixed; want 258, 24865, 3`);
  return next;
}

export default function draw(pen: Pen, ink: Ink) {
  const next = census(machine());
  const fill = word(CARPET, CARPET);
  const area = pen.frame(0.07);
  const gap = area.w * 0.025;
  const tile = Math.floor((area.w - (TILES - 1) * gap) / TILES / SIDE) * SIDE;
  const span = TILES * tile + (TILES - 1) * gap;
  const left = area.x + (area.w - span) / 2;
  const top = area.y + (area.h - span) / 2;
  let x = 1;
  for (let t = 0; t < TILES * TILES; t++) {
    const at = (want: (p: number) => boolean) =>
      Array.from({ length: SIDE }, (_, i) => Array.from({ length: SIDE }, (_, j) => want(SIDE * i + j)));
    const on = (p: number) => ((x >> p) & 1) === 1;
    const side = (p: number) => ((fill >> p) & 1) === 1;
    const grid = new Grid(
      frame(left + (t % TILES) * (tile + gap), top + Math.floor(t / TILES) * (tile + gap), tile, tile),
      SIDE,
      SIDE,
      0,
    );
    grid.carpet(
      pen,
      at((p) => side(p) && !on(p)),
      ink.line,
    );
    grid.carpet(
      pen,
      at((p) => side(p) && on(p)),
      ink.blue,
    );
    grid.carpet(
      pen,
      at((p) => !side(p) && on(p)),
      ink.orange,
    );
    x = next[x];
  }
}
