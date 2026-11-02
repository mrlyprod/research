import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const LEVEL = 5;
const SIDE = 1 << LEVEL;
const CELLS = SIDE * SIDE;
const ROUNDS = 8;
const TILES = 4;

// CONE

type Spread = (i: number, j: number) => [number, number];

const shift: Spread = (i, j) => [i, j + 1];
const hyper: Spread = (i, j) => [2 * i + j + 1, 3 * i + 2 * j];

function board() {
  const fill = Array.from({ length: CELLS }, (_, p) => ((p >> LEVEL) + (p & (SIDE - 1))) % 2 === 0);
  const reads: [number, number[]][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!fill[a]) continue;
    const r: number[] = [];
    for (let di = -1; di <= 1; di++)
      for (let dj = -1; dj <= 1; dj++) {
        const q = ((((a >> LEVEL) + di + SIDE) % SIDE) << LEVEL) | (((a & (SIDE - 1)) + dj + SIDE) % SIDE);
        if (q !== a && !fill[q]) r.push(q);
      }
    reads.push([a, r]);
  }
  return { fill, reads };
}

function cone(spread: Spread) {
  const { fill, reads } = board();
  const map = Array.from({ length: CELLS }, (_, p) => {
    const [i, j] = spread(p >> LEVEL, p & (SIDE - 1));
    return ((i % SIDE) << LEVEL) | (j % SIDE);
  });
  if (new Set(map).size !== CELLS) throw new Error("research-random-mixing: a spread that is not a bijection");
  const out: boolean[][] = [];
  let set = new Array<boolean>(CELLS).fill(false);
  set[0] = true;
  for (let t = 0; t < ROUNDS; t++) {
    const grown = set.slice();
    for (const [a, r] of reads) if (r.some((q) => set[q])) grown[a] = true;
    set = new Array<boolean>(CELLS).fill(false);
    for (let q = 0; q < CELLS; q++) if (grown[q]) set[map[q]] = true;
    out.push(set);
  }
  return { fill, out };
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const slow = cone(shift);
  const fast = cone(hyper);
  if (!fast.out[ROUNDS - 1].every((b) => b))
    throw new Error(`research-random-mixing: the hyperbolic spread does not cover the torus in ${ROUNDS} rounds`);
  if (slow.out[ROUNDS - 1].filter((b) => b).length > (2 * ROUNDS + 1) ** 2)
    throw new Error("research-random-mixing: the shift leaves its light cone");
  const area = pen.frame(0.07);
  const gap = area.w * 0.03;
  const tile = Math.floor((area.w - (TILES - 1) * gap) / TILES / SIDE) * SIDE;
  const span = TILES * tile + (TILES - 1) * gap;
  const left = area.x + (area.w - span) / 2;
  const top = area.y + (area.h - span) / 2;
  for (let t = 0; t < TILES * TILES; t++) {
    const run = t < 2 * TILES ? slow : fast;
    const set = run.out[t % ROUNDS];
    const at = (want: (p: number) => boolean) =>
      Array.from({ length: SIDE }, (_, i) => Array.from({ length: SIDE }, (_, j) => want(SIDE * i + j)));
    const grid = new Grid(
      frame(left + (t % TILES) * (tile + gap), top + Math.floor(t / TILES) * (tile + gap), tile, tile),
      SIDE,
      SIDE,
      0,
    );
    grid.carpet(
      pen,
      at((p) => !set[p]),
      ink.line,
    );
    grid.carpet(
      pen,
      at((p) => set[p] && run.fill[p]),
      ink.blue,
    );
    grid.carpet(
      pen,
      at((p) => set[p] && !run.fill[p]),
      ink.orange,
    );
  }
}
