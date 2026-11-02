import type { Ink, Pen } from "@mrlyprod/draw";

const SIDE = 4;
const CELLS = SIDE * SIDE;
const STATES = 1 << CELLS;
const CARPET = 7;
const FILL = 4959;
const STAMP = 1030;
const RULE = [2, 3, 5, 7].reduce((r, v) => r | (1 << v), 0);
const CYCLES = 12;
const LONGEST = 40739;

// ROUND

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

const ones = (v: number) => {
  let c = 0;
  for (; v; v &= v - 1) c++;
  return c;
};

function round() {
  const fill = word(CARPET, CARPET);
  if (fill !== FILL) throw new Error(`paper-design-rounds: the level-2 carpet is ${fill}, want ${FILL}`);
  const reads: [number, number][] = [];
  for (let a = 0; a < CELLS; a++) {
    if (!((fill >> a) & 1)) continue;
    let mask = 0;
    for (let di = -1; di <= 1; di++)
      for (let dj = -1; dj <= 1; dj++) {
        if (!di && !dj) continue;
        mask |= 1 << ((((a >> 2) + di + SIDE) % SIDE) * SIDE + (((a & 3) + dj + SIDE) % SIDE));
      }
    reads.push([a, mask & ~fill & (STATES - 1)]);
  }
  const slide = Array.from({ length: CELLS }, (_, p) => (((p >> 2) + 1) % SIDE) * SIDE + (((p & 3) + 1) % SIDE));
  const next = new Uint16Array(STATES);
  for (let x = 0; x < STATES; x++) {
    let y = x;
    for (const [a, mask] of reads) y ^= ((RULE >> ones(x & mask)) & 1) << a;
    let z = 0;
    for (let p = 0; p < CELLS; p++) z |= ((y >> p) & 1) << slide[p];
    next[x] = z ^ STAMP;
  }
  const hit = new Uint8Array(STATES);
  for (const y of next) if (hit[y]++) throw new Error(`paper-design-rounds: the round is not a bijection at ${y}`);
  return next;
}

// CYCLES

function prime(p: number) {
  if (p < 2) return false;
  for (let d = 2; d * d <= p; d++) if (p % d === 0) return false;
  return true;
}

function cycles(next: Uint16Array) {
  const seen = new Uint8Array(STATES);
  const lengths: number[] = [];
  for (let s = 0; s < STATES; s++) {
    if (seen[s]) continue;
    let len = 0;
    for (let x = s; !seen[x]; x = next[x]) {
      seen[x] = 1;
      len++;
    }
    lengths.push(len);
  }
  lengths.sort((a, b) => b - a);
  const [top] = lengths;
  if (lengths.length !== CYCLES || top !== LONGEST)
    throw new Error(`paper-design-rounds: ${lengths.length} cycles, longest ${top}; want ${CYCLES}, ${LONGEST}`);
  if (!(2 * top > STATES && top + 3 <= STATES && prime(top)))
    throw new Error(`paper-design-rounds: the longest cycle ${top} is not a prime in (M/2, M - 3]`);
  if ((STATES - lengths.length) % 2) throw new Error("paper-design-rounds: the round is odd");
  return lengths;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const lengths = cycles(round());
  const area = pen.frame(0.07);
  const [cx, cy] = area.center();
  let tail = STATES;
  const sides = lengths.map((len) => {
    const side = area.w * Math.sqrt(tail / STATES);
    tail -= len;
    return side;
  });
  sides.push(0);
  lengths.forEach((_, k) => {
    const outer = sides[k];
    const inner = sides[k + 1];
    const gap = Math.min(area.w * 0.004, (outer - inner) * 0.18);
    const color = k ? ink.fade(ink.orange, k % 2 ? 1 : 0.5) : ink.blue;
    pen.rect(cx - outer / 2, cy - outer / 2, outer, outer, color);
    if (inner > 0) pen.rect(cx - inner / 2 - gap, cy - inner / 2 - gap, inner + 2 * gap, inner + 2 * gap, ink.ground);
  });
}
