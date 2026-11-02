import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 32;
const COUNT = SIDE * SIDE;
const NAME = "research-random-readings";

let memo: { thue: boolean[][]; rudin: boolean[][] } | undefined;

function ones(m: number) {
  let c = 0;
  for (let x = m; x > 0; x >>>= 1) c += x & 1;
  return c;
}

const thue = (m: number) => (ones(m) & 1) === 1;
const rudin = (m: number) => (ones(m & (m >>> 1)) & 1) === 1;

function pair(read: (n: number) => boolean, at: (n: number) => number) {
  const out: boolean[][] = Array.from({ length: SIDE }, () => new Array<boolean>(2 * SIDE).fill(false));
  for (let n = 0; n < COUNT; n++) {
    const row = Math.floor(n / SIDE);
    const col = n % SIDE;
    out[row][col] = read(at(n));
    out[row][SIDE + col] = read(at(2 * n));
  }
  return out;
}

function count(mask: boolean[][], test: (row: boolean[], col: number) => boolean) {
  let c = 0;
  for (const row of mask) for (let col = 0; col < SIDE; col++) if (test(row, col)) c++;
  return c;
}

function masks() {
  if (memo) return memo;
  if (!thue(7) || thue(9) || !rudin(3) || rudin(5) || !rudin(11) || rudin(27))
    throw new Error(`${NAME}: the digit rules broke`);
  const t = pair(thue, (n) => n * n);
  const r = pair(rudin, (n) => n * n + n);
  const tones = count(t, (row, col) => row[col]);
  const tapart = count(t, (row, col) => row[col] !== row[SIDE + col]);
  const rones = count(r, (row, col) => row[col]);
  const rapart = count(r, (row, col) => row[col] !== row[SIDE + col]);
  if (tapart !== 0) throw new Error(`${NAME}: t(4n^2) differs from t(n^2) at ${tapart} cells`);
  if ([tones, rones, rapart].join() !== "570,502,529")
    throw new Error(`${NAME}: counts ${tones},${rones},${rapart}, want 570,502,529`);
  memo = { thue: t, rudin: r };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { thue: t, rudin: r } = masks();
  const area = pen.frame(0.07);
  const gutter = Math.max(2, Math.round(area.w * 0.03));
  const cell = Math.floor((Math.min(area.w, area.h) - gutter) / (2 * SIDE));
  const block = cell * SIDE;
  const span = 2 * block + gutter;
  const left = area.x + (area.w - span) / 2;
  const top = area.y + (area.h - span) / 2;
  const split = (mask: boolean[][], half: number) => mask.map((row) => row.slice(half * SIDE, (half + 1) * SIDE));
  for (const [mask, y, color] of [
    [t, top, ink.blue],
    [r, top + block + gutter, ink.orange],
  ] as const) {
    new Grid(frame(left, y, block, block), SIDE, SIDE, 0).carpet(pen, split(mask, 0), color);
    new Grid(frame(left + block + gutter, y, block, block), SIDE, SIDE, 0).carpet(pen, split(mask, 1), color);
  }
}
