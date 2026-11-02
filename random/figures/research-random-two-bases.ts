import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 64;
const LAGS = SIDE * SIDE;
const NAME = "research-random-two-bases";

let memo: { binary: boolean[][]; ternary: boolean[][]; both: boolean[][] } | undefined;

function thue(h: number) {
  const bits = h === 0 ? 0 : 32 - Math.clz32(h);
  let x = 3;
  let y = -1;
  for (let i = bits - 1; i >= 0; i--) {
    const s = -(x + y);
    if ((h >> i) & 1) {
      x = s;
      y *= 2;
    } else {
      x *= 2;
      y = s;
    }
  }
  return { x, bits };
}

function twos(h: number) {
  let m = 0;
  for (;;) {
    if (h === 0) return m;
    if (h === 1) return -1;
    const r = h % 3;
    h = r === 1 ? (h - 1) / 3 + 1 : (h - r) / 3;
    if (r) m++;
  }
}

function masks() {
  if (memo) return memo;
  for (const [h, x, bits] of [
    [1, -2, 1],
    [2, -4, 2],
    [3, 4, 2],
    [5, 0, 3],
  ]) {
    const got = thue(h);
    if (got.x !== x || got.bits !== bits) throw new Error(`${NAME}: C_t(${h}) is ${got.x}/(3 2^${got.bits})`);
  }
  if (twos(2) !== 1 || twos(4) !== 2 || twos(1) !== -1) throw new Error(`${NAME}: the ternary recursion broke`);
  const grid = () => Array.from({ length: SIDE }, () => new Array<boolean>(SIDE).fill(false));
  const binary = grid();
  const ternary = grid();
  const both = grid();
  let counts = [0, 0, 0];
  for (let h = 0; h < LAGS; h++) {
    const { x, bits } = thue(h);
    const m = twos(h);
    const two = 16 * Math.abs(x) >= 3 * 2 ** bits;
    const three = m >= 0 && m <= 3;
    const row = Math.floor(h / SIDE);
    const col = h % SIDE;
    if (two && three) both[row][col] = true;
    else if (two) binary[row][col] = true;
    else if (three) ternary[row][col] = true;
    counts = [counts[0] + +two, counts[1] + +three, counts[2] + +(two && three)];
  }
  if (counts.join() !== "436,214,63") throw new Error(`${NAME}: counts ${counts}, want 436,214,63`);
  memo = { binary, ternary, both };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { binary, ternary, both } = masks();
  const area = pen.frame(0.07);
  const cell = Math.floor(Math.min(area.w, area.h) / SIDE);
  const span = cell * SIDE;
  const box = frame(area.x + (area.w - span) / 2, area.y + (area.h - span) / 2, span, span);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, binary, ink.blue);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, ternary, ink.orange);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, both, ink.fg);
}
