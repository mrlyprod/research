import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";

const SIDE = 64;
const CELLS = SIDE * SIDE;
const NAME = "research-random-squares";

let memo: { order: boolean[][]; noise: boolean[][] } | undefined;

function thue(n: number) {
  let p = 0;
  for (let x = n; x > 0; x = Math.floor(x / 2)) p ^= x % 2;
  return p;
}

function blocks(bits: number[]) {
  const seen = new Set<number>();
  for (let n = 0; n + 8 <= bits.length; n++) {
    let w = 0;
    for (let i = 0; i < 8; i++) w |= bits[n + i] << i;
    seen.add(w);
  }
  return seen.size;
}

function masks() {
  if (memo) return memo;
  const along = Array.from({ length: CELLS + 7 }, (_, n) => thue(n));
  const squares = Array.from({ length: CELLS + 7 }, (_, n) => thue(n * n));
  const seen = [blocks(along), blocks(squares)];
  if (seen.join() !== "22,256") throw new Error(`${NAME}: 8-blocks ${seen}, want 22,256`);
  const bias = squares.slice(0, CELLS).reduce((s, b) => s + 1 - 2 * b, 0);
  if (bias !== -260) throw new Error(`${NAME}: S(4096) is ${bias}, want -260`);
  const grid = () => Array.from({ length: SIDE }, () => new Array<boolean>(SIDE).fill(false));
  const order = grid();
  const noise = grid();
  for (let i = 0; i < SIDE; i++) {
    for (let j = 0; j < SIDE; j++) {
      const n = SIDE * i + j;
      if (j > i) order[i][j] = along[n] === 1;
      else if (j < i) noise[i][j] = squares[n] === 1;
    }
  }
  memo = { order, noise };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { order, noise } = masks();
  const area = pen.frame(0.07);
  const cell = Math.floor(Math.min(area.w, area.h) / SIDE);
  const span = cell * SIDE;
  const box = frame(area.x + (area.w - span) / 2, area.y + (area.h - span) / 2, span, span);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, order, ink.blue);
  new Grid(box, SIDE, SIDE, 0).carpet(pen, noise, ink.orange);
}
