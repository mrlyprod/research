import { frame as box, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LEVEL = 6;
const SIDE = 1 << LEVEL;
const STRIP = 5;
const GAP = 2;
const LUCAS = 11;
const ANF = "00100110001100100110001101110011";

function primes() {
  const flags = new Array<boolean>(SIDE).fill(true);
  flags[0] = flags[1] = false;
  for (let p = 2; p * p < SIDE; p++) {
    if (flags[p]) for (let m = p * p; m < SIDE; m += p) flags[m] = false;
  }
  return flags;
}

function lucas() {
  const cells = math.two.create(LUCAS, 2, LEVEL, 0, 2);
  if (cells.shape[0] !== SIDE || cells.shape[1] !== SIDE)
    throw new Error(`research-random-mess: level ${LEVEL} is ${cells.shape}, want ${SIDE} square`);
  const mask: boolean[][] = [];
  for (let row = 0; row < SIDE; row++) {
    mask.push([]);
    for (let col = 0; col < SIDE; col++) {
      const on = cells.types[row * SIDE + col] !== 0;
      if (on !== ((row & col) === row))
        throw new Error(
          `research-random-mess: cell ${row},${col} is not the binary domination of the row by the column`,
        );
      mask[row].push(on);
    }
  }
  return mask;
}

function parity(mask: boolean[][], chosen: boolean[]) {
  const out: boolean[] = [];
  for (let col = 0; col < SIDE; col++) {
    let odd = false;
    for (let row = 0; row < SIDE; row++) if (chosen[row] && mask[row][col]) odd = !odd;
    out.push(odd);
  }
  const head = out
    .slice(0, ANF.length)
    .map((bit) => (bit ? "1" : "0"))
    .join("");
  if (head !== ANF) throw new Error(`research-random-mess: the prime parities begin ${head}, want ${ANF}`);
  return out;
}

export default function draw(pen: Pen, ink: Ink) {
  const mask = lucas();
  const chosen = primes();
  const result = parity(mask, chosen);
  const area = pen.frame(0.07);
  const rows = SIDE + GAP + STRIP;
  const cell = Math.floor(Math.min(area.w / SIDE, area.h / rows));
  const left = area.x + (area.w - cell * SIDE) / 2;
  const top = area.y + (area.h - cell * rows) / 2;
  const grid = new Grid(box(left, top, cell * SIDE, cell * rows), SIDE, rows, 0);
  const quiet = ink.fade(ink.blue, 0.55);
  for (let row = 0; row < SIDE; row++) {
    for (let col = 0; col < SIDE; col++) {
      if (mask[row][col]) grid.fill(pen, col, row, chosen[row] ? ink.orange : quiet);
    }
  }
  for (let col = 0; col < SIDE; col++) {
    if (!result[col]) continue;
    for (let band = 0; band < STRIP; band++) grid.fill(pen, col, SIDE + GAP + band, ink.green);
  }
}
