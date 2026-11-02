import { Grid, type Ink, type Pen } from "@mrlyprod/draw";

const BASE = 19;
const DIGITS = [0, 1, 3, 12, 14, 15];
const SIDE = BASE;

function design(): number[] {
  const out: number[] = [];
  for (const high of DIGITS) for (const low of DIGITS) out.push(high * BASE + low);
  return out.sort((a, b) => a - b);
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const grid = new Grid(area, SIDE, SIDE, 0.12);
  const cells = design();
  const member = new Set(cells);
  const shadow = new Set<number>();
  for (const x of cells) for (const z of cells) if (x < z && (x + z) % 2 === 0) shadow.add((x + z) / 2);
  for (const n of shadow) {
    if (member.has(n)) throw new Error(`research-progressions: ${n} is a midpoint inside the design`);
    grid.fill(pen, n % BASE, SIDE - 1 - Math.floor(n / BASE), ink.orange);
  }
  for (const n of cells) grid.fill(pen, n % BASE, SIDE - 1 - Math.floor(n / BASE), ink.blue);
  if (cells.length !== DIGITS.length ** 2)
    throw new Error(`research-progressions: ${cells.length} cells, want ${DIGITS.length ** 2}`);
  if (!(shadow.size > cells.length && shadow.size < SIDE * SIDE - cells.length))
    throw new Error(`research-progressions: ${shadow.size} midpoints`);
}
