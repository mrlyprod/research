import type { Ink, Pen, Point } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const BASE = 2;
const CODE = 7;
const LEVEL = 2;
const DIGITS = 3;
const CELLS = 9;
const DRAWS = CELLS + LEVEL - 1;
const TAIL = 24000;
const SEED = 3;
const MARGIN = 0.08;
const DOT = 2;
const THICK = 2.2;
const KNOT = 5.5;

// THE SEQUENCE

function bruijn(m: number, k: number): number[] {
  const a = new Array<number>(m * k).fill(0);
  const out: number[] = [];
  const walk = (t: number, p: number) => {
    if (t > k) {
      if (k % p === 0) out.push(...a.slice(1, p + 1));
      return;
    }
    a[t] = a[t - p]!;
    walk(t + 1, p);
    for (let j = a[t - p]! + 1; j < m; j++) {
      a[t] = j;
      walk(t + 1, t);
    }
  };
  walk(1, 1);
  return out;
}

// THE GAME

function play(digits: Point[], draws: number[], start: Point): Point[] {
  const path: Point[] = [start];
  for (const d of draws) {
    const [x, y] = path[path.length - 1]!;
    path.push([(x + digits[d]![0]) / BASE, (y + digits[d]![1]) / BASE]);
  }
  return path;
}

const cellOf = ([x, y]: Point) => `${Math.floor(x * BASE ** LEVEL)},${Math.floor(y * BASE ** LEVEL)}`;

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const design = math.radix.tile(BASE, CODE);
  if (design.size() !== DIGITS) throw new Error(`research-chaos-game: ${design.size()} digits, want ${DIGITS}`);
  if (design.fill(LEVEL) !== String(CELLS))
    throw new Error(`research-chaos-game: fill ${design.fill(LEVEL)}, want ${CELLS}`);
  const digits = design.digits().map(([x, y]): Point => [Number(x), Number(y)]);
  const level = new Set(design.plane(LEVEL).map((p) => cellOf(p)));
  if (level.size !== CELLS) throw new Error(`research-chaos-game: ${level.size} cells, want ${CELLS}`);

  const word = bruijn(DIGITS, LEVEL);
  if (word.length !== CELLS) throw new Error(`research-chaos-game: de Bruijn word of ${word.length}, want ${CELLS}`);
  const lap = [...word, ...word.slice(0, LEVEL - 1)];
  const path = play(digits, lap, [0, 0]);
  if (path.length !== DRAWS + 1) throw new Error(`research-chaos-game: ${path.length} points, want ${DRAWS + 1}`);
  const painted = path.slice(LEVEL).map(cellOf);
  if (new Set(painted).size !== CELLS || !painted.every((c) => level.has(c)))
    throw new Error(
      `research-chaos-game: the lap paints ${new Set(painted).size} cells, want the ${CELLS} of the design`,
    );

  const rng = new math.Rng(SEED);
  const draws = Array.from({ length: TAIL }, () => rng.below(DIGITS));
  const tail = play(digits, draws, path[path.length - 1]!).slice(1);
  if (!tail.every((p) => level.has(cellOf(p)))) throw new Error("research-chaos-game: the game left the design");

  const area = pen.frame(MARGIN).square();
  const side = BASE ** LEVEL;
  const step = area.w / side;
  const at = ([x, y]: Point): Point => [area.x + x * area.w, area.y + (1 - y) * area.h];

  pen.rect(area.x, area.y, area.w, area.h, ink.panel);
  painted.forEach((c, i) => {
    const [cx, cy] = c.split(",").map(Number) as [number, number];
    const shade = ink.fade(ink.blue, 0.12 + (0.6 * i) / (CELLS - 1));
    pen.rect(area.x + cx * step, area.y + (side - 1 - cy) * step, step, step, shade);
  });
  for (const p of tail) {
    const [x, y] = at(p);
    pen.rect(x - DOT / 2, y - DOT / 2, DOT, DOT, ink.fg);
  }
  pen.polyline(path.map(at), THICK, ink.orange);
  for (const p of path.slice(1)) {
    const [x, y] = at(p);
    pen.disc(x, y, KNOT, ink.orange);
  }
}
