import type { Ink, Pen, Point } from "@mrlyprod/draw";
import { iso } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const SIDE = 3;
const HALF = SIDE / 2;

type Mark = { depth: number; face: iso.Face } | { depth: number; edges: [Point, Point][] };

function node(i: number, j: number, k: number) {
  return iso.project((i - HALF) / HALF, (j - HALF) / HALF, (k - HALF) / HALF);
}

function cage(c: number[]) {
  const [x, y, z] = c;
  const corner = (b: number) => node(x + (b & 1), y + ((b >> 1) & 1), z + ((b >> 2) & 1));
  return iso.EDGES.map(([a, b]): [Point, Point] => [corner(a), corner(b)]);
}

let memo: { marks: Mark[]; faces: iso.Face[] } | undefined;

function marks() {
  if (memo) return memo;
  const solid = math.three.carpet(SIDE, 1);
  if (solid.shape.length !== 3 || solid.shape.some((n) => n !== SIDE))
    throw new Error(`paper-divisor-avatars: the solid is ${solid.shape}, want ${SIDE} by ${SIDE} by ${SIDE}`);
  const filled = math.three.fills(solid);
  if (filled !== 20) throw new Error(`paper-divisor-avatars: the solid fills ${filled}, want 20`);
  const drilled: number[][] = [];
  for (let i = 0; i < SIDE; i++)
    for (let j = 0; j < SIDE; j++)
      for (let k = 0; k < SIDE; k++) if (solid.types[(i * SIDE + j) * SIDE + k] === 0) drilled.push([i, j, k]);
  if (drilled.length !== 7) throw new Error(`paper-divisor-avatars: ${drilled.length} drilled corners, want 7`);
  const faces = iso.faces(math.three.quads(solid));
  const list: Mark[] = faces.map((face) => ({ depth: face.depth, face }));
  for (const cell of drilled) {
    const mid = cell.reduce((total, v) => total + (v + 0.5), 0);
    list.push({ depth: (mid - 3 * HALF) / HALF, edges: cage(cell) });
  }
  list.sort((a, b) => a.depth - b.depth);
  memo = { marks: list, faces };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const { marks: list, faces } = marks();
  const put = iso.fit(faces, frame);
  const gold = [ink.yellow, ink.mix(ink.yellow, ink.ground, 0.4), ink.mix(ink.yellow, ink.ground, 0.66)];
  const unit = frame.h / (SIDE * 2);
  for (const mark of list) {
    if ("face" in mark) {
      const screen = mark.face.quad.map(put);
      pen.polygon(screen, gold[mark.face.tone]);
      pen.polyline([...screen, screen[0]], unit / 70, ink.line);
    } else {
      for (const [a, b] of mark.edges) pen.segment(put(a), put(b), unit / 34, ink.dim);
    }
  }
}
