import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-baselines.json" with { type: "json" };

const MARGIN = 0.08;
const GAP = 0.14;
const FADE = 0.62;

export const units = {};

type Point = [number, number];

const key = ([a, b]: Point) => `${a},${b}`;

function peel(m: number, points: Point[]): Map<string, number> {
  const half = (m + 1) / 2;
  const round = new Map<string, number>();
  let live = points.slice();
  for (let r = 0; live.length; r++) {
    const mids = new Set<string>();
    for (let i = 0; i < live.length; i++)
      for (let j = i + 1; j < live.length; j++)
        mids.add(key([((live[i][0] + live[j][0]) * half) % m, ((live[i][1] + live[j][1]) * half) % m]));
    const out = live.filter((p) => !mids.has(key(p)));
    if (!out.length) throw new Error(`research-baselines: ${live.length} points do not peel`);
    for (const p of out) round.set(key(p), r);
    live = live.filter((p) => mids.has(key(p)));
  }
  return round;
}

export default function draw(pen: Pen, ink: Ink) {
  const { m, rounds, points } = census as { m: number; rounds: number; points: Point[] };
  if (points.length !== 19) throw new Error(`research-baselines: ${points.length} points, want 19`);
  const round = peel(m, points);
  const depth = Math.max(...round.values()) + 1;
  if (depth !== rounds) throw new Error(`research-baselines: ${depth} peel rounds, want ${rounds}`);
  const grid = new Grid(pen.frame(MARGIN), m, m, GAP);
  const empty = ink.mix(ink.line, ink.ground, 0.45);
  for (let a = 0; a < m; a++)
    for (let b = 0; b < m; b++) {
      const r = round.get(key([a, b]));
      grid.fill(pen, a, m - 1 - b, r === undefined ? empty : ink.mix(ink.blue, ink.ground, (FADE * r) / (depth - 1)));
    }
}
