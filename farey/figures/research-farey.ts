import type { Ink, Pen } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";

export const units = { num };

const ORDER = 60;

type Stack = { x: number; y: number; weight: number }[];

let memo: Stack | undefined;

function stack() {
  if (memo) return memo;
  const nodes = num.lattice.farey(ORDER);
  if (nodes.length !== 1103) throw new Error(`research-farey: ${nodes.length} nodes, want 1103`);
  const grids = new Map<number, number>();
  const grid = (a: number, b: number) => {
    const key = a * (ORDER + 1) + b;
    let found = grids.get(key);
    if (found === undefined) {
      found = num.factor.lcm(a, b);
      grids.set(key, found);
    }
    return found;
  };
  const found: Stack = [];
  for (const across of nodes) {
    for (const down of nodes) {
      const lcm = grid(across.den, down.den);
      if (lcm > ORDER) continue;
      found.push({ x: across.num / across.den, y: down.num / down.den, weight: Math.floor(ORDER / lcm) });
    }
  }
  if (found.length !== 63261) throw new Error(`research-farey: ${found.length} stacked, want 63261`);
  memo = found.sort((a, b) => a.weight - b.weight);
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const ramp = new ink.Ramp([ink.ground, ink.blue, ink.yellow]);
  const frame = pen.frame(0.08);
  for (const { x, y, weight } of stack()) {
    const share = weight / ORDER;
    const [px, py] = frame.at(x, 1 - y);
    pen.disc(px, py, 1 + 3 * share ** 0.7, ramp.at(share ** 0.3));
  }
}
