import type { Ink, Pen } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";
import census from "./research-stack.json" with { type: "json" };

const RING = "Gaussian";
const TOP = 46;
const WRAPPED = 917;

export const units = { num };

type Node = [number, number, number];

let memo: { top: number; lit: Map<number, number> } | undefined;

function classes(bound: number) {
  return num.gauss.classes(RING, bound).length;
}

function lights(bound: number, nodes: Node[]) {
  if (memo) return memo;
  const top = classes(bound);
  if (top !== TOP) throw new Error(`research-stack: ${top} classes, want ${TOP}`);
  const lit = new Map<number, number>();
  for (const [, , n] of nodes) if (!lit.has(n)) lit.set(n, classes(Math.floor(bound / n)));
  memo = { top, lit };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { bound, nodes } = census as { bound: number; nodes: Node[] };
  if (nodes.some(([x, y, n]) => !(n > 0 && x >= 0 && x < n && y >= 0 && y < n)))
    throw new Error("research-stack: a node is off its residue square");
  const { top, lit } = lights(bound, nodes);
  const stack: [number, number, number][] = [];
  for (const [xn, yn, n] of nodes) {
    const [x, y, light] = [xn / n, yn / n, lit.get(n) as number];
    stack.push([x, y, light]);
    if (xn === 0) stack.push([1, y, light]);
    if (yn === 0) stack.push([x, 1, light]);
    if (xn === 0 && yn === 0) stack.push([1, 1, light]);
  }
  if (stack.length !== WRAPPED) throw new Error(`research-stack: ${stack.length} wrapped nodes, want ${WRAPPED}`);
  stack.sort((a, b) => a[2] - b[2]);
  const ramp = new ink.Ramp([ink.ground, ink.blue, ink.yellow]);
  const frame = pen.frame(0.08);
  for (const [x, y, light] of stack) {
    const share = light / top;
    const [px, py] = frame.at(x, 1 - y);
    pen.disc(px, py, 3.5 + 15 * share ** 0.7, ramp.at(share ** 0.3));
  }
}
