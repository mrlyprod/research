import { type Color, field, type Ink, type Pen, type Point } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";

export const units = { num };

const TOP = 2048;
const DEEP = 32;
const CIRCLES = 2448;
const RESTING = 323;
const THICK = 1.7;
const HAIR = 1.2;
const RULE = 2;
const OVER = 0.035;

type Circle = { k: number; x: number; y: number };
type Packing = { seed: Circle[]; grown: Circle[]; resting: Circle[]; depth: Map<number, number> };

let memo: Packing | undefined;

function packing() {
  if (memo) return memo;
  const root = num.apollonian.root("strip");
  const seed = root.map((c) => ({ k: Number(c.k), x: Number(c.x), y: Number(c.y) }));
  const fords = root.slice(2).every((c) => num.apollonian.is_ford(c) && c.k === 2n);
  for (const c of root) c.free();
  const grown = num.apollonian.grow("strip", TOP);
  const touches = num.apollonian.touches(grown);
  const phi = num.factor.totients(DEEP);
  let nodes = -1;
  for (let n = 1; n <= DEEP; n++) nodes += Number(phi[n]);
  const resting = grown.circles.filter((c) => c.y === 1);
  if (grown.broken !== 0)
    throw new Error(`research-apollonian: ${grown.broken} quadruples break the Descartes form, want none`);
  if (grown.circles.length !== CIRCLES)
    throw new Error(`research-apollonian: ${grown.circles.length} circles, want ${CIRCLES}`);
  if (resting.length !== RESTING) throw new Error(`research-apollonian: ${resting.length} resting, want ${RESTING}`);
  if (resting.length !== nodes)
    throw new Error(`research-apollonian: ${resting.length} resting, want ${nodes} Farey nodes`);
  if (!grown.circles.every((c) => c.k > 0 && c.x > 0 && c.x < c.k))
    throw new Error("research-apollonian: a circle leaves the open period");
  if (touches.length !== resting.length)
    throw new Error(`research-apollonian: ${touches.length} touches, want ${resting.length}`);
  const depth = new Map<number, number>();
  let deepest = 1;
  for (const { num: a, den: b, k } of touches) {
    if (num.factor.gcd(a, b) !== "1") throw new Error(`research-apollonian: ${a}/${b} is not reduced`);
    if (!(0 < a && a < b && b <= DEEP)) throw new Error(`research-apollonian: ${a}/${b} is outside the stack`);
    depth.set(k, b);
    deepest = Math.max(deepest, b);
  }
  if (deepest !== DEEP) throw new Error(`research-apollonian: deepest ${deepest}, want ${DEEP}`);
  for (const c of resting) {
    const circle = num.apollonian.Circle.from(c);
    const ford = num.apollonian.is_ford(circle);
    circle.free();
    if (!ford) throw new Error(`research-apollonian: ${c.k},${c.x},${c.y} is not a Ford circle`);
    if (!depth.has(c.k)) throw new Error(`research-apollonian: no touch at curvature ${c.k}`);
  }
  if (!fords || seed[2].x !== 0 || seed[3].x !== 2)
    throw new Error("research-apollonian: the seed circles are not the Ford circles over 0 and 1");
  depth.set(2, 1);
  memo = { seed, grown: grown.circles, resting, depth };
  return memo;
}

function mark(pen: Pen, at: Point, r: number, thick: number, span: [number, number], color: Color) {
  const half = thick / 2;
  if (at[0] - r - half >= span[0] && at[0] + r + half <= span[1]) {
    pen.ring(at[0], at[1], r, thick, color);
    return;
  }
  const reach = r + half + 1;
  const x0 = Math.floor(Math.max(at[0] - reach, span[0] - 1, 0));
  const x1 = Math.min(Math.ceil(Math.max(Math.min(at[0] + reach, span[1] + 1), 0)), pen.width);
  const y0 = Math.floor(Math.max(at[1] - reach, 0));
  const y1 = Math.min(Math.ceil(Math.max(at[1] + reach, 0)), pen.height);
  const cut = field.patch(x0, y0, Math.max(x1 - x0, 0), Math.max(y1 - y0, 0));
  for (let py = y0; py < y1; py++) {
    for (let px = x0; px < x1; px++) {
      const fx = px + 0.5;
      const fy = py + 0.5;
      const ring = Math.abs(Math.sqrt((fx - at[0]) ** 2 + (fy - at[1]) ** 2) - r) - half;
      const clip = Math.max(span[0] - fx, fx - span[1]);
      cut.blend(px, py, color, 0.5 - Math.max(ring, clip));
    }
  }
  cut.paint(pen);
}

export default function draw(pen: Pen, ink: Ink) {
  const { seed, grown, resting, depth } = packing();
  const frame = pen.frame(0.08);
  const span: [number, number] = [frame.x, frame.x + frame.w];
  const over = frame.w * OVER;
  const place = (c: Circle): [Point, number] => [
    [frame.x + (frame.w * c.x) / c.k, frame.y + frame.h * (1 - c.y / c.k)],
    frame.w / c.k,
  ];
  for (const foot of [frame.y, frame.y + frame.h])
    pen.segment([span[0] - over, foot], [span[1] + over, foot], RULE, ink.dim);
  const rest = ink.fade(ink.dim, 0.75);
  for (const c of grown) {
    if (c.y === 1) continue;
    const [at, r] = place(c);
    mark(pen, at, r, HAIR, span, rest);
  }
  const ramp = ink.Ramp.tone(ink.blue, ink.fg);
  const reach = Math.log(DEEP);
  for (const c of [...resting, seed[2], seed[3]]) {
    const [at, r] = place(c);
    mark(pen, at, r, THICK, span, ramp.at(Math.log(depth.get(c.k)!) / reach));
  }
}
