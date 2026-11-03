import { field, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./research-spirograph.json" with { type: "json" };

type Wheel = { bounds: math.spirograph.Disc; shape: math.spirograph.Cover; points: Float32Array };

const CODE = 495;
const SIDE = 3;
const LEVEL = 1;
const BASE = 3;
const RING = 7;
const WHEEL = 3;
const REACH = 0.9;
const SAMPLES = 8000;
const CURVES = 8;
const CROSSINGS = 1288;
const RASTER = 1024;
const MARGIN = 0.08;
const SUB = 4;
const WASH = 0.3;
const THIN = 1.4;
const HAIR = 1.1;
const DOT = 2.6;

export const units = { math };

let memo: Wheel | undefined;

function wheel(): Wheel {
  if (memo) return memo;
  const tile = math.two.create(CODE, SIDE, LEVEL, 0, BASE);
  const path = math.spirograph.track("in", RING, WHEEL, 4, 1);
  const pens = math.spirograph.pencils(tile.types, tile.shape[1], tile.shape[0], "fill", REACH, 0, 1);
  const curves = math.spirograph.distinct(path, pens, true);
  if (curves !== CURVES) throw new Error(`research-spirograph: ${curves} curves, want ${CURVES}`);
  const law = math.spirograph.nodes(path, pens, true);
  if (law !== BigInt(CROSSINGS))
    throw new Error(`research-spirograph: the law counts ${law} crossings, want ${CROSSINGS}`);
  const points = math.spirograph.trace(path, pens, SAMPLES);
  if (points.length !== 2 * CURVES * SAMPLES)
    throw new Error(`research-spirograph: trace of ${points.length} numbers, want ${2 * CURVES * SAMPLES}`);
  const shape = math.spirograph.cover(path, pens, true, 2, RASTER);
  if (shape.side !== RASTER || shape.mask.length !== RASTER * RASTER)
    throw new Error(`research-spirograph: cover of side ${shape.side} and ${shape.mask.length} cells, want ${RASTER}`);
  memo = { bounds: math.spirograph.disc(path, pens), shape, points };
  return memo;
}

function wash(pen: Pen, ink: Ink, shape: math.spirograph.Cover, corner: Point, span: number) {
  const n = shape.side;
  const cell = span / n;
  const x0 = Math.floor(Math.max(corner[0], 0));
  const y0 = Math.floor(Math.max(corner[1], 0));
  const x1 = Math.min(Math.ceil(Math.max(corner[0] + span, 0)), pen.width);
  const y1 = Math.min(Math.ceil(Math.max(corner[1] + span, 0)), pen.height);
  const tint = field.patch(x0, y0, x1 - x0, y1 - y0);
  for (let py = y0; py < y1; py++) {
    for (let px = x0; px < x1; px++) {
      let hit = 0;
      for (let sy = 0; sy < SUB; sy++) {
        for (let sx = 0; sx < SUB; sx++) {
          const u = (px + (sx + 0.5) / SUB - corner[0]) / cell;
          const v = (py + (sy + 0.5) / SUB - corner[1]) / cell;
          if (u < 0 || v < 0) continue;
          const column = Math.trunc(u);
          const row = Math.trunc(v);
          if (column < n && row < n && shape.mask[row * n + column] === 3) hit++;
        }
      }
      if (hit > 0) tint.blend(px, py, ink.dim, (WASH * hit) / (SUB * SUB));
    }
  }
  tint.paint(pen);
}

export default function draw(pen: Pen, ink: Ink) {
  const { bounds, shape, points } = wheel();
  if (census.nodes.length !== CROSSINGS || census.nodes.some((node) => node.length !== 2))
    throw new Error(`research-spirograph: ${census.nodes.length} nodes, want ${CROSSINGS} pairs`);
  const area = pen.frame(MARGIN);
  const [cx, cy] = area.center();
  const scale = area.radius() / RING;
  const at = (x: number, y: number): Point => [cx + (x - bounds.x) * scale, cy - (y - bounds.y) * scale];
  const span = 2 * bounds.radius * scale;
  wash(pen, ink, shape, [cx - span / 2, cy - span / 2], span);
  pen.ring(cx, cy, RING * scale, HAIR, ink.line);
  for (let k = 0; k < CURVES; k++) {
    const curve = Array.from({ length: SAMPLES }, (_, i): Point => {
      const j = 2 * (k * SAMPLES + i);
      return at(points[j], points[j + 1]);
    });
    pen.polyline(curve, THIN, ink.blue);
  }
  plot.dots(
    pen,
    census.nodes.map(([x, y]) => at(x, y)),
    DOT,
    ink.orange,
  );
}
