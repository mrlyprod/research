import type { Ink, Pen, Point } from "@mrlyprod/draw";
import { plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const LEVEL = 6;

let memo: { values: Float64Array; ones: number } | undefined;

// FACTS

function spectrum() {
  if (memo) return memo;
  const cell = math.two.create(7, 2, LEVEL, 0, 2);
  const network = math.graph.core_graph({ shape: cell.shape, data: cell.types });
  const values = math.spectrum.laplacian_spectrum(network, true);
  const ones = math.spectrum.multiplicity(values, 1, 1e-9);
  if (values.length !== 729) throw new Error(`research-complexity: ${values.length} eigenvalues, want 729`);
  if (ones !== 243) throw new Error(`research-complexity: ${ones} ones, want 243`);
  memo = { values, ones };
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { values, ones } = spectrum();
  const frame = pen.frame(0.08);
  plot.axis(pen, frame, ink.line);
  const total = values.length;
  const at = (lambda: number, rank: number): Point => [
    frame.x + (frame.w * lambda) / 2,
    frame.y + frame.h * (1 - rank / total),
  ];
  const steps: Point[] = [at(0, 0)];
  for (const [index, value] of values.entries()) steps.push(at(value, index), at(value, index + 1));
  steps.push(at(2, total));
  pen.polyline(steps, 3, ink.blue);
  const below = values.reduce((n, v) => n + (v < 1 - 1e-9 ? 1 : 0), 0);
  pen.segment(at(1, below), at(1, below + ones), 7, ink.yellow);
}
