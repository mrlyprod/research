import type { Ink, Pen, Point } from "@mrlyprod/draw";
import { plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";

export const units = { num, math };

const LENGTH = 16;
const STEPS = 2 ** LENGTH;
const DRAWN = 4096;
const REACH = Math.sqrt(STEPS);

let memo: { low: Float64Array; high: Float64Array } | undefined;

// FACTS

function walk(digits: number[]) {
  const values = math.design.elements(4, digits, LENGTH);
  if (values.length !== STEPS - 1)
    throw new Error(`research-mobius: ${values.length} elements of ${digits}, want ${STEPS - 1}`);
  const mu = [num.factor.mobius(0), ...Array.from(values, (v) => num.factor.mobius(Number(v)))];
  return Float64Array.from([0, ...Array.from(math.design.meter(mu), Number)]);
}

function walks() {
  if (memo) return memo;
  const low = walk([0, 1]);
  const high = walk([0, 2]);
  for (let i = 0; i <= STEPS; i++)
    if (high[i] !== -low[i])
      throw new Error(`research-mobius: step ${i} walks ${low[i]} and ${high[i]}, want a mirror`);
  const peak = low.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  if (!(peak < REACH)) throw new Error(`research-mobius: peak ${peak}, want under ${REACH}`);
  memo = { low, high };
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { low, high } = walks();
  const frame = pen.frame(0.08);
  const at = (index: number, meter: number): Point => [
    frame.x + (frame.w * index) / STEPS,
    frame.y + frame.h * (0.5 - meter / (2 * REACH)),
  ];
  pen.segment(at(0, 0), at(STEPS, 0), 1.6, ink.line);
  plot.axis(pen, frame, ink.line);
  for (const sign of [1, -1]) {
    const envelope = Array.from({ length: DRAWN + 1 }, (_, k): Point => {
      const index = (STEPS * k) / DRAWN;
      return at(index, sign * Math.sqrt(index));
    });
    pen.polyline(envelope, 2, ink.fade(ink.dim, 0.75));
  }
  const stride = STEPS / DRAWN;
  for (const [steps, color] of [
    [low, ink.blue],
    [high, ink.orange],
  ] as const) {
    pen.polyline(
      Array.from({ length: DRAWN + 1 }, (_, k) => at(k * stride, steps[k * stride])),
      2.4,
      color,
    );
  }
}
