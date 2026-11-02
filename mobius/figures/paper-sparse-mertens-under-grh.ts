import type { Ink, Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./paper-sparse-mertens-under-grh.json" with { type: "json" };

export const units = { math };

const BASE = 7;
const DROP = 3;
const LEVEL = 3;

function masses(kept: number[], span: number, half: number) {
  const weights = math.dissection.weights(BASE, kept, LEVEL);
  if (weights.length !== span)
    throw new Error(`paper-sparse-mertens-under-grh: ${weights.length} weights, want ${span}`);
  return Array.from({ length: span }, (_, i) => weights[(((i - half) % span) + span) % span]);
}

export default function draw(pen: Pen, ink: Ink) {
  const { bound } = census as { bound: number };
  const kept = Array.from({ length: BASE }, (_, d) => d).filter((d) => d !== DROP);
  const span = BASE ** LEVEL;
  const peak = kept.length ** LEVEL;
  const half = (span - 1) / 2;
  const mass = masses(kept, span, half);

  if (kept.length !== BASE - 1)
    throw new Error(`paper-sparse-mertens-under-grh: ${kept.length} digits kept, want ${BASE - 1}`);
  if (!(Math.abs(mass[half] - peak) < 1e-9))
    throw new Error(`paper-sparse-mertens-under-grh: centre ${mass[half]}, want ${peak}`);
  const energy = mass.reduce((sum, v) => sum + v * v, 0);
  if (!(Math.abs(energy - span * peak) < 1e-6 * span * peak))
    throw new Error(`paper-sparse-mertens-under-grh: energy ${energy}, want ${span * peak}`);
  const total = mass.reduce((sum, v) => sum + v, 0);
  if (!(total < bound)) throw new Error(`paper-sparse-mertens-under-grh: total ${total} is not under the proved bound`);
  if (!(total >= span)) throw new Error(`paper-sparse-mertens-under-grh: total ${total} is under ${span}`);

  const frame = pen.frame(0.07);
  const ramp = ink.Ramp.tone(ink.blue, ink.yellow);
  const slot = frame.w / mass.length;
  const pad = slot * 0.18;
  const axis = frame.y + frame.h / 2;
  mass.forEach((value, i) => {
    const reach = (frame.h / 2) * (value / peak) ** (1 / LEVEL);
    pen.rect(frame.x + i * slot + pad, axis - reach, slot - 2 * pad, 2 * reach, ramp.at((2 * reach) / frame.h));
  });
}
