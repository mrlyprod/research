import type { Ink, Pen } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";

export const units = { num };

const N = 100;

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const step = frame.w / N;
  const lit = step * 0.52;
  const hid = step * 0.24;
  let visible = 0;
  for (let a = 1; a <= N; a++) {
    for (let b = 1; b <= N; b++) {
      const x = frame.x + (a - 0.5) * step;
      const y = frame.y + frame.h - (b - 0.5) * step;
      if (num.factor.coprime(a, b)) {
        visible++;
        pen.rect(x - lit / 2, y - lit / 2, lit, lit, ink.yellow);
      } else {
        pen.rect(x - hid / 2, y - hid / 2, hid, hid, ink.fade(ink.dim, 0.9));
      }
    }
  }
  if (BigInt(visible) !== num.lattice.coprime_pairs(N))
    throw new Error(`research-pi: ${visible} visible, want ${num.lattice.coprime_pairs(N)}`);
}
