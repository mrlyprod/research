import type { Ink, Pen } from "@mrlyprod/draw";
import * as num from "@mrlyprod/num";

export const units = { num };

const RADIUS = 30;
const VISIBLE = 1668;

export default function draw(pen: Pen, ink: Ink) {
  const frame = pen.frame(0.08);
  const ring = "Eisenstein";
  const [cx, cy] = frame.center();
  const unit = frame.w / (2 * RADIUS + 0.8);
  let total = 0;
  let visible = 0;
  for (let a = -RADIUS; a <= RADIUS; a++) {
    for (let b = -RADIUS; b <= RADIUS; b++) {
      if (num.gauss.Ring.reach(ring, a, b) > BigInt(RADIUS)) continue;
      total++;
      const [u, v] = num.gauss.Ring.place(ring, a, b);
      const [x, y] = [cx + u * unit, cy - v * unit];
      if (a === 0 && b === 0) {
        pen.disc(x, y, unit * 0.2, ink.fade(ink.dim, 0.9));
      } else if (num.factor.gcd(Math.abs(a), Math.abs(b)) === "1") {
        visible++;
        pen.disc(x, y, unit * 0.3, ink.blue);
      } else {
        pen.disc(x, y, unit * 0.13, ink.fade(ink.dim, 0.8));
      }
    }
  }
  if (total !== num.gauss.Ring.count(ring, RADIUS))
    throw new Error(`research-bases: ${total} points, want ${num.gauss.Ring.count(ring, RADIUS)}`);
  if (visible !== VISIBLE) throw new Error(`research-bases: ${visible} visible, want ${VISIBLE}`);
}
