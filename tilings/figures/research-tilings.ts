import * as core from "@mrlyprod/core";
import { frame as box, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-tilings.json" with { type: "json" };

const LEVELS = 3;

export const units = { core };

export default function draw(pen: Pen, ink: Ink) {
  const owners: number[][] = census.owners;
  if (owners.length !== LEVELS) throw new Error(`research-tilings: ${owners.length} levels, want ${LEVELS}`);
  const frame = pen.frame(0.08);
  const unit = frame.w / 30;
  const [tint, blue, shadow] = core.colors.shades(core.colors.BLUE());
  const tones = [ink.yellow, tint, blue, shadow];
  let corner = 0;
  owners.forEach((owner, n) => {
    const side = 1 << (n + 2);
    if (owner.length !== side * side)
      throw new Error(`research-tilings: level ${n + 1} holds ${owner.length} cells, want ${side * side}`);
    const size = side * unit;
    const grid = new Grid(box(frame.x + corner, frame.y + corner, size, size), side, side, 0.14);
    owner.forEach((k, r) => {
      if (!(k >= 0 && k < 4)) throw new Error(`research-tilings: cell ${r} of level ${n + 1} owned by ${k}`);
      grid.fill(pen, r % side, Math.floor(r / side), tones[k]);
    });
    corner += size + unit;
  });
}
