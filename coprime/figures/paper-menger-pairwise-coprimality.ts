import { frame, Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

export default function draw(pen: Pen, ink: Ink) {
  const sponge = math.three.carpet(3, 1);
  if (sponge.shape.join() !== "3,3,3")
    throw new Error(`paper-menger-pairwise-coprimality: sponge ${sponge.shape}, want 3,3,3`);
  const types = sponge.types;
  const area = pen.frame(0.08);
  const step = area.w / 3;
  const side = step - 10;
  let kept = 0;
  let near = 0;
  for (let c = 0; c < 3; c++) {
    const panel = frame(area.x + c * step, area.y + c * step, side, side);
    const grid = new Grid(panel, 3, 3, 0.07);
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const [x, y, w, h] = grid.cell(b, a);
        if (types[(a * 3 + b) * 3 + c] === 0) {
          pen.polyline(
            [
              [x, y],
              [x + w, y],
              [x + w, y + h],
              [x, y + h],
              [x, y],
            ],
            4,
            ink.line,
          );
          continue;
        }
        kept++;
        const zeros = [a, b, c].filter((value) => value === 0).length;
        if (zeros <= 1) {
          near++;
          pen.rect(x, y, w, h, ink.green);
        } else {
          pen.rect(x, y, w, h, ink.dim);
        }
      }
    }
  }
  if (kept !== 20 || near !== 13)
    throw new Error(`paper-menger-pairwise-coprimality: kept ${kept} and near ${near}, want 20 and 13`);
}
