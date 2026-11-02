import { type Color, type Frame, frame, type Ink, type Pen } from "@mrlyprod/draw";
import census from "./research-multilayers.json" with { type: "json" };

const LEVELS = 12;
const BINS = 128;

export const units = {};

// PANEL

function panel(pen: Pen, ink: Ink, box: Frame, grid: number[][], hue: Color) {
  pen.rect(box.x, box.y, box.w, box.h, ink.line);
  const gap = 3;
  box.rows(LEVELS).forEach((row, k) => {
    const w = row.w / BINS;
    grid[k].forEach((open, i) => {
      const c = ink.mix(ink.panel, hue, open);
      pen.rect(row.x + i * w, row.y + gap / 2, w + 0.5, row.h - gap, c);
    });
  });
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { light, fading } = census;
  for (const [name, grid] of [
    ["light", light],
    ["fading", fading],
  ] as const) {
    if (
      grid.length !== LEVELS ||
      grid.some((row) => row.length !== BINS || row.some((open) => !(open >= 0 && open <= 1 + 1e-9)))
    )
      throw new Error(`research-multilayers: ${name} is not ${LEVELS} rows of ${BINS} shares of 0 to 1`);
  }
  const area = pen.frame(0.08);
  const side = area.w * 0.47;
  panel(pen, ink, frame(area.x, area.y, side, side), light, ink.blue);
  panel(pen, ink, frame(area.x + area.w - side, area.y + area.h - side, side, side), fading, ink.orange);
}
