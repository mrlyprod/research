import type { Ink, Pen } from "@mrlyprod/draw";
import census from "./research-automata-encyclopedia.json" with { type: "json" };

export const units = {};

const FAINT = 0.35;
const GAP = 0.12;

type Census = { masks: string[]; rules: string[]; grid: number[][][]; counted: number };

export default function draw(pen: Pen, ink: Ink) {
  const { masks, rules, grid, counted } = census as Census;
  if (grid.length !== rules.length || grid.length === 0)
    throw new Error(`research-automata-encyclopedia: ${grid.length} grid rows for ${rules.length} rules`);
  let seen = 0;
  for (const row of grid) {
    if (row.length !== masks.length)
      throw new Error(`research-automata-encyclopedia: ${row.length} cells for ${masks.length} masks`);
    for (const cell of row) seen += cell[0] + cell[1] + cell[2];
  }
  if (seen !== counted) throw new Error(`research-automata-encyclopedia: ${seen} counted, census says ${counted}`);
  const area = pen.frame(0.08);
  const pitch_x = area.w / masks.length;
  const pitch_y = area.h / rules.length;
  const gap_x = pitch_x * GAP;
  const gap_y = pitch_y * GAP;
  const faint = ink.fade(ink.dim, FAINT);
  grid.forEach((row, r) => {
    row.forEach((cell, c) => {
      const total = cell[0] + cell[1] + cell[2];
      if (total === 0) return;
      const x = area.x + c * pitch_x + gap_x / 2;
      const y = area.y + r * pitch_y + gap_y / 2;
      const w = pitch_x - gap_x;
      const h = pitch_y - gap_y;
      const parts = [
        [cell[2] / total, ink.orange],
        [cell[1] / total, ink.blue],
        [cell[0] / total, faint],
      ] as const;
      let at = x;
      for (const [share, color] of parts) {
        const width = w * share;
        if (width > 0) pen.rect(at, y, width, h, color);
        at += width;
      }
    });
  });
}
