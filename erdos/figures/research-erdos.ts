import type { Color, Ink, Pen } from "@mrlyprod/draw";
import census from "./research-erdos.json" with { type: "json" };

const TERMS = 29;
const WIDTH = 860;
const PITCH = 29;
const GAP = 3;
const FLOOR = 0.15;
const REACH = 0.6;

export const units = {};

export default function draw(pen: Pen, ink: Ink) {
  const { counts, widths } = census;
  if (counts.length !== TERMS || widths.length !== TERMS)
    throw new Error(`research-erdos: ${counts.length} rows of counts and ${widths.length} of widths, want ${TERMS}`);
  counts.forEach((half, n) => {
    if (half.length !== WIDTH / 2 || widths[n].length !== WIDTH / 2)
      throw new Error(
        `research-erdos: row ${n} holds ${half.length} counts and ${widths[n].length} widths, want ${WIDTH / 2}`,
      );
    if (half.some((hit, c) => hit < 0 || hit > widths[n][c] || widths[n][c] < 1))
      throw new Error(`research-erdos: row ${n} holds a count past its width`);
  });
  const tone = (share: number): Color | null =>
    share >= 1 ? ink.blue : share > 0 ? ink.mix(ink.ground, ink.dim, FLOOR + REACH * share) : null;
  const same = (a: Color | null, b: Color | null) =>
    a === b || (a !== null && b !== null && a.every((v, i) => v === b[i]));
  const x = Math.trunc((pen.width - WIDTH) / 2);
  const top = Math.trunc((pen.height - (TERMS * PITCH - GAP)) / 2);
  counts.forEach((half, n) => {
    const left = half.map((hit, c) => hit / widths[n][c]);
    const row = [...left, ...left.toReversed()];
    const tones = row.map(tone);
    let start = 0;
    for (let c = 1; c <= row.length; c++) {
      if (c < row.length && same(tones[c], tones[start])) continue;
      const color = tones[start];
      if (color) pen.rect(x + start, top + n * PITCH, c - start, PITCH - GAP, color);
      start = c;
    }
  });
}
