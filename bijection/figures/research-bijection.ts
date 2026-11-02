import type { Ink, Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const CENSUS = 16;
const STAMP_COLS = 4;
const STAMP_ROWS = 2;

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const universe = math.bang.bang(3);
  if (universe.total !== 256) throw new Error(`research-bijection: ${universe.total} codes, want 256`);
  if (universe.distinct() !== 22) throw new Error(`research-bijection: ${universe.distinct()} orbits, want 22`);
  if (math.bang.baseq.distinct_designs(2, 3) !== "22")
    throw new Error(`research-bijection: ${math.bang.baseq.distinct_designs(2, 3)} distinct designs, want 22`);
  if (math.bang.corners(3).length !== 8)
    throw new Error(`research-bijection: ${math.bang.corners(3).length} corners, want 8`);

  const pitch = area.w / CENSUS;
  const bit = (pitch * 0.84) / STAMP_COLS;
  const band = bit * (STAMP_ROWS + 1);
  const top = area.y + (area.h - band * CENSUS) / 2;
  const pad = bit * 0.1;
  let gold = 0;
  for (let code = 0; code < 256; code++) {
    const design = universe.design(code);
    let color = ink.blue;
    if (design.canonical) {
      gold++;
      color = ink.yellow;
    }
    const ox = area.x + (code % CENSUS) * pitch + (pitch - bit * 4) / 2;
    const oy = top + Math.floor(code / CENSUS) * band + bit / 2;
    for (let slot = 0; slot < 8; slot++) {
      if (((code >> slot) & 1) === 0) continue;
      const x = ox + (slot % STAMP_COLS) * bit;
      const y = oy + Math.floor(slot / STAMP_COLS) * bit;
      pen.rect(x + pad, y + pad, bit - 2 * pad, bit - 2 * pad, color);
    }
  }
  if (gold !== 22) throw new Error(`research-bijection: ${gold} canonical, want 22`);
}
