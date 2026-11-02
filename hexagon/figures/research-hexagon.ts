import type { Ink, Pen } from "@mrlyprod/draw";
import { hex } from "@mrlyprod/draw";
import census from "./research-hexagon.json" with { type: "json" };

export const units = {};

const NONE = 65535;

type Sampled = { shades: Float64Array; codes: Uint16Array; starts: number[] };

let memo: Sampled | undefined;

// FIELD

const unpack = (text: string) => new DataView(Uint8Array.from(atob(text), (c) => c.charCodeAt(0)).buffer);

function sampled(): Sampled {
  if (memo) return memo;
  const { mesh, reach } = census;
  const table = unpack(census.shades);
  const grid = unpack(census.field);
  if (!(reach > 0)) throw new Error(`research-hexagon: reach ${reach}, want above zero`);
  if (table.byteLength % 8 !== 0)
    throw new Error(`research-hexagon: ${table.byteLength} shade bytes, want a multiple of 8`);
  if (grid.byteLength !== 2 * hex.count(mesh))
    throw new Error(`research-hexagon: ${grid.byteLength} field bytes, want ${2 * hex.count(mesh)}`);
  const shades = new Float64Array(table.byteLength / 8);
  for (let i = 0; i < shades.length; i++) shades[i] = table.getFloat64(8 * i, true);
  const codes = new Uint16Array(grid.byteLength / 2);
  for (let i = 0; i < codes.length; i++) {
    codes[i] = grid.getUint16(2 * i, true);
    if (codes[i] !== NONE && codes[i] >= shades.length)
      throw new Error(`research-hexagon: code ${codes[i]} past ${shades.length} shades`);
  }
  const starts: number[] = [];
  let at = 0;
  for (let row = 0; row < 2 * mesh; row++) {
    starts.push(at);
    at += hex.row_len(mesh, row);
  }
  memo = { shades, codes, starts };
  return memo;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const { mesh, background, reach } = census;
  const { shades, codes, starts } = sampled();
  const ramp = new ink.Ramp([ink.blue, ink.panel, ink.yellow]);
  const frame = pen.frame(0.08);
  hex.hexagon(pen, frame, mesh, 0, (row, col) => {
    const code = codes[starts[row] + col];
    return code === NONE ? null : ramp.at(0.5 + (shades[code] - background) / (2 * reach));
  });
}
