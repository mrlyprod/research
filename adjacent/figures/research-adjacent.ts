import { Grid, type Ink, type Pen } from "@mrlyprod/draw";

const B = 4;
const L = 6;
const SIDE = B ** (L / 2);

function adjacent(digits: number[]): boolean {
  if (digits[L - 1] === 0) return false;
  for (let i = 0; i + 1 < L; i++) if (digits[i] === digits[i + 1]) return false;
  return true;
}

function missing(digits: number[]): boolean {
  return digits.every((d) => d !== 0);
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const grid = new Grid(area, SIDE, SIDE, 0.14);
  let inR = 0;
  let inF = 0;
  let both = 0;
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      const digits: number[] = new Array(L);
      for (let k = 0; k < L / 2; k++) {
        digits[L - 1 - 2 * k] = Math.floor(x / B ** k) % B;
        digits[L - 2 - 2 * k] = Math.floor(y / B ** k) % B;
      }
      const r = adjacent(digits);
      const f = missing(digits);
      if (r) inR++;
      if (f) inF++;
      if (r && f) both++;
      const color = r && f ? ink.blue : r ? ink.orange : f ? ink.dim : null;
      if (color) grid.fill(pen, x, SIDE - 1 - y, color);
    }
  }
  const fill = B - 1;
  if (inR !== fill ** L) throw new Error(`research-adjacent: ${inR} in R, want ${fill ** L}`);
  if (inF !== fill ** L) throw new Error(`research-adjacent: ${inF} missing 0, want ${fill ** L}`);
  if (both !== fill * (fill - 1) ** (L - 1))
    throw new Error(`research-adjacent: ${both} shared, want ${fill * (fill - 1) ** (L - 1)}`);
}
