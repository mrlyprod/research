import { Grid, type Ink, type Pen } from "@mrlyprod/draw";

const P = 11;
const SPAN = Math.floor((2 * P - 1) / 3) + 1;
const SIDE = SPAN * SPAN;

function kummer(digits: number[]): [boolean, number] {
  let state = 0;
  for (const d of digits) {
    if (3 * d + state <= P - 1) state = 0;
    else if (2 * d + state >= P && 3 * d + state <= 2 * P - 1) state = 1;
    else return [false, state];
  }
  return [true, state];
}

function rowZero(digits: number[]): boolean {
  return digits.every((d) => d <= Math.floor((P - 1) / 3) || (d >= (P + 1) / 2 && d <= Math.floor((2 * P - 1) / 3)));
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const grid = new Grid(area, SIDE, SIDE, 0.14);
  let inK = 0;
  let inF = 0;
  let endOne = 0;
  let both = 0;
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      const digits = [x % SPAN, y % SPAN, Math.floor(x / SPAN), Math.floor(y / SPAN)];
      const [k, state] = kummer(digits);
      const f = rowZero(digits);
      if (k) inK++;
      if (f) inF++;
      if (k && state === 1) endOne++;
      if (k && f) both++;
      const color = k && f ? ink.blue : k ? ink.orange : f ? ink.dim : null;
      if (color) grid.fill(pen, x, SIDE - 1 - y, color);
    }
  }
  const fill = (P + 1) / 2;
  if (inK !== fill ** 4) throw new Error(`research-kummer: ${inK} in K, want ${fill ** 4}`);
  if (inF !== fill ** 4) throw new Error(`research-kummer: ${inF} in F*, want ${fill ** 4}`);
  if (endOne !== fill ** 4 / 3) throw new Error(`research-kummer: ${endOne} end in state 1, want ${fill ** 4 / 3}`);
  if (!(both > 0 && both < inK)) throw new Error(`research-kummer: ${both} shared, want some but not all`);
}
