import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

export const units = { math };

const CODES = 15;
const LENGTH = 32;
const SLACK = 0.25;

type Verdict = { row: number; col: number; meets: boolean; constant: boolean };

let memo: Verdict[] | undefined;

function alternating(a: number, b: number) {
  return Array.from({ length: LENGTH }, (_, place) =>
    math.bang.MagicLayer.new(new math.name.Bang(place % 2 === 0 ? a : b, 2, 2), 2),
  );
}

function verdict(a: number, b: number): [boolean, boolean] {
  const letters = alternating(a, b);
  const rates = math.bang.word.rates(letters);
  const [component, fill] = rates[rates.length - 1];
  const constant = math.bang.word.constant_functional(letters);
  return [Math.abs(component - fill) < SLACK, Math.abs(component - constant) < SLACK];
}

function verdicts() {
  if (memo) return memo;
  memo = [];
  for (let row = 0; row < CODES; row++) {
    for (let col = row + 1; col < CODES; col++) {
      const [meets, constant] = verdict(row + 1, col + 1);
      memo.push({ row, col, meets, constant });
    }
  }
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const grid = new Grid(pen.frame(0.08), CODES, CODES, 0.12);
  let alphabets = 0;
  let ceiling = 0;
  let exact = 0;
  for (const { row, col, meets, constant } of verdicts()) {
    alphabets++;
    const [x, y, w, h] = grid.cell(col, row);
    if (meets) {
      ceiling++;
      pen.rect(x, y, w, h, ink.green);
    } else {
      pen.rect(x, y, w, h, ink.orange);
    }
    if (constant) {
      exact++;
      const pip = w * 0.34;
      pen.rect(x + (w - pip) / 2, y + (h - pip) / 2, pip, pip, ink.yellow);
    }
  }
  if (alphabets !== 105 || ceiling !== 89 || exact !== 27)
    throw new Error(
      `paper-component-exponent-of-kronecker-words: ${alphabets} alphabets, ${ceiling} at the ceiling, ${exact} exact, want 105, 89, 27`,
    );
}
