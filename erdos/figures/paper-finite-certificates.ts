import { Grid, type Ink, type Pen } from "@mrlyprod/draw";

export const units = {};

const SET = [302, 409, 447, 459, 465, 466, 474];
const TOP = 474;
const SIDE = 78;
const FAINT = 0.25;
const NAME = "paper-finite-certificates";

type Board = { pad: number; span: number; ternary: Uint8Array; subset: Uint8Array };

let memo: Board | undefined;

function board() {
  if (memo) return memo;
  if (SET.length !== 7 || Math.max(...SET) !== TOP || SET.some((a, i) => i > 0 && a <= SET[i - 1]))
    throw new Error(`${NAME}: the set is not seven increasing elements up to ${TOP}`);
  const total = SET.reduce((s, a) => s + a, 0);
  const span = 2 * total + 1;
  const ternary = new Uint8Array(span);
  let sums = [0];
  for (const a of SET) sums = sums.flatMap((s) => [s, s + a, s + 2 * a]);
  for (const s of sums) {
    if (ternary[s]) throw new Error(`${NAME}: the ternary sum ${s} repeats`);
    ternary[s] = 1;
  }
  if (sums.length !== 3 ** SET.length) throw new Error(`${NAME}: ${sums.length} ternary sums, want ${3 ** SET.length}`);
  const subset = new Uint8Array(span);
  let plain = [0];
  for (const a of SET) plain = plain.flatMap((s) => [s, s + a]);
  for (const s of plain) subset[s] = 1;
  const listed = [...subset.keys()].filter((s) => subset[s]);
  if (listed.length !== 2 ** SET.length)
    throw new Error(`${NAME}: ${listed.length} subset sums, want ${2 ** SET.length}`);
  for (let i = 0; i < listed.length; i++)
    for (let j = i + 1; j < listed.length; j++) {
      const twice = listed[i] + listed[j];
      if (twice % 2 === 0 && subset[twice / 2])
        throw new Error(`${NAME}: ${listed[i]}, ${twice / 2}, ${listed[j]} is a 3-term progression`);
    }
  for (const s of listed) if (!ternary[s]) throw new Error(`${NAME}: the subset sum ${s} is no ternary sum`);
  const pad = Math.floor((SIDE * SIDE - span) / 2);
  if (pad < 0 || SIDE * SIDE - span >= SIDE) throw new Error(`${NAME}: ${span} integers do not fill a ${SIDE}-square`);
  if (sums.some((s) => s + pad < 0 || s + pad >= SIDE * SIDE))
    throw new Error(`${NAME}: a ternary sum falls outside the ${SIDE}-square`);
  memo = { pad, span, ternary, subset };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { pad, span, ternary, subset } = board();
  const grid = new Grid(pen.frame(0.08), SIDE, SIDE, 0.14);
  const faint = ink.mix(ink.ground, ink.dim, FAINT);
  for (let at = 0; at < SIDE * SIDE; at++) {
    const x = at - pad;
    const inside = x >= 0 && x < span;
    const color = inside && subset[x] ? ink.orange : inside && ternary[x] ? ink.blue : faint;
    grid.fill(pen, at % SIDE, SIDE - 1 - Math.floor(at / SIDE), color);
  }
}
