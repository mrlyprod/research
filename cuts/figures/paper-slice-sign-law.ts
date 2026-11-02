import type { Ink, Pen } from "@mrlyprod/draw";

export const units = {};

const LOW = 6;
const HIGH = 40;
const TOP = 91;
const TERMS = 40;

// THE CENSUS EXCESS

function binomial(n: number, k: number) {
  let out = 1n;
  for (let i = 1; i <= k; i++) out = (out * BigInt(n - k + i)) / BigInt(i);
  return out;
}

function digits(dim: number) {
  const even = Array.from({ length: dim + 1 }, (_, k) => binomial(dim, k));
  const odd = Array.from({ length: dim }, (_, k) => BigInt(dim) * binomial(dim - 1, k));
  const out = Array<bigint>(2 * dim + 1).fill(0n);
  for (const [k, v] of even.entries()) out[2 * k] = v;
  for (const [k, v] of odd.entries()) out[2 * k + 1] = v;
  return out;
}

function excess(dim: number, top: number) {
  const P = digits(dim);
  const at = (s: number) => (s >= 0 && s < P.length ? P[s] : 0n);
  const half = (dim - 1) >> 1;
  const weight = (c: number, k: number) => (k === 0 ? at(c + dim) : at(c + dim - 3 * k) + at(c + dim + 3 * k));
  let row = Array.from({ length: half + 1 }, (_, c) => BigInt((c % 3 === 0 ? 3 : 0) - 1));
  const signs = [row[0] < 0n];
  let stop = 0;
  for (let level = 1; level <= top; level++) {
    row = row.map((_, c) => row.reduce((sum, r, k) => sum + weight(c, k) * r, 0n));
    signs.push(row[0] < 0n);
    if (!stop && row.every((r) => r > 0n)) stop = level;
  }
  return { signs, stop };
}

// THE CROSSING

function crossing(dim: number) {
  let logR = 0;
  let shift = 0;
  for (let i = 2; i < TERMS; i++) {
    const a = Math.cos(Math.PI / 3 ** i);
    const b = Math.cos((2 * Math.PI) / 3 ** i);
    logR += Math.log(a / b);
    shift += Math.log((dim - 2 * a) / (dim - 2)) - Math.log((dim + 2 * b) / (dim + 2));
  }
  return ((dim - 1) * logR + shift) / Math.log((dim + 2) / (dim - 2));
}

function rows() {
  const out = [];
  for (let dim = LOW; dim <= HIGH; dim += 2) {
    const { signs, stop } = excess(dim, TOP);
    const last = signs.lastIndexOf(true);
    const star = crossing(dim);
    const floor = Math.floor(star);
    const below = floor % 2 ? floor : floor - 1;
    if (last !== below) throw new Error(`paper-slice-sign-law: last negative level ${last} at ${dim}, want ${below}`);
    if (!stop || stop - last < 1 || stop - last > 2)
      throw new Error(`paper-slice-sign-law: stopping level ${stop} at ${dim}, want ${last} + 1 or + 2`);
    if (signs.some((neg, level) => neg !== (level % 2 === 1 && level <= last)))
      throw new Error(`paper-slice-sign-law: a broken run at ${dim}`);
    out.push({ dim, last, star });
  }
  const forty = out.find((one) => one.dim === 40);
  if (forty?.last !== 87 || Math.abs(forty.star - 88.7719) > 1e-3)
    throw new Error(`paper-slice-sign-law: ${forty?.last} and ${forty?.star} at 40, want 87 and 88.7719`);
  return out;
}

// DRAW

export default function draw(pen: Pen, ink: Ink) {
  const table = rows();
  const frame = pen.frame(0.06);
  const cols = (TOP + 1) / 2;
  const slot = frame.w / cols;
  const band = frame.h / table.length;
  const gap = slot * 0.16;
  const place = (level: number) => frame.x + ((level - 1) / 2 + 0.5) * slot;
  table.forEach(({ last, star }, i) => {
    const top = frame.y + i * band + band * 0.14;
    const height = band * 0.72;
    for (let level = 1; level <= last; level += 2)
      pen.rect(place(level) - slot / 2 + gap / 2, top, slot - gap, height, ink.blue);
    pen.rect(place(star) - 1.5, frame.y + i * band + band * 0.04, 3, band * 0.92, ink.yellow);
  });
}
