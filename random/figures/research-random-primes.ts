import type { Ink, Pen, Point } from "@mrlyprod/draw";

const NAME = "research-random-primes";
const TOP = 22;
const LOW = 14;
const STEPS = 64;
const THETA = Math.log2(3) / 2;
const PRIMES = 295947;
const UNITS = -177147;

let memo: { us: number[]; primes: number[]; model: number[]; three: number[] } | undefined;

function parity(n: number) {
  let p = 0;
  for (let x = n; x > 0; x &= x - 1) p ^= 1;
  return p;
}

function series() {
  if (memo) return memo;
  const size = 2 ** TOP;
  const composite = new Uint8Array(size >> 1);
  for (let p = 3; p * p < size; p += 2) {
    if (composite[p >> 1]) continue;
    for (let m = p * p; m < size; m += 2 * p) composite[m >> 1] = 1;
  }
  const isPrime = (n: number) => n === 2 || (n > 2 && n % 2 === 1 && !composite[n >> 1]);
  const powers: [number, number][] = [];
  for (let p = 2; p * p < size; p++) {
    if (!isPrime(p)) continue;
    for (let m = p * p; m < size; m *= p) powers.push([m, Math.log(p)]);
  }
  powers.sort((a, b) => a[0] - b[0]);
  const us = Array.from({ length: (TOP - LOW) * STEPS + 1 }, (_, i) => LOW + i / STEPS);
  const marks = us.map((u) => Math.round(2 ** u));
  const primes: number[] = [];
  const model: number[] = [];
  const three: number[] = [];
  let [psi, m30, m6, c3, count, next, k] = [0, 0, 0, 2, 0, 0, 0];
  for (let n = 1; n <= size; n++) {
    while (next < marks.length && marks[next] === n) {
      const scale = n ** THETA;
      primes.push(-psi / scale);
      model.push((-15 / 4) * (m30 / scale));
      three.push(c3 / 2 / scale);
      next++;
    }
    if (n === size) break;
    const t = 1 - 2 * parity(n);
    if (isPrime(n)) {
      psi += t * Math.log(n);
      count++;
    }
    while (k < powers.length && powers[k][0] === n) psi += t * powers[k++][1];
    if (n % 2 && n % 3) {
      m6 += t;
      if (n % 5) m30 += t;
    }
    c3 += n % 3 ? -t : 2 * t;
  }
  if (count !== PRIMES) throw new Error(`${NAME}: pi(2^${TOP}) is ${count}, want ${PRIMES}`);
  if (3 * m6 !== UNITS) throw new Error(`${NAME}: the q | 6 part is ${3 * m6}, want ${UNITS}`);
  if (-c3 / 2 !== UNITS) throw new Error(`${NAME}: T_3 is ${-c3 / 2}, want ${UNITS}`);
  memo = { us, primes, model, three };
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const { us, primes, model, three } = series();
  const area = pen.frame(0.07);
  const all = [...model, ...three, ...primes];
  const lo = Math.min(...all) - 0.04;
  const hi = Math.max(...all) + 0.04;
  const at = (u: number, g: number): Point => [
    area.x + (area.w * (u - LOW)) / (TOP - LOW),
    area.y + area.h * (1 - (g - lo) / (hi - lo)),
  ];
  const line = area.w / 260;
  pen.polyline(
    us.map((u, i) => at(u, three[i])),
    line * 0.6,
    ink.dim,
  );
  pen.polyline(
    us.map((u, i) => at(u, model[i])),
    line,
    ink.blue,
  );
  for (let i = 0; i < us.length; i++) {
    const [x, y] = at(us[i], primes[i]);
    pen.disc(x, y, line * 0.9, ink.orange);
  }
}
