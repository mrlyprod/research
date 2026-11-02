import { Grid, type Ink, type Pen } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import * as num from "@mrlyprod/num";

export const units = { math, num };

const SIDE = 32;
const STRONG = "286747428434478734107892712789";
const SPONGE = "251620868451255";
const EXPONENTS = ["1.7268", "0.8340", "1.2618"];
const CENSUS = [0, 60, 1434, 32268, 721524, 15141288];
const GAPS = ["0.251621", "0.101621", "0.072371", "0.049946", "0.026145", "0.015038"];
const RATIOS = ["0.404", "0.712", "0.690", "0.523", "0.575"];
const PREDICTED = ["0.222232", "0.234019"];
const RATES: [number, string][] = [
  [2, "0.600"],
  [5, "0.447"],
  [7, "0.200"],
  [11, "0.374"],
  [13, "0.532"],
  [17, "0.338"],
  [19, "0.411"],
];

type Design = [string, number, number, number[][], [number, number], string];

let checked = false;

function fail(what: string): never {
  throw new Error(`paper-coprime-density-above-dimension-one: ${what}`);
}

function gcd(a: number, b: number) {
  while (b) [a, b] = [b, a % b];
  return a;
}

function popcount(word: number) {
  word -= (word >>> 1) & 0x55555555;
  word = (word & 0x33333333) + ((word >>> 2) & 0x33333333);
  return (((word + (word >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
}

function cube(base: number, dim: number) {
  let out: number[][] = [[]];
  for (let i = 0; i < dim; i++) out = out.flatMap((v) => Array.from({ length: base }, (_, d) => [...v, d]));
  return out;
}

function primes(n: number) {
  const out: number[] = [];
  for (let p = 2; p <= n; p++) if (n % p === 0 && out.every((q) => p % q !== 0)) out.push(p);
  return out;
}

function designs(): Design[] {
  const carpet = cube(3, 2).filter((v) => v.filter((d) => d === 1).length <= 1);
  const sponge = cube(3, 3).filter((v) => v.filter((d) => d === 1).length <= 1);
  return [
    [
      "gasket",
      2,
      2,
      [
        [0, 0],
        [1, 0],
        [0, 1],
      ],
      [2, 3],
      "0.540379646092",
    ],
    [
      "or-triangle",
      2,
      2,
      [
        [1, 0],
        [0, 1],
        [1, 1],
      ],
      [1, 1],
      "0.810569469138",
    ],
    ["carpet", 3, 2, carpet, [7, 8], "0.598428240887"],
    [
      "plus",
      3,
      2,
      [
        [1, 0],
        [0, 1],
        [1, 1],
        [2, 1],
        [1, 2],
      ],
      [1, 1],
      "0.683917989585",
    ],
    ["sponge", 3, 3, sponge, [19, 20], "0.820708619488"],
    [
      "base-6 sample",
      6,
      2,
      [
        [0, 0],
        [0, 3],
        [1, 1],
        [1, 2],
        [2, 0],
        [2, 1],
        [4, 0],
        [5, 5],
      ],
      [1, 2],
      "0.455945326390",
    ],
    [
      "base-6 pair",
      6,
      2,
      [
        [0, 0],
        [0, 1],
        [1, 0],
        [1, 1],
        [2, 3],
        [3, 2],
        [4, 5],
        [5, 4],
      ],
      [7, 8],
      "0.797904321183",
    ],
    [
      "dust",
      3,
      2,
      [
        [0, 0],
        [0, 2],
        [2, 0],
        [2, 2],
      ],
      [3, 4],
      "0.512938492189",
    ],
  ];
}

function brackets() {
  const zeta = (dim: number) => (dim === 2 ? num.series.BASEL() : num.series.APERY());
  for (const [name, base, dim, corners, [top, bottom], printed] of designs()) {
    const fill = corners.length;
    const good = corners.filter((v) => gcd(v.reduce(gcd, 0), base) === 1).length;
    const radical = primes(base);
    let mobius = 0;
    for (let s = 0; s < 1 << radical.length; s++) {
      const e = radical.filter((_, i) => (s >> i) & 1).reduce((a, p) => a * p, 1);
      const sign = radical.filter((_, i) => (s >> i) & 1).length % 2 ? -1 : 1;
      mobius += sign * corners.filter((v) => v.every((d) => d % e === 0)).length;
    }
    if (good !== mobius || good * bottom !== top * fill)
      fail(`${name} bracket ${good}/${fill} by gcd and ${mobius}/${fill} by Mobius, want ${top}/${bottom}`);
    const delta = radical.reduce((acc, p) => acc / (1 - p ** -dim), good / fill) / zeta(dim);
    const shown = (Math.floor(delta * 1e12) / 1e12).toFixed(12);
    if (shown !== printed) fail(`${name} density ${delta}, printed ${printed}`);
  }
  const sample = designs()[5][3];
  const counts = [2, 3, 6].map((e) => sample.filter((v) => v.every((d) => d % e === 0)).length);
  if (counts.join() !== "3,2,1" || (1 - 3 / 8) * (1 - 2 / 8) !== 0.46875) fail(`base-6 corner counts ${counts}`);
  let xs = [0];
  let ys = [0];
  for (let level = 1; level <= 6; level++) {
    const nx: number[] = [];
    const ny: number[] = [];
    for (const [a, b] of sample)
      for (let i = 0; i < xs.length; i++) {
        nx.push(a + 6 * xs[i]);
        ny.push(b + 6 * ys[i]);
      }
    [xs, ys] = [nx, ny];
    const prime = xs.filter((x, i) => gcd(gcd(x, ys[i]), 6) === 1).length;
    if (2 * prime !== 8 ** level) fail(`base-6 sample level ${level}: ${prime} of ${8 ** level} prime to 6`);
  }
}

function pairs(level: number) {
  const n = 3 ** level;
  const words = Math.ceil(n / 32);
  const masks = new Int32Array(n);
  for (let v = 0; v < n; v++)
    for (let j = 0, r = v; j < level; j++, r = Math.floor(r / 3)) if (r % 3 === 1) masks[v] |= 1 << j;
  const coprime = new Uint32Array(n * words);
  for (let x = 0; x < n; x++)
    for (let z = 0; z < n; z++) if (gcd(x, z) === 1) coprime[x * words + (z >> 5)] |= 1 << (z & 31);
  const open = new Uint32Array((1 << level) * words);
  const size = new Int32Array(1 << level);
  const odd = new Int32Array(1 << level);
  for (let u = 0; u < 1 << level; u++)
    for (let z = 0; z < n; z++)
      if (!(masks[z] & u)) {
        open[u * words + (z >> 5)] |= 1 << (z & 31);
        size[u]++;
        odd[u] += z & 1;
      }
  let points = 0;
  let pairwise = 0;
  let parity = 0;
  for (let x = 0; x < n; x++)
    for (let y = 0; y < n; y++) {
      if (masks[x] & masks[y]) continue;
      const u = masks[x] | masks[y];
      points += size[u];
      parity += x & y & 1 ? size[u] : (x ^ y) & 1 ? odd[u] : 0;
      if (!((coprime[x * words + (y >> 5)] >>> (y & 31)) & 1)) continue;
      for (let w = 0; w < words; w++)
        pairwise += popcount(coprime[x * words + w] & coprime[y * words + w] & open[u * words + w]);
    }
  return { points, pairwise, parity };
}

function sponge() {
  const digits = cube(3, 3);
  const kept = digits.filter((v) => v.filter((d) => d === 1).length <= 1);
  const zeros = (v: number[]) => v.filter((d) => d === 0).length;
  const fibre = (axes: number[]) =>
    Math.max(...cube(3, axes.length).map((w) => kept.filter((v) => axes.every((a, i) => v[a] === w[i])).length));
  const near = kept.filter((v) => zeros(v) <= 1).length;
  const full = digits.filter((v) => zeros(v) <= 1).length;
  const kappa = [
    [0, 1],
    [0, 2],
    [1, 2],
  ].map(fibre);
  const single = [[0], [1], [2]].map(fibre);
  if (kept.length !== 20 || near !== 13 || full !== 20 || kappa.join() !== "3,3,3" || single.join() !== "8,8,8")
    fail(`sponge digits ${kept.length}, ${near}, ${full}, fibres ${kappa} and ${single}`);
  const exponents = [20 / 3, 20 / 8, 4].map((x) => (Math.floor((Math.log(x) / Math.log(3)) * 1e4) / 1e4).toFixed(4));
  if (exponents.join() !== EXPONENTS.join()) fail(`exponents ${exponents}, printed ${EXPONENTS}`);
  const product = ((BigInt(STRONG) * 351n) / 400n).toString().padStart(30, "0");
  if (!product.startsWith(SPONGE)) fail(`351/400 of the strongly carefree digits ${product}, printed ${SPONGE}`);
  const limit = Number(`0.${product}`);
  let last = 0;
  for (let level = 1; level <= 6; level++) {
    const { points, pairwise, parity } = pairs(level);
    const exact = (2 * 20 ** level - 3 * 12 ** level + (-4) ** level) / 4;
    if (points !== 20 ** level || pairwise !== CENSUS[level - 1] || parity !== exact)
      fail(
        `sponge level ${level}: ${points} points, ${pairwise} pairwise coprime, ${parity} at most one even against ${exact}`,
      );
    const gap = limit - pairwise / 20 ** level;
    if (gap.toFixed(6) !== GAPS[level - 1]) fail(`sponge gap ${gap} at level ${level}`);
    if (level > 1 && (gap / last).toFixed(3) !== RATIOS[level - 2])
      fail(`sponge gap ratio ${gap / last} at level ${level}`);
    if (level >= 5) {
      const beta = 0.5 - 0.75 * 0.6 ** level + 0.25 * (-0.2) ** level;
      if ((2 * beta * limit).toFixed(6) !== PREDICTED[level - 5])
        fail(`parity prediction ${2 * beta * limit} at level ${level}`);
    }
    last = gap;
  }
  for (const [p, printed] of RATES) {
    let order = 1;
    while (3 ** order % p !== 1 % p) order++;
    const size = Array.from({ length: p ** 3 }, (_, t) => {
      const [a, b, c] = [t % p, Math.floor(t / p) % p, Math.floor(t / p ** 2)];
      let re = 0;
      let im = 0;
      for (const [x, y, z] of kept) {
        const angle = (2 * Math.PI * ((a * x + b * y + c * z) % p)) / p;
        re += Math.cos(angle);
        im += Math.sin(angle);
      }
      return Math.hypot(re, im) / 20;
    });
    let best = 0;
    for (let t = 1; t < p ** 3; t++) {
      let product = 1;
      for (let j = 0, u = t; j < order; j++) {
        product *= size[u];
        const [a, b, c] = [u % p, Math.floor(u / p) % p, Math.floor(u / p ** 2)];
        u = ((3 * a) % p) + p * ((3 * b) % p) + p * p * ((3 * c) % p);
      }
      best = Math.max(best, product ** (1 / order));
    }
    if (best.toFixed(3) !== printed || best > 0.6 + 1e-12 || (p !== 2 && best > 0.6 - 1e-3))
      fail(`character rate ${best} at ${p}, printed ${printed}`);
  }
}

function numbers() {
  if (checked) return;
  brackets();
  sponge();
  checked = true;
}

export default function draw(pen: Pen, ink: Ink) {
  numbers();
  const gasket = math.two.from_corners(
    [
      [0, 0],
      [1, 0],
      [0, 1],
    ],
    2,
    5,
    0,
    2,
  );
  if (gasket.shape.join() !== `${SIDE},${SIDE}`) fail(`gasket ${gasket.shape}, want 32,32`);
  const types = gasket.types;
  const grid = new Grid(pen.frame(0.08), SIDE, SIDE, 0.1);
  let filled = 0;
  let visible = 0;
  for (let y = 0; y < SIDE; y++) {
    for (let x = 0; x < SIDE; x++) {
      if (types[y * SIDE + x] === 0) continue;
      filled++;
      const coprime = gcd(x, y) === 1;
      if (coprime) visible++;
      grid.fill(pen, x, SIDE - 1 - y, coprime ? ink.yellow : ink.dim);
    }
  }
  if (filled !== 243 || visible !== 122) fail(`${filled} filled and ${visible} coprime, want 243 and 122`);
}
