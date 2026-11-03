import { frame as box, type Color, type Frame, type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";

const RE_WIDTH = 1;
const IM_REACH = 60;
const RHO = [0.72079, 28.60568];
const BOX_RE = [0.72074, 0.72084];
const BOX_IM = [28.60563, 28.60573];
const BANDS = [
  [22.01, 24.01],
  [56, 58],
];
const BAND_RE = 3.02;
const OVER = 22;

export const units = { math };

function comb(design: math.ladder.Design): Point[] {
  const poles: Point[] = [];
  for (let j = 0; ; j++) {
    const pole = design.pole(0, j);
    const point: Point = [pole.re, pole.im];
    if (point[1] > IM_REACH) return poles;
    poles.push(point);
  }
}

function check(alpha: number, period: number, poles: Point[]) {
  if (Math.floor(IM_REACH / period) !== 10)
    throw new Error(`paper-design-dirichlet-inverse: ${Math.floor(IM_REACH / period)} periods in reach, want 10`);
  if (poles.length !== 11) throw new Error(`paper-design-dirichlet-inverse: ${poles.length} poles, want 11`);
  if (BANDS.length !== 2) throw new Error(`paper-design-dirichlet-inverse: ${BANDS.length} bands, want 2`);
  if (!poles.every((p) => p[0] === alpha)) throw new Error("paper-design-dirichlet-inverse: a pole off the abscissa");
  if (!(RHO[0] > alpha && RHO[0] < alpha + RE_WIDTH))
    throw new Error(`paper-design-dirichlet-inverse: rho ${RHO[0]} off the strip`);
  if (!(RHO[0] > BOX_RE[0] && RHO[0] < BOX_RE[1]))
    throw new Error(`paper-design-dirichlet-inverse: rho ${RHO[0]} outside its box`);
  if (!(RHO[1] > BOX_IM[0] && RHO[1] < BOX_IM[1]))
    throw new Error(`paper-design-dirichlet-inverse: rho ${RHO[1]} outside its box`);
  if (!BANDS.every((b) => b[1] > b[0] && b[0] > 0 && b[1] < IM_REACH))
    throw new Error("paper-design-dirichlet-inverse: a band off the reach");
  if (!(BAND_RE > RE_WIDTH)) throw new Error(`paper-design-dirichlet-inverse: band reach ${BAND_RE} under ${RE_WIDTH}`);
  if (!(Math.abs(alpha - 0.6309297536) < 1e-9))
    throw new Error(`paper-design-dirichlet-inverse: abscissa ${alpha}, want 0.6309297536`);
  if (!(Math.abs(period - 5.7192017348) < 1e-9))
    throw new Error(`paper-design-dirichlet-inverse: period ${period}, want 5.7192017348`);
}

function at(rect: Frame, abscissa: number, re: number, im: number): Point {
  return [rect.x + (rect.w * (re - abscissa)) / RE_WIDTH, rect.y + rect.h * (1 - im / IM_REACH)];
}

function strokeBox(pen: Pen, rect: Frame, thick: number, color: Color) {
  const pts: Point[] = [
    [rect.x, rect.y],
    [rect.x + rect.w, rect.y],
    [rect.x + rect.w, rect.y + rect.h],
    [rect.x, rect.y + rect.h],
    [rect.x, rect.y],
  ];
  pen.polyline(pts, thick, color);
}

function halfPlane(pen: Pen, ink: Ink, rect: Frame) {
  pen.rect(rect.x, rect.y, rect.w, rect.h, ink.panel);
  pen.segment([rect.x - OVER, rect.y + rect.h], [rect.x + rect.w, rect.y + rect.h], 1.6, ink.line);
  pen.segment([rect.x, rect.y - OVER], [rect.x, rect.y + rect.h + OVER], 2.8, ink.fade(ink.indigo, 0.5));
}

function band(pen: Pen, ink: Ink, rect: Frame, abscissa: number, span: number[]) {
  const low = at(rect, abscissa, abscissa, span[0])[1];
  const high = at(rect, abscissa, abscissa, span[1])[1];
  const run = rect.w + OVER;
  pen.rect(rect.x, high, run, low - high, ink.fade(ink.yellow, 0.16));
  pen.segment([rect.x, high], [rect.x + run, high], 1.8, ink.fade(ink.yellow, 0.5));
  pen.segment([rect.x, low], [rect.x + run, low], 1.8, ink.fade(ink.yellow, 0.5));
}

export default function draw(pen: Pen, ink: Ink) {
  const ladder = new math.ladder.Design(3, [0, 1]);
  const alpha = ladder.abscissa();
  const period = ladder.period();
  const poles = comb(ladder);
  ladder.free();
  check(alpha, period, poles);

  const work = pen.frame(0.08).inset(OVER + 2);
  const pw = work.w * 0.455;
  const ph = work.h * 0.865;
  const design = box(work.x, work.y, pw, ph);
  const control = box(work.x + work.w - pw, work.y + work.h - ph, pw, ph);

  halfPlane(pen, ink, control);
  plot.dots(pen, [at(control, 1, 1, 0)], 6, ink.fade(ink.indigo, 0.92));

  halfPlane(pen, ink, design);
  for (const span of BANDS) band(pen, ink, design, alpha, span);
  plot.dots(
    pen,
    poles.map((p) => at(design, alpha, p[0], p[1])),
    6,
    ink.fade(ink.indigo, 0.92),
  );

  const emblem = at(design, alpha, RHO[0], RHO[1]);
  const half = design.w * 0.057;
  strokeBox(pen, box(emblem[0] - half, emblem[1] - half, 2 * half, 2 * half), 2, ink.fade(ink.yellow, 0.85));
  plot.dots(pen, [emblem], 12, ink.yellow);
}
