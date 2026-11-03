import { type Ink, type Pen, type Point, plot } from "@mrlyprod/draw";
import * as math from "@mrlyprod/math";
import census from "./research-dimensions.json" with { type: "json" };

const RE_LO = -1;
const RE_HI = 1;
const IM_REACH = 40;
const CONTROL = 21;
const LATTICE = 13;

export const units = { math };

function lattice(): { real: number; poles: Point[] } {
  const design = new math.ladder.Design(3, [0, 1]);
  const real = design.abscissa();
  const poles: Point[] = [];
  for (let m = -6; m <= 6; m++) {
    const pole = design.pole(0, m);
    poles.push([pole.re, pole.im]);
  }
  design.free();
  if (poles.length !== LATTICE) throw new Error(`research-dimensions: ${poles.length} lattice poles, want ${LATTICE}`);
  return { real, poles };
}

function control(): Point[] {
  const poles = census.control as Point[];
  if (poles.length !== CONTROL) throw new Error(`research-dimensions: ${poles.length} control poles, want ${CONTROL}`);
  if (poles.some((p) => p.length !== 2)) throw new Error("research-dimensions: a control pole that is not a pair");
  return poles;
}

export default function draw(pen: Pen, ink: Ink) {
  const { real, poles } = lattice();
  const frame = pen.frame(0.08);
  plot.axis(pen, frame, ink.line);
  const at = (re: number, im: number): Point => [
    frame.x + (frame.w * (re - RE_LO)) / (RE_HI - RE_LO),
    frame.y + frame.h * (1 - (im + IM_REACH) / (2 * IM_REACH)),
  ];
  pen.segment(at(0, -IM_REACH), at(0, IM_REACH), 1.6, ink.line);
  pen.segment(at(RE_LO, 0), at(RE_HI, 0), 1.6, ink.line);
  pen.segment(at(real, -IM_REACH), at(real, IM_REACH), 1.8, ink.fade(ink.blue, 0.35));
  for (const [re, im] of control()) pen.disc(...at(re, im), 7, ink.orange);
  for (const [re, im] of poles) pen.disc(...at(re, im), 9, ink.blue);
}
