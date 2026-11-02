import type { Color, Ink, Pen } from "@mrlyprod/draw";
import * as life from "@mrlyprod/life";

export const units = { life };

const CENSUS = 16;
const STEPS = 31;
const WINDOW = 63;
const CROP = 48;
const FAINT = 0.5;

type Rule = { degree: number; surjective: boolean; diagram: ArrayLike<number> };

let memo: Rule[] | undefined;

function rules() {
  if (memo) return memo;
  const found: Rule[] = [];
  for (let rule = 0; rule < 256; rule++) {
    const diagram = life.single_seed(rule, STEPS);
    if (diagram.shape.length !== 2 || diagram.shape[0] !== STEPS + 1 || diagram.shape[1] !== WINDOW)
      throw new Error(`research-automata: rule ${rule} diagram ${diagram.shape}, want ${STEPS + 1},${WINDOW}`);
    found.push({ degree: life.rule_degree(rule), surjective: life.surjective(rule), diagram: diagram.data });
  }
  const affine = found.filter((r) => r.degree === 1).length;
  const surjective = found.filter((r) => r.surjective).length;
  if (affine !== 14) throw new Error(`research-automata: ${affine} affine rules, want 14`);
  if (surjective !== 30) throw new Error(`research-automata: ${surjective} surjective rules, want 30`);
  if (surjective - affine !== 16)
    throw new Error(`research-automata: ${surjective - affine} surjective rules beyond affine, want 16`);
  memo = found;
  return memo;
}

export default function draw(pen: Pen, ink: Ink) {
  const area = pen.frame(0.08);
  const rows = STEPS + 1;
  const pitch = area.w / CENSUS;
  const rise = (pitch - rows) / 2;
  const left = Math.floor((WINDOW - CROP) / 2);
  const faint = ink.fade(ink.dim, FAINT);
  rules().forEach((rule, n) => {
    const color: Color = rule.degree === 1 ? ink.yellow : rule.surjective ? ink.blue : faint;
    const ox = Math.round(area.x + (n % CENSUS) * pitch);
    const oy = Math.round(area.y + Math.floor(n / CENSUS) * pitch + rise);
    const colors = new Uint8Array(rows * CROP * 4);
    for (let t = 0; t < rows; t++) {
      for (let c = 0; c < CROP; c++) {
        if (rule.diagram[t * WINDOW + left + c] !== 0) colors.set(color, (t * CROP + c) * 4);
      }
    }
    pen.image(ox, oy, CROP, rows, { shape: [rows, CROP], colors });
  });
}
