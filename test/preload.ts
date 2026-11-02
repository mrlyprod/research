import { existsSync, readFileSync } from "node:fs";

const UNITS = ["core", "font", "life", "math", "num", "tile"];

const wasm = (unit: string) =>
  new URL(`../node_modules/@mrlyprod/${unit}/pkg/mrly_${unit}_js_bg.wasm`, import.meta.url);

const lost = UNITS.filter((unit) => !existsSync(wasm(unit)));
if (lost.length)
  throw new Error(
    `test: no wasm under node_modules/@mrlyprod/{${lost.join(",")}}/pkg; build the siblings' js/pkg first`,
  );

for (const unit of UNITS) {
  const door = await import(`@mrlyprod/${unit}`);
  door.initSync({ module: readFileSync(wasm(unit)) });
}
