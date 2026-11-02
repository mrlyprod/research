import { afterAll, expect, test } from "bun:test";
import { existsSync, readdirSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, relative } from "node:path";

const REPO = join(import.meta.dir, "..");

const NAME = /^[a-z0-9]+(-[a-z0-9]+)*$/;

const out = join(tmpdir(), `figures-${process.pid}`);

afterAll(() => rmSync(out, { recursive: true, force: true }));

type Figure = { name: string; file: string };

function roster(): Figure[] {
  const found: Figure[] = [];
  for (const topic of readdirSync(REPO, { withFileTypes: true })) {
    const dir = join(REPO, topic.name, "figures");
    if (!topic.isDirectory() || topic.name.startsWith(".") || !existsSync(dir)) continue;
    for (const name of readdirSync(dir)
      .filter((one) => one.endsWith(".ts"))
      .sort())
      found.push({ name: name.slice(0, -3), file: join(dir, name) });
  }
  const seen = new Map<string, string>();
  for (const one of found) {
    const rel = relative(REPO, one.file);
    if (!NAME.test(one.name))
      throw new Error(
        `${rel} draws the figure ${one.name}, and a figure name is lowercase words and digits joined by hyphens`,
      );
    if (seen.has(one.name)) throw new Error(`two figures take the name ${one.name}: ${seen.get(one.name)} and ${rel}`);
    seen.set(one.name, rel);
  }
  return found.sort((a, b) => a.name.localeCompare(b.name));
}

test("every figure loads: default is a draw function and every unit it names is defined", async () => {
  const list = roster();
  const built = await Bun.build({
    entrypoints: list.map((one) => one.file),
    root: REPO,
    outdir: out,
    target: "bun",
  });
  expect(built.logs.filter((log) => log.level === "error").map(String)).toEqual([]);
  const bad: string[] = [];
  for (const one of list) {
    try {
      const figure = await import(join(out, relative(REPO, one.file).replace(/\.ts$/, ".js")));
      if (typeof figure.default !== "function") bad.push(`${one.name}: no default draw`);
      for (const [unit, value] of Object.entries(figure.units ?? {}))
        if (value === undefined) bad.push(`${one.name}: units.${unit} is undefined`);
    } catch (error) {
      bad.push(`${one.name}: ${error instanceof Error ? error.message : String(error)}`);
    }
  }
  expect(list.length).toBeGreaterThan(0);
  expect(bad).toEqual([]);
}, 60_000);
