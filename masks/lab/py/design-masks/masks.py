import argparse
import gc
import glob
import itertools
import json
import math
import os
import resource
import subprocess
import time

os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.5"
os.environ["PYTORCH_MPS_LOW_WATERMARK_RATIO"] = "0.4"

import numpy as np
import torch
import torch.nn.functional as F
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
from scipy.stats import t as tdist

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
BIN = os.path.join(REPO, "target", "release")
FIGURE = os.path.join(REPO, "research", "masks", "figures", "research-masks.json")
S = 81
FAMILIES = ["Fractal", "Magic", "Special", "Mosaic", "General"]
SEEDS = [1, 2, 3]
DIGITS = [127, 239]
LOW, HIGH = 0.03, 0.97
LAGS = [1, 2, 3, 5, 8, 13, 21, 34, 55]
VOCAB, WIDTH, HEADS, LAYERS = 32, 64, 4, 2
PAIRS_KV, QUERY = 20, 40

# MASKS


def mrly(fn, *args):
    out = subprocess.run(
        [os.path.join(BIN, "mrly-" + fn.split(".")[0]), fn, json.dumps(list(args))],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(out.stdout)


def grid(cell):
    t = cell["types"]
    return np.array(t["data"]["U8"], dtype=np.uint8).reshape(t["shape"])


def fit(a):
    a = np.tile(a, (-(-S // a.shape[0]), -(-S // a.shape[1])))
    top, left = (a.shape[0] - S) // 2, (a.shape[1] - S) // 2
    return a[top : top + S, left : left + S]


def draw():
    out = []
    for family in FAMILIES:
        for seed in SEEDS:
            config = {
                "groups": [family],
                "catalog": "Universe",
                "min_size": S,
                "max_size": S,
                "parity": "Odds",
                "invert": None,
            }
            recipe = mrly("tile.build.create_2d", config, seed)
            out.append(
                {
                    "name": f"{family.lower()}-{seed}",
                    "family": family,
                    "recipe": recipe,
                    "shape": [S, S],
                    "mask": grid(mrly("tile.build.build_2d", recipe)["cell"]),
                }
            )
    for seed in SEEDS:
        hextile = mrly("tile.build.random_tile_6d", S, seed)
        cell = grid(mrly("tile.build.build_6d", hextile)["cell"]["cell"])
        out.append(
            {
                "name": f"hex-{seed}",
                "family": "Hex",
                "recipe": hextile,
                "shape": list(cell.shape),
                "mask": fit(((cell == 1) | (cell >= 3)).astype(np.uint8)),
            }
        )
    for code in DIGITS:
        recipe = {"code": code, "number": 3, "level": 4, "rotation": 0, "base": 3}
        out.append(
            {
                "name": f"digit-{code}",
                "family": "Digit",
                "recipe": recipe,
                "shape": [S, S],
                "mask": grid(mrly("math.two.create", code, 3, 4, 0, 3)["cell"]),
            }
        )
    for d in out:
        assert d["mask"].shape == (S, S) and set(np.unique(d["mask"])) <= {0, 1}
    return out


def masks():
    path = os.path.join(DATA_DIR, "masks.json")
    if not os.path.exists(path):
        os.makedirs(DATA_DIR, exist_ok=True)
        rows = [
            {**d, "mask": ["".join(map(str, r)) for r in d["mask"]]} for d in draw()
        ]
        with open(path, "w") as f:
            json.dump(rows, f)
    with open(path) as f:
        rows = json.load(f)
    for d in rows:
        d["mask"] = np.array([[c == "1" for c in r] for r in d["mask"]])
        d["density"] = float(d["mask"].mean())
        d["kept"] = LOW <= d["density"] <= HIGH
    return rows


# STRUCTURE

EYE = np.eye(S, dtype=bool)
LOWER = np.tril(np.ones((S, S), dtype=bool))
PAIRS = int(LOWER.sum())


def allowed(m):
    return LOWER & (m | EYE)


def hops(a, k):
    r = a.copy()
    for _ in range(k - 1):
        r = (r.astype(np.int64) @ a.astype(np.int64)) > 0
    return r


def fiedler(a):
    a = (a | a.T) & ~EYE
    lap = np.diag(a.sum(1)) - a.astype(float)
    return max(0.0, float(np.sort(np.linalg.eigvalsh(lap))[1]))


def band(w):
    i, j = np.indices((S, S))
    return (j <= i) & (i - j <= w)


def window(count):
    return next(w for w in range(S) if band(w).sum() >= count)


def twin(count, index, seed):
    rng = np.random.default_rng([index, seed])
    lower = np.flatnonzero(np.tril(np.ones((S, S), dtype=bool), -1))
    m = EYE.copy()
    m.flat[rng.choice(lower, count - S, replace=False)] = True
    return m


def describe(m):
    c = allowed(m)
    return {
        "effective": c.sum() / PAIRS,
        "reach1": (m | EYE).mean(),
        "reach2": hops(m | EYE, 2).mean(),
        "causal2": hops(c, 2)[LOWER].mean(),
        "fiedler": fiedler(m),
        "cfiedler": fiedler(c),
    }


def short(recipe):
    if "tile" in recipe:
        t = recipe["tile"]
        return f"{recipe['projection']} {t['group']} {t['sources'][0]['design']} n{t['numbers'][0]} L{t['levels'][0]} r{t['rotations'][0]}{' inv' if t['invert'] else ''} {t['width']}x{t['height']}"
    if "code" in recipe:
        return f"code {recipe['code']} base {recipe['base']} L{recipe['level']}"
    codes = ",".join(s["code"] for s in recipe["sources"])
    return f"codes {codes} n{','.join(map(str, recipe['numbers']))} L{','.join(map(str, recipe['levels']))} r{','.join(map(str, recipe['rotations']))} f{recipe['factor']}{' inv' if recipe['invert'] else ''}{' flip' if recipe['flip'] else ''}"


def census(args):
    print(
        f"side {S}; on = 1 in the 2d cell; hex on = fill or a cube face, tiled and center-cropped to {S} x {S}"
    )
    print(
        f"{'mask':<12} {'recipe':<44} {'dens':>6} {'eff':>6} {'r1':>6} {'r2':>6} {'c2':>6} {'l2':>7} {'cl2':>7} | {'w':>3} {'wdens':>6} {'wc2':>6} {'wl2':>7} | {'rdens':>6} {'rc2':>6} {'rl2':>7}"
    )
    for i, d in enumerate(masks()):
        m = d["mask"]
        if not d["kept"]:
            print(
                f"{d['name']:<12} {short(d['recipe']):<44} {d['density']:6.4f} skipped: density outside [{LOW}, {HIGH}]"
            )
            continue
        s = describe(m)
        count = int(allowed(m).sum())
        w = window(count)
        b = band(w)
        r = twin(count, i, SEEDS[0])
        print(
            f"{d['name']:<12} {short(d['recipe']):<44} {d['density']:6.4f} {s['effective']:6.4f} {s['reach1']:6.4f} {s['reach2']:6.4f} {s['causal2']:6.4f} {s['fiedler']:7.3f} {s['cfiedler']:7.3f} | {w:3d} {b.sum() / PAIRS:6.4f} {hops(b, 2)[LOWER].mean():6.4f} {fiedler(b):7.3f} | {r.sum() / PAIRS:6.4f} {hops(r, 2)[LOWER].mean():6.4f} {fiedler(r):7.3f}"
        )
    for d in masks():
        if d["family"] == "Digit":
            tile = d["mask"][:3, :3]
            power = tile
            for _ in range(d["recipe"]["level"] - 1):
                power = np.kron(power, tile)
            assert np.array_equal(power, d["mask"])
            print(
                f"{d['name']}: the mask is the {d['recipe']['level']}th Kronecker power of its 3 x 3 corner tile {''.join(str(int(x)) for x in tile.ravel())}, fill {int(tile.sum())}"
            )


# TASKS


def copy_batch(g, batch):
    lag = torch.randint(0, len(LAGS), (batch,), generator=g)
    x = torch.randint(0, VOCAB - len(LAGS), (batch, S), generator=g)
    x[:, 0] = VOCAB - len(LAGS) + lag
    src = torch.arange(S)[None, :] - torch.tensor(LAGS)[lag][:, None]
    return x, x.gather(1, src.clamp(min=0)), src >= 0, lag


def recall_batch(g, batch):
    keys = torch.rand(batch, VOCAB, generator=g).argsort(1)[:, :PAIRS_KV]
    values = torch.randint(0, VOCAB, (batch, PAIRS_KV), generator=g)
    pick = torch.randint(0, PAIRS_KV, (batch, S - QUERY), generator=g)
    x = torch.zeros(batch, S, dtype=torch.long)
    y = torch.zeros(batch, S, dtype=torch.long)
    x[:, 0:QUERY:2], x[:, 1:QUERY:2] = keys, values
    x[:, QUERY:], y[:, QUERY:] = keys.gather(1, pick), values.gather(1, pick)
    on = torch.zeros(batch, S, dtype=torch.bool)
    on[:, QUERY:] = True
    return x, y, on, torch.zeros(batch, dtype=torch.long)


TASKS = {"copy": copy_batch, "recall": recall_batch}

# MODEL


def init(seed, vocab=VOCAB):
    g = torch.Generator().manual_seed(seed)
    n = lambda *shape, std=0.02: torch.randn(*shape, generator=g) * std
    out = {
        "tok": n(vocab, WIDTH),
        "pos": n(S, WIDTH),
        "lnf_g": torch.ones(WIDTH),
        "lnf_b": torch.zeros(WIDTH),
        "head": n(WIDTH, vocab),
        "head_b": torch.zeros(vocab),
    }
    for i in range(LAYERS):
        out |= {
            f"{i}.ln1_g": torch.ones(WIDTH),
            f"{i}.ln1_b": torch.zeros(WIDTH),
            f"{i}.qkv": n(WIDTH, 3 * WIDTH),
            f"{i}.qkv_b": torch.zeros(3 * WIDTH),
            f"{i}.out": n(WIDTH, WIDTH, std=0.02 / math.sqrt(2 * LAYERS)),
            f"{i}.out_b": torch.zeros(WIDTH),
            f"{i}.ln2_g": torch.ones(WIDTH),
            f"{i}.ln2_b": torch.zeros(WIDTH),
            f"{i}.fc1": n(WIDTH, 4 * WIDTH),
            f"{i}.fc1_b": torch.zeros(4 * WIDTH),
            f"{i}.fc2": n(4 * WIDTH, WIDTH, std=0.02 / math.sqrt(2 * LAYERS)),
            f"{i}.fc2_b": torch.zeros(WIDTH),
        }
    return out


def ensemble(seeds, device, vocab=VOCAB):
    inits = {s: init(s, vocab) for s in set(seeds)}
    return {
        k: torch.nn.Parameter(torch.stack([inits[s][k] for s in seeds]).to(device))
        for k in inits[seeds[0]]
    }


def norm(h, g, b):
    return F.layer_norm(h, h.shape[-1:]) * g[:, None, :] + b[:, None, :]


def forward(p, x, bias):
    e, batch, t = x.shape
    h = (
        p["tok"][torch.arange(e, device=x.device)[:, None, None], x]
        + p["pos"][:, None, :t]
    )
    h = h.reshape(e, batch * t, WIDTH)
    for i in range(LAYERS):
        a = torch.baddbmm(
            p[f"{i}.qkv_b"][:, None, :],
            norm(h, p[f"{i}.ln1_g"], p[f"{i}.ln1_b"]),
            p[f"{i}.qkv"],
        )
        q, k, v = a.view(e, batch, t, 3, HEADS, WIDTH // HEADS).permute(
            3, 0, 1, 4, 2, 5
        )
        w = torch.softmax(
            q @ k.transpose(-1, -2) / math.sqrt(WIDTH // HEADS) + bias, dim=-1
        )
        o = (w @ v).permute(0, 1, 3, 2, 4).reshape(e, batch * t, WIDTH)
        h = h + torch.baddbmm(p[f"{i}.out_b"][:, None, :], o, p[f"{i}.out"])
        f = F.gelu(
            torch.baddbmm(
                p[f"{i}.fc1_b"][:, None, :],
                norm(h, p[f"{i}.ln2_g"], p[f"{i}.ln2_b"]),
                p[f"{i}.fc1"],
            )
        )
        h = h + torch.baddbmm(p[f"{i}.fc2_b"][:, None, :], f, p[f"{i}.fc2"])
    return torch.baddbmm(
        p["head_b"][:, None, :], norm(h, p["lnf_g"], p["lnf_b"]), p["head"]
    ).view(e, batch, t, -1)


def losses(logits, y, on):
    ce = F.cross_entropy(
        logits.reshape(-1, logits.shape[-1]), y.reshape(-1), reduction="none"
    ).view(y.shape)
    return (ce * on).sum((1, 2)) / on.sum((1, 2))


def schedule(opt, steps, warm=50):
    return torch.optim.lr_scheduler.LambdaLR(
        opt,
        lambda k: (
            (k + 1) / warm
            if k < warm
            else 0.5 * (1 + math.cos(math.pi * (k - warm) / max(1, steps - warm)))
        ),
    )


def device():
    return torch.device("mps" if torch.backends.mps.is_available() else "cpu")


def batches(task, seeds, step, batch, offset):
    rows = [
        TASKS[task](torch.Generator().manual_seed(1000003 * (s + offset) + step), batch)
        for s in seeds
    ]
    return [torch.stack([r[i] for r in rows]) for i in range(4)]


def train(task, conditions, args):
    dev = device()
    seeds = [s for _ in conditions for s in SEEDS]
    which = [i for i in range(len(conditions)) for _ in SEEDS]
    bias = torch.stack(
        [
            torch.where(torch.from_numpy(c[s - 1]), 0.0, float("-inf"))
            for c in conditions
            for s in SEEDS
        ]
    ).to(dev)[:, None, None]
    p = ensemble(seeds, dev)
    opt = torch.optim.AdamW(p.values(), lr=args.lr)
    sched = schedule(opt, args.steps)
    tail = torch.zeros(len(seeds), device=dev)
    for step in range(args.steps):
        x, y, on, _ = (t.to(dev) for t in batches(task, SEEDS, step, args.batch, 0))
        pick = torch.tensor([s - 1 for s in seeds], device=dev)
        loss = losses(forward(p, x[pick], bias), y[pick], on[pick])
        opt.zero_grad(set_to_none=True)
        loss.sum().backward()
        opt.step()
        sched.step()
        if step >= args.steps - 20:
            tail += loss.detach()
    tail = (tail / min(20, args.steps)).tolist()
    hits = torch.zeros(len(seeds), len(LAGS), device=dev)
    totals = torch.zeros(len(seeds), len(LAGS), device=dev)
    with torch.no_grad():
        for k in range(args.evals):
            x, y, on, lag = (
                t.to(dev) for t in batches(task, SEEDS, k, args.batch, 7919)
            )
            pick = torch.tensor([s - 1 for s in seeds], device=dev)
            right = (forward(p, x[pick], bias).argmax(-1) == y[pick]) & on[pick]
            onehot = F.one_hot(lag[pick], len(LAGS)).float()
            hits += torch.einsum("ebt,ebl->el", right.float(), onehot)
            totals += torch.einsum("ebt,ebl->el", on[pick].float(), onehot)
    acc = (hits.sum(1) / totals.sum(1)).tolist()
    bylag = (hits / totals.clamp(min=1)).tolist()
    del p, opt, sched, bias, hits, totals
    gc.collect()
    if dev.type == "mps":
        torch.mps.empty_cache()
    return [
        {
            "condition": which[i],
            "seed": seeds[i],
            "loss": tail[i],
            "acc": acc[i],
            "bylag": bylag[i],
        }
        for i in range(len(seeds))
    ]


# SWEEP


def conditions():
    out, seen = (
        [
            {
                "name": "dense",
                "kind": "dense",
                "w": S,
                "effective": 1.0,
                "masks": [LOWER] * 3,
            }
        ],
        set(),
    )
    for i, d in enumerate(masks()):
        if not d["kept"]:
            continue
        c = allowed(d["mask"])
        count = int(c.sum())
        w = window(count)
        out.append(
            {
                "name": d["name"],
                "kind": "design",
                "w": w,
                "effective": count / PAIRS,
                "masks": [c] * 3,
            }
        )
        out.append(
            {
                "name": f"random~{d['name']}",
                "kind": "random",
                "w": w,
                "effective": count / PAIRS,
                "masks": [twin(count, i, s) for s in SEEDS],
            }
        )
        if w not in seen:
            seen.add(w)
            out.append(
                {
                    "name": f"window-{w}",
                    "kind": "window",
                    "w": w,
                    "effective": float(band(w).sum() / PAIRS),
                    "masks": [band(w)] * 3,
                }
            )
    return out


def results_path():
    return os.path.join(DATA_DIR, "results.json")


def load():
    if os.path.exists(results_path()):
        with open(results_path()) as f:
            return json.load(f)
    return None


def sweep(args):
    conds = conditions()
    res = load()
    if args.fresh or not res or res["steps"] != args.steps:
        res = {
            "steps": args.steps,
            "batch": args.batch,
            "lr": args.lr,
            "evals": args.evals,
            "seeds": SEEDS,
            "runs": {"copy": {}, "recall": {}},
            "wall": {},
        }
    start = time.time()
    for task in TASKS:
        todo = [c for c in conds if c["name"] not in res["runs"][task]]
        for k in range(0, len(todo), args.group):
            group = todo[k : k + args.group]
            t = time.time()
            for r in train(task, [c["masks"] for c in group], args):
                res["runs"][task].setdefault(group[r["condition"]]["name"], []).append(
                    {key: r[key] for key in ("seed", "loss", "acc", "bylag")}
                )
            res["wall"][f"{task}:{k // args.group}"] = time.time() - t
            res["peak_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
            with open(results_path(), "w") as f:
                json.dump(res, f)
            print(
                f"{task} group {k // args.group + 1} of {-(-len(todo) // args.group)}: {len(group)} conditions x {len(SEEDS)} seeds, {time.time() - t:.1f} s, peak {res['peak_gb']:.2f} GB",
                flush=True,
            )
    print(
        f"sweep {time.time() - start:.0f} s this call, {sum(res['wall'].values()):.0f} s in all groups, {args.steps} steps, device {device().type}"
    )
    table(args)


def stats(runs, key):
    v = np.array([r[key] for r in sorted(runs, key=lambda r: r["seed"])])
    return v, v.mean(), v.std(ddof=1)


def sign_p(won, lost):
    n = won + lost
    return min(1.0, 2 * sum(math.comb(n, k) for k in range(min(won, lost) + 1)) / 2**n)


def table(args):
    res = load()
    conds = conditions()
    meta = {c["name"]: c for c in conds}
    designs = [c for c in conds if c["kind"] == "design"]
    pairs = lambda name: round(meta[name]["effective"] * PAIRS)
    print(
        f"{res['steps']} steps, batch {res['batch']}, lr {res['lr']}, seeds {res['seeds']}; loss = mean cross-entropy of the last 20 steps; acc on {res['evals']} fresh batches; sd over seeds, ddof 1; chance loss copy ln 23 = {math.log(23):.3f}, recall ln 32 = {math.log(32):.3f}"
    )
    print(
        f"pairs = allowed causal pairs of {PAIRS}; a window holds at or above its design's pairs, excess beside it"
    )
    print(
        f"{'condition':<20} {'eff':>6} {'pairs':>5} {'excess':>7} | {'copy loss per seed':<20} {'mean':>6} {'sd':>6} {'acc':>6} | {'recall loss per seed':<20} {'mean':>6} {'sd':>6} {'acc':>6}"
    )

    def line(name, ref=None):
        out = f"{name:<20} {meta[name]['effective']:6.4f} {pairs(name):5d} {f'{100 * (pairs(name) / ref - 1):+.1f}%' if ref else '':>7}"
        for task in TASKS:
            v, mu, sd = stats(res["runs"][task][name], "loss")
            _, acc, _ = stats(res["runs"][task][name], "acc")
            out += f" | {' '.join(f'{x:6.3f}' for x in v):<20} {mu:6.3f} {sd:6.3f} {acc:6.3f}"
        print(out)

    line("dense")
    score = {(t, b): [0, 0] for t in TASKS for b in ("random", "window", "dense")}
    paired, excess = [], []
    for c in designs:
        w = c["w"]
        excess.append((pairs(f"window-{w}") / pairs(c["name"]) - 1, c["name"]))
        print()
        line(c["name"])
        line(f"random~{c['name']}")
        line(f"window-{w}", pairs(c["name"]))
        row = []
        for task in TASKS:
            mine, mu, _ = stats(res["runs"][task][c["name"]], "loss")
            for kind, other in (
                ("random", f"random~{c['name']}"),
                ("window", f"window-{w}"),
                ("dense", "dense"),
            ):
                theirs, nu, _ = stats(res["runs"][task][other], "loss")
                score[(task, kind)][0 if mu < nu else 1] += 1
                row.append(f"{task} vs {kind} {int((mine < theirs).sum())}/3")
        paired.append(f"{c['name']:<12} " + ", ".join(row))
    print()
    print("paired seeds where the design loss is lower:")
    for r in paired:
        print(f"  {r}")
    print(
        "design mean loss below baseline mean loss, masks won:lost; the window at or above the design's density"
    )
    for (task, kind), (won, lost) in score.items():
        print(f"  {task} vs {kind}: {won}:{lost}")
    excess.sort()
    print(
        f"window pair excess over the design: {100 * excess[0][0]:.1f}% ({excess[0][1]}) to {100 * excess[-1][0]:.1f}% ({excess[-1][1]}), "
        + ", ".join(f"{n} {100 * e:.1f}%" for e, n in excess)
    )
    print(
        f"design minus random twin, loss mean over seeds per mask, over {len(designs)} masks:"
    )
    for task in TASKS:
        d = np.array([stats(res["runs"][task][c["name"]], "loss")[0] for c in designs])
        t = np.array(
            [
                stats(res["runs"][task][f"random~{c['name']}"], "loss")[0]
                for c in designs
            ]
        )
        g = (d - t).mean(1)
        won, lost = int((d < t).all(1).sum()), int((d > t).all(1).sum())
        print(
            f"  {task}: mean {g.mean():+.3f} nats, sd {g.std(ddof=1):.3f}, t {g.mean() / (g.std(ddof=1) / math.sqrt(len(g))):.2f}; lower on all 3 seeds {won}, higher on all 3 {lost}, two-sided sign test p {sign_p(won, lost):.3g}"
        )


def coverage(c):
    r = hops(c, 2)
    return [
        [
            float(np.mean([c[t, t - n] for t in range(n, S)])),
            float(np.mean([c[t, t - n] and c[t, 0] for t in range(n, S)])),
            float(np.mean([r[t, t - n] and r[t, 0] for t in range(n, S)])),
        ]
        for n in LAGS
    ]


def corr(runs, center=False):
    a = np.concatenate(
        [(x - x.mean(1, keepdims=True) if center else x).ravel() for x, _ in runs]
    )
    v = np.concatenate(
        [(y - y.mean(1, keepdims=True) if center else y).ravel() for _, y in runs]
    )
    return float(np.corrcoef(a, v)[0, 1])


def lags(args):
    res = load()
    conds = conditions()
    rows, seen = {}, set()
    print(
        "copy accuracy by lag n over the offset share, the share of queries t >= n that see t - n in one hop, means over seeds"
    )
    print(f"{'condition':<20} " + " ".join(f"{n:>9}" for n in LAGS))
    for c in conds:
        acc = np.array(
            [
                r["bylag"]
                for r in sorted(res["runs"]["copy"][c["name"]], key=lambda r: r["seed"])
            ]
        )
        cov = np.array([coverage(m) for m in c["masks"]])
        rows[c["name"]] = (c, acc, cov)
        print(
            f"{c['name']:<20} "
            + " ".join(
                f"{a:4.2f}/{v:4.2f}" for a, v in zip(acc.mean(0), cov[:, :, 0].mean(0))
            )
        )
        if c["kind"] == "design":
            seen.add((c["name"], float(c["masks"][0][:, 0].mean())))
    cells = np.array(
        [
            (a, *v)
            for _, acc, cov in rows.values()
            for a, v in zip(acc.ravel(), cov.reshape(-1, 3))
        ]
    )
    acc, *covs = cells.T
    print(
        f"{len(cells)} cells = {len(conds)} conditions x {len(SEEDS)} seeds x {len(LAGS)} lags, not independent; chance {1 / 23:.3f}; correlation of accuracy with"
    )
    for name, v in zip(
        (
            "the offset share",
            "the share seeing both t - n and 0 in one hop",
            "the share reaching both t - n and 0 within 2 hops",
        ),
        covs,
    ):
        print(f"  {name}: {np.corrcoef(v, acc)[0, 1]:.3f}")
    v = covs[0]
    for name, sel in (
        ("offset share 0", v == 0),
        ("offset share in (0, 0.5)", (v > 0) & (v < 0.5)),
        ("offset share in [0.5, 1)", (v >= 0.5) & (v < 1)),
        ("offset share 1", v == 1),
    ):
        print(
            f"  {name}: {int(sel.sum())} cells, accuracy mean {acc[sel].mean():.3f}, max {acc[sel].max():.3f}"
        )
    runs = {
        k: [(a[i : i + 1], cv[i : i + 1, :, 0]) for i in range(len(SEEDS))]
        for k, (_, a, cv) in rows.items()
    }
    every = [r for k in runs for r in runs[k]]
    rng = np.random.default_rng(0)
    names = list(runs)
    boot = [
        corr(
            [r for i in rng.integers(0, len(names), len(names)) for r in runs[names[i]]]
        )
        for _ in range(2000)
    ]
    print("offset share against accuracy:")
    print(
        f"  pooled {corr(every):.3f}; within run, each run's lag profile centered, {corr(every, True):.3f}; by-condition cluster bootstrap, 2000 draws, 95% [{np.percentile(boot, 2.5):.3f}, {np.percentile(boot, 97.5):.3f}]"
    )
    for name, kinds in (
        ("without windows", ("dense", "design", "random")),
        ("designs alone", ("design",)),
        ("twins alone", ("random",)),
    ):
        print(
            f"  {name}: {corr([r for k in runs for r in runs[k] if rows[k][0]['kind'] in kinds]):.3f}"
        )
    bound = 1 + (len(LAGS) - 1) / 23
    blind = [k for k, (_, a, _) in rows.items() if a.mean(0).sum() <= bound]
    print(
        f"blind default: a model that cannot read the lag token scores about 1 at one favoured lag and 1/23 at the rest, summed accuracy {bound:.3f}"
    )
    print(
        f"  {len(blind)} of {len(rows)} conditions sum at or under it, {sum(rows[k][0]['kind'] == 'design' for k in blind)} of {sum(c['kind'] == 'design' for c in conds)} designs: "
        + ", ".join(blind)
    )
    print(
        f"  above it: all kinds {corr([r for k in runs if k not in blind for r in runs[k]]):.3f}, designs and twins {corr([r for k in runs if k not in blind and rows[k][0]['kind'] in ('design', 'random') for r in runs[k]]):.3f}"
    )
    gaps = [
        (
            c["name"],
            n,
            rows[c["name"]][2][:, j, 0].mean()
            - rows[f"random~{c['name']}"][2][:, j, 0].mean(),
            rows[c["name"]][1][:, j].mean()
            - rows[f"random~{c['name']}"][1][:, j].mean(),
        )
        for c in conds
        if c["kind"] == "design"
        for j, n in enumerate(LAGS)
    ]
    go, ga = np.array([g[2] for g in gaps]), np.array([g[3] for g in gaps])
    nz, big = go != 0, np.abs(go) >= 0.2
    print(
        "design minus twin per (design, lag), means over seeds: offset gap against accuracy gap"
    )
    print(
        f"  correlation {np.corrcoef(go, ga)[0, 1]:.3f} over {len(gaps)} cells; same sign {int((np.sign(go[nz]) == np.sign(ga[nz])).sum())} of {int(nz.sum())} with offset gap nonzero, {int((np.sign(go[big]) == np.sign(ga[big])).sum())} of {int(big.sum())} with offset gap at least 0.2"
    )
    clear = [
        g
        for g in gaps
        if abs(g[2]) >= 0.1 and abs(g[3]) >= 0.03 and np.sign(g[2]) != np.sign(g[3])
    ]
    print(
        f"  {len(clear)} clear violations, offset gap at least 0.1 and accuracy gap of the other sign at least 0.03: "
        + ", ".join(f"{n} lag {l} ({o:+.2f}, {a:+.2f})" for n, l, o, a in clear)
    )
    print(
        "share of queries that see position 0, the lag token: "
        + ", ".join(f"{n} {z:.3f}" for n, z in sorted(seen, key=lambda x: (x[1], x[0])))
    )


def figure(args):
    out, seen = {"side": S, "masks": []}, set()
    for d in masks():
        if d["kept"] and d["family"] not in seen:
            seen.add(d["family"])
            out["masks"].append(
                {
                    "family": d["family"],
                    "name": d["name"],
                    "recipe": d["recipe"],
                    "density": round(d["density"], 4),
                    "effective": round(float(allowed(d["mask"]).sum() / PAIRS), 4),
                    "rows": d["mask"].astype(int).tolist(),
                }
            )
    with open(FIGURE, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print(
        f"{len(out['masks'])} masks, {os.path.getsize(FIGURE)} bytes: "
        + ", ".join(m["name"] for m in out["masks"])
    )


# OFFSETS

BASES = [2, 3, 4, 5, 9]
DRAWS = 200


def digits(n, b):
    out = []
    while n:
        n, d = divmod(n, b)
        out.append(d)
    return out


def digit_sets(b):
    top = 4 if b == 9 else b - 1
    return [
        (0, *c) for k in range(1, top) for c in itertools.combinations(range(1, b), k)
    ]


def offset_design(b, f):
    return [n for n in range(S) if set(digits(n, b)) <= set(f)]


def dname(b, f):
    return f"d{b}-{''.join(map(str, f))}"


def toeplitz(offsets):
    i, j = np.indices((S, S))
    return np.isin(i - j, offsets)


def weight(offsets):
    return int(sum(S - n for n in offsets))


def reach(offsets):
    v = np.zeros(S, dtype=np.int64)
    v[list(offsets)] = 1
    return np.flatnonzero(np.convolve(v, v)[:S])


def first_gap(offsets):
    r = set(reach(offsets).tolist())
    return next((n for n in range(S) if n not in r), S)


def strided(l):
    return sorted(set(range(l)) | set(range(0, S, l)))


def profile(offsets):
    r = reach(offsets)
    m = toeplitz(offsets)
    assert int(m.sum()) == weight(offsets) and int(hops(m, 2)[LOWER].sum()) == weight(r)
    return {
        "count": len(offsets),
        "pairs": weight(offsets),
        "cover1": len(offsets) / S,
        "cover2": len(r) / S,
        "causal2": weight(r) / PAIRS,
        "gap": first_gap(offsets),
    }


def matched(count, pairs, rng, draws):
    total = count * S - pairs
    ways = np.zeros((S, count, total + 1))
    ways[0, 0, 0] = 1
    for i in range(1, S):
        ways[i] = ways[i - 1]
        if i <= total:
            ways[i, 1:, i:] += ways[i - 1, :-1, : total + 1 - i]
    j, t = np.full(draws, count - 1), np.full(draws, total)
    pick = np.zeros((draws, S), dtype=bool)
    pick[:, 0] = True
    for i in range(S - 1, 0, -1):
        live = (j > 0) & (t >= i)
        p = np.where(
            live,
            ways[i - 1, np.maximum(j - 1, 0), np.maximum(t - i, 0)] / ways[i, j, t],
            0.0,
        )
        take = rng.random(draws) < p
        pick[take, i] = True
        j, t = j - take, t - i * take
    assert (j == 0).all() and (t == 0).all()
    return [np.flatnonzero(r).tolist() for r in pick]


def covers(sets):
    c = np.array([len(reach(o)) / S for o in sets])
    return float(c.mean()), float(np.percentile(c, 5)), int((c == 1).sum())


def nearest_strided(count):
    return min(
        range(1, S), key=lambda l: (abs(len(strided(l)) - count), weight(strided(l)))
    )


def offset_rows():
    rows, seen = [], {}
    for b in BASES:
        for f in digit_sets(b):
            d = offset_design(b, f)
            fs = {x + y for x in f for y in f}
            ok = set(range(b)) <= fs
            g = S if ok else min(set(range(b)) - fs)
            if tuple(d) in seen:
                seen[tuple(d)]["aliases"].append(dname(b, f))
                seen[tuple(d)]["criterion"].append([ok, g])
                continue
            row = {
                "name": dname(b, f),
                "base": b,
                "digits": list(f),
                "aliases": [],
                "criterion": [[ok, g]],
                "offsets": d,
                **profile(d),
            }
            seen[tuple(d)] = row
            rows.append(row)
    for k, row in enumerate(rows):
        count, pairs = row["count"], row["pairs"]
        rng = np.random.default_rng([7, k])
        row["random"] = covers(
            [
                [0, *rng.choice(np.arange(1, S), count - 1, replace=False)]
                for _ in range(DRAWS)
            ]
        )
        row["twin"] = covers(matched(count, pairs, rng, DRAWS))
        l = nearest_strided(count)
        row["strided"] = {"l": l, **profile(strided(l))}
        row["window"] = profile(list(range(count)))
    return rows


def offsets(args):
    rows = offset_rows()
    strides = [{"l": l, **profile(strided(l))} for l in range(1, S)]
    named = sum(1 + len(r["aliases"]) for r in rows)
    print(
        f"side {S}, {PAIRS} causal pairs; Toeplitz mask: key j allowed for query i when i - j is an offset; D(F) = the n < {S} whose base-b digits all lie in F; coverage = share of offsets 0..{S - 1} in D (1 hop) or D + D (2 hops); causal2 = share of causal pairs within 2 hops; gap = least offset outside D + D, so queries 0..gap-1 reach every earlier position in 2 hops and query gap does not"
    )
    print(
        f"bases {BASES}: {named} digit sets F with 0 in F and 2 <= |F| < b (|F| <= 4 at b = 9; none at b = 2), {len(rows)} distinct offset sets after de-duplication"
    )
    print(
        f"baselines at the design's offset count k: strided {{0..l-1}} + lZ with l nearest in count; random = {{0}} + k-1 uniform offsets, {DRAWS} draws; twin = {{0}} + k-1 offsets with the same pair count, uniform by exact counting, {DRAWS} draws; window = offsets 0..k-1; mean, 5th percentile and full-coverage count of 2-hop coverage"
    )
    print(
        f"{'design':<10} {'aliases':<10} {'k':>3} {'pairs':>5} {'dens':>6} {'c1':>6} {'c2':>6} {'causal2':>7} {'gap':>3} | {'l':>2} {'k':>3} {'pairs':>5} {'c2':>5} | {'rand':>5} {'p5':>5} {'full':>4} | {'twin':>5} {'p5':>5} {'full':>4} | {'win':>5} {'pairs':>5}"
    )
    for r in rows:
        s, w = r["strided"], r["window"]
        print(
            f"{r['name']:<10} {','.join(r['aliases']):<10} {r['count']:3d} {r['pairs']:5d} {r['pairs'] / PAIRS:6.4f} {r['cover1']:6.4f} {r['cover2']:6.4f} {r['causal2']:7.4f} {r['gap']:3d} | {s['l']:2d} {s['count']:3d} {s['pairs']:5d} {s['cover2']:5.3f} | {r['random'][0]:5.3f} {r['random'][1]:5.3f} {r['random'][2]:4d} | {r['twin'][0]:5.3f} {r['twin'][1]:5.3f} {r['twin'][2]:4d} | {w['cover2']:5.3f} {w['pairs']:5d}"
        )
    agree = sum(
        (c[0] == (r["cover2"] == 1.0)) and (c[1] == r["gap"])
        for r in rows
        for c in r["criterion"]
    )
    print(
        f"digit criterion, full 2-hop coverage iff F + F contains 0..b-1 and otherwise gap = least digit outside F + F: agrees on {agree} of {named} digit sets"
    )
    print(
        f"matrix check: every design's 2-hop causal pairs equal the offset count of D + D, by the {S} x {S} boolean square"
    )
    front, best = [], -1.0
    for r in sorted(rows, key=lambda r: (r["pairs"], -r["cover2"])):
        if r["cover2"] > best:
            front.append(r)
            best = r["cover2"]
    print(
        "Pareto front of the digit designs, pairs against 2-hop coverage: "
        + ", ".join(f"{r['name']} {r['pairs']}/{r['cover2']:.4f}" for r in front)
    )
    full = sorted((r for r in rows if r["cover2"] == 1.0), key=lambda r: r["pairs"])
    print(
        f"full 2-hop coverage, {len(full)} of {len(rows)} distinct designs, lowest pairs first:"
    )
    for r in full:
        print(
            f"  {r['name']}{' = ' + ', '.join(r['aliases']) if r['aliases'] else ''}: k {r['count']}, pairs {r['pairs']}/{PAIRS} = {r['pairs'] / PAIRS:.4f}, offsets {' '.join(map(str, r['offsets'])) if r['count'] <= 24 else '...'}"
        )
    sfull = sorted((s for s in strides if s["cover2"] == 1.0), key=lambda s: s["pairs"])
    print(
        f"strided, every l in 1..{S - 1} has full 2-hop coverage ({len(sfull)} of {len(strides)}), lowest pairs first: "
        + ", ".join(f"l {s['l']} k {s['count']} pairs {s['pairs']}" for s in sfull[:6])
    )
    fewest = min(strides, key=lambda s: (s["count"], s["pairs"]))
    print(
        f"fewest strided offsets: l {fewest['l']}, k {fewest['count']}, pairs {fewest['pairs']}; fewest design offsets with full coverage: {min(full, key=lambda r: (r['count'], r['pairs']))['name']}, k {min(r['count'] for r in full)}"
    )
    print(
        f"lowest density with full 2-hop coverage: design {full[0]['name']} {full[0]['pairs']}/{PAIRS} = {full[0]['pairs'] / PAIRS:.6f} against strided l {sfull[0]['l']} {sfull[0]['pairs']}/{PAIRS} = {sfull[0]['pairs'] / PAIRS:.6f}, difference {full[0]['pairs'] - sfull[0]['pairs']:+d} pairs"
    )
    print(
        f"counting bound: k offsets give at most k(k+1)/2 sums, so full 2-hop coverage of {S} offsets needs k >= {next(k for k in range(S) if k * (k + 1) // 2 >= S)}"
    )
    for d in masks():
        c = allowed(d["mask"])
        offs = [n for n in range(S) if all(c[t, t - n] for t in range(n, S))]
        if d["kept"] and np.array_equal(c, toeplitz(offs)):
            print(
                f"{d['name']} from the 2d census is the Toeplitz mask of {len(offs)} offsets {' '.join(map(str, offs[:12]))} ..., {int(c.sum())} pairs, 2-hop coverage {len(reach(offs)) / S:.4f}"
            )
    with open(os.path.join(DATA_DIR, "offsets.json"), "w") as f:
        json.dump({"side": S, "pairs": PAIRS, "rows": rows, "strided": strides}, f)


def basis(args):
    sums = np.array([(a, b) for a in range(S) for b in range(a, S - a)])
    m, (a, b), k = len(sums), sums.T, np.arange(len(sums))
    rows = np.r_[a + b, S + k, S + k, S + m + k, S + m + k]
    cols = np.r_[S + k, S + k, a, S + k, b]
    vals = np.r_[np.ones(m), np.ones(m), -np.ones(m), np.ones(m), -np.ones(m)]
    cover = LinearConstraint(
        coo_matrix((vals, (rows, cols)), shape=(S + 2 * m, S + m)).tocsr(),
        np.r_[np.ones(S), np.full(2 * m, -np.inf)],
        np.r_[np.full(S, np.inf), np.zeros(2 * m)],
    )
    cost = np.r_[S - np.arange(S), np.zeros(m)]
    lp = milp(cost, constraints=cover, integrality=np.zeros(S + m), bounds=Bounds(0, 1))
    t = time.time()
    res = milp(
        cost,
        constraints=cover,
        integrality=np.r_[np.ones(S), np.zeros(m)],
        bounds=Bounds(0, 1),
        options={"node_limit": args.nodes},
    )
    offs = np.flatnonzero(res.x[:S] > 0.5).tolist()
    found = profile(offs)
    assert found["cover2"] == 1.0 and found["pairs"] == round(res.fun)
    print(
        f"integer program over all offset sets at side {S}: minimise pairs subject to D + D containing 0..{S - 1}; {S} binary offsets, {m} pair variables, HiGHS, node limit {args.nodes}, {res.mip_node_count} nodes, {time.time() - t:.0f} s, status: {res.message}"
    )
    print(
        f"best found: k {found['count']}, pairs {found['pairs']}/{PAIRS} = {found['pairs'] / PAIRS:.4f}, offsets {' '.join(map(str, offs))}; 2-hop coverage rechecked by convolution and by the {S} x {S} boolean square: {found['cover2']:.4f}"
    )
    design, stride = (
        weight(offset_design(3, (0, 1))),
        min(weight(strided(l)) for l in range(1, S)),
    )
    print(
        f"LP relaxation {lp.fun:.1f}; dual bound {res.mip_dual_bound:.1f}, read from the solver, not re-derived: no offset set with full 2-hop coverage has fewer pairs if HiGHS is right; the best found is {design - found['pairs']} pairs below d3-01 ({design}) and {stride - found['pairs']} below the best strided mask ({stride})"
    )
    with open(os.path.join(DATA_DIR, "basis.json"), "w") as f:
        json.dump(
            {
                "nodes": args.nodes,
                "offsets": offs,
                "pairs": found["pairs"],
                "bound": float(res.mip_dual_bound),
                "message": res.message,
            },
            f,
        )


# LANGUAGE


def corpus():
    with open(os.path.join(DATA_DIR, "tinyshakespeare.txt")) as f:
        s = f.read()
    chars = sorted(set(s))
    lookup = {c: i for i, c in enumerate(chars)}
    ids = torch.tensor([lookup[c] for c in s], dtype=torch.long)
    cut = int(0.9 * len(ids))
    return ids[:cut], ids[cut:], len(chars)


def lm_batch(ids, seed, step, batch, offset):
    g = torch.Generator().manual_seed(1000003 * (seed + offset) + step)
    start = torch.randint(0, len(ids) - S, (batch,), generator=g)
    chunk = ids[start[:, None] + torch.arange(S + 1)]
    return chunk[:, :-1], chunk[:, 1:]


def lm_stack(ids, step, args, offset, dev):
    rows = [lm_batch(ids, s, step, args.batch, offset) for s in args.seeds]
    return torch.stack([r[0] for r in rows]).to(dev), torch.stack(
        [r[1] for r in rows]
    ).to(dev)


def lm_train(conditions, args, data):
    train_ids, val_ids, vocab = data
    dev = device()
    n = len(args.seeds)
    seeds = [s for _ in conditions for s in args.seeds]
    which = [i for i in range(len(conditions)) for _ in args.seeds]
    bias = torch.stack(
        [
            torch.where(torch.from_numpy(c[k]), 0.0, float("-inf"))
            for c in conditions
            for k in range(n)
        ]
    ).to(dev)[:, None, None]
    p = ensemble(seeds, dev, vocab)
    opt = torch.optim.AdamW(p.values(), lr=args.lr)
    sched = schedule(opt, args.steps)
    pick = torch.tensor([k for _ in conditions for k in range(n)], device=dev)
    every = torch.ones(len(seeds), args.batch, S, dtype=torch.bool, device=dev)
    tail = torch.zeros(len(seeds), device=dev)
    val = torch.zeros(len(seeds), device=dev)
    for step in range(args.steps):
        x, y = lm_stack(train_ids, step, args, 0, dev)
        loss = losses(forward(p, x[pick], bias), y[pick], every)
        opt.zero_grad(set_to_none=True)
        loss.sum().backward()
        opt.step()
        sched.step()
        if step >= args.steps - 20:
            tail += loss.detach()
            vx, vy = lm_stack(val_ids, step, args, 7919, dev)
            with torch.no_grad():
                val += losses(forward(p, vx[pick], bias), vy[pick], every)
    tail, val = (
        (tail / min(20, args.steps)).tolist(),
        (val / min(20, args.steps)).tolist(),
    )
    del p, opt, sched, bias, every
    gc.collect()
    if dev.type == "mps":
        torch.mps.empty_cache()
    return [
        {"condition": which[i], "seed": seeds[i], "train": tail[i], "val": val[i]}
        for i in range(len(seeds))
    ]


def lm_conditions(seeds):
    rows = offset_rows()
    e = len(seeds)
    full = sorted((r for r in rows if r["cover2"] == 1.0), key=lambda r: r["pairs"])[:3]
    l = min(
        range(1, S),
        key=lambda l: (abs(weight(strided(l)) - full[0]["pairs"]), weight(strided(l))),
    )
    toe = [(r["name"], "design", r["offsets"]) for r in full] + [
        (f"strided-{l}", "strided", strided(l)),
        ("d3-02", "cantor", offset_design(3, (0, 2))),
    ]
    out, seen = (
        [
            {
                "name": "dense",
                "kind": "dense",
                "w": S,
                "pairs": PAIRS,
                "offsets": None,
                "masks": [LOWER] * e,
            }
        ],
        set(),
    )

    def windowed(w):
        if w not in seen:
            seen.add(w)
            out.append(
                {
                    "name": f"window-{w}",
                    "kind": "window",
                    "w": w,
                    "pairs": int(band(w).sum()),
                    "offsets": [list(range(w + 1))] * e,
                    "masks": [band(w)] * e,
                }
            )

    for k, (name, kind, offs) in enumerate(toe):
        pairs = weight(offs)
        w = window(pairs)
        twins = [
            matched(len(offs), pairs, np.random.default_rng([11, k, s]), 1)[0]
            for s in seeds
        ]
        out.append(
            {
                "name": name,
                "kind": kind,
                "w": w,
                "pairs": pairs,
                "offsets": [offs] * e,
                "masks": [toeplitz(offs)] * e,
            }
        )
        out.append(
            {
                "name": f"random~{name}",
                "kind": "random",
                "w": w,
                "pairs": pairs,
                "offsets": twins,
                "masks": [toeplitz(o) for o in twins],
            }
        )
        windowed(w)
    for i, d in enumerate(masks()):
        if d["name"] == "general-2":
            c = allowed(d["mask"])
            count = int(c.sum())
            w = window(count)
            out.append(
                {
                    "name": "general-2",
                    "kind": "2d",
                    "w": w,
                    "pairs": count,
                    "offsets": None,
                    "masks": [c] * e,
                }
            )
            out.append(
                {
                    "name": "random~general-2",
                    "kind": "random",
                    "w": w,
                    "pairs": count,
                    "offsets": None,
                    "masks": [twin(count, i, s) for s in seeds],
                }
            )
            windowed(w)
    for c in out:
        assert all(int(m.sum()) == c["pairs"] for m in c["masks"])
    return out


def lm_path(args):
    if args.steps == 600 and args.seeds == SEEDS and not args.only:
        return os.path.join(DATA_DIR, "shakespeare.json")
    return os.path.join(
        DATA_DIR, f"shakespeare-{args.steps}-{'-'.join(map(str, args.seeds))}.json"
    )


def shakespeare(args):
    conds = [
        c for c in lm_conditions(args.seeds) if not args.only or c["name"] in args.only
    ]
    assert not args.only or len(conds) == len(args.only)
    data = corpus()
    res = None
    if os.path.exists(lm_path(args)) and not args.fresh:
        with open(lm_path(args)) as f:
            res = json.load(f)
    if not res or res["steps"] != args.steps or res["seeds"] != args.seeds:
        res = {
            "steps": args.steps,
            "batch": args.batch,
            "lr": args.lr,
            "seeds": args.seeds,
            "vocab": data[2],
            "chars": [len(data[0]), len(data[1])],
            "conditions": {},
            "runs": {},
            "wall": {},
        }
    for c in conds:
        res["conditions"][c["name"]] = {
            k: c[k] for k in ("kind", "w", "pairs", "offsets")
        }
    todo = [c for c in conds if c["name"] not in res["runs"]]
    start = time.time()
    for k in range(0, len(todo), args.group):
        group = todo[k : k + args.group]
        t = time.time()
        for r in lm_train([c["masks"] for c in group], args, data):
            res["runs"].setdefault(group[r["condition"]]["name"], []).append(
                {key: r[key] for key in ("seed", "train", "val")}
            )
        res["wall"][str(len(res["wall"]))] = time.time() - t
        res["peak_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
        with open(lm_path(args), "w") as f:
            json.dump(res, f)
        print(
            f"group {k // args.group + 1} of {-(-len(todo) // args.group)}: {len(group)} conditions x {len(args.seeds)} seeds, {time.time() - t:.1f} s, peak {res['peak_gb']:.2f} GB",
            flush=True,
        )
    print(
        f"shakespeare {time.time() - start:.0f} s this call, {sum(res['wall'].values()):.0f} s in all groups, {args.steps} steps, device {device().type}"
    )
    offsets_table(args)


def seed_vals(runs, name, seeds):
    by = {r["seed"]: r["val"] for r in runs[name]}
    return np.array([by[s] for s in seeds])


def paired(res, a, b):
    return seed_vals(res["runs"], a, res["seeds"]) - seed_vals(
        res["runs"], b, res["seeds"]
    )


def tstat(d):
    sd = d.std(ddof=1)
    t = d.mean() / (sd / math.sqrt(len(d)))
    return f"mean {d.mean():+.4f}, sd {sd:.4f}, t {t:+.2f}, df {len(d) - 1}, two-sided p {2 * tdist.sf(abs(t), len(d) - 1):.3g}, lower on {int((d < 0).sum())} of {len(d)}, per seed {' '.join(f'{x:+.4f}' for x in d)}"


def lm_runs(steps):
    runs = {}
    for path in sorted(glob.glob(os.path.join(DATA_DIR, "shakespeare*.json"))):
        with open(path) as f:
            r = json.load(f)
        if r["steps"] == steps:
            for name, rs in r["runs"].items():
                have = {x["seed"] for x in runs.setdefault(name, [])}
                runs[name] += [x for x in rs if x["seed"] not in have]
    return runs


def offsets_table(args):
    with open(lm_path(args)) as f:
        res = json.load(f)
    meta = res["conditions"]
    done = [n for n in meta if n in res["runs"]]
    ln2 = math.log(2)
    print(
        f"tiny Shakespeare, {res['vocab']} characters, train {res['chars'][0]} and validation {res['chars'][1]} characters (last 10%); {res['steps']} steps, batch {res['batch']}, lr {res['lr']}, seeds {res['seeds']}; val = validation cross-entropy, one fresh batch at each of the last 20 steps, averaged; train = the last 20 training losses; sd over seeds, ddof 1; bits = nats / ln 2; uniform {math.log(res['vocab']):.3f} nats"
    )
    print(
        f"{'condition':<18} {'kind':<8} {'pairs':>5} {'dens':>6} | {'val nats per seed':<23} {'mean':>6} {'sd':>6} | {'bits':>6} {'sd':>6} | {'train':>6} | per seed: 2-hop coverage[offsets 1..4 held] offsets held in 0..8"
    )
    for name in done:
        c = meta[name]
        v = seed_vals(res["runs"], name, res["seeds"])
        tr = np.mean([r["train"] for r in res["runs"][name]])
        extra = (
            " ".join(
                f"{len(reach(o)) / S:.3f}[{','.join(str(n) for n in o if 1 <= n <= 4) or '-'}]{sum(n <= 8 for n in o)}"
                for o in c["offsets"]
            )
            if c["offsets"]
            else ""
        )
        print(
            f"{name:<18} {c['kind']:<8} {c['pairs']:5d} {c['pairs'] / PAIRS:6.4f} | {' '.join(f'{x:7.4f}' for x in v):<23} {v.mean():6.4f} {v.std(ddof=1):6.4f} | {v.mean() / ln2:6.4f} {v.std(ddof=1) / ln2:6.4f} | {tr:6.4f} | {extra}"
        )
    print(
        "ranked by mean validation nats: "
        + ", ".join(
            f"{n} {seed_vals(res['runs'], n, res['seeds']).mean():.4f}"
            for n in sorted(
                done, key=lambda n: seed_vals(res["runs"], n, res["seeds"]).mean()
            )
        )
    )
    margins = {1: [], 2: [], 3: []}
    for name in done:
        c = meta[name]
        if c["offsets"] and c["kind"] != "window" and f"window-{c['w']}" in done:
            for o, g in zip(c["offsets"], paired(res, name, f"window-{c['w']}")):
                margins[min(3, next(n for n in range(1, S) if n not in o))].append(g)
    if margins[1]:
        print(
            "offset-mask runs by the least positive offset missing, minus the same seed's window at or above their pairs:"
        )
        for m, g in margins.items():
            print(
                f"  {'missing ' + str(m) if m < 3 else 'holding 1 and 2'}: {len(g)} runs, from {min(g):+.4f} to {max(g):+.4f}, mean {np.mean(g):+.4f}"
            )
    if "random~general-2" in done:
        i = next(i for i, d in enumerate(masks()) if d["name"] == "general-2")
        count = meta["general-2"]["pairs"]
        sub = [
            float(np.mean([m[t, t - 1] for t in range(1, S)]))
            for m in (twin(count, i, s) for s in res["seeds"])
        ]
        print(
            f"random~general-2 holds offset 1 in {', '.join(f'{x:.4f}' for x in sub)} of the rows t >= 1 per seed, expected ({count} - {S})/{PAIRS - S} = {(count - S) / (PAIRS - S):.4f}; general-2 in all"
        )
    designs = [n for n in done if meta[n]["kind"] == "design"]
    stride = next(n for n in done if meta[n]["kind"] == "strided")
    print("paired differences in validation nats per seed, a minus b, n = seeds:")
    for a in [n for n in [*designs, "d3-02", stride, "general-2"] if n in done]:
        for b in (
            f"random~{a}",
            stride if meta[a]["kind"] in ("design", "cantor") else None,
            f"window-{meta[a]['w']}",
            "dense",
        ):
            if b in done:
                excess = (
                    f"; window pairs {100 * (meta[b]['pairs'] / meta[a]['pairs'] - 1):+.1f}% over {a}"
                    if b.startswith("window")
                    else ""
                )
                print(f"  {a} - {b}: {tstat(paired(res, a, b))}{excess}")
    for b in ("twin", "strided"):
        if all((f"random~{a}" if b == "twin" else stride) in done for a in designs):
            d = np.mean(
                [
                    paired(res, a, f"random~{a}" if b == "twin" else stride)
                    for a in designs
                ],
                axis=0,
            )
            print(
                f"  mean over {', '.join(designs)} of design - {'its twin' if b == 'twin' else stride}, per seed: {tstat(d)}"
            )
    for a in sorted(
        (n for n in done if meta[n]["kind"] == "window"), key=lambda n: meta[n]["w"]
    ):
        if "dense" in done:
            print(f"  {a} - dense: {tstat(paired(res, a, 'dense'))}")
    near = f"window-{meta[stride]['w']}"
    if designs and near in done:
        print(
            f"  {designs[0]} - {near}: {tstat(paired(res, designs[0], near))}; window pairs {100 * (meta[near]['pairs'] / meta[designs[0]]['pairs'] - 1):+.1f}% over {designs[0]}"
        )
    runs = lm_runs(res["steps"])
    common = sorted(
        set.intersection(*({r["seed"] for r in runs[n]} for n in [*designs, stride]))
    )
    print(
        f"design minus {stride} over every cached {res['steps']}-step seed, {len(common)} seeds {common}:"
    )
    for a in designs:
        print(
            f"  {a} - {stride}: {tstat(seed_vals(runs, a, common) - seed_vals(runs, stride, common))}"
        )
    d = np.mean(
        [seed_vals(runs, a, common) - seed_vals(runs, stride, common) for a in designs],
        axis=0,
    )
    print(f"  mean over {', '.join(designs)} - {stride}, per seed: {tstat(d)}")
    for a, b in itertools.combinations(designs, 2):
        print(
            f"  {b} - {a}: {tstat(seed_vals(runs, b, common) - seed_vals(runs, a, common))}"
        )


# MAIN


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "verb",
        choices=[
            "census",
            "sweep",
            "table",
            "lags",
            "figure",
            "offsets",
            "basis",
            "shakespeare",
            "offsets_table",
        ],
    )
    p.add_argument("--steps", type=int, default=None)
    p.add_argument("--batch", type=int, default=64)
    p.add_argument("--lr", type=float, default=3e-3)
    p.add_argument("--evals", type=int, default=4)
    p.add_argument("--group", type=int, default=6)
    p.add_argument("--fresh", action="store_true")
    p.add_argument("--nodes", type=int, default=3000)
    p.add_argument("--seeds", default="1,2,3")
    p.add_argument("--only", default="")
    args = p.parse_args()
    args.seeds = [int(x) for x in args.seeds.split(",")]
    args.only = [x for x in args.only.split(",") if x]
    args.steps = args.steps or (
        600 if args.verb in ("shakespeare", "offsets_table") else 300
    )
    {
        "census": census,
        "sweep": sweep,
        "table": table,
        "lags": lags,
        "figure": figure,
        "offsets": offsets,
        "basis": basis,
        "shakespeare": shakespeare,
        "offsets_table": offsets_table,
    }[args.verb](args)


if __name__ == "__main__":
    main()
