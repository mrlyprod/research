import argparse
import gc
import json
import math
import os
import resource
import sys
import time

os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.5"
os.environ["PYTORCH_MPS_LOW_WATERMARK_RATIO"] = "0.4"

import numpy as np
import torch
import torch.nn.functional as F
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
RESULTS = os.path.join(DATA_DIR, "results.json")
C = 16
LAGS = [129, 143, 144, 151, 168, 214, 254, 255]
GRID = [1.0, 0.3, 0.1, 0.03, 0.01]
DEPTH, LENGTH, START, BATCH, PROBES, HIT = 8, 512, 255, 32, 8, 0.5

# PATH COUNTS


def product(q, b, L):
    r = np.array([1.0])
    for j in range(L):
        s = b**j
        out = np.zeros(len(r) + s * (len(q) - 1))
        for k, w in enumerate(q):
            out[k * s : k * s + len(r)] += w * r
        r = out
    return r


def popcount(n):
    return np.array([int(m).bit_count() for m in n], dtype=float)


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not math.isfinite(float(x)) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    return x


def load():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {}


def save(key, value):
    res = load()
    res[key] = clean(value)
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(RESULTS, "w") as f:
        json.dump(res, f, indent=1)


def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


# FIELD


def arms():
    out = [
        {
            "id": "pure",
            "name": "pure K=3",
            "block": "pure",
            "K": 3,
            "L": 8,
            "gain": 1.0,
            "c": None,
        }
    ]
    for c in [2.0, 0.5, 0.125]:
        out.append(
            {
                "id": f"c{c:g}",
                "name": f"residual c={c:g}",
                "block": "res",
                "K": 2,
                "L": 10,
                "gain": c,
                "c": c,
            }
        )
    out.append(
        {
            "id": "flat",
            "name": "flat",
            "block": "flat",
            "K": 2,
            "L": 10,
            "gain": 1.0,
            "c": None,
        }
    )
    return out


def q_of(arm):
    if arm["block"] == "pure":
        return [arm["gain"]] * arm["K"]
    if arm["block"] == "res":
        return [1 + arm["gain"], arm["gain"]]
    return [1.0, arm["gain"]]


def field_draw(arm, relu, gen):
    K, L = arm["K"], arm["L"]
    T = 2**L
    std = ((2 if relu else 1) * arm["gain"] / C) ** 0.5
    x = torch.randn(1, C, T, generator=gen, dtype=torch.float64).requires_grad_()
    u = torch.randn(C, generator=gen, dtype=torch.float64)
    u = u / u.norm()
    h = x
    for j in range(L):
        d = 2**j
        W = torch.randn(C, C, K, generator=gen, dtype=torch.float64) * std
        if arm["block"] == "flat":
            W[:, :, K - 1] = 0
        z = F.conv1d(F.pad(h, ((K - 1) * d, 0)), W, dilation=d)
        a = torch.relu(z) if relu else z
        h = a if arm["block"] == "pure" else h + a
    (g,) = torch.autograd.grad(h[0, :, -1] @ u, x)
    return (g[0] ** 2).sum(0).flip(0).numpy()


def lsq(cols, y):
    A = np.column_stack([np.ones(len(y))] + cols)
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = A @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return pred, coef, r2


def corr(a, b):
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return math.nan
    return float(np.corrcoef(a, b)[0, 1])


def field_stats(E, se, r, c):
    n = np.arange(len(E), dtype=float)
    logE, logr = np.log(E), np.log(r)
    varying = np.ptp(logr) > 1e-9
    z = np.abs(E - r) / np.maximum(se, 1e-300)
    out = {
        "lags": len(E),
        "rel": float(np.median(np.abs(E - r) / r)),
        "corr": corr(logE, logr),
        "z3": float(np.mean(z > 3)),
        "zmax": float(z.max()),
    }
    out["slope_r"] = float(np.polyfit(logr, logE, 1)[0]) if varying else math.nan
    out["scale"] = float(np.exp(np.median(logE - logr)))
    x = n / n.max()
    pred, _, _ = lsq([x, x**2], logE)
    out["gauss_corr"], out["gauss_rel"] = (
        corr(logE, pred),
        float(np.median(np.abs(E - np.exp(pred)) / np.exp(pred))),
    )
    pred, _, _ = lsq([np.log(n + 1)], logE)
    out["dist_corr"], out["dist_rel"] = (
        corr(logE, pred),
        float(np.median(np.abs(E - np.exp(pred)) / np.exp(pred))),
    )
    p = popcount(n)
    _, coef, r2 = lsq([p], logE)
    out["pop_slope"], out["pop_r2"] = float(coef[1]), float(r2)
    out["pop_law"] = (
        math.log(c / (1 + c)) if c else (0.0 if np.ptp(r) == 0 else math.nan)
    )
    ratio = E / r
    out["ratio_by_r"] = (
        {f"{v:g}": float(ratio[r == v].mean()) for v in sorted(set(r.tolist()))}
        if len(set(r.tolist())) <= 40
        else {}
    )
    out["ratio_by_pop"] = {
        str(int(v)): float(ratio[p == v].mean()) for v in sorted(set(p.tolist()))
    }
    return out


def field(a):
    key = f"draws={a.draws} seed={a.seed}"
    store = load().get("field", {})
    if "rows" in store:
        store = {}
    run = store.get(key, {"draws": a.draws, "seed": a.seed, "channels": C, "rows": {}})
    for i, arm in enumerate(arms()):
        if a.arms and arm["id"] not in a.arms:
            continue
        r = product(q_of(arm), 2, arm["L"])[: 2 ** arm["L"]]
        for relu in [False, True]:
            t0 = time.time()
            gen = torch.Generator().manual_seed(a.seed + 2 * i + relu)
            acc = np.zeros(2 ** arm["L"])
            sq = np.zeros(2 ** arm["L"])
            for _ in range(a.draws):
                g = field_draw(arm, relu, gen)
                acc += g
                sq += g * g
            E = acc / a.draws
            se = np.sqrt(np.maximum(sq / a.draws - E**2, 0) / a.draws)
            row = {
                "arm": arm["name"],
                "act": "relu" if relu else "linear",
                "K": arm["K"],
                "L": arm["L"],
                "c": arm["c"],
                "draws": a.draws,
                "seed": a.seed,
            }
            row.update(field_stats(E, se, r, arm["c"]))
            row["E"], row["se"], row["wall"] = E.tolist(), se.tolist(), time.time() - t0
            run["rows"][f"{arm['name']} {row['act']}"] = row
            gc.collect()
    run["peak_gb"] = peak_gb()
    store[key] = run
    save("field", store)
    table_field(key, load()["field"][key])


def num(x, f=".4f"):
    return "-" if x is None else format(x, f)


def table_field(key, res):
    print(
        f"field {key}: C = {res['channels']}, base 2, energy = mean over draws of the squared input gradient summed over channels, |u| = 1, per-lag E and SE stored"
    )
    for w in res["rows"].values():
        tag = f"field {key}: {w['arm']} {w['act']}"
        print(
            f"{tag} L={w['L']} lags={w['lags']}: r_Q median rel err {num(w['rel'])}, corr(log E, log r) {num(w['corr'])}, slope on log r {num(w['slope_r'])}, median E/r {num(w['scale'])}, share |z|>3 {num(w['z3'])}, max |z| {num(w['zmax'], '.2f')}, {w['wall']:.1f} s"
        )
        print(
            f"{tag} rivals: gaussian in n corr {num(w['gauss_corr'])} median rel err {num(w['gauss_rel'])}; log distance corr {num(w['dist_corr'])} median rel err {num(w['dist_rel'])}"
        )
        law = w["pop_law"]
        gap = (
            "-"
            if law in (None, 0.0) or w["pop_slope"] is None
            else format(w["pop_slope"] / law - 1, "+.4f")
        )
        print(
            f"{tag} popcount slope {num(w['pop_slope'])} against log(c/(1+c)) {num(law)}, relative gap {gap}, R^2 {num(w['pop_r2'])}"
        )
        if w["ratio_by_r"]:
            print(
                f"{tag} mean E/r by r: "
                + ", ".join(f"{k}: {v:.3f}" for k, v in w["ratio_by_r"].items())
            )
            high = [v for k, v in w["ratio_by_r"].items() if float(k) >= 24]
            if high and "1" in w["ratio_by_r"]:
                print(
                    f"{tag} mean E/r at r = 1 {w['ratio_by_r']['1']:.3f}, over r >= 24 from {min(high):.3f} to {max(high):.3f}"
                )
        if w["c"] is not None or w["arm"] == "flat":
            print(
                f"{tag} mean E/r by popcount: "
                + ", ".join(f"{k}: {v:.3f}" for k, v in w["ratio_by_pop"].items())
            )
    print(f"field {key}: peak resident {res['peak_gb']:.3f} GB")


def spread(store):
    groups = {}
    for run in store.values():
        for name, w in run["rows"].items():
            groups.setdefault((run["draws"], name), []).append(w)
    for (draws, name), ws in groups.items():
        if len(ws) < 2:
            continue
        seeds = [w["seed"] for w in ws]
        print(
            f"field spread: {name}, {draws} draws, seeds {seeds}: slope on log r {span(ws, 'slope_r')}, median rel err {span(ws, 'rel')}, corr {span(ws, 'corr')}, share |z|>3 {span(ws, 'z3')}, max |z| {span(ws, 'zmax')}, popcount slope {span(ws, 'pop_slope')}, R^2 {span(ws, 'pop_r2')}"
        )


def span(ws, k):
    vals = [w[k] for w in ws if w[k] is not None]
    return f"{min(vals):.4f} to {max(vals):.4f}" if vals else "-"


# LEARN


def device(a):
    return torch.device(a.device)


def init(K, gen, dev):
    params = []
    for j in range(DEPTH):
        cin = 1 if j == 0 else C
        params.append(torch.randn(K * cin, C, generator=gen) * (2 / (K * cin)) ** 0.5)
    params.append(torch.randn(C, 1, generator=gen) * (1 / C) ** 0.5)
    return [p.to(dev).requires_grad_() for p in params]


def forward(params, x, K):
    h = x
    T = x.shape[1]
    for j, W in enumerate(params[:-1]):
        d = 2**j
        hp = F.pad(h, (0, 0, (K - 1) * d, 0))
        taps = torch.cat(
            [hp[:, (K - 1 - k) * d : (K - 1 - k) * d + T] for k in range(K)], dim=-1
        )
        h = F.gelu(taps @ W)
    return h @ params[-1]


def target(x):
    return sum(F.pad(x, (0, 0, n, 0))[:, : x.shape[1]] for n in LAGS)


def filt(params, probe, K):
    p = probe.clone().requires_grad_()
    (g,) = torch.autograd.grad(forward(params, p, K)[:, -1, 0].sum(), p)
    return g[:, :, 0].mean(0).flip(0).cpu().numpy()


def run(K, opt, lr, seed, steps, probe, dev, every, read=True):
    params = init(K, torch.Generator().manual_seed(1000 + seed), dev)
    data = torch.Generator().manual_seed(2000 + seed)
    o = torch.optim.Adam(params, lr) if opt == "adam" else torch.optim.SGD(params, lr)
    hits = [None] * len(LAGS)
    losses, fs = [], []
    status = "capped"
    for step in range(steps + 1):
        if read and step % every == 0:
            f = filt(params, probe, K)
            fs.append([float(f[n]) for n in LAGS])
            for i, n in enumerate(LAGS):
                if hits[i] is None and f[n] >= HIT:
                    hits[i] = step
            if all(h is not None for h in hits):
                status = "done"
                break
        if step == steps:
            break
        x = torch.randn(BATCH, LENGTH, 1, generator=data).to(dev)
        loss = ((forward(params, x, K) - target(x))[:, START:] ** 2).mean()
        o.zero_grad()
        loss.backward()
        o.step()
        losses.append(loss.item())
        if not math.isfinite(losses[-1]):
            status = "diverged"
            break
    del params, o
    return {"hits": hits, "losses": losses, "f": fs, "steps": step, "status": status}


def stable(K, probe, dev, tune):
    seeds, steps = ([99], 100) if tune == "registered" else ([97, 98, 99], 300)
    tried = []
    for lr in GRID:
        ok = True
        for seed in seeds:
            L = np.array(
                run(K, "sgd", lr, seed, steps, probe, dev, 1, read=False)["losses"]
            )
            good = (
                len(L) == steps
                and np.all(np.isfinite(L))
                and L.max() < 10 * L[0]
                and L[-25:].mean() < L[:25].mean()
            )
            tried.append(
                {
                    "lr": lr,
                    "seed": seed,
                    "ok": bool(good),
                    "first": float(L[:25].mean()),
                    "last": float(L[-25:].mean()) if len(L) >= 25 else math.nan,
                }
            )
            ok = ok and good
            gc.collect()
        if ok:
            return lr, tried
    return GRID[-1], tried


def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    return corr(ra, rb)


def counts(K):
    own = product([1.0] * K, 2, DEPTH)
    stern = product([1.0] * 3, 2, DEPTH)
    return np.array([own[n] for n in LAGS]), np.array([stern[n] for n in LAGS])


def score(times, inv, rng, shuffles):
    times = np.array(times, dtype=float)
    med = np.median(times, axis=0)
    rho = spearman(med, inv)
    null = (
        np.array([spearman(med, rng.permutation(inv)) for _ in range(shuffles)])
        if math.isfinite(rho)
        else np.array([math.nan])
    )
    pooled = spearman(times.ravel(), np.tile(inv, times.shape[0]))
    return {
        "rho": rho,
        "p": float(np.mean(null >= rho - 1e-12)) if math.isfinite(rho) else math.nan,
        "null95": float(np.nanquantile(null, 0.95)) if math.isfinite(rho) else math.nan,
        "pooled": pooled,
    }


def analyse(K, runs, shuffles):
    own, stern = counts(K)
    times = [[math.inf if h is None else h for h in r["hits"]] for r in runs]
    med = np.median(np.array(times, dtype=float), axis=0)
    rng = np.random.default_rng(5)
    n = np.array(LAGS, dtype=float)
    out = {"median": med.tolist(), "own_r": own.tolist(), "stern_r": stern.tolist()}
    out["own"] = score(times, 1 / own, rng, shuffles)
    out["stern"] = score(times, 1 / stern, rng, shuffles)
    out["distance"] = score(times, n, rng, shuffles)
    out["popcount"] = score(times, popcount(n), rng, shuffles)
    out["capped"] = int(sum(h is None for r in runs for h in r["hits"]))
    out["diverged"] = int(sum(r["status"] == "diverged" for r in runs))
    return out


def learn(a):
    torch.set_num_threads(a.threads)
    dev = device(a)
    probe = torch.randn(
        PROBES, LENGTH, 1, generator=torch.Generator().manual_seed(7)
    ).to(dev)
    res = {
        "device": dev.type,
        "threads": a.threads,
        "seeds": a.seeds,
        "cap": a.cap,
        "every": a.every,
        "tune": a.tune,
        "stacks": {},
    }
    t0 = time.time()
    for K in [3, 2, 4]:
        tk = time.time()
        lr, tried = stable(K, probe, dev, a.tune)
        entry = {"sgd_lr": lr, "tried": tried, "opts": {}}
        for opt, rate in [("sgd", lr), ("adam", a.adam)]:
            runs = []
            for seed in range(a.seeds):
                ts = time.time()
                out = run(K, opt, rate, seed, a.cap, probe, dev, a.every)
                out["wall"] = time.time() - ts
                out["losses"] = out["losses"][::25]
                out["f"] = out["f"][:: max(1, 25 // a.every)]
                runs.append(out)
                print(
                    f"learn {a.key}: K={K} {opt} lr {rate:g} seed {seed}: hits {out['hits']} {out['status']} at step {out['steps']}, {out['wall']:.1f} s"
                )
                sys.stdout.flush()
                gc.collect()
                if dev.type == "mps":
                    torch.mps.empty_cache()
            entry["opts"][opt] = {
                "lr": rate,
                "runs": runs,
                "stats": analyse(K, runs, a.shuffles),
            }
        entry["wall"] = time.time() - tk
        res["stacks"][str(K)] = entry
        res["wall"] = time.time() - t0
        res["peak_gb"] = peak_gb()
        save("learn", {**load().get("learn", {}), a.key: res})
    table_learn(a.key, load()["learn"][a.key])


def table_learn(key, res):
    print(
        f"learn {key}: GELU stack, b = 2, L = {DEPTH}, C = {C}, lags {LAGS}, batch {BATCH}, length {LENGTH}, read every {res['every']} steps at {PROBES} probes, hit f(n) >= {HIT}, cap {res['cap']}, seeds {res['seeds']}, sgd lr rule {res['tune']}, device {res['device']}, {res['threads']} threads"
    )
    _, stern = counts(3)
    n = np.array(LAGS, dtype=float)
    print(
        f"learn {key}: lag set r = s(n+1) {[round(v) for v in stern]}, popcount {[int(v) for v in popcount(n)]}, Spearman of r with n {spearman(stern, n):.3f}, with popcount {spearman(stern, popcount(n)):.3f}"
    )
    for K, entry in res["stacks"].items():
        tried = ", ".join(
            f"{t['lr']:g} seed {t['seed']} {'ok' if t['ok'] else 'no'} ({num(t['first'], '.3f')} -> {num(t['last'], '.3f')})"
            for t in entry["tried"]
        )
        print(f"learn {key}: K={K} sgd lr probe: {tried}; chosen {entry['sgd_lr']:g}")
        own, _ = counts(int(K))
        pred = 1 / own if np.ptp(own) > 0 else 1 / stern
        print(
            f"learn {key}: K={K} scored against {'1/r own' if np.ptp(own) > 0 else '1/s(n+1)'}, whose Spearman with n on the lag set is {spearman(pred, n):.3f}"
        )
        for opt, block in entry["opts"].items():
            s = block["stats"]
            per = [
                spearman([math.inf if h is None else h for h in r["hits"]], pred)
                for r in block["runs"]
            ]
            print(
                f"learn {key}: K={K} {opt} per-seed Spearman(time, predictor) {[num(v, '.3f') for v in per]}"
            )
            print(
                f"learn {key}: K={K} {opt} adjacent pairs, higher s(n+1) first / tie / reversed / unscored over runs: {order(block['runs'], stern)}"
            )
            for seed, r in enumerate(block["runs"]):
                print(
                    f"learn {key}: K={K} {opt} seed {seed}: hits {r['hits']}, {r['status']} at step {r['steps']}, final loss {num(r['losses'][-1] if r['losses'] else None, '.4f')}, {r['wall']:.1f} s"
                )
            print(
                f"learn {key}: K={K} {opt} lr {block['lr']:g}: median learn time {[None if m is None else round(m) for m in s['median']]}, own r {[round(v) for v in s['own_r']]}, stern r {[round(v) for v in s['stern_r']]}, capped {s['capped']} cells, diverged {s['diverged']} runs"
            )
            for name, label in [
                ("own", "1/r own"),
                ("stern", "1/s(n+1)"),
                ("distance", "n"),
                ("popcount", "popcount"),
            ]:
                v = s[name]
                print(
                    f"learn {key}: K={K} {opt} Spearman(median time, {label}) {num(v['rho'], '.3f')}, shuffle share >= {num(v['p'], '.4f')}, shuffle 95% {num(v['null95'], '.3f')}, pooled over cells {num(v['pooled'], '.3f')}"
                )
            med = dict(zip(LAGS, s["median"]))
            pairs = [(143, 144), (254, 255), (129, 151), (168, 214)]
            print(
                f"learn {key}: K={K} {opt} pairs: "
                + ", ".join(
                    f"{p}/{q} {num(med[p], '.0f')}/{num(med[q], '.0f')}"
                    for p, q in pairs
                )
            )
        print(f"learn {key}: K={K} {entry['wall']:.0f} seconds")
    contrast(key, res)
    print(
        f"learn {key}: {res['wall']:.0f} seconds, peak resident {res['peak_gb']:.3f} GB"
    )


def order(runs, stern):
    tally = [0, 0, 0, 0]
    for r in runs:
        for p, q in [(143, 144), (254, 255)]:
            i, j = LAGS.index(p), LAGS.index(q)
            hi, lo = (i, j) if stern[i] > stern[j] else (j, i)
            a, b = r["hits"][hi], r["hits"][lo]
            tally[
                3 if a is None or b is None else 0 if a < b else 1 if a == b else 2
            ] += 1
    return "/".join(str(t) for t in tally)


def contrast(key, res):
    st = res["stacks"]
    if "3" not in st or "2" not in st:
        return
    _, stern = counts(3)
    rng = np.random.default_rng(5)
    for opt in ["sgd", "adam"]:
        m3, m2 = (
            np.array(
                [
                    math.inf if m is None else m
                    for m in st[K]["opts"][opt]["stats"]["median"]
                ]
            )
            for K in ["3", "2"]
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = m3 / m2
        if not np.all(np.isfinite(ratio)):
            print(
                f"learn {key}: post hoc {opt} K=3 over K=2 median time: capped cells, no contrast"
            )
            continue
        v = score([ratio], 1 / stern, rng, 10000)
        print(
            f"learn {key}: post hoc {opt} K=3 over K=2 median time {[round(float(x), 3) for x in ratio]}, Spearman with 1/s(n+1) {num(v['rho'], '.3f')}, shuffle share >= {num(v['p'], '.4f')}, shuffle 95% {num(v['null95'], '.3f')}"
        )
        k = np.array(LAGS) != 255
        inv = 1 / stern
        print(
            f"learn {key}: post hoc {opt} without lag 255, Spearman with 1/s(n+1): K=3 median {spearman(m3[k], inv[k]):.3f}, K=2 median {spearman(m2[k], inv[k]):.3f}, ratio {spearman(ratio[k], inv[k]):.3f}; K=2 median on all 8 lags {spearman(m2, inv):.3f}"
        )


# MAIN


def table(a):
    res = load()
    store = res.get("field", {})
    if "rows" not in store:
        for key, run in store.items():
            table_field(key, run)
        spread(store)
    for key, block in res.get("learn", {}).items():
        table_learn(key, block)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("verb", choices=["field", "learn", "table", "all"])
    p.add_argument("--draws", type=int, default=200)
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--arms", nargs="*", default=[])
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--cap", type=int, default=3000)
    p.add_argument("--adam", type=float, default=1e-3)
    p.add_argument("--shuffles", type=int, default=10000)
    p.add_argument("--device", default="cpu")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--every", type=int, default=25)
    p.add_argument("--tune", choices=["registered", "strict"], default="registered")
    p.add_argument("--key", default="registered")
    a = p.parse_args()
    for verb in ["field", "learn"] if a.verb == "all" else [a.verb]:
        t = time.time()
        globals()[verb](a)
        print(f"{verb}: {time.time() - t:.1f} seconds")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
