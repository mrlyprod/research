use crate::engine::{Frame, Run, Stop};
use crate::grid::Grid;
use mrly_core::{json, Json};
use mrly_life::{churn, entropy};
use mrly_core::Result;

pub fn growth(boxes: &[u32]) -> f64 {
    let points: Vec<(f64, f64)> = boxes
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, &b)| b > 0)
        .map(|(g, &b)| ((g as f64).ln(), f64::from(b).ln()))
        .collect();
    if points.len() < 2 {
        return 0.0;
    }
    let n = points.len() as f64;
    let (sx, sy) = points
        .iter()
        .fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
    let (mx, my) = (sx / n, sy / n);
    let (num, den) = points.iter().fold((0.0, 0.0), |(p, q), (x, y)| {
        (p + (x - mx) * (y - my), q + (x - mx) * (x - mx))
    });
    if den == 0.0 {
        0.0
    } else {
        num / den
    }
}

fn frame_json(grid: &Grid) -> Json {
    json!({ "h": grid.h, "w": grid.w, "rows": grid.hex_rows() })
}

fn shape_json(grid: &Grid) -> Json {
    let (proper, tilings) = grid.cuts();
    json!({
        "h": grid.h,
        "w": grid.w,
        "pop": grid.population(),
        "symmetry": grid.symmetry(),
        "cuts": proper,
        "tilings": tilings,
    })
}

fn pair_churn(prev: &Frame, last: &Frame) -> Result<f64> {
    if prev.grid.h == 0 && last.grid.h == 0 {
        return Ok(0.0);
    }
    let r0 = prev.origin.0.min(last.origin.0);
    let c0 = prev.origin.1.min(last.origin.1);
    let r1 = (prev.origin.0 + prev.grid.h as i64).max(last.origin.0 + last.grid.h as i64);
    let c1 = (prev.origin.1 + prev.grid.w as i64).max(last.origin.1 + last.grid.w as i64);
    let (h, w) = ((r1 - r0) as usize, (c1 - c0) as usize);
    let place = |f: &Frame| {
        f.grid
            .embed(h, w, (f.origin.0 - r0) as usize, (f.origin.1 - c0) as usize)
            .cell()
    };
    Ok(churn(&[place(prev)?, place(last)?]))
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

fn sample(curve: &[u32]) -> Json {
    let first: Vec<u32> = curve.iter().take(16).copied().collect();
    let every: Vec<u32> = curve
        .iter()
        .enumerate()
        .skip(16)
        .filter(|(i, _)| i.is_multiple_of(16))
        .map(|(_, &v)| v)
        .collect();
    let last = curve.last().copied().unwrap_or(0);
    let (max_gen, max) =
        curve.iter().enumerate().fold(
            (0usize, 0u32),
            |acc, (g, &v)| if v > acc.1 { (g, v) } else { acc },
        );
    json!([first, every, last, [max, max_gen]])
}

pub fn short(run: &Run) -> Result<Json> {
    let crop = run.last.grid.crop();
    let (cuts, tilings) = crop.cuts();
    let entropy = if crop.h == 0 {
        0
    } else {
        entropy(&crop.cell()?)
    };
    let churn = match &run.prev {
        Some(prev) => pair_churn(prev, &run.last)?,
        None => 0.0,
    };
    let visited = run.heat.counts.iter().filter(|&&v| v > 0).count();
    let peak = run.heat.counts.iter().copied().max().unwrap_or(0);
    let mut out = json!({
        "stop": run.stop.kind(),
        "gen": run.gen,
    });
    let map = out.as_object_mut().unwrap();
    if let Some(period) = run.stop.period() {
        map.insert("period".into(), json!(period));
    }
    if let Stop::Mover {
        dr, dc, wall_by, ..
    } = &run.stop
    {
        map.insert("d".into(), json!([dr, dc]));
        if let Some(by) = wall_by {
            map.insert("wall_by".into(), json!(by));
        }
    }
    if let Some(settled) = run.settled {
        map.insert("settled".into(), json!(settled));
    }
    map.insert("growth".into(), json!(round3(growth(&run.boxes))));
    map.insert("entropy".into(), json!(entropy));
    map.insert("churn".into(), json!(round3(churn)));
    map.insert("symmetry".into(), json!(crop.symmetry()));
    if !cuts.is_empty() {
        map.insert("cuts".into(), json!(cuts));
    }
    if !tilings.is_empty() {
        map.insert("tilings".into(), json!(tilings));
    }
    map.insert("heat".into(), json!([visited, peak]));
    if run.flips {
        map.insert("flips".into(), json!(true));
    }
    if run.background != 0 {
        map.insert("background".into(), json!(run.background));
    }
    map.insert("pop".into(), sample(&run.pop));
    map.insert("box".into(), sample(&run.boxes));
    Ok(out)
}

pub fn readings(run: &Run) -> Result<Json> {
    let crop = run.last.grid.crop();
    let stop = match &run.stop {
        Stop::Mover {
            period,
            dr,
            dc,
            wall_by,
        } => json!({
            "kind": "mover", "gen": run.gen, "period": period, "dr": dr, "dc": dc, "wall_by": wall_by
        }),
        other => json!({ "kind": other.kind(), "gen": run.gen, "period": other.period() }),
    };
    let churn = match &run.prev {
        Some(prev) => pair_churn(prev, &run.last)?,
        None => 0.0,
    };
    let entropy = if crop.h == 0 {
        0
    } else {
        entropy(&crop.cell()?)
    };
    let level = run.heat.level_set();
    let visited = run.heat.counts.iter().filter(|&&v| v > 0).count();
    let peak = run.heat.counts.iter().copied().max().unwrap_or(0);
    Ok(json!({
        "mode": run.mode.token(),
        "stop": stop,
        "settled": run.settled,
        "gens": run.gen,
        "pop": run.pop,
        "box": run.boxes,
        "growth": growth(&run.boxes),
        "entropy": entropy,
        "churn": churn,
        "flips": run.flips,
        "background": run.background,
        "frame": frame_json(&crop),
        "origin": [run.last.origin.0, run.last.origin.1],
        "shape": shape_json(&crop),
        "heat": {
            "visited": visited,
            "peak": peak,
            "level": shape_json(&level.crop()),
        },
    }))
}
