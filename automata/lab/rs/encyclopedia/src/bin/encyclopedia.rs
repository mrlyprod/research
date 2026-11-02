use encyclopedia::cache::{append, chunk_file, done_keys, free_percent, line};
use encyclopedia::{control, laws, library};
use encyclopedia::{readings, run, run_named, short, Grammar, MaskSpec, Mode, Recipe, Setup};
use std::collections::{BTreeMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Instant;

struct Options {
    words: Vec<String>,
    x: usize,
    cap: usize,
    threads: usize,
    per: usize,
    minutes: f64,
    chunk: usize,
    floor: f64,
    full: bool,
    oeis: Option<String>,
}

fn options() -> Options {
    let mut out = Options {
        words: Vec::new(),
        x: 1024,
        cap: 301,
        threads: 4,
        per: 8,
        minutes: 20.0,
        chunk: 50000,
        floor: 20.0,
        full: false,
        oeis: None,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let value = args.get(i + 1).and_then(|v| v.parse::<usize>().ok());
        match (args[i].as_str(), value) {
            ("--x", Some(v)) => out.x = v,
            ("--cap", Some(v)) => out.cap = v,
            ("--threads", Some(v)) => out.threads = v.max(1),
            ("--per", Some(v)) => out.per = v.max(1),
            ("--minutes", Some(v)) => out.minutes = v as f64,
            ("--chunk", Some(v)) => out.chunk = v.max(1),
            ("--floor", Some(v)) => out.floor = v as f64,
            ("--full", _) => {
                out.full = true;
                i += 1;
                continue;
            }
            ("--oeis", _) => {
                out.oeis = args.get(i + 1).cloned();
            }
            (word, _) => {
                out.words.push(word.to_string());
                i += 1;
                continue;
            }
        }
        i += 2;
    }
    out
}

fn number(words: &[String], at: usize) -> usize {
    words
        .get(at)
        .and_then(|w| w.parse::<usize>().ok())
        .unwrap_or_else(|| panic!("argument {at} is a number"))
}

fn pool<I: Sync, T: Send>(
    recipes: Vec<I>,
    threads: usize,
    work: impl Fn(&I) -> T + Sync,
    mut sink: impl FnMut(usize, &I, T),
) {
    let next = AtomicUsize::new(0);
    let (tx, rx) = mpsc::channel();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next, recipes, work) = (&next, &recipes, &work);
            scope.spawn(move || loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= recipes.len() {
                    break;
                }
                let out = work(&recipes[i]);
                if tx.send((i, out)).is_err() {
                    break;
                }
            });
        }
        drop(tx);
        for (i, out) in rx {
            sink(i, &recipes[i], out);
        }
    });
}

fn class(recipe: &Recipe) -> String {
    match recipe.mode {
        Mode::Island => format!("island x{}", recipe.tiling),
        Mode::Wallpaper => "wallpaper".to_string(),
    }
}

fn count(grammar: &Grammar) {
    let seeds = grammar.seeds().unwrap();
    let masks = grammar.masks().unwrap();
    let rules = grammar.rules();
    let mut classes: BTreeMap<String, usize> = BTreeMap::new();
    let mut copies = 0;
    let total = grammar
        .walk(0, usize::MAX, |_, recipe| {
            *classes.entry(class(&recipe)).or_default() += 1;
            copies += usize::from(matches!(recipe.mask, MaskSpec::Copy { .. }));
        })
        .unwrap();
    println!(
        "seeds {} (distinct designs over sides {:?} levels {:?})",
        seeds.len(),
        grammar.sides,
        grammar.levels
    );
    println!(
        "masks {} (distinct popped designs over sides {:?} levels {:?}) plus copy and copyinv where they differ",
        masks.len(),
        grammar.mask_sides,
        grammar.mask_levels
    );
    println!(
        "rules {} ({} sources x {} rolls, life {}, {} pairs x {} rolls)",
        rules.len(),
        grammar.sources.len(),
        grammar.rolls.len(),
        grammar.life,
        grammar.pairs.len(),
        grammar.rolls.len()
    );
    for (name, n) in &classes {
        println!("{name}: {n}");
    }
    println!("copy-mask recipes: {copies}");
    println!("recipes: {total} (x {} cap {})", grammar.x, grammar.cap);
}

struct Done {
    rows: usize,
    bytes: usize,
    thread_seconds: f64,
    last: Option<(usize, String)>,
}

fn chunk(
    grammar: &Grammar,
    from: usize,
    to: usize,
    threads: usize,
    done: &HashSet<String>,
) -> Done {
    let mut recipes: Vec<(usize, Recipe)> = Vec::new();
    grammar
        .walk(from, to, |i, recipe| {
            if !done.contains(&recipe.name()) {
                recipes.push((i, recipe));
            }
        })
        .unwrap();
    let path = chunk_file(from, to);
    println!(
        "chunk {from}..{to}: {} to run, {} cached, writing {}",
        recipes.len(),
        done.len(),
        path.display()
    );
    let started = Instant::now();
    let mut batch = Vec::new();
    let mut out = Done {
        rows: 0,
        bytes: 0,
        thread_seconds: 0.0,
        last: None,
    };
    pool(
        recipes,
        threads,
        |(i, recipe)| {
            let clock = Instant::now();
            let setup = Setup::of(recipe).unwrap();
            let text = line(&recipe.name(), *i, &short(&run(&setup)).unwrap());
            (text, clock.elapsed().as_secs_f64())
        },
        |_, (i, recipe), (text, seconds)| {
            out.bytes += text.len() + 1;
            batch.push(text);
            out.thread_seconds += seconds;
            if out.last.as_ref().is_none_or(|(j, _)| *i > *j) {
                out.last = Some((*i, recipe.name()));
            }
            if batch.len() >= 64 {
                append(&path, &batch).unwrap();
                out.rows += batch.len();
                batch.clear();
            }
        },
    );
    append(&path, &batch).unwrap();
    out.rows += batch.len();
    println!(
        "wrote {} rows, {} bytes, in {:.1} s wall, {:.1} thread-seconds",
        out.rows,
        out.bytes,
        started.elapsed().as_secs_f64(),
        out.thread_seconds
    );
    out
}

fn big(grammar: &Grammar, o: &Options) {
    let total = grammar.count().unwrap();
    let clock = Instant::now();
    let budget = o.minutes * 60.0;
    let cached = done_keys();
    let (mut done, mut bytes, mut thread_seconds, mut from) = (0usize, 0usize, 0.0f64, 0usize);
    let mut cut: Option<(usize, String)> = None;
    while from < total {
        let free = free_percent();
        if free < o.floor {
            println!(
                "disk {free:.1}% free is under the floor {:.0}%: stopping",
                o.floor
            );
            break;
        }
        if clock.elapsed().as_secs_f64() > budget {
            println!(
                "clock {:.0} s is past the budget {budget:.0} s: no new chunk",
                clock.elapsed().as_secs_f64()
            );
            break;
        }
        let to = (from + o.chunk).min(total);
        println!(
            "disk {free:.1}% free, clock {:.0} s",
            clock.elapsed().as_secs_f64()
        );
        let out = chunk(grammar, from, to, o.threads, &cached);
        done += out.rows;
        bytes += out.bytes;
        thread_seconds += out.thread_seconds;
        if let Some(last) = out.last {
            cut = Some(last);
        }
        from = to;
    }
    let (index, name) = cut.unwrap_or((0, String::new()));
    println!(
        "cut: index {index} of {total}, last recipe {name:?}, {done} run now, {} cached before, {bytes} bytes ({:.0} a recipe), {thread_seconds:.0} thread-seconds, {:.0} s wall",
        cached.len(),
        bytes as f64 / done.max(1) as f64,
        clock.elapsed().as_secs_f64()
    );
}

fn timing(grammar: &Grammar, threads: usize, per: usize) {
    let mut totals: BTreeMap<String, usize> = BTreeMap::new();
    let total = grammar
        .walk(0, usize::MAX, |_, recipe| {
            *totals.entry(class(&recipe)).or_default() += 1;
        })
        .unwrap();
    println!("recipes: {total}");
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    let mut picked: Vec<Recipe> = Vec::new();
    grammar
        .walk(0, usize::MAX, |_, recipe| {
            let name = class(&recipe);
            let stride = (totals[&name] / per).max(1);
            let k = seen.entry(name).or_default();
            if k.is_multiple_of(stride)
                && picked.iter().filter(|r| class(r) == class(&recipe)).count() < per
            {
                picked.push(recipe);
            }
            *k += 1;
        })
        .unwrap();
    for x in [64usize, 256, 1024] {
        let recipes: Vec<Recipe> = picked
            .iter()
            .map(|r| Recipe {
                x,
                cap: grammar.cap,
                ..r.clone()
            })
            .collect();
        let mut rows: BTreeMap<String, Vec<(f64, &'static str)>> = BTreeMap::new();
        let started = Instant::now();
        pool(
            recipes,
            threads,
            |recipe| {
                let clock = Instant::now();
                let setup = Setup::of(recipe).unwrap();
                let out = run(&setup);
                (clock.elapsed().as_secs_f64() * 1000.0, out.stop.kind())
            },
            |_, recipe, (ms, kind)| rows.entry(class(recipe)).or_default().push((ms, kind)),
        );
        let wall = started.elapsed().as_secs_f64();
        for (name, list) in &rows {
            let mut ms: Vec<f64> = list.iter().map(|(t, _)| *t).collect();
            ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mean = ms.iter().sum::<f64>() / ms.len() as f64;
            let mut kinds: BTreeMap<&str, usize> = BTreeMap::new();
            for (_, k) in list {
                *kinds.entry(k).or_default() += 1;
            }
            println!(
                "x {x} cap {} {name}: n {} mean {mean:.1} ms median {:.1} ms max {:.1} ms, {} recipes in class, stops {:?}",
                grammar.cap,
                ms.len(),
                ms[ms.len() / 2],
                ms[ms.len() - 1],
                totals[name],
                kinds
            );
        }
        println!(
            "x {x}: {} recipes on {threads} threads in {wall:.1} s wall",
            picked.len()
        );
    }
}

fn mini() {
    let s = control::summary().unwrap();
    println!(
        "mini control on canvas {} for {} generations: {} runs",
        control::CANVAS,
        control::GENERATIONS,
        s.runs
    );
    for (name, counts) in &s.table {
        println!("{name}: {counts:?}");
    }
    println!(
        "numbers {:?} claimed {:?}",
        s.numbers,
        control::CLAIMED_NUMBERS
    );
    println!(
        "old designs {:?} claimed {:?}",
        s.old_designs,
        control::CLAIMED_DESIGNS
    );
    println!("all designs {:?}", s.all_designs);
    println!("life {:?} claimed {:?}", s.life, control::CLAIMED_LIFE);
}

fn main() {
    let o = options();
    let grammar = Grammar::standard(o.x, o.cap);
    match o.words.first().map(String::as_str) {
        Some("count") => count(&grammar),
        Some("list") => {
            let (from, to) = (number(&o.words, 1), number(&o.words, 2));
            for recipe in grammar.recipes(from, to).unwrap() {
                println!("{}", recipe.name());
            }
        }
        Some("run") => {
            let (recipe, out) = run_named(&o.words[1..].join(" ")).unwrap();
            let readings = if o.full {
                readings(&out).unwrap()
            } else {
                short(&out).unwrap()
            };
            println!("{}", line(&recipe.name(), 0, &readings));
        }
        Some("chunk") => {
            chunk(
                &grammar,
                number(&o.words, 1),
                number(&o.words, 2),
                o.threads,
                &done_keys(),
            );
        }
        Some("big") => big(&grammar, &o),
        Some("laws") => laws::report(&grammar, o.oeis.as_deref().map(std::path::Path::new)).unwrap(),
        Some("library") => library::write().unwrap(),
        Some("time") => timing(&grammar, o.threads, o.per),
        Some("mini") => mini(),
        _ => println!("commands: count | list FROM TO | run NAME [--full] | chunk FROM TO | big | laws | library | time | mini; options --x --cap --threads --per --minutes --chunk --floor --oeis PATH"),
    }
}
