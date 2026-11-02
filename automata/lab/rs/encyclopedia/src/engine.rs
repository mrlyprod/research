use crate::grid::Grid;
use crate::mask::{Kernel, MaskSpec};
use crate::recipe::{Mode, Recipe};
use mrly_core::Result;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Rule {
    pub budget: usize,
    pub birth: Vec<usize>,
    pub survive: Vec<usize>,
    tables: [Vec<u8>; 2],
    next_background: [u8; 2],
}

impl Rule {
    pub fn new(birth: &[usize], survive: &[usize], budget: usize) -> Rule {
        let n = budget;
        let b: Vec<bool> = (0..=n).map(|k| birth.contains(&k)).collect();
        let s: Vec<bool> = (0..=n).map(|k| survive.contains(&k)).collect();
        let next_background = [u8::from(b[0]), u8::from(s[n])];
        let mut tables = [vec![0u8; 2 * (n + 1)], vec![0u8; 2 * (n + 1)]];
        for bg in 0..2usize {
            for y in 0..2usize {
                for ny in 0..=n {
                    let x = y ^ bg;
                    let nx = if bg == 0 { ny } else { n - ny };
                    let fired = if x == 1 { s[nx] } else { b[nx] };
                    tables[bg][y * (n + 1) + ny] = u8::from(fired) ^ next_background[bg];
                }
            }
        }
        Rule {
            budget,
            birth: birth.to_vec(),
            survive: survive.to_vec(),
            tables,
            next_background,
        }
    }

    pub fn flips(&self) -> bool {
        self.next_background[0] == 1
    }

    pub fn next_background(&self, background: u8) -> u8 {
        self.next_background[usize::from(background)]
    }

    pub fn apply(&self, background: u8, cells: &[u8], counts: &[i32]) -> Vec<u8> {
        let table = &self.tables[usize::from(background)];
        let width = self.budget + 1;
        cells
            .iter()
            .zip(counts)
            .map(|(&y, &n)| {
                debug_assert!((0..=self.budget as i32).contains(&n));
                table[usize::from(y) * width + n as usize]
            })
            .collect()
    }
}

pub fn step(grid: &Grid, kernel: &Kernel, rule: &Rule, background: u8, wrap: bool) -> Grid {
    let counts = kernel.count(grid, wrap);
    Grid::of(grid.h, grid.w, rule.apply(background, &grid.bits, &counts))
}

#[derive(Clone, Debug)]
pub struct Setup {
    pub board: Grid,
    pub kernel: Kernel,
    pub rule: Rule,
    pub mode: Mode,
    pub x: usize,
    pub cap: usize,
}

impl Setup {
    pub fn of(recipe: &Recipe) -> Result<Setup> {
        let board = recipe
            .seed
            .build(recipe.side, recipe.level)?
            .tile(recipe.tiling);
        let kernel = match recipe.mask {
            MaskSpec::Design { shape, side, level } => Kernel::design(shape, side, level)?,
            MaskSpec::Copy { inverted } => Kernel::copy(
                recipe.seed,
                recipe.side,
                recipe.level,
                recipe.tiling,
                inverted,
            )?,
        };
        let birth = recipe.birth.values(kernel.budget)?;
        let survive = recipe.survive.values(kernel.budget)?;
        let rule = Rule::new(&birth, &survive, kernel.budget);
        Ok(Setup {
            board,
            kernel,
            rule,
            mode: recipe.mode,
            x: recipe.x,
            cap: recipe.cap,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Stop {
    Dead,
    Still,
    Loop {
        period: usize,
    },
    Mover {
        period: usize,
        dr: i64,
        dc: i64,
        wall_by: Option<usize>,
    },
    Wall,
    Limit,
}

impl Stop {
    pub fn kind(&self) -> &'static str {
        match self {
            Stop::Dead => "dead",
            Stop::Still => "still",
            Stop::Loop { .. } => "loop",
            Stop::Mover { .. } => "mover",
            Stop::Wall => "wall",
            Stop::Limit => "limit",
        }
    }

    pub fn period(&self) -> Option<usize> {
        match self {
            Stop::Dead | Stop::Still => Some(1),
            Stop::Loop { period } | Stop::Mover { period, .. } => Some(*period),
            Stop::Wall | Stop::Limit => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Frame {
    pub grid: Grid,
    pub origin: (i64, i64),
}

impl Frame {
    pub fn cropped(grid: &Grid, origin: (i64, i64)) -> Frame {
        match grid.bbox() {
            Some((r0, c0, r1, c1)) => Frame {
                grid: grid.window(r0, c0, r1 - r0 + 1, c1 - c0 + 1),
                origin: (origin.0 + r0 as i64, origin.1 + c0 as i64),
            },
            None => Frame {
                grid: Grid::new(0, 0),
                origin,
            },
        }
    }

    pub fn side(&self) -> u32 {
        self.grid.h.max(self.grid.w) as u32
    }
}

#[derive(Clone, Debug)]
pub struct Heat {
    pub h: usize,
    pub w: usize,
    pub origin: (i64, i64),
    pub counts: Vec<u32>,
}

impl Heat {
    fn new(h: usize, w: usize, origin: (i64, i64)) -> Heat {
        Heat {
            h,
            w,
            origin,
            counts: vec![0; h * w],
        }
    }

    fn add(&mut self, frame: &Frame) {
        let (r0, c0) = (
            frame.origin.0 - self.origin.0,
            frame.origin.1 - self.origin.1,
        );
        for r in 0..frame.grid.h {
            for c in 0..frame.grid.w {
                if frame.grid.get(r, c) != 0 {
                    let (hr, hc) = (r0 + r as i64, c0 + c as i64);
                    if (0..self.h as i64).contains(&hr) && (0..self.w as i64).contains(&hc) {
                        self.counts[hr as usize * self.w + hc as usize] += 1;
                    }
                }
            }
        }
    }

    pub fn level_set(&self) -> Grid {
        let mut visited: Vec<u32> = self.counts.iter().copied().filter(|&v| v > 0).collect();
        if visited.is_empty() {
            return Grid::new(0, 0);
        }
        visited.sort_unstable();
        let median = visited[(visited.len() - 1) / 2];
        Grid::of(
            self.h,
            self.w,
            self.counts.iter().map(|&v| u8::from(v >= median)).collect(),
        )
    }
}

#[derive(Clone, Debug)]
pub struct Run {
    pub mode: Mode,
    pub stop: Stop,
    pub gen: usize,
    pub settled: Option<usize>,
    pub pop: Vec<u32>,
    pub boxes: Vec<u32>,
    pub last: Frame,
    pub prev: Option<Frame>,
    pub heat: Heat,
    pub background: u8,
    pub flips: bool,
}

struct Walls {
    lo: i64,
    hi: i64,
    reach: i64,
}

impl Walls {
    fn hit(&self, frame: &Frame) -> bool {
        if frame.grid.h == 0 {
            return false;
        }
        let (r0, c0) = frame.origin;
        let (r1, c1) = (r0 + frame.grid.h as i64 - 1, c0 + frame.grid.w as i64 - 1);
        r0 - self.lo < self.reach
            || c0 - self.lo < self.reach
            || self.hi - r1 < self.reach
            || self.hi - c1 < self.reach
    }

    fn wall_by(
        &self,
        frame: &Frame,
        dr: i64,
        dc: i64,
        period: usize,
        gen: usize,
        x: usize,
    ) -> Option<usize> {
        let (r0, c0) = frame.origin;
        let (r1, c1) = (r0 + frame.grid.h as i64 - 1, c0 + frame.grid.w as i64 - 1);
        let mut k = u64::MAX;
        let mut room = |space: i64, speed: i64| {
            if speed > 0 {
                k = k.min((space / speed) as u64 + 1);
            }
        };
        room(self.hi - self.reach - r1, dr);
        room(r0 - self.lo - self.reach, -dr);
        room(self.hi - self.reach - c1, dc);
        room(c0 - self.lo - self.reach, -dc);
        if k == u64::MAX {
            return None;
        }
        let at = gen as u64 + k * period as u64;
        (at <= x as u64).then_some(at as usize)
    }
}

pub fn run(setup: &Setup) -> Run {
    match setup.mode {
        Mode::Island => island(setup),
        Mode::Wallpaper => wallpaper(setup),
    }
}

fn island(setup: &Setup) -> Run {
    let reach = setup.kernel.reach;
    let side = setup.board.h;
    let cap = setup.cap.max(side);
    let lo = -(((cap - side) / 2) as i64);
    let walls = Walls {
        lo,
        hi: lo + cap as i64 - 1,
        reach: reach as i64,
    };
    let mut frame = Frame::cropped(&setup.board, (0, 0));
    let mut heat = Heat::new(cap, cap, (lo, lo));
    let mut history: HashMap<u128, (usize, i64, i64)> = HashMap::new();
    let (mut pop, mut boxes) = (Vec::new(), Vec::new());
    let mut background = 0u8;
    let mut prev = None;
    let mut gen = 0;
    let record = |frame: &Frame, pop: &mut Vec<u32>, boxes: &mut Vec<u32>, heat: &mut Heat| {
        pop.push(frame.grid.population() as u32);
        boxes.push(frame.side());
        heat.add(frame);
    };
    record(&frame, &mut pop, &mut boxes, &mut heat);
    history.insert(
        frame.grid.hash() ^ u128::from(background),
        (0, frame.origin.0, frame.origin.1),
    );
    let stop = if walls.hit(&frame) {
        Stop::Wall
    } else {
        loop {
            if gen >= setup.x {
                break Stop::Limit;
            }
            let (grown, counts) = setup.kernel.count_grown(&frame.grid);
            let next = Grid::of(
                grown.h,
                grown.w,
                setup.rule.apply(background, &grown.bits, &counts),
            );
            let next = Frame::cropped(
                &next,
                (frame.origin.0 - reach as i64, frame.origin.1 - reach as i64),
            );
            background = setup.rule.next_background(background);
            gen += 1;
            prev = Some(std::mem::replace(&mut frame, next));
            record(&frame, &mut pop, &mut boxes, &mut heat);
            if frame.grid.h == 0 {
                break Stop::Dead;
            }
            if walls.hit(&frame) {
                break Stop::Wall;
            }
            let key = frame.grid.hash() ^ u128::from(background);
            if let Some(&(seen, r0, c0)) = history.get(&key) {
                let period = gen - seen;
                let (dr, dc) = (frame.origin.0 - r0, frame.origin.1 - c0);
                break if (dr, dc) == (0, 0) {
                    if period == 1 {
                        Stop::Still
                    } else {
                        Stop::Loop { period }
                    }
                } else {
                    let wall_by = walls.wall_by(&frame, dr, dc, period, gen, setup.x);
                    Stop::Mover {
                        period,
                        dr,
                        dc,
                        wall_by,
                    }
                };
            }
            history.insert(key, (gen, frame.origin.0, frame.origin.1));
        }
    };
    let settled = settled_at(&stop, gen);
    Run {
        mode: Mode::Island,
        stop,
        gen,
        settled,
        pop,
        boxes,
        last: frame,
        prev,
        heat,
        background,
        flips: setup.rule.flips(),
    }
}

fn wallpaper(setup: &Setup) -> Run {
    let n = setup.board.h;
    let mut grid = setup.board.clone();
    let mut heat = Heat::new(n, n, (0, 0));
    let mut history: HashMap<u128, usize> = HashMap::new();
    let (mut pop, mut boxes) = (Vec::new(), Vec::new());
    let mut background = 0u8;
    let mut prev = None;
    let mut gen = 0;
    let whole = |grid: &Grid| Frame {
        grid: grid.clone(),
        origin: (0, 0),
    };
    let record = |grid: &Grid, pop: &mut Vec<u32>, boxes: &mut Vec<u32>, heat: &mut Heat| {
        pop.push(grid.population() as u32);
        boxes.push(grid.crop().h.max(grid.crop().w) as u32);
        heat.add(&Frame {
            grid: grid.clone(),
            origin: (0, 0),
        });
    };
    record(&grid, &mut pop, &mut boxes, &mut heat);
    history.insert(grid.hash() ^ u128::from(background), 0);
    let stop = loop {
        if gen >= setup.x {
            break Stop::Limit;
        }
        let next = step(&grid, &setup.kernel, &setup.rule, background, true);
        background = setup.rule.next_background(background);
        gen += 1;
        prev = Some(whole(&grid));
        grid = next;
        record(&grid, &mut pop, &mut boxes, &mut heat);
        if grid.is_empty() {
            break Stop::Dead;
        }
        let key = grid.hash() ^ u128::from(background);
        if let Some(&seen) = history.get(&key) {
            let period = gen - seen;
            break if period == 1 {
                Stop::Still
            } else {
                Stop::Loop { period }
            };
        }
        history.insert(key, gen);
    };
    let settled = settled_at(&stop, gen);
    Run {
        mode: Mode::Wallpaper,
        stop,
        gen,
        settled,
        pop,
        boxes,
        last: whole(&grid),
        prev,
        heat,
        background,
        flips: setup.rule.flips(),
    }
}

fn settled_at(stop: &Stop, gen: usize) -> Option<usize> {
    match stop {
        Stop::Dead => Some(gen),
        Stop::Still | Stop::Loop { .. } | Stop::Mover { .. } => {
            Some(gen - stop.period().unwrap_or(0))
        }
        Stop::Wall | Stop::Limit => None,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fate {
    Dead,
    Alive,
    Loop,
    Timeout,
}

pub fn torus(
    board: &Grid,
    kernel: &Kernel,
    rule: &Rule,
    canvas: usize,
    generations: usize,
) -> (Fate, usize) {
    assert!(
        !rule.flips(),
        "the torus control runs rules without birth on zero"
    );
    let mut current = board.centred(canvas);
    let mut history: HashMap<Vec<u8>, usize> = HashMap::new();
    history.insert(current.bits.clone(), 0);
    for i in 1..generations {
        let next = step(&current, kernel, rule, 0, true);
        if next == current {
            return (
                if next.is_empty() {
                    Fate::Dead
                } else {
                    Fate::Alive
                },
                1,
            );
        }
        if let Some(&seen) = history.get(&next.bits) {
            return (Fate::Loop, i - seen);
        }
        history.insert(next.bits.clone(), i);
        current = next;
    }
    (Fate::Timeout, 0)
}
