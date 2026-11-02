use crate::engine::{run, step, torus, Rule, Setup, Stop};
use crate::grid::Grid;
use crate::mask::{Factor, Kernel};
use crate::readings::readings;
use crate::recipe::{Mode, Recipe};
use crate::{run_named, Grammar, Shape};
use mrly_core::Rng;
use mrly_life::{next_grid, Boundary};

fn moore() -> Kernel {
    Kernel::from_grid(&Grid::from_rows(&["111", "101", "111"])).unwrap()
}

fn life(board: Grid, cap: usize, x: usize) -> Setup {
    Setup {
        board,
        kernel: moore(),
        rule: Rule::new(&[3], &[2, 3], 8),
        mode: Mode::Island,
        x,
        cap,
    }
}

fn noise(h: usize, w: usize, seed: u64, density: f64) -> Grid {
    let mut rng = Rng::new(seed);
    Grid::of(
        h,
        w,
        (0..h * w).map(|_| u8::from(rng.chance(density))).collect(),
    )
}

fn reference(grid: &Grid, birth: &[usize], survive: &[usize], mask: &Grid, wrap: bool) -> Grid {
    let boundary = if wrap {
        Boundary::Wrap
    } else {
        Boundary::Constant
    };
    let next = next_grid(
        &grid.cell().unwrap(),
        birth,
        survive,
        &mask.tensor().unwrap(),
        boundary,
    )
    .unwrap();
    Grid::from_cell(&next).unwrap()
}

fn pentomino() -> Grid {
    Grid::from_rows(&[".##", "##.", ".#."])
}

#[test]
fn the_blinker_loops_with_period_two_from_generation_zero() {
    let out = run(&life(Grid::from_rows(&["...", "###", "..."]), 31, 64));
    assert_eq!(out.stop, Stop::Loop { period: 2 });
    assert_eq!(out.gen, 2);
    assert_eq!(out.settled, Some(0));
    assert_eq!(out.pop, vec![3, 3, 3]);
    assert_eq!(out.last.grid, Grid::from_rows(&["###"]));
    let json = readings(&out).unwrap();
    assert_eq!(json["stop"]["period"], 2);
    assert_eq!(json["shape"]["symmetry"], 4);
}

#[test]
fn the_glider_moves_one_cell_down_and_right_every_four_generations() {
    let out = run(&life(Grid::from_rows(&[".#.", "..#", "###"]), 41, 1024));
    assert_eq!(
        out.stop,
        Stop::Mover {
            period: 4,
            dr: 1,
            dc: 1,
            wall_by: Some(76),
        }
    );
    assert_eq!(out.gen, 4);
    assert_eq!(out.last.origin, (1, 1));
    assert_eq!(out.pop, vec![5; 5]);
}

#[test]
fn the_r_pentomino_reaches_a_close_wall_and_settles_at_116_cells_by_generation_1103() {
    let close = run(&life(pentomino(), 41, 2000));
    assert_eq!(close.stop, Stop::Wall);
    assert!(close.gen < 100, "wall at {}", close.gen);
    let wide = run(&life(pentomino(), 901, 1200));
    assert_eq!(wide.stop, Stop::Limit);
    assert_eq!(wide.pop[1103], 116);
    assert!(wide.pop[1103..].iter().all(|&p| p == 116));
    assert!(wide.pop[1000] != 116 || wide.pop[500] != 116);
}

#[test]
fn an_island_agrees_with_a_much_larger_constant_canvas_frame_for_frame() {
    let big = 301;
    let mask = moore().dense();
    let mut frame = pentomino().centred(big);
    let mut checks = 0;
    for gen in 1..=80 {
        frame = reference(&frame, &[3], &[2, 3], &mask, false);
        if [1, 5, 17, 40, 63, 80].contains(&gen) {
            let out = run(&life(pentomino(), 201, gen));
            assert_eq!(out.stop, Stop::Limit);
            let origin = (big as i64 - 3) / 2;
            let placed = out.last.grid.embed(
                big,
                big,
                (origin + out.last.origin.0) as usize,
                (origin + out.last.origin.1) as usize,
            );
            assert_eq!(placed, frame, "generation {gen}");
            checks += 1;
        }
    }
    assert_eq!(checks, 6);
}

#[test]
fn a_wallpaper_is_the_seed_torus_under_a_mask_wider_than_the_torus() {
    let seed = Shape::parse("carpet").unwrap().build(3, 2).unwrap();
    let kernel = Kernel::design(Shape::parse("net").unwrap(), 5, 2).unwrap();
    assert_eq!(kernel.side, 25);
    let (birth, survive) = (vec![3, 5, 7], vec![2, 3, 4, 9]);
    let setup = Setup {
        board: seed.clone(),
        kernel: kernel.clone(),
        rule: Rule::new(&birth, &survive, kernel.budget),
        mode: Mode::Wallpaper,
        x: 40,
        cap: 0,
    };
    let out = run(&setup);
    let mut frame = seed;
    let mask = kernel.dense();
    for _ in 0..out.gen {
        frame = reference(&frame, &birth, &survive, &mask, true);
    }
    assert_eq!(out.last.grid, frame);
    assert!(out.gen <= 40);
}

#[test]
fn the_sublattice_lemma_holds_on_the_plane_and_on_the_torus() {
    let b = Grid::from_rows(&["111", "101", "111"]);
    let centre = Grid::from_rows(&["000", "010", "000"]);
    let mut rng = Rng::new(3);
    let mut checked = 0;
    for trial in 0..12 {
        let p = noise(5, 5, 100 + trial, 0.5);
        let budget = p.popped().population();
        if budget == 0 {
            continue;
        }
        let mask = Kernel::of(
            vec![(1, vec![Factor::of(&p), Factor::of(&centre)])],
            p.get(2, 2),
        )
        .unwrap();
        let quotient = Kernel::from_grid(&p.popped()).unwrap();
        assert_eq!(mask.budget, budget);
        let birth: Vec<usize> = (1..=budget).filter(|_| rng.chance(0.4)).collect();
        let survive: Vec<usize> = (0..=budget).filter(|_| rng.chance(0.4)).collect();
        let rule = Rule::new(&birth, &survive, budget);
        let a = noise(7, 7, 200 + trial, 0.5);
        let (grown, counts) = mask.count_grown(&a.kron(&b));
        let full = Grid::of(grown.h, grown.w, rule.apply(0, &grown.bits, &counts)).crop();
        let (grown_a, counts_a) = quotient.count_grown(&a);
        let small = Grid::of(
            grown_a.h,
            grown_a.w,
            rule.apply(0, &grown_a.bits, &counts_a),
        )
        .crop();
        assert_eq!(full, small.kron(&b), "plane trial {trial}");
        let a3 = noise(3, 3, 300 + trial, 0.6);
        let full = step(&a3.kron(&b), &mask, &rule, 0, true);
        let small = step(&a3, &quotient, &rule, 0, true);
        assert_eq!(full, small.kron(&b), "torus trial {trial}");
        checked += 1;
    }
    assert!(checked >= 10);
}

#[test]
fn a_zero_in_birth_runs_as_the_rule_pair_on_a_flipping_background() {
    let seed = Grid::from_rows(&[".##", "##.", ".#."]);
    let mask = moore().dense();
    for (birth, survive) in [
        (vec![0, 3], vec![2, 3]),
        (vec![0, 1], vec![8]),
        (vec![0, 2], vec![0, 8]),
    ] {
        let rule = Rule::new(&birth, &survive, 8);
        assert!(rule.flips());
        let x = 12;
        let setup = Setup {
            board: seed.clone(),
            kernel: moore(),
            rule: rule.clone(),
            mode: Mode::Island,
            x,
            cap: 101,
        };
        let out = run(&setup);
        assert_eq!(
            out.gen, x,
            "{birth:?}/{survive:?} stopped as {:?}",
            out.stop
        );
        let big = 101;
        let mut frame = seed.centred(big);
        for _ in 0..x {
            frame = reference(&frame, &birth, &survive, &mask, true);
        }
        let origin = (big as i64 - 3) / 2;
        let placed = out.last.grid.embed(
            big,
            big,
            (origin + out.last.origin.0) as usize,
            (origin + out.last.origin.1) as usize,
        );
        let unflipped = Grid::of(
            big,
            big,
            placed.bits.iter().map(|&y| y ^ out.background).collect(),
        );
        assert_eq!(unflipped, frame, "{birth:?}/{survive:?}");
        let paper = Setup {
            board: seed.tile(3),
            kernel: moore(),
            rule: rule.clone(),
            mode: Mode::Wallpaper,
            x,
            cap: 0,
        };
        let out = run(&paper);
        let mut frame = seed.tile(3);
        for _ in 0..out.gen {
            frame = reference(&frame, &birth, &survive, &mask, true);
        }
        let unflipped = Grid::of(
            9,
            9,
            out.last
                .grid
                .bits
                .iter()
                .map(|&y| y ^ out.background)
                .collect(),
        );
        assert_eq!(unflipped, frame, "wallpaper {birth:?}/{survive:?}");
    }
}

#[test]
fn the_torus_control_matches_the_engine_on_a_padded_canvas() {
    let board = Grid::from_rows(&["111", "101", "111"]).tile(3);
    let kernel = moore();
    let rule = Rule::new(&[3], &[2, 3], 8);
    let mut frames = vec![board.centred(27)];
    let mut want = ("Timeout".to_string(), 0);
    for i in 1..64 {
        let next = reference(&frames[i - 1], &[3], &[2, 3], &kernel.dense(), true);
        if next == frames[i - 1] {
            want = (
                if next.is_empty() { "Dead" } else { "Alive" }.to_string(),
                1,
            );
            break;
        }
        if let Some(seen) = frames.iter().position(|f| *f == next) {
            want = ("Loop".to_string(), i - seen);
            break;
        }
        frames.push(next);
    }
    let (fate, period) = torus(&board, &kernel, &rule, 27, 64);
    assert_eq!((format!("{fate:?}"), period), want);
    assert_ne!(want.0, "Timeout");
}

#[test]
fn one_recipe_of_every_class_runs_from_its_name_and_reads() {
    let mut grammar = Grammar::standard(32, 61);
    grammar.sides = vec![3];
    grammar.levels = vec![1];
    grammar.mask_sides = vec![3, 5];
    grammar.mask_levels = vec![1];
    let all = grammar.recipes(0, usize::MAX).unwrap();
    let mut seen = std::collections::BTreeSet::new();
    for recipe in &all {
        let class = (
            recipe.mode.token(),
            recipe.tiling,
            matches!(recipe.mask, crate::MaskSpec::Copy { .. }),
        );
        if !seen.insert(class) {
            continue;
        }
        let (parsed, out) = run_named(&recipe.name()).unwrap();
        assert_eq!(parsed, *recipe);
        let json = readings(&out).unwrap();
        assert_eq!(json["mode"], recipe.mode.token());
        assert!(json["pop"].as_array().unwrap().len() == out.gen + 1);
    }
    assert!(seen.len() >= 6, "{seen:?}");
    let copy =
        Recipe::parse("carpet.3.1 x3 island mask:copyinv birth:3 survive:2,3 x:64 cap:61").unwrap();
    let setup = Setup::of(&copy).unwrap();
    assert_eq!(setup.kernel.side, 9);
    assert_eq!(setup.kernel.budget, 8);
}
