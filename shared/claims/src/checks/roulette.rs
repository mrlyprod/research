use mrly_math::roulette::{nodes, spread, Nodes};
use mrly_math::spirograph::{pencils, track, Pencil};

const TOL: f64 = 4e-4;

const SAMPLES: usize = 4000;

const CARPET: [u8; 9] = [1, 1, 1, 1, 0, 1, 1, 1, 1];

const TWO: [u8; 9] = [1, 0, 0, 0, 0, 0, 0, 1, 0];

// THE SEATS

fn pin(target: f64) -> Vec<Pencil> {
    let reach = target * 1.5f64.hypot(0.5);
    pencils(&[1, 0, 0], 3, 1, "fill", reach, 0.0, 1).expect("the pin seats")
}

fn tile(types: &[u8], target: f64) -> Vec<Pencil> {
    let reach = target * 1.5f64.hypot(1.5) / 2f64.sqrt();
    pencils(types, 3, 3, "fill", reach, 0.0, 1).expect("the tile seats")
}

fn count(
    kind: &str,
    a: usize,
    b: usize,
    seats: &[Pencil],
    samples: usize,
) -> Result<Nodes, String> {
    let path = track(kind, a, b, 4, 1).expect("the track");
    let seats = spread(&path, seats, true);
    let edge = if kind == "in" {
        (1.0f64).min((a - b) as f64 / b as f64)
    } else {
        1.0
    };
    for seat in &seats {
        let reach = seat.x.hypot(seat.y);
        if reach <= 0.0 || reach >= edge {
            return Err(format!(
                "{kind} {a}/{b} seats a pencil at {reach}, outside the law's 0 to {edge}"
            ));
        }
    }
    Ok(nodes(&path, &seats, samples, TOL).expect("the nodes"))
}

// THE LAWS

pub fn one_curve_of_a_circle_roulette() -> Result<(), String> {
    for (kind, a, b, target) in [
        ("in", 5, 2, 0.5),
        ("in", 7, 3, 0.9),
        ("in", 11, 4, 0.5),
        ("out", 5, 2, 0.9),
        ("out", 7, 3, 0.5),
    ] {
        let selved = count(kind, a, b, &pin(target), SAMPLES)?.selved();
        let want = a * (b - 1);
        if selved != want {
            return Err(format!(
                "{kind} {a}/{b} at a seat of {target} crosses itself {selved} times against {want}"
            ));
        }
    }
    Ok(())
}

pub fn two_distinct_curves_of_one_wheel() -> Result<(), String> {
    let count = count("in", 7, 4, &tile(&TWO, 0.6), 6000)?;
    if (count.selved(), count.paired()) != (42, 56) {
        return Err(format!(
            "the two curves of 7/4 cross {} times against 56 and themselves {} against 42",
            count.paired(),
            count.selved()
        ));
    }
    Ok(())
}

pub fn a_whole_design_carries_one_node() -> Result<(), String> {
    let count = count("in", 7, 3, &tile(&CARPET, 0.4), SAMPLES)?;
    let (a, b, k) = (7, 3, 8);
    let want = (k * a * (b - 1), a * b * k * (k - 1));
    if (count.selved(), count.paired()) != want {
        return Err(format!(
            "the carpet on 7/3 carries {} and {} against {} and {}",
            count.selved(),
            count.paired(),
            want.0,
            want.1
        ));
    }
    Ok(())
}

pub fn a_roulette_cuts_the_plane_into() -> Result<(), String> {
    let cells: [(usize, usize, Vec<Pencil>, usize); 5] = [
        (3, 1, pin(0.5), 2),
        (5, 2, pin(0.5), 7),
        (7, 3, pin(0.5), 16),
        (3, 1, tile(&TWO, 0.5), 8),
        (5, 2, tile(&TWO, 0.5), 32),
    ];
    for (a, b, seats, want) in cells {
        let count = count("in", a, b, &seats, SAMPLES)?;
        let faces = count.branches - count.points + 2;
        if faces != want || count.crowded > 0 {
            return Err(format!(
                "in {a}/{b} cuts the plane into {faces} regions against {want}, {} nodes crowded",
                count.crowded
            ));
        }
    }
    Ok(())
}

// THE ENDS

pub fn the_loop_threshold_is_not_where() -> Result<(), String> {
    let path = track("in", 7, 6, 4, 1).expect("the track");
    let selved = nodes(&path, &pin(0.9), 6000, TOL)
        .expect("the nodes")
        .selved();
    if selved != 7 {
        return Err(format!(
            "7/6 at a seat of 0.9 crosses itself {selved} times against 7, the law saying 35"
        ));
    }
    let path = track("in", 7, 4, 4, 1).expect("the track");
    let count = nodes(&path, &spread(&path, &tile(&TWO, 0.9), true), 6000, TOL).expect("the nodes");
    if (count.paired(), count.selves.clone()) != (42, vec![21, 21]) {
        return Err(format!(
            "7/4 at seats of 0.9 and 0.636 crosses {} times against 42, the law saying 56, with selves {:?}",
            count.paired(),
            count.selves
        ));
    }
    Ok(())
}

pub fn the_self_law_ends_at_the() -> Result<(), String> {
    let path = track("in", 6, 5, 4, 1).expect("the track");
    let edge = 1.0 / 5.0;
    for (level, want) in [(0.9, 24), (1.02, 18), (1.3, 6)] {
        let selved = nodes(&path, &pin(level * edge), 6000, TOL)
            .expect("the nodes")
            .selved();
        if selved != want {
            return Err(format!(
                "6/5 at a level of {level} crosses itself {selved} times against {want}"
            ));
        }
    }
    Ok(())
}
