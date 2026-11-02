use std::f64::consts::PI;

// SETUP

const GAMMA: f64 = 0.577_215_664_901_532_9;
const GUARD: f64 = 1e-12;
const SCI_GUARD: f64 = 1e-10;
const BASES: [u64; 10] = [
    1_000,
    2_000,
    3_000,
    3_689,
    3_690,
    5_000,
    10_000,
    100_000,
    1_000_000,
    1_000_000_000,
];
const SWEEP_LO: u64 = 3_690;
const SWEEP_HI: u64 = 100_000;
const CORO: [u64; 3] = [10_000, 100_000, 1_000_000];
const COST: [u64; 7] = [
    3_689,
    3_690,
    5_000,
    10_000,
    100_000,
    1_000_000,
    1_000_000_000,
];
const WALL: [u64; 2] = [3_689, 3_690];
const SHARP_HI: u64 = 4_000_000;
const SHARP_COST_HI: u64 = 100_000;

// KERNEL CONSTANTS

fn harmonic_bound(n: f64) -> f64 {
    n.ln() + GAMMA + 0.5 / n
}

fn phi(q: f64) -> f64 {
    let n = ((q - 2.0) / 2.0).ceil();
    (4.0 / PI) * q + (2.0 * q / PI) * harmonic_bound(n) + (1.0 - 2.0 / PI) * (q - 2.0) + 0.727
}

fn pb(q: f64, m: f64) -> f64 {
    m.sqrt() + phi(q) / q
}

fn alpha(q: f64, m: f64) -> f64 {
    (q - m).ln() / q.ln()
}

fn c_exp(q: f64, m: f64) -> f64 {
    pb(q, m).ln() / q.ln()
}

fn delta(q: f64, m: f64) -> f64 {
    let a = alpha(q, m);
    (a - 0.75 - c_exp(q, m)) / a
}

fn delta_stable(q: f64, m: f64) -> f64 {
    (gap(q, m) / pb(q, m)).ln_1p() / (q.ln() * alpha(q, m))
}

fn defect(q: f64, m: f64) -> f64 {
    m / (2.0 * (q - m) * q.ln())
}

fn psi(q: f64) -> f64 {
    let n = ((q - 2.0) / 2.0).ceil();
    (2.0 * q / PI) * harmonic_bound(n) + (1.0 - 2.0 / PI) * q
}

fn psi_chord(q: f64) -> f64 {
    let p = (q / 2.0).floor();
    let h = harmonic_bound(p - 1.0);
    if (q as u64).is_multiple_of(2) {
        (q / PI) * (2.0 * h - 1.0 + 1.0 / p) + (1.0 - 2.0 / PI) * q / 2.0
    } else {
        (q / PI) * (2.0 * h - 1.0 + 2.0 / p) + (1.0 - 2.0 / PI) * (q / 2.0 + 1.0 / (2.0 * q))
    }
}

fn sec_off(q: f64, extreme: bool) -> f64 {
    let c = if extreme {
        (q - 1.0) / 2.0
    } else if (q as u64) % 2 == 1 {
        0.0
    } else {
        0.5
    };
    1.0 / (PI * c / q).cos()
}

#[derive(Clone, Copy)]
enum Cst {
    Step3,
    Sharp(bool),
    Chord(bool),
}

fn onestep(q: f64, c: Cst) -> f64 {
    match c {
        Cst::Step3 => 1.0 + phi(q) / q,
        Cst::Sharp(e) => 4.0 / PI + psi(q) / q + 0.5 - sec_off(q, e) / (2.0 * q),
        Cst::Chord(e) => 4.0 / PI + psi_chord(q) / q + 0.5 - sec_off(q, e) / (2.0 * q),
    }
}

fn onestep_floor(c: Cst) -> u64 {
    match c {
        Cst::Step3 => 3,
        Cst::Sharp(_) => 17,
        Cst::Chord(_) => 36,
    }
}

fn onestep_delta(q: f64, c: Cst) -> f64 {
    let a = alpha(q, 1.0);
    (a - 0.75 - onestep(q, c).ln() / q.ln()) / a
}

fn onestep_wall(c: Cst, hi: u64) -> (u64, bool) {
    let (mut first, mut held) = (0u64, 0u64);
    for q in onestep_floor(c)..=hi {
        let qf = q as f64;
        if (qf - 1.0) * qf.powf(-0.75) > onestep(qf, c) {
            held += 1;
            if first == 0 {
                first = q;
            }
        }
    }
    (first, first > 0 && held == hi - first + 1)
}

fn onestep_cross(c: Cst, hi: u64) -> (u64, bool, f64, f64) {
    let (mut first, mut held) = (0u64, 0u64);
    for q in onestep_floor(c)..=hi {
        let qf = q as f64;
        if onestep_delta(qf, c) > defect(qf, 1.0) {
            held += 1;
            if first == 0 {
                first = q;
            }
        }
    }
    let qf = first as f64;
    (
        first,
        first > 0 && held == hi - first + 1,
        onestep_delta(qf, c),
        defect(qf, 1.0),
    )
}

fn closes(q: f64, m: f64) -> bool {
    0.75 + c_exp(q, m) < alpha(q, m)
}

fn gap_b(q: f64, b: f64, m: f64) -> f64 {
    (q - m) * q.powf(-b) - pb(q, m)
}

fn gap(q: f64, m: f64) -> f64 {
    gap_b(q, 0.75, m)
}

// SAFE ROUNDING

fn ceil_units(x: f64, digits: u32) -> i64 {
    let s = 10f64.powi(digits as i32);
    ((x + GUARD) * s).ceil() as i64
}

fn floor_units(x: f64, digits: u32) -> i64 {
    let s = 10f64.powi(digits as i32);
    ((x - GUARD) * s).floor() as i64
}

fn fixed(units: i64, digits: u32) -> String {
    let s = 10i64.pow(digits);
    let sign = if units < 0 { "-" } else { "" };
    let a = units.abs();
    format!(
        "{}{}.{:0width$}",
        sign,
        a / s,
        a % s,
        width = digits as usize
    )
}

fn sci(x: f64, decimals: u32, up: bool) -> String {
    let neg = x < 0.0;
    let a = x.abs();
    assert!(a > 0.0 && a.is_finite(), "sci needs a nonzero finite value");
    let mag_up = up != neg;
    let b = if mag_up {
        a * (1.0 + SCI_GUARD)
    } else {
        a * (1.0 - SCI_GUARD)
    };
    let mut e = b.log10().floor() as i32;
    let mut u = {
        let s = 10f64.powi(decimals as i32 - e);
        if mag_up {
            (b * s).ceil() as i64
        } else {
            (b * s).floor() as i64
        }
    };
    let lo = 10i64.pow(decimals);
    if u >= 10 * lo {
        u /= 10;
        e += 1;
    }
    if u < lo {
        u *= 10;
        e -= 1;
    }
    format!("{}{}e{}", if neg { "-" } else { "" }, fixed(u, decimals), e)
}

fn label(q: u64) -> String {
    let mut p = q;
    let mut e = 0u32;
    while p.is_multiple_of(10) {
        p /= 10;
        e += 1;
    }
    if p == 1 && e >= 4 {
        format!("10^{e}")
    } else {
        q.to_string()
    }
}

// TABLE ROWS

struct Row {
    q: u64,
    m: u64,
    alpha: i64,
    c: i64,
    delta: i64,
    closes: bool,
}

fn row(q: u64, m: u64) -> Row {
    let (qf, mf) = (q as f64, m as f64);
    Row {
        q,
        m,
        alpha: floor_units(alpha(qf, mf), 6),
        c: ceil_units(c_exp(qf, mf), 5),
        delta: floor_units(delta(qf, mf), 5),
        closes: closes(qf, mf),
    }
}

fn render(r: &Row) -> String {
    format!(
        "| {} | {} | {} | {} | {} |",
        label(r.q),
        fixed(r.alpha, 6),
        fixed(r.c, 5),
        fixed(r.delta, 5),
        if r.closes { "yes" } else { "no" }
    )
}

fn render_cost(q: u64, m: u64) -> String {
    let (qf, mf) = (q as f64, m as f64);
    let d = delta_stable(qf, mf);
    let f = defect(qf, mf);
    format!(
        "| {} | {} | {} | {} | {} |",
        label(q),
        m,
        sci(d, 5, false),
        sci(f, 5, true),
        if d > f { "yes" } else { "no" }
    )
}

fn render_sharp(name: &str, digits: &str, c: Cst) -> String {
    let (w, w_up) = onestep_wall(c, SHARP_HI);
    let (x, x_up, d, f) = onestep_cross(c, SHARP_COST_HI);
    format!(
        "| {} | {} | {} | {} | {} | {} | {} |",
        name,
        digits,
        label(w),
        label(x),
        sci(d, 5, false),
        sci(f, 5, true),
        if w_up && x_up { "yes" } else { "no" }
    )
}

fn render_coro(r: &Row) -> String {
    format!(
        "| {} | {} | {} | {} | {} |",
        label(r.q),
        r.m,
        fixed(r.alpha, 6),
        fixed(r.c, 5),
        fixed(r.delta, 5)
    )
}

// CHECKS

fn gap_sweep(lo: u64, hi: u64, m: u64) -> (u64, f64, f64) {
    let mf = m as f64;
    let mut prev = gap(lo as f64, mf);
    let mut worst = f64::INFINITY;
    let mut arg = lo;
    let mut min_gap = prev;
    for q in (lo + 1)..=hi {
        let g = gap(q as f64, mf);
        let step = g - prev;
        if step < worst {
            worst = step;
            arg = q - 1;
        }
        if g < min_gap {
            min_gap = g;
        }
        prev = g;
    }
    (arg, worst, min_gap)
}

fn crossover(lo: u64, hi: u64) -> (u64, u64) {
    let mut first = 0u64;
    let mut down = 0u64;
    let mut prev = f64::NEG_INFINITY;
    for q in lo..=hi {
        let qf = q as f64;
        let d = delta_stable(qf, 1.0) - defect(qf, 1.0);
        if first == 0 && d > 0.0 {
            first = q;
        }
        if q > lo && d <= prev {
            down += 1;
        }
        prev = d;
    }
    (first, down)
}

fn max_excluded(q: u64) -> u64 {
    let qf = q as f64;
    let mut m = 0u64;
    while gap(qf, (m + 1) as f64) > 0.0 {
        m += 1;
    }
    m
}

// LADDER

const LADDER: [(i64, i64); 10] = [
    (1, 2),
    (13, 25),
    (11, 20),
    (4, 7),
    (3, 5),
    (2, 3),
    (3, 4),
    (4, 5),
    (9, 10),
    (19, 20),
];
const MONO_C: f64 = 1.291;
const MONO_LO: f64 = 40.0;
const EXACT_ULPS: f64 = 1024.0;
const DYADIC_TOP: f64 = 9_007_199_254_740_992.0;
const CORO_Q: f64 = 1e7;

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.abs()
}

#[derive(Clone, Copy)]
struct Rat {
    n: i64,
    d: i64,
}

impl Rat {
    fn new(n: i64, d: i64) -> Rat {
        let g = gcd(n, d);
        Rat { n: n / g, d: d / g }
    }

    fn val(self) -> f64 {
        self.n as f64 / self.d as f64
    }

    fn lt(self, o: Rat) -> bool {
        self.n * o.d < o.n * self.d
    }

    fn same(self, o: Rat) -> bool {
        self.n * o.d == o.n * self.d
    }

    fn show(self) -> String {
        if self.d == 1 {
            self.n.to_string()
        } else {
            format!("{}/{}", self.n, self.d)
        }
    }
}

fn zhang(a: Rat) -> Option<Rat> {
    if a.lt(Rat::new(1, 2)) || Rat::new(4, 7).lt(a) {
        return None;
    }
    Some(Rat::new(
        a.n * (8 * a.d - 7 * a.n),
        2 * a.d * (2 * a.d - a.n),
    ))
}

fn baker_harman(a: Rat) -> Rat {
    if a.lt(Rat::new(11, 20)) {
        Rat::new(4 * a.n + a.d, 4 * a.d)
    } else if a.lt(Rat::new(3, 5)) {
        Rat::new(4, 5)
    } else {
        Rat::new(a.n + a.d, 2 * a.d)
    }
}

fn ladder_b(a: Rat) -> (Rat, &'static str) {
    let h = baker_harman(a);
    match zhang(a) {
        None => (h, "BH"),
        Some(z) if z.lt(h) => (z, "Zhang"),
        Some(z) if h.lt(z) => (h, "BH"),
        Some(z) => (z, "both"),
    }
}

fn mono_ok(q: f64, b: f64) -> bool {
    (1.0 - b) * (q - 2.0) * (q + 1.0).powf(-b) >= MONO_C
}

fn mono_floor(b: f64) -> f64 {
    let mut hi = MONO_LO;
    while !mono_ok(hi, b) {
        hi *= 2.0;
    }
    let mut lo = MONO_LO;
    for _ in 0..400 {
        let mid = lo + (hi - lo) / 2.0;
        if mid <= lo || mid >= hi {
            break;
        }
        if mono_ok(mid, b) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi.ceil()
}

fn pb_low(q: f64) -> f64 {
    1.0 + 4.0 / PI
        + (2.0 / PI) * (((q - 2.0) / 2.0).ln() + GAMMA)
        + (1.0 - 2.0 / PI) * (q - 2.0) / q
}

fn u_bound(q: f64, b: f64) -> f64 {
    q.powf(1.0 - b) - pb_low(q)
}

fn ladder_wall(b: f64) -> f64 {
    let floor = mono_floor(b).max(SWEEP_LO as f64);
    if gap_b(floor, b, 1.0) > 0.0 {
        return floor;
    }
    let mut hi = floor;
    while gap_b(hi, b, 1.0) <= 0.0 {
        hi *= 2.0;
    }
    let mut lo = floor;
    if hi <= DYADIC_TOP {
        while hi - lo > 1.0 {
            let mid = (lo + (hi - lo) / 2.0).floor();
            if gap_b(mid, b, 1.0) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        return hi;
    }
    for _ in 0..400 {
        let mid = lo + (hi - lo) / 2.0;
        if mid <= lo || mid >= hi {
            break;
        }
        if gap_b(mid, b, 1.0) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

fn ladder_exact(q: f64, b: f64) -> bool {
    if q >= DYADIC_TOP {
        return false;
    }
    let noise = EXACT_ULPS * f64::EPSILON * ((q - 1.0) * q.powf(-b)).max(pb(q, 1.0));
    gap_b(q, b, 1.0) > noise && gap_b(q - 1.0, b, 1.0) < -noise
}

fn integer_or_bound(q: f64, exact: bool) -> String {
    if exact {
        format!("{}", q as u64)
    } else {
        format!("<= {}", sci(q, 5, true))
    }
}

fn max_excluded_b(q: f64, b: f64) -> u64 {
    let mut m = 0u64;
    while gap_b(q, b, (m + 1) as f64) > 0.0 {
        m += 1;
    }
    m
}

fn render_ladder(a: Rat) -> String {
    let (bq, src) = ladder_b(a);
    let b = bq.val();
    let floor = mono_floor(b);
    let q0 = ladder_wall(b);
    format!(
        "| {} | {} | {} | {} | {} |",
        a.show(),
        bq.show(),
        src,
        integer_or_bound(q0, ladder_exact(q0, b)),
        integer_or_bound(floor, floor < DYADIC_TOP)
    )
}

fn ladder_trend() -> String {
    let mut parts: Vec<String> = Vec::new();
    for (n, d) in LADDER {
        let b = ladder_b(Rat::new(n, d)).0.val();
        parts.push(fixed(ceil_units(ladder_wall(b).log10(), 2), 2));
    }
    parts.join(" ")
}

fn render_ladder_coro(a: Rat) -> String {
    let (bq, src) = ladder_b(a);
    format!(
        "| {} | {} | {} | {} |",
        a.show(),
        bq.show(),
        src,
        max_excluded_b(CORO_Q, bq.val())
    )
}

// MAIN

fn main() {
    println!("mertens-numerology");
    println!(
        "guard {GUARD:e}; alpha_q truncated down at 6 digits, c_q rounded up at 5, delta_q rounded down at 5"
    );
    println!();
    println!("TABLE, one excluded digit");
    println!("| `q` | `alpha_q` | `c_q` (proved, up) | `delta_q` (down) | closes |");
    println!("|---|---|---|---|---|");
    for q in BASES {
        let r = row(q, 1);
        assert!(c_exp(q as f64, 1.0) >= 0.0, "l^1 floor broken at q = {q}");
        assert_eq!(
            r.closes,
            gap(q as f64, 1.0) > 0.0,
            "tests disagree at q = {q}"
        );
        println!("{}", render(&r));
    }
    println!();
    println!("WALL, one excluded digit");
    assert!(!closes(3689.0, 1.0), "3689 must not close");
    assert!(closes(3690.0, 1.0), "3690 must close");
    println!(
        "- closes fails at `q = 3689`, holds at `q = 3690`; gap `<= {}` and `>= {}`",
        fixed(ceil_units(gap(3689.0, 1.0), 6), 6),
        fixed(floor_units(gap(3690.0, 1.0), 6), 6)
    );
    let (arg, worst, min_gap) = gap_sweep(SWEEP_LO, SWEEP_HI, 1);
    assert!(worst > 0.0, "gap not stepwise increasing");
    assert!(min_gap > 0.0, "gap not positive on the sweep");
    println!(
        "- gap `(q-1) q^(-3/4) - PB_q(1)` steps up at all {} steps of `{SWEEP_LO}..{SWEEP_HI}`; smallest step `>= {}` at `q = {arg}`",
        SWEEP_HI - SWEEP_LO,
        fixed(floor_units(worst, 8), 8)
    );
    println!();
    println!("MARGIN AT THE WALL, one excluded digit, 10 significant digits");
    println!(
        "- relative guard {SCI_GUARD:e}; `delta_q` from the cancellation-free form `ln(1 + gap_q/PB_q) / (alpha_q ln q)`"
    );
    for q in WALL {
        let qf = q as f64;
        let up = !closes(qf, 1.0);
        let rel = if up { "<=" } else { ">=" };
        println!(
            "- `q = {q}`: `delta_q {rel} {}`, `gap_q(1) {rel} {}`",
            sci(delta_stable(qf, 1.0), 9, up),
            sci(gap(qf, 1.0), 9, up)
        );
    }
    println!("- the fifth-digit table columns above cannot display either sign");
    println!();
    println!("COST-OUT, `delta_q` against the defect exponent `m/(2(q-m) ln q)` of the level-`x^(alpha/2)` bound");
    println!("| `q` | `m` | `delta_q` (down) | defect (up) | saving beats defect |");
    println!("|---|---|---|---|---|");
    for q in COST {
        println!("{}", render_cost(q, 1));
    }
    for q in CORO {
        println!("{}", render_cost(q, max_excluded(q)));
    }
    let (first, down) = crossover(SWEEP_LO, SWEEP_HI);
    assert!(first > 0, "no crossover in the scan");
    println!(
        "- least `q` in `{SWEEP_LO}..{SWEEP_HI}` with `delta_q > 1/(2(q-1) ln q)`: `q = {first}`"
    );
    println!(
        "- the difference `delta_q - 1/(2(q-1) ln q)` rises at every one of the {} steps of that scan ({down} exceptions); monotonicity beyond the scan is not proved",
        SWEEP_HI - SWEEP_LO
    );
    println!();
    println!("SHARPENED COST-OUT, the GRH rung re-costed at every one-step constant proved here");
    println!(
        "`1 + Phi_q/q` is step 3; `Psi_q` is the sharpening with the phase identity, `q >= 17`;"
    );
    println!("`Psi'_q` is the same with the chord kernel bound in place of `Phi_q`, `q >= 36`");
    println!("`H` is the harmonic upper bound `ln n + gamma + 1/(2n)` in both, and the `36` is that convention:");
    println!("with the harmonic number itself the monotone step's floor `Psi'_q >= (1 + pi) q/2` first holds at `q = 37`");
    println!("the wall is the least `q` with `(q-1) q^(-3/4)` above the constant, the crossing the least `q` where");
    println!("`delta_q` exceeds the defect `1/(2(q-1) ln q)`; no wall is inherited and no crossing is inherited");
    println!("| constant | excluded digit | wall `q_0(1/2)` | crossing `q` | `delta_q` there (down) | defect (up) | both up-sets |");
    println!("|---|---|---|---|---|---|---|");
    for (name, digits, c) in [
        ("`1 + Phi_q/q`", "any", Cst::Step3),
        ("`Psi_q`", "any", Cst::Sharp(false)),
        ("`Psi_q`", "`0` or `q-1`", Cst::Sharp(true)),
        ("`Psi'_q`", "any", Cst::Chord(false)),
        ("`Psi'_q`", "`0` or `q-1`", Cst::Chord(true)),
    ] {
        println!("{}", render_sharp(name, digits, c));
    }
    println!(
        "- each row is scanned exhaustively from its own floor, `q >= 3`, `17` and `36`, to `{SHARP_HI}` for the wall and `{SHARP_COST_HI}` for the crossing"
    );
    println!("- every crossing sits within five steps of its own wall, so a lower wall costs out at once and is never inherited");
    println!("- the rungs above `a = 1/2` are not costed here and stay owed");
    println!();
    println!("M-COROLLARY, largest `m` with `PB_q(m) < (q-m) q^(-3/4)`");
    println!("| `q` | max `m` | `alpha_q` | `c_q` (proved, up) | `delta_q` (down) |");
    println!("|---|---|---|---|---|");
    for q in CORO {
        let m = max_excluded(q);
        let r = row(q, m);
        assert!(r.closes, "corollary maximum must close at q = {q}");
        assert!(
            !closes(q as f64, (m + 1) as f64),
            "maximum not maximal at q = {q}"
        );
        println!("{}", render_coro(&r));
    }
    println!();
    println!("LADDER, one excluded digit, `q_0(a)` = least `q >= 3` with `(q-1) q^(-b(a)) - PB_q(1) > 0`");
    println!("- `b(a)` is the smaller of Baker-Harman and Zhang where both apply; `Q(b)` is the proved monotone floor");
    println!("| `a` | `b(a)` | source | `q_0(a)` | `Q(b)` |");
    println!("|---|---|---|---|---|");
    for (n, d) in LADDER {
        let a = Rat::new(n, d);
        let b = ladder_b(a).0.val();
        let floor = mono_floor(b);
        let q0 = ladder_wall(b);
        assert!(
            floor < q0,
            "monotone floor above the wall at a = {}",
            a.show()
        );
        for q in 3..SWEEP_LO {
            assert!(
                gap_b(q as f64, b, 1.0) < 0.0,
                "close below 3690 at a = {}",
                a.show()
            );
        }
        if floor > SWEEP_LO as f64 {
            assert!(
                u_bound(SWEEP_LO as f64, b) < 0.0,
                "floor range unclear at a = {}",
                a.show()
            );
            assert!(
                u_bound(floor, b) < 0.0,
                "floor range unclear at a = {}",
                a.show()
            );
        }
        println!("{}", render_ladder(a));
    }
    assert!(
        ladder_b(Rat::new(1, 2)).0.same(Rat::new(3, 4)),
        "the GRH rung must read b = 3/4"
    );
    assert_eq!(
        ladder_wall(0.75),
        SWEEP_LO as f64,
        "the GRH rung must be the wall"
    );
    println!(
        "- `log10 q_0(a)` across the rungs, rounded up: {}",
        ladder_trend()
    );
    println!("- `gap_q(a, 1)` steps up at every `q >= Q(b)`, and `Q(b) < q_0(a)` at every rung, so it steps up from `q_0(a)` on");
    println!(
        "- no `q < {SWEEP_LO}` closes at any printed rung; `<=` marks a rung whose wall is past `2^53` or whose neighbouring gap steps fall under the {EXACT_ULPS}-ulp noise floor"
    );
    println!();
    println!("LADDER M-COROLLARY, largest `m` with `PB_q(m) < (q-m) q^(-b(a))` at `q = 10^7`");
    println!("| `a` | `b(a)` | source | max `m` |");
    println!("|---|---|---|---|");
    for (n, d) in LADDER {
        let a = Rat::new(n, d);
        let b = ladder_b(a).0.val();
        if ladder_wall(b) >= CORO_Q {
            continue;
        }
        let m = max_excluded_b(CORO_Q, b);
        assert!(
            m >= 1,
            "rung past its wall must admit m = 1 at a = {}",
            a.show()
        );
        assert!(
            gap_b(CORO_Q, b, (m + 1) as f64) <= 0.0,
            "maximum not maximal at a = {}",
            a.show()
        );
        println!("{}", render_ladder_coro(a));
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    fn kernel_sum(q: u64, t: f64) -> f64 {
        let qf = q as f64;
        let num = (PI * t).sin().abs();
        let mut s = 0.0;
        for r in 0..q {
            let x = (t + r as f64) / qf;
            let d = if x <= 0.5 { x } else { 1.0 - x };
            let den = (PI * d).sin();
            s += if den < 1e-12 { qf } else { num / den };
        }
        s
    }

    #[test]
    fn table_rows_pinned() {
        let got: Vec<String> = BASES.iter().map(|&q| render(&row(q, 1))).collect();
        let want = [
            "| 1000 | 0.999855 | 0.28087 | -0.03102 | no |",
            "| 2000 | 0.999934 | 0.26335 | -0.01342 | no |",
            "| 3000 | 0.999958 | 0.25430 | -0.00434 | no |",
            "| 3689 | 0.999966 | 0.24997 | -0.00001 | no |",
            "| 3690 | 0.999967 | 0.24997 | 0.00000 | yes |",
            "| 5000 | 0.999976 | 0.24393 | 0.00605 | yes |",
            "| 10^4 | 0.999989 | 0.23141 | 0.01858 | yes |",
            "| 10^5 | 0.999999 | 0.19906 | 0.05094 | yes |",
            "| 10^6 | 0.999999 | 0.17589 | 0.07411 | yes |",
            "| 10^9 | 0.999999 | 0.13305 | 0.11695 | yes |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn wall_sign_flip() {
        assert!(!closes(3689.0, 1.0));
        assert!(closes(3690.0, 1.0));
        assert!(gap(3689.0, 1.0) < 0.0);
        assert!(gap(3690.0, 1.0) > 0.0);
        for q in 3..3689u64 {
            assert!(
                !closes(q as f64, 1.0),
                "unexpected close below the wall at {q}"
            );
        }
    }

    #[test]
    fn gap_steps_up_and_stays_positive() {
        let (arg, worst, min_gap) = gap_sweep(SWEEP_LO, SWEEP_HI, 1);
        assert!(worst > 0.0);
        assert!(min_gap > 0.0);
        assert_eq!(arg, SWEEP_HI - 2);
        assert_eq!(fixed(floor_units(worst, 8), 8), "0.00003172");
        for q in SWEEP_LO..5_000 {
            assert!(
                gap((q + 1) as f64, 1.0) > gap(q as f64, 1.0),
                "step down at {q}"
            );
        }
    }

    #[test]
    fn closes_matches_gap_sign() {
        for q in 3..20_000u64 {
            assert_eq!(
                closes(q as f64, 1.0),
                gap(q as f64, 1.0) > 0.0,
                "at q = {q}"
            );
        }
    }

    #[test]
    fn corollary_maxima() {
        assert_eq!(max_excluded(10_000), 6);
        assert_eq!(max_excluded(100_000), 78);
        assert_eq!(max_excluded(1_000_000), 451);
        let got: Vec<String> = CORO
            .iter()
            .map(|&q| render_coro(&row(q, max_excluded(q))))
            .collect();
        let want = [
            "| 10^4 | 6 | 0.999934 | 0.24865 | 0.00129 |",
            "| 10^5 | 78 | 0.999932 | 0.24972 | 0.00022 |",
            "| 10^6 | 451 | 0.999967 | 0.24994 | 0.00003 |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn margin_at_the_wall_pinned() {
        let got: Vec<String> = WALL
            .iter()
            .flat_map(|&q| {
                let qf = q as f64;
                let up = !closes(qf, 1.0);
                [sci(delta_stable(qf, 1.0), 9, up), sci(gap(qf, 1.0), 9, up)]
            })
            .collect();
        let want = [
            "-2.395807653e-6",
            "-1.533059397e-4",
            "5.863425182e-6",
            "3.752213034e-4",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn margin_bounds_are_safe() {
        for q in WALL {
            let qf = q as f64;
            let up = !closes(qf, 1.0);
            for (v, s) in [
                (delta_stable(qf, 1.0), sci(delta_stable(qf, 1.0), 9, up)),
                (gap(qf, 1.0), sci(gap(qf, 1.0), 9, up)),
            ] {
                let p: f64 = s.parse().unwrap();
                if up {
                    assert!(p >= v, "not an upper bound at q = {q}: {s}");
                } else {
                    assert!(p <= v, "not a lower bound at q = {q}: {s}");
                }
                assert!(
                    (p - v).abs() <= 1e-8 * v.abs(),
                    "bound too loose at q = {q}"
                );
            }
        }
    }

    #[test]
    fn delta_stable_matches_delta() {
        for q in BASES.iter().chain(COST.iter()) {
            let qf = *q as f64;
            let (a, b) = (delta(qf, 1.0), delta_stable(qf, 1.0));
            assert!(
                (a - b).abs() <= 1e-9 * b.abs(),
                "delta forms disagree at q = {q}"
            );
        }
        for q in CORO {
            let (qf, mf) = (q as f64, max_excluded(q) as f64);
            let (a, b) = (delta(qf, mf), delta_stable(qf, mf));
            assert!(
                (a - b).abs() <= 1e-9 * b.abs(),
                "delta forms disagree at q = {q}"
            );
        }
    }

    #[test]
    fn cost_out_rows_pinned() {
        let mut got: Vec<String> = COST.iter().map(|&q| render_cost(q, 1)).collect();
        got.extend(CORO.iter().map(|&q| render_cost(q, max_excluded(q))));
        let want = [
            "| 3689 | 1 | -2.39581e-6 | 1.65072e-5 | no |",
            "| 3690 | 1 | 5.86342e-6 | 1.65022e-5 | no |",
            "| 5000 | 1 | 6.05211e-3 | 1.17434e-5 | yes |",
            "| 10^4 | 1 | 1.85809e-2 | 5.42923e-6 | yes |",
            "| 10^5 | 1 | 5.09409e-2 | 4.34299e-7 | yes |",
            "| 10^6 | 1 | 7.41160e-2 | 3.61913e-8 | yes |",
            "| 10^9 | 1 | 1.16951e-1 | 2.41275e-11 | yes |",
            "| 10^4 | 6 | 1.29265e-3 | 3.25917e-5 | yes |",
            "| 10^5 | 78 | 2.20253e-4 | 3.39015e-5 | yes |",
            "| 10^6 | 451 | 3.14081e-5 | 1.63296e-5 | yes |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn sharpened_cost_out_is_pinned() {
        let got: Vec<String> = [
            ("`1 + Phi_q/q`", "any", Cst::Step3),
            ("`Psi_q`", "any", Cst::Sharp(false)),
            ("`Psi_q`", "`0` or `q-1`", Cst::Sharp(true)),
            ("`Psi'_q`", "any", Cst::Chord(false)),
            ("`Psi'_q`", "`0` or `q-1`", Cst::Chord(true)),
        ]
        .iter()
        .map(|(n, d, c)| render_sharp(n, d, *c))
        .collect();
        let want = [
            "| `1 + Phi_q/q` | any | 3690 | 3692 | 1.69819e-5 | 1.64921e-5 | yes |",
            "| `Psi_q` | any | 2446 | 2450 | 3.90025e-5 | 2.61622e-5 | yes |",
            "| `Psi_q` | `0` or `q-1` | 1812 | 1815 | 3.85795e-5 | 3.67324e-5 | yes |",
            "| `Psi'_q` | any | 1499 | 1502 | 4.80703e-5 | 4.55409e-5 | yes |",
            "| `Psi'_q` | `0` or `q-1` | 1032 | 1036 | 8.24697e-5 | 6.95785e-5 | yes |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn step3_row_reproduces_the_ladder() {
        let (w, up) = onestep_wall(Cst::Step3, SHARP_HI);
        assert_eq!((w, up), (3_690, true));
        assert_eq!(onestep_cross(Cst::Step3, SHARP_COST_HI).0, 3_692);
        for q in [3_690u64, 10_000, 100_000] {
            let qf = q as f64;
            assert!((onestep(qf, Cst::Step3) - pb(qf, 1.0)).abs() < 1e-12 * pb(qf, 1.0));
        }
    }

    #[test]
    fn the_chord_constant_saves_half_a_base() {
        for q in (36u64..4_000).step_by(1) {
            let qf = q as f64;
            let want = if q % 2 == 0 {
                psi(qf) - qf / 2.0 + 2.0 / PI
            } else {
                psi_chord(qf)
            };
            assert!(
                (psi_chord(qf) - want).abs() <= 1e-9 * psi(qf),
                "chord identity fails at q = {q}"
            );
            assert!(psi_chord(qf) < psi(qf), "chord not sharper at q = {q}");
            assert!(
                psi_chord(qf) >= (1.0 + PI) * qf / 2.0,
                "monotone floor fails at q = {q}"
            );
        }
        assert!(psi_chord(35.0) < (1.0 + PI) * 35.0 / 2.0);
    }

    #[test]
    fn every_sharpened_wall_lies_under_its_predecessor() {
        let walls: Vec<u64> = [
            Cst::Step3,
            Cst::Sharp(false),
            Cst::Sharp(true),
            Cst::Chord(false),
            Cst::Chord(true),
        ]
        .iter()
        .map(|&c| onestep_wall(c, SHARP_HI).0)
        .collect();
        for w in walls.windows(2) {
            assert!(w[1] < w[0], "walls not decreasing: {walls:?}");
        }
        for (i, &c) in [
            Cst::Step3,
            Cst::Sharp(false),
            Cst::Sharp(true),
            Cst::Chord(false),
            Cst::Chord(true),
        ]
        .iter()
        .enumerate()
        {
            let (x, up, _, _) = onestep_cross(c, SHARP_COST_HI);
            assert!(up, "crossing not an up-set at row {i}");
            assert!(
                x >= walls[i] && x - walls[i] <= 5,
                "crossing far from wall at row {i}"
            );
        }
    }

    #[test]
    fn crossover_is_pinned() {
        let (first, down) = crossover(SWEEP_LO, SWEEP_HI);
        assert_eq!(first, 3_692);
        assert_eq!(down, 0);
        for q in [3_690u64, 3_691] {
            let qf = q as f64;
            assert!(
                delta_stable(qf, 1.0) < defect(qf, 1.0),
                "premature crossing at {q}"
            );
        }
        assert!(delta_stable(3_692.0, 1.0) > defect(3_692.0, 1.0));
    }

    #[test]
    fn sci_is_directional() {
        for x in [1.0e-6, -1.0e-6, 9.999999999e-3, 1.23456789e7, -7.5e-11] {
            for d in [3u32, 5, 9] {
                let up: f64 = sci(x, d, true).parse().unwrap();
                let dn: f64 = sci(x, d, false).parse().unwrap();
                assert!(up >= x, "sci up failed on {x} at {d}");
                assert!(dn <= x, "sci down failed on {x} at {d}");
                let w = 10f64.powi(-(d as i32)) * x.abs() * 20.0;
                assert!(up - dn <= w, "sci band too wide on {x} at {d}");
            }
        }
    }

    #[test]
    fn rounding_is_directional() {
        for q in BASES {
            let (qf, mf) = (q as f64, 1.0);
            let r = row(q, 1);
            assert!(fixed(r.c, 5).parse::<f64>().unwrap() >= c_exp(qf, mf) - GUARD);
            assert!(fixed(r.delta, 5).parse::<f64>().unwrap() <= delta(qf, mf) + GUARD);
            assert!(fixed(r.alpha, 6).parse::<f64>().unwrap() <= alpha(qf, mf) + GUARD);
            assert!(c_exp(qf, mf) - fixed(r.c, 5).parse::<f64>().unwrap() > -1e-5);
        }
    }

    #[test]
    fn l1_floor_holds() {
        for q in 3..5_000u64 {
            assert!(c_exp(q as f64, 1.0) >= 0.0, "negative c_q at {q}");
            assert!(pb(q as f64, 1.0) >= 1.0, "PB below the l^1 floor at {q}");
        }
    }

    #[test]
    fn the_chord_floor_carries_its_harmonic_convention() {
        let exact = |n: f64| (1..=(n as u64)).map(|j| 1.0 / j as f64).sum::<f64>();
        let psi_exact = |q: f64| {
            let p = (q / 2.0).floor();
            let h = exact(p - 1.0);
            if (q as u64).is_multiple_of(2) {
                (q / PI) * (2.0 * h - 1.0 + 1.0 / p) + (1.0 - 2.0 / PI) * q / 2.0
            } else {
                (q / PI) * (2.0 * h - 1.0 + 2.0 / p)
                    + (1.0 - 2.0 / PI) * (q / 2.0 + 1.0 / (2.0 * q))
            }
        };
        assert!(psi_chord(36.0) >= (1.0 + PI) * 18.0);
        assert!(psi_exact(36.0) < (1.0 + PI) * 18.0);
        for q in 37u64..4_000 {
            let qf = q as f64;
            assert!(
                psi_exact(qf) >= (1.0 + PI) * qf / 2.0,
                "harmonic-number floor fails at q = {q}"
            );
        }
    }

    #[test]
    fn harmonic_bound_dominates() {
        let mut h = 0.0;
        for n in 1..=2_000u64 {
            h += 1.0 / n as f64;
            assert!(harmonic_bound(n as f64) >= h, "harmonic bound fails at {n}");
        }
    }

    #[test]
    fn kernel_bound_dominates_sampled_sup() {
        for q in [50u64, 101, 200] {
            let bound = phi(q as f64);
            let mut worst: f64 = 0.0;
            for i in 0..4_001u64 {
                let t = i as f64 / 4_001.0;
                worst = worst.max(kernel_sum(q, t));
            }
            assert!(
                worst <= bound,
                "kernel bound fails at q = {q}: {worst} > {bound}"
            );
            assert!(
                worst > 0.8 * bound,
                "kernel bound absurdly loose at q = {q}"
            );
        }
    }

    #[test]
    fn ladder_rows_pinned() {
        let got: Vec<String> = LADDER
            .iter()
            .map(|&(n, d)| render_ladder(Rat::new(n, d)))
            .collect();
        let want = [
            "| 1/2 | 3/4 | both | 3690 | 723 |",
            "| 13/25 | 1417/1850 | Zhang | 8578 | 1486 |",
            "| 11/20 | 913/1160 | Zhang | 33547 | 4754 |",
            "| 4/7 | 4/5 | both | 92317 | 11221 |",
            "| 3/5 | 4/5 | BH | 92317 | 11221 |",
            "| 2/3 | 5/6 | BH | 3107080 | 216023 |",
            "| 3/4 | 7/8 | BH | 6939524168 | 129458304 |",
            "| 4/5 | 9/10 | BH | <= 3.09358e13 | 128606353005 |",
            "| 9/10 | 19/20 | BH | <= 3.23663e34 | <= 1.73431e28 |",
            "| 19/20 | 39/40 | BH | <= 9.24614e83 | <= 3.30712e68 |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn ladder_trend_pinned() {
        assert_eq!(
            ladder_trend(),
            "3.57 3.94 4.53 4.97 4.97 6.50 9.85 13.50 34.52 83.97"
        );
    }

    #[test]
    fn last_down_step_below_the_floor() {
        assert!((4.0 / PI) / 662.0 > 0.25 * 664f64.powf(-0.75));
        let mut last = 0u64;
        for q in 3..SWEEP_LO {
            if gap((q + 1) as f64, 1.0) <= gap(q as f64, 1.0) {
                last = q;
            }
        }
        assert_eq!(last, 662);
    }

    #[test]
    fn general_b_floor_bound() {
        let k = 2.0 / PI;
        let slack = 2.0 * 723f64.powf(-0.75);
        assert!(slack < 0.015, "step-A slack above 0.015: {slack}");
        let konst = slack - 1.0 - 4.0 / PI - k * GAMMA + k * (1448.0f64 / 721.0).ln()
            - (1.0 - k) * (721.0 / 723.0);
        assert!(
            konst < -2.544,
            "general-b constant weaker than -2.544: {konst}"
        );
        assert!(4.0 * k - 2.544 > 0.0, "bracket not increasing on (0, 1/4]");
        let bracket = 1.291 - k * (1.291 * 4.0f64).ln() - 2.544 / 4.0;
        assert!(bracket < -0.39014, "bracket above -0.39014: {bracket}");
        assert!(
            4.0 * bracket < -1.56,
            "floor gap bound weaker than -1.56: {bracket}"
        );
        let over =
            3691f64.powf(-0.75) + k * ((3690.0f64 / 3689.0).ln() + 1.0 / 3689.0) + 0.727 / 3691.0;
        assert!(over < 0.004, "majorant overshoot above 0.004: {over}");
        assert!(u_bound(3690.0, 1417.0 / 1850.0) < -0.95);
        let mut b = 0.75;
        while b < 0.98 {
            let floor = mono_floor(b);
            assert!(
                gap_b(floor, b, 1.0) < -1.56,
                "gap at the floor above the bound at b = {b}"
            );
            if floor > SWEEP_LO as f64 {
                assert!(
                    u_bound(SWEEP_LO as f64, b) < 0.0,
                    "majorant positive at b = {b}"
                );
                assert!(
                    u_bound(floor, b) < 0.0,
                    "majorant positive at the floor, b = {b}"
                );
            }
            b += 0.005;
        }
    }

    #[test]
    fn ladder_grh_rung_is_the_wall() {
        let (b, src) = ladder_b(Rat::new(1, 2));
        assert!(b.same(Rat::new(3, 4)));
        assert_eq!(src, "both");
        assert_eq!(ladder_wall(0.75), 3_690.0);
        assert!(gap(3_690.0, 1.0) > 0.0 && gap(3_689.0, 1.0) < 0.0);
    }

    #[test]
    fn zhang_beats_baker_harman_strictly_inside() {
        assert!(zhang(Rat::new(1, 2))
            .unwrap()
            .same(baker_harman(Rat::new(1, 2))));
        assert!(zhang(Rat::new(4, 7))
            .unwrap()
            .same(baker_harman(Rat::new(4, 7))));
        let (lo, hi) = (Rat::new(1, 2), Rat::new(4, 7));
        let mut seen = 0u64;
        for d in 2..=200i64 {
            for n in 1..d {
                let a = Rat::new(n, d);
                if !lo.lt(a) || !a.lt(hi) {
                    continue;
                }
                let z = zhang(a).unwrap();
                assert!(
                    z.lt(baker_harman(a)),
                    "Zhang not smaller at a = {}",
                    a.show()
                );
                assert_eq!(ladder_b(a).1, "Zhang");
                seen += 1;
            }
        }
        assert!(seen > 400, "too few interior rationals tested");
        assert!(zhang(Rat::new(3, 5)).is_none());
    }

    #[test]
    fn ladder_floor_forces_steps_up() {
        for (n, d) in LADDER {
            let a = Rat::new(n, d);
            let b = ladder_b(a).0.val();
            let floor = mono_floor(b);
            assert!(
                mono_ok(floor, b),
                "floor fails its own test at a = {}",
                a.show()
            );
            if floor < DYADIC_TOP {
                assert!(
                    !mono_ok(floor - 1.0, b),
                    "floor not least at a = {}",
                    a.show()
                );
            }
            assert!(
                floor < ladder_wall(b),
                "floor above the wall at a = {}",
                a.show()
            );
            let lo = floor.max(3_690.0);
            if lo < 4.0e6 {
                let mut q = lo;
                while q < (lo + 20_000.0).min(4.0e6) {
                    assert!(
                        gap_b(q + 1.0, b, 1.0) > gap_b(q, b, 1.0),
                        "step down at q = {q}"
                    );
                    q += 1.0;
                }
            }
        }
    }

    #[test]
    fn pb_step_below_the_constant() {
        let mut worst = 0.0f64;
        for q in 40..400_000u64 {
            let qf = q as f64;
            worst = worst.max((pb(qf + 1.0, 1.0) - pb(qf, 1.0)) * (qf - 2.0));
        }
        assert!(
            worst < MONO_C,
            "PB step above the monotone constant: {worst}"
        );
        assert!(worst > 1.27, "PB step bound absurdly loose: {worst}");
    }

    #[test]
    fn pb_low_is_a_lower_bound() {
        for q in 3..20_000u64 {
            let qf = q as f64;
            assert!(pb_low(qf) <= pb(qf, 1.0), "pb_low above pb at q = {q}");
            assert!(
                gap_b(qf, 0.9, 1.0) <= u_bound(qf, 0.9),
                "u_bound below gap at q = {q}"
            );
        }
    }

    #[test]
    fn ladder_below_the_floor_is_clear() {
        for (n, d) in LADDER {
            let a = Rat::new(n, d);
            let b = ladder_b(a).0.val();
            for q in 3..3_690u64 {
                assert!(
                    gap_b(q as f64, b, 1.0) < 0.0,
                    "close below the wall at a = {}",
                    a.show()
                );
            }
            let floor = mono_floor(b);
            if floor > 3_690.0 {
                assert!(u_bound(3_690.0, b) < 0.0);
                assert!(u_bound(floor, b) < 0.0);
                assert!(u_bound(3_690.0, b).max(u_bound(floor, b)) < -1.5);
            }
        }
    }

    #[test]
    fn ladder_wall_matches_exhaustive_scan() {
        for (n, d) in LADDER {
            let a = Rat::new(n, d);
            let b = ladder_b(a).0.val();
            let q0 = ladder_wall(b);
            if q0 > 4.0e6 {
                continue;
            }
            let mut found = 0.0f64;
            let mut q = 3.0;
            while q <= q0 {
                if gap_b(q, b, 1.0) > 0.0 {
                    found = q;
                    break;
                }
                q += 1.0;
            }
            assert_eq!(
                found,
                q0,
                "scan disagrees with bisection at a = {}",
                a.show()
            );
        }
    }

    #[test]
    fn ladder_m_corollary_pinned() {
        let got: Vec<String> = LADDER
            .iter()
            .map(|&(n, d)| Rat::new(n, d))
            .filter(|a| ladder_wall(ladder_b(*a).0.val()) < CORO_Q)
            .map(render_ladder_coro)
            .collect();
        let want = [
            "| 1/2 | 3/4 | both | 1971 |",
            "| 13/25 | 1417/1850 | Zhang | 1002 |",
            "| 11/20 | 913/1160 | Zhang | 365 |",
            "| 4/7 | 4/5 | both | 176 |",
            "| 3/5 | 4/5 | BH | 176 |",
            "| 2/3 | 5/6 | BH | 8 |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn ladder_wall_rises_with_a() {
        let mut prev = 0.0f64;
        for (n, d) in LADDER {
            let b = ladder_b(Rat::new(n, d)).0.val();
            let q0 = ladder_wall(b);
            assert!(q0 >= prev, "wall not rising at {n}/{d}");
            prev = q0;
        }
        assert!(prev > 1.0e83);
    }

    #[test]
    fn asymptotic_shape() {
        let q = 1e12;
        let c = c_exp(q, 1.0);
        let model = (q.ln().ln() + (2.0 / PI).ln()) / q.ln();
        assert!((c - model).abs() < 0.01);
        assert!(delta(q, 1.0) > 0.13 && delta(q, 1.0) < 0.25);
    }
}
