use crate::menergy::{bits_of, box_row, boxes, column, exps, qpow, random_column, Bits, Rng};

// SIEVES

pub fn spf_table(lim: usize) -> Vec<u32> {
    let mut spf = vec![0u32; lim + 1];
    let mut i = 2usize;
    while i <= lim {
        if spf[i] == 0 {
            let mut j = i;
            while j <= lim {
                if spf[j] == 0 {
                    spf[j] = i as u32;
                }
                j += i;
            }
        }
        i += 1;
    }
    spf
}

pub fn mobius(lim: usize, spf: &[u32]) -> Vec<i8> {
    let mut mu = vec![0i8; lim + 1];
    if lim >= 1 {
        mu[1] = 1;
    }
    for n in 2..=lim {
        let p = spf[n] as usize;
        let m = n / p;
        mu[n] = if m.is_multiple_of(p) { 0 } else { -mu[m] };
    }
    mu
}

fn liouville(lim: usize, spf: &[u32]) -> Vec<i8> {
    let mut lam = vec![0i8; lim + 1];
    if lim >= 1 {
        lam[1] = 1;
    }
    for n in 2..=lim {
        lam[n] = -lam[n / spf[n] as usize];
    }
    lam
}

fn sign_draw(lim: usize, seed: u64) -> Vec<i8> {
    let mut rng = Rng::new(seed);
    (0..=lim)
        .map(|_| if rng.sign() > 0.0 { 1i8 } else { -1i8 })
        .collect()
}

pub fn coefficients(lim: usize) -> (Vec<&'static str>, Vec<Vec<i8>>) {
    let spf = spf_table(lim);
    let mu = mobius(lim, &spf);
    let names = vec![
        "1", "mu", "lambda", "rand1", "rand2", "rand3", "sfr1", "sfr2", "sfr3",
    ];
    let masked = |seed: u64| -> Vec<i8> {
        let d = sign_draw(lim, seed);
        (0..=lim).map(|n| d[n] * mu[n] * mu[n]).collect()
    };
    let vals = vec![
        vec![1i8; lim + 1],
        mu.clone(),
        liouville(lim, &spf),
        sign_draw(lim, 0x51_6e_ed_01),
        sign_draw(lim, 0x51_6e_ed_02),
        sign_draw(lim, 0x51_6e_ed_03),
        masked(0x5f_6e_ed_01),
        masked(0x5f_6e_ed_02),
        masked(0x5f_6e_ed_03),
    ];
    (names, vals)
}

// THE ROW

pub struct SignedRow {
    pub m: u64,
    pub n: u64,
    pub r: u128,
    pub e: u128,
    pub num: Vec<i128>,
    pub diag: Vec<f64>,
    pub rms: Option<f64>,
    pub csloss: Vec<f64>,
}

pub fn signed_row(bits: &Bits, mm: u64, nn: u64, g: u128, x: u128, coefs: &[Vec<i8>]) -> SignedRow {
    let nb = coefs.len();
    let mut cprime = vec![0i64; nn as usize];
    let mut qq = vec![0i128; nb];
    let mut cc = vec![0i128; nb];
    let mut u = vec![0i64; nb];
    let mut r: u128 = 0;
    let mut e: u128 = 0;
    let dl = g as f64 / x as f64;
    let bsum: Vec<f64> = coefs
        .iter()
        .map(|b| (nn..2 * nn).map(|l| b[l as usize] as f64).sum())
        .collect();
    let mut l1 = vec![0.0f64; nb];
    let mut sq = vec![0.0f64; nb];
    for m in mm..2 * mm {
        u.iter_mut().for_each(|v| *v = 0);
        let mut c = 0i64;
        for l in nn..2 * nn {
            if bits.get((m * l) as usize) {
                c += 1;
                cprime[(l - nn) as usize] += 1;
                for (j, b) in coefs.iter().enumerate() {
                    u[j] += b[l as usize] as i64;
                }
            }
        }
        r += c as u128;
        e += (c as u128) * (c as u128);
        for j in 0..nb {
            qq[j] += (u[j] as i128) * (u[j] as i128);
            cc[j] += u[j] as i128;
            let v = u[j] as f64 - dl * bsum[j];
            l1[j] += v.abs();
            sq[j] += v * v;
        }
    }
    let csloss: Vec<f64> = (0..nb)
        .map(|j| {
            let bd = (mm as f64 * sq[j]).sqrt();
            if bd > 0.0 {
                l1[j] / bd
            } else {
                0.0
            }
        })
        .collect();
    let xi = x as i128;
    let gi = g as i128;
    let mi = mm as i128;
    let delta = g as f64 / x as f64;
    let mut num = Vec::with_capacity(nb);
    let mut diag = Vec::with_capacity(nb);
    for (j, b) in coefs.iter().enumerate() {
        let mut bs = 0i128;
        let mut s2 = 0i128;
        let mut pp = 0i128;
        for l in nn..2 * nn {
            let v = b[l as usize] as i128;
            bs += v;
            s2 += v * v;
            pp += v * v * cprime[(l - nn) as usize] as i128;
        }
        num.push(
            xi * xi * (qq[j] - pp) - 2 * gi * xi * (bs * cc[j] - pp)
                + mi * gi * gi * (bs * bs - s2),
        );
        diag.push(
            pp as f64 * (1.0 - delta) * (1.0 - delta) + (mi * s2 - pp) as f64 * delta * delta,
        );
    }
    SignedRow {
        m: mm,
        n: nn,
        r,
        e,
        num,
        diag,
        rms: None,
        csloss,
    }
}

const RMSCAP: u64 = 2048;

pub fn pair_rms(bits: &Bits, mm: u64, nn: u64, g: u128, x: u128) -> Option<f64> {
    if mm > RMSCAP {
        return None;
    }
    let ms = mm as usize;
    let mut tt = vec![0u32; ms * ms];
    let mut cm = vec![0u32; ms];
    let mut cprime = vec![0u32; nn as usize];
    let mut hits: Vec<u32> = Vec::with_capacity(ms);
    for l in nn..2 * nn {
        hits.clear();
        for m in mm..2 * mm {
            if bits.get((m * l) as usize) {
                hits.push((m - mm) as u32);
            }
        }
        cprime[(l - nn) as usize] = hits.len() as u32;
        for &a in hits.iter() {
            cm[a as usize] += 1;
            for &b in hits.iter() {
                tt[a as usize * ms + b as usize] += 1;
            }
        }
    }
    let delta = g as f64 / x as f64;
    let nf = nn as f64;
    let mut all = 0.0f64;
    for a in 0..ms {
        for b in 0..ms {
            let gg =
                tt[a * ms + b] as f64 - delta * (cm[a] as f64 + cm[b] as f64) + nf * delta * delta;
            all += gg * gg;
        }
    }
    let mut dg = 0.0f64;
    for &c in cprime.iter() {
        let w = c as f64 * (1.0 - delta) * (1.0 - delta) + (mm as f64 - c as f64) * delta * delta;
        dg += w * w;
    }
    Some((2.0 * (all - dg)).max(0.0).sqrt())
}

pub fn engineered(
    bits: &Bits,
    mm: u64,
    nn: u64,
    g: u128,
    x: u128,
    rounds: usize,
) -> (f64, f64, usize) {
    let delta = g as f64 / x as f64;
    let ms = mm as usize;
    let ns = nn as usize;
    let mut b = vec![1.0f64; ns];
    let mut v = vec![0.0f64; ms];
    let mut cprime = vec![0u32; ns];
    for (i, m) in (mm..2 * mm).enumerate() {
        let mut acc = 0.0f64;
        for (j, l) in (nn..2 * nn).enumerate() {
            if bits.get((m * l) as usize) {
                acc += 1.0 - delta;
                cprime[j] += 1;
            } else {
                acc -= delta;
            }
        }
        v[i] = acc;
    }
    let w: Vec<f64> = cprime
        .iter()
        .map(|&c| c as f64 * (1.0 - delta) * (1.0 - delta) + (mm as f64 - c as f64) * delta * delta)
        .collect();
    let mut flips = 0usize;
    let mut grad = vec![0.0f64; ns];
    for _ in 0..rounds {
        grad.iter_mut().for_each(|z| *z = 0.0);
        for (i, m) in (mm..2 * mm).enumerate() {
            for (j, l) in (nn..2 * nn).enumerate() {
                let ps = if bits.get((m * l) as usize) {
                    1.0 - delta
                } else {
                    -delta
                };
                grad[j] += v[i] * ps;
            }
        }
        let mut best = 0usize;
        let mut gain = 0.0f64;
        for j in 0..ns {
            let d = -4.0 * b[j] * grad[j] + 4.0 * w[j];
            if d < gain {
                gain = d;
                best = j;
            }
        }
        if gain >= 0.0 {
            break;
        }
        let l = nn + best as u64;
        for (i, m) in (mm..2 * mm).enumerate() {
            let ps = if bits.get((m * l) as usize) {
                1.0 - delta
            } else {
                -delta
            };
            v[i] -= 2.0 * b[best] * ps;
        }
        b[best] = -b[best];
        flips += 1;
    }
    let form: f64 = v.iter().map(|z| z * z).sum();
    let dg: f64 = w.iter().sum();
    (form, dg, flips)
}

fn round_div(num: i128, den: i128) -> i128 {
    if num >= 0 {
        (num + den / 2) / den
    } else {
        -((-num + den / 2) / den)
    }
}

// THE SWEEP

pub struct SignedSweep {
    pub top: SignedRow,
    pub l5: Option<SignedRow>,
    pub seen: usize,
    pub live: usize,
}

pub fn sweep(bits: &Bits, x: u128, g: u128, coefs: &[Vec<i8>]) -> Option<SignedSweep> {
    let all = boxes(x);
    let edge = (x as f64).powf(0.4);
    let mut top: Option<SignedRow> = None;
    let mut l5: Option<SignedRow> = None;
    let mut live = 0usize;
    for &(mm, nn) in all.iter() {
        let row = signed_row(bits, mm, nn, g, x, coefs);
        if row.r == 0 {
            continue;
        }
        live += 1;
        let score = row.num[0].abs();
        if top.as_ref().map(|b| score > b.num[0].abs()).unwrap_or(true) {
            top = Some(signed_row(bits, mm, nn, g, x, coefs));
        }
        if mm as f64 >= edge
            && nn as f64 >= edge
            && l5.as_ref().map(|b| score > b.num[0].abs()).unwrap_or(true)
        {
            l5 = Some(row);
        }
    }
    top.map(|t| SignedSweep {
        top: t,
        l5,
        seen: all.len(),
        live,
    })
}

// STUDY

pub struct Cell {
    pub q: u64,
    pub digits: Vec<u64>,
    pub label: &'static str,
    pub depths: [usize; 2],
}

fn ex(q: u64, e: u64) -> Vec<u64> {
    (0..q).filter(|&f| f != e).collect()
}

pub fn cells() -> Vec<Cell> {
    vec![
        Cell {
            q: 3,
            digits: vec![0, 1],
            label: "01",
            depths: [12, 14],
        },
        Cell {
            q: 4,
            digits: vec![0, 1, 2],
            label: "012",
            depths: [8, 9],
        },
        Cell {
            q: 5,
            digits: vec![0, 1, 2, 3],
            label: "0123",
            depths: [7, 8],
        },
        Cell {
            q: 10,
            digits: ex(10, 7),
            label: "ex7",
            depths: [5, 6],
        },
    ]
}

const PAIRCAP: usize = 8_000_000;

pub struct Tally {
    pub cells: usize,
    pub rand: [usize; 3],
    pub sf: [usize; 3],
    pub mu_under_lam: usize,
    pub musq: (f64, f64),
    pub lamsq: (f64, f64),
    pub under_r: usize,
    pub one_over_r: usize,
    pub fd: (f64, f64),
    pub murms: (f64, f64),
    pub sfrms: (f64, f64),
    pub cs: (f64, f64),
}

fn place(v: f64, lo: f64, hi: f64) -> usize {
    if v < lo {
        0
    } else if v > hi {
        2
    } else {
        1
    }
}

fn emit(
    cell: &Cell,
    l: usize,
    kind: &str,
    tag: &str,
    sw: &SignedSweep,
    row: &SignedRow,
    tally: &mut Tally,
) {
    let (q, label, x) = (cell.q, cell.label, qpow(cell.q, l));
    let den = (x as i128) * (x as i128);
    let sig: Vec<i128> = row.num.iter().map(|&n| round_div(n, den)).collect();
    let val = |n: i128| (n as f64) / (den as f64);
    let base = row.num[0].abs() as f64;
    let ratio = |n: i128| {
        if base == 0.0 {
            0.0
        } else {
            (n.abs() as f64) / base
        }
    };
    let band = |a: usize, b: usize, c: usize| -> (f64, f64) {
        let v = [ratio(row.num[a]), ratio(row.num[b]), ratio(row.num[c])];
        (
            v.iter().cloned().fold(f64::INFINITY, f64::min),
            v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    let (rlo, rhi) = band(3, 4, 5);
    let (slo, shi) = band(6, 7, 8);
    let mut fdlo = f64::INFINITY;
    let mut fdhi = f64::NEG_INFINITY;
    for j in 1..row.num.len() {
        let f = (val(row.num[j]) + row.diag[j]) / row.diag[j];
        fdlo = fdlo.min(f);
        fdhi = fdhi.max(f);
    }
    let (rtxt, mrms, lrms) = match row.rms {
        Some(v) if v > 0.0 => (
            format!("{:.4e}", v),
            format!("{:.4}", val(row.num[1]).abs() / v),
            format!("{:.4}", val(row.num[2]).abs() / v),
        ),
        _ => ("-".to_string(), "-".to_string(), "-".to_string()),
    };
    let root = val(row.num[0]).abs().sqrt();
    let sq = |n: i128| {
        if root == 0.0 {
            0.0
        } else {
            val(n).abs() / root
        }
    };
    let sqband = [sq(row.num[6]), sq(row.num[7]), sq(row.num[8])];
    let sqlo = sqband.iter().cloned().fold(f64::INFINITY, f64::min);
    let sqhi = sqband.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if kind == "digit" {
        tally.cells += 1;
        tally.rand[place(ratio(row.num[1]), rlo, rhi)] += 1;
        tally.sf[place(ratio(row.num[1]), slo, shi)] += 1;
        if row.num[1].abs() < row.num[2].abs() {
            tally.mu_under_lam += 1;
        }
        tally.musq.0 = tally.musq.0.min(sq(row.num[1]));
        tally.musq.1 = tally.musq.1.max(sq(row.num[1]));
        tally.lamsq.0 = tally.lamsq.0.min(sq(row.num[2]));
        tally.lamsq.1 = tally.lamsq.1.max(sq(row.num[2]));
        let rr = row.r as f64;
        if val(row.num[1]).abs() < rr && val(row.num[2]).abs() < rr {
            tally.under_r += 1;
        }
        if val(row.num[0]).abs() > rr {
            tally.one_over_r += 1;
        }
        for j in 1..row.num.len() {
            let f = (val(row.num[j]) + row.diag[j]) / row.diag[j];
            tally.fd.0 = tally.fd.0.min(f);
            tally.fd.1 = tally.fd.1.max(f);
        }
        tally.cs.0 = tally.cs.0.min(row.csloss[1]);
        tally.cs.1 = tally.cs.1.max(row.csloss[1]);
        if let Some(rv) = row.rms {
            if rv > 0.0 {
                tally.murms.0 = tally.murms.0.min(val(row.num[1]).abs() / rv);
                tally.murms.1 = tally.murms.1.max(val(row.num[1]).abs() / rv);
                for j in 6..9 {
                    tally.sfrms.0 = tally.sfrms.0.min(val(row.num[j]).abs() / rv);
                    tally.sfrms.1 = tally.sfrms.1.max(val(row.num[j]).abs() / rv);
                }
            }
        }
    }
    println!(
        "| {} | {} | {} | {} | {} | {} | {}/{} | {} | {} | {} | {} | {} | {} | {:.6} | {:.6} | {:.6} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} | {} | {:.4} | {:.4} |",
        q,
        label,
        l,
        x,
        kind,
        tag,
        sw.live,
        sw.seen,
        row.m,
        row.n,
        row.r,
        sig[0],
        sig[1],
        sig[2],
        exps(val(row.num[0]).abs(), x),
        exps(val(row.num[1]).abs(), x),
        exps(val(row.num[2]).abs(), x),
        ratio(row.num[1]),
        ratio(row.num[2]),
        rlo,
        rhi,
        slo,
        shi,
        sq(row.num[1]),
        sq(row.num[2]),
        sqlo,
        sqhi,
        rtxt,
        mrms,
        lrms,
        fdlo,
        fdhi
    );
}

pub fn run() {
    let cs = cells();
    let mut lim = 16usize;
    for c in cs.iter() {
        for &l in c.depths.iter() {
            let x = qpow(c.q, l);
            for &(_, nn) in boxes(x).iter() {
                lim = lim.max(2 * nn as usize);
            }
        }
    }
    let (_, coefs) = coefficients(lim);
    println!("menergy signed");
    println!("| q | F | L | x | kind | box | live/seen | M | N | R | Sigma_1 | Sigma_mu | Sigma_lam | exp_1 | exp_mu | exp_lam | mu/1 | lam/1 | rand lo | rand hi | sf lo | sf hi | mu/sq | lam/sq | sf sq lo | sf sq hi | rms | mu/rms | lam/rms | fd lo | fd hi |");
    let mut tally = Tally {
        cells: 0,
        rand: [0; 3],
        sf: [0; 3],
        mu_under_lam: 0,
        musq: (f64::INFINITY, f64::NEG_INFINITY),
        lamsq: (f64::INFINITY, f64::NEG_INFINITY),
        under_r: 0,
        one_over_r: 0,
        fd: (f64::INFINITY, f64::NEG_INFINITY),
        murms: (f64::INFINITY, f64::NEG_INFINITY),
        sfrms: (f64::INFINITY, f64::NEG_INFINITY),
        cs: (f64::INFINITY, f64::NEG_INFINITY),
    };
    for c in cs.iter() {
        for (d, &l) in c.depths.iter().enumerate() {
            let x = qpow(c.q, l);
            let vals = column(c.q, &c.digits, l);
            let g = vals.len() as u128;
            let bits = bits_of(&vals, x);
            let mut sw = match sweep(&bits, x, g, &coefs) {
                Some(s) => s,
                None => continue,
            };
            sw.top.rms = pair_rms(&bits, sw.top.m, sw.top.n, g, x);
            if let Some(r) = sw.l5.as_mut() {
                r.rms = pair_rms(&bits, r.m, r.n, g, x);
            }
            let delta = g as f64 / x as f64;
            for (tag, row) in [Some(("top", &sw.top)), sw.l5.as_ref().map(|r| ("l5", r))]
                .into_iter()
                .flatten()
            {
                if let Some(chk) = box_row(&bits, row.m, row.n, delta, PAIRCAP) {
                    assert_eq!(chk.r, row.r, "R disagrees with the census routine");
                    assert_eq!(chk.e, row.e, "E_x(M,N) disagrees with the census routine");
                }
                emit(c, l, "digit", tag, &sw, row, &mut tally);
            }
            if d + 1 == c.depths.len() {
                let mut rng = Rng::new(0x51_6e_ed_00 + c.q * 131 + l as u64);
                let rv = random_column(x, vals.len(), &mut rng);
                let rbits = bits_of(&rv, x);
                if let Some(mut rs) = sweep(&rbits, x, g, &coefs) {
                    rs.top.rms = pair_rms(&rbits, rs.top.m, rs.top.n, g, x);
                    if let Some(r) = rs.l5.as_mut() {
                        r.rms = pair_rms(&rbits, r.m, r.n, g, x);
                    }
                    let tag = if rs.l5.is_some() { "l5" } else { "top" };
                    let row = rs.l5.as_ref().unwrap_or(&rs.top);
                    emit(c, l, "random", tag, &rs, row, &mut tally);
                }
            }
        }
    }
    println!("menergy signed summary");
    println!("| digit cells | mu below rand | mu inside rand | mu above rand | mu below sf | mu inside sf | mu above sf | mu under lam | mu/sq lo | mu/sq hi | lam/sq lo | lam/sq hi | signed under R | unsigned over R | form/diag lo | form/diag hi | mu/rms lo | mu/rms hi | sf/rms lo | sf/rms hi | mu cs lo | mu cs hi |");
    println!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} |",
        tally.cells,
        tally.rand[0],
        tally.rand[1],
        tally.rand[2],
        tally.sf[0],
        tally.sf[1],
        tally.sf[2],
        tally.mu_under_lam,
        tally.musq.0,
        tally.musq.1,
        tally.lamsq.0,
        tally.lamsq.1,
        tally.under_r,
        tally.one_over_r,
        tally.fd.0,
        tally.fd.1,
        tally.murms.0,
        tally.murms.1,
        tally.sfrms.0,
        tally.sfrms.1,
        tally.cs.0,
        tally.cs.1
    );
    println!("menergy signed engineered");
    println!("| q | F | L | box | M | N | R | flips | form | diag | form/diag | bound | floor | bound/floor |");
    for c in cs.iter() {
        let l = c.depths[0];
        let x = qpow(c.q, l);
        let vals = column(c.q, &c.digits, l);
        let g = vals.len() as u128;
        let delta = g as f64 / x as f64;
        let bits = bits_of(&vals, x);
        let sw = match sweep(&bits, x, g, &coefs) {
            Some(s) => s,
            None => continue,
        };
        for (tag, row) in [Some(("top", &sw.top)), sw.l5.as_ref().map(|r| ("l5", r))]
            .into_iter()
            .flatten()
        {
            let (form, dg, flips) = engineered(&bits, row.m, row.n, g, x, 4000);
            let bound = (row.m as f64 * form).max(0.0).sqrt();
            let floor = (1.0 - delta) * (row.m as f64 * row.r as f64).sqrt();
            println!(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {:.4e} | {:.4e} | {:.4} | {:.4e} | {:.4e} | {:.4} |",
                c.q, c.label, l, tag, row.m, row.n, row.r, flips, form, dg, form / dg, bound, floor, bound / floor
            );
        }
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sieves_are_right() {
        let spf = spf_table(20);
        let mu = mobius(20, &spf);
        let lam = liouville(20, &spf);
        assert_eq!(&mu[1..=12], &[1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0][..]);
        assert_eq!(
            &lam[1..=12],
            &[1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1][..]
        );
    }

    fn brute_num(bits: &Bits, mm: u64, nn: u64, g: u128, x: u128, b: &[i8]) -> i128 {
        let xi = x as i128;
        let gi = g as i128;
        let mut total = 0i128;
        for l1 in nn..2 * nn {
            for l2 in nn..2 * nn {
                if l1 == l2 {
                    continue;
                }
                let mut w = 0i128;
                for m in mm..2 * mm {
                    let p1 = if bits.get((m * l1) as usize) { xi } else { 0 } - gi;
                    let p2 = if bits.get((m * l2) as usize) { xi } else { 0 } - gi;
                    w += p1 * p2;
                }
                total += b[l1 as usize] as i128 * b[l2 as usize] as i128 * w;
            }
        }
        total
    }

    #[test]
    fn signed_row_matches_the_definition() {
        let (_, coefs) = coefficients(256);
        for (q, digits, l, mm, nn) in [
            (3u64, vec![0u64, 1], 8usize, 8u64, 16u64),
            (3, vec![1, 2], 8, 16, 16),
            (4, vec![0, 1, 2], 6, 8, 32),
            (5, vec![0, 1, 2, 3], 5, 16, 32),
        ] {
            let x = qpow(q, l);
            let vals = column(q, &digits, l);
            let g = vals.len() as u128;
            let bits = bits_of(&vals, x);
            let row = signed_row(&bits, mm, nn, g, x, &coefs);
            for (j, b) in coefs.iter().enumerate() {
                assert_eq!(
                    row.num[j],
                    brute_num(&bits, mm, nn, g, x, b),
                    "coefficient {j} at q={q} L={l}"
                );
            }
        }
    }

    #[test]
    fn the_unsigned_row_matches_the_census() {
        let (_, coefs) = coefficients(1024);
        for (q, digits, l, mm, nn) in [
            (3u64, vec![0u64, 1], 12usize, 256u64, 512u64),
            (4, vec![0, 1, 2], 9, 256u64, 256u64),
            (5, vec![0, 1, 2, 3], 8, 256u64, 256u64),
        ] {
            let x = qpow(q, l);
            let vals = column(q, &digits, l);
            let g = vals.len() as u128;
            let delta = g as f64 / x as f64;
            let bits = bits_of(&vals, x);
            let row = signed_row(&bits, mm, nn, g, x, &coefs);
            let chk = box_row(&bits, mm, nn, delta, PAIRCAP).unwrap();
            assert_eq!(chk.r, row.r);
            assert_eq!(chk.e, row.e);
        }
    }

    #[test]
    fn the_random_column_kills_the_unsigned_sum() {
        let (_, coefs) = coefficients(1024);
        let (q, digits, l, mm, nn) = (3u64, vec![0u64, 1], 12usize, 128u64, 128u64);
        let x = qpow(q, l);
        let vals = column(q, &digits, l);
        let g = vals.len() as u128;
        let bits = bits_of(&vals, x);
        let digit = signed_row(&bits, mm, nn, g, x, &coefs);
        let mut rng = Rng::new(0x7e57);
        let rv = random_column(x, vals.len(), &mut rng);
        let rbits = bits_of(&rv, x);
        let rand = signed_row(&rbits, mm, nn, g, x, &coefs);
        assert!(
            digit.num[0] > 20 * rand.num[0].abs(),
            "digit {} against random {}",
            digit.num[0],
            rand.num[0]
        );
    }
}
