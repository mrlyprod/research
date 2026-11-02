use crate::menergy::{bits_of, boxes, column, qpow, Bits};
use crate::signed::{cells, mobius, spf_table, sweep};

// THE PIECES

pub fn iroot(n: u128, k: u32) -> u64 {
    let mut r = (n as f64).powf(1.0 / k as f64) as u64;
    while r > 0 && (r as u128).pow(k) > n {
        r -= 1;
    }
    while ((r + 1) as u128).pow(k) <= n {
        r += 1;
    }
    r
}

pub fn pieces(mu: &[i8], lim: usize, u: u64) -> (Vec<f64>, Vec<f64>) {
    let mut b = vec![0.0f64; lim + 1];
    let mut c = vec![0.0f64; lim + 1];
    let cap = (u as usize).min(lim);
    for (d, &sign) in mu.iter().enumerate().take(cap + 1).skip(1) {
        if sign == 0 {
            continue;
        }
        let m = sign as f64;
        let ld = (d as f64).ln();
        let mut l = d;
        while l <= lim {
            b[l] += m;
            c[l] += m * ld;
            l += d;
        }
    }
    let blog: Vec<f64> = (0..=lim)
        .map(|l| {
            if l == 0 {
                0.0
            } else {
                (l as f64).ln() * b[l] - c[l]
            }
        })
        .collect();
    (b, blog)
}

pub fn signs(b: &[f64]) -> Vec<f64> {
    b.iter()
        .map(|&v| {
            if v > 0.0 {
                1.0
            } else if v < 0.0 {
                -1.0
            } else {
                0.0
            }
        })
        .collect()
}

// THE ROW

pub struct VaughanRow {
    pub r: u128,
    pub form: Vec<f64>,
    pub diag: Vec<f64>,
    pub bmax: Vec<f64>,
    pub nz: Vec<usize>,
}

pub fn vaughan_row(
    bits: &Bits,
    mm: u64,
    nn: u64,
    g: u128,
    x: u128,
    coefs: &[Vec<f64>],
) -> VaughanRow {
    let nb = coefs.len();
    let delta = g as f64 / x as f64;
    let bsum: Vec<f64> = coefs
        .iter()
        .map(|b| (nn..2 * nn).map(|l| b[l as usize]).sum())
        .collect();
    let mut cprime = vec![0u32; nn as usize];
    let mut form = vec![0.0f64; nb];
    let mut u = vec![0.0f64; nb];
    let mut r: u128 = 0;
    for m in mm..2 * mm {
        u.iter_mut().for_each(|v| *v = 0.0);
        for l in nn..2 * nn {
            if bits.get((m * l) as usize) {
                r += 1;
                cprime[(l - nn) as usize] += 1;
                for (j, b) in coefs.iter().enumerate() {
                    u[j] += b[l as usize];
                }
            }
        }
        for j in 0..nb {
            let v = u[j] - delta * bsum[j];
            form[j] += v * v;
        }
    }
    let mut diag = vec![0.0f64; nb];
    let mut bmax = vec![0.0f64; nb];
    let mut nz = vec![0usize; nb];
    for (j, b) in coefs.iter().enumerate() {
        let mut d = 0.0f64;
        for l in nn..2 * nn {
            let v = b[l as usize];
            if v != 0.0 {
                nz[j] += 1;
            }
            if v.abs() > bmax[j] {
                bmax[j] = v.abs();
            }
            let c = cprime[(l - nn) as usize] as f64;
            d += v * v * ((1.0 - delta) * (1.0 - delta) * c + delta * delta * (mm as f64 - c));
        }
        diag[j] = d;
    }
    VaughanRow {
        r,
        form,
        diag,
        bmax,
        nz,
    }
}

// THE STUDY

const BAND: (f64, f64) = (0.2343, 2.4727);

pub fn run() {
    println!("menergy signed vaughan");
    println!("| q | F | L | box | live/seen | M | N | R | U3 | U5 | 2N/U5 | nz3 | nz5 | bmax3 | bmax5 | fd3 | fd5 | fdlog3 | fdlog5 | fdsgn | w5 | bf5 | raw5 |");
    let mut count = 0usize;
    let mut inband = 0usize;
    let mut under = 0usize;
    let mut inrange = 0usize;
    let mut boxrows = 0usize;
    let mut fdlo = f64::INFINITY;
    let mut fdhi = f64::NEG_INFINITY;
    let mut bflo = f64::INFINITY;
    let mut bfhi = f64::NEG_INFINITY;
    let mut rawlo = f64::INFINITY;
    let mut rawhi = f64::NEG_INFINITY;
    let mut nzlo = f64::INFINITY;
    let mut nzhi = f64::NEG_INFINITY;
    let mut fd5 = (f64::INFINITY, f64::NEG_INFINITY);
    let mut fd5l5 = (f64::INFINITY, f64::NEG_INFINITY);
    let mut bfl5 = (f64::INFINITY, f64::NEG_INFINITY);
    let mut pml5 = 0usize;
    let mut piece = [0usize; 5];
    let mut edge = f64::INFINITY;
    let mut msmall = 0usize;
    let mut nzall = (f64::INFINITY, f64::NEG_INFINITY);
    let mut bmaxall = 0.0f64;
    let mut bmaxhi = 0.0f64;
    let mut plusminus = 0usize;
    let mut wlo = f64::INFINITY;
    let mut whi = f64::NEG_INFINITY;
    for c in cells().iter() {
        for &l in c.depths.iter() {
            let x = qpow(c.q, l);
            let mut lim = 16usize;
            for &(_, nn) in boxes(x).iter() {
                lim = lim.max(2 * nn as usize - 1);
            }
            let vals = column(c.q, &c.digits, l);
            let g = vals.len() as u128;
            let bits = bits_of(&vals, x);
            let ones = vec![vec![1i8; lim + 1]];
            let sw = match sweep(&bits, x, g, &ones) {
                Some(s) => s,
                None => continue,
            };
            let spf = spf_table(lim);
            let mu = mobius(lim, &spf);
            let u3 = iroot(x, 3);
            let u5 = iroot(x * x, 5);
            let (b3, b3l) = pieces(&mu, lim, u3);
            let (b5, b5l) = pieces(&mu, lim, u5);
            let sg = signs(&b5);
            let coefs = vec![b3, b5, b3l, b5l, sg];
            for (tag, row) in [Some(("top", &sw.top)), sw.l5.as_ref().map(|r| ("l5", r))]
                .into_iter()
                .flatten()
            {
                let v = vaughan_row(&bits, row.m, row.n, g, x, &coefs);
                assert_eq!(v.r, row.r, "R disagrees with the signed sweep");
                let delta = g as f64 / x as f64;
                let floor = (1.0 - delta) * (row.m as f64 * v.r as f64).sqrt();
                let fd: Vec<f64> = (0..coefs.len())
                    .map(|j| {
                        if v.diag[j] > 0.0 {
                            v.form[j] / v.diag[j]
                        } else {
                            0.0
                        }
                    })
                    .collect();
                let raw = (row.m as f64 * v.form[1]).sqrt() / floor;
                let bf = if v.bmax[1] > 0.0 {
                    raw / v.bmax[1]
                } else {
                    0.0
                };
                let w = v.diag[1]
                    / (v.bmax[1] * v.bmax[1] * (1.0 - delta) * (1.0 - delta) * v.r as f64);
                assert!(
                    (bf - (fd[1] * w).sqrt()).abs() < 1e-9 * (1.0 + bf),
                    "bound/floor is not (form/diag times the diagonal share)^(1/2)"
                );
                let share = |j: usize| v.nz[j] as f64 / row.n as f64;
                boxrows += 1;
                if row.n > u5 {
                    inrange += 1;
                }
                if row.m > u5 {
                    msmall += 1;
                }
                for j in 0..coefs.len() {
                    nzall.0 = nzall.0.min(share(j));
                    nzall.1 = nzall.1.max(share(j));
                    bmaxall = bmaxall.max(v.bmax[j]);
                }
                for j in 0..coefs.len() {
                    count += 1;
                    fdlo = fdlo.min(fd[j]);
                    fdhi = fdhi.max(fd[j]);
                    if fd[j] >= BAND.0 && fd[j] <= BAND.1 {
                        inband += 1;
                        piece[j] += 1;
                    }
                    edge = edge.min((fd[j] - BAND.0).abs().min((fd[j] - BAND.1).abs()));
                    if fd[j] < 0.1 {
                        under += 1;
                    }
                }
                bflo = bflo.min(bf);
                bfhi = bfhi.max(bf);
                rawlo = rawlo.min(raw);
                rawhi = rawhi.max(raw);
                nzlo = nzlo.min(share(0).min(share(1)));
                nzhi = nzhi.max(share(0).max(share(1)));
                fd5.0 = fd5.0.min(fd[1]);
                fd5.1 = fd5.1.max(fd[1]);
                if tag == "l5" {
                    fd5l5.0 = fd5l5.0.min(fd[1]);
                    fd5l5.1 = fd5l5.1.max(fd[1]);
                    bfl5.0 = bfl5.0.min(bf);
                    bfl5.1 = bfl5.1.max(bf);
                    if v.bmax[1] == 1.0 {
                        pml5 += 1;
                    }
                }
                bmaxhi = bmaxhi.max(v.bmax[0].max(v.bmax[1]));
                wlo = wlo.min(w);
                whi = whi.max(w);
                if v.bmax[1] == 1.0 {
                    plusminus += 1;
                }
                println!(
                    "| {} | {} | {} | {} | {}/{} | {} | {} | {} | {} | {} | {:.2} | {:.4} | {:.4} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} |",
                    c.q,
                    c.label,
                    l,
                    tag,
                    sw.live,
                    sw.seen,
                    row.m,
                    row.n,
                    v.r,
                    u3,
                    u5,
                    2.0 * row.n as f64 / u5 as f64,
                    share(0),
                    share(1),
                    v.bmax[0],
                    v.bmax[1],
                    fd[0],
                    fd[1],
                    fd[2],
                    fd[3],
                    fd[4],
                    w,
                    bf,
                    raw
                );
            }
        }
    }
    println!("menergy signed vaughan summary");
    println!("| box rows | N > U5 | M > U5 | coefficient cells | fd lo | fd hi | inside band | under 0.1 | bf lo | bf hi | raw lo | raw hi | nz lo | nz hi | bmax hi | bmax5 = 1 | at l5 | w lo | w hi | fd5 lo | fd5 hi | fd5 l5 lo | fd5 l5 hi | bf l5 lo | bf l5 hi | in band by piece | band edge | nz all | bmax all |");
    println!(
        "| {} | {} | {} | {} | {:.4} | {:.4} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} {} {} {} {} | {:.4} | {:.4} {:.4} | {:.4} |",
        boxrows, inrange, msmall, count, fdlo, fdhi, inband, under, bflo, bfhi, rawlo, rawhi, nzlo, nzhi,
        bmaxhi, plusminus, pml5, wlo, whi, fd5.0, fd5.1, fd5l5.0, fd5l5.1, bfl5.0, bfl5.1, piece[0], piece[1], piece[2],
        piece[3], piece[4], edge, nzall.0, nzall.1, bmaxall
    );
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_vaughan_pieces_match_a_direct_divisor_sum() {
        let lim = 1200usize;
        let spf = spf_table(lim);
        let mu = mobius(lim, &spf);
        for u in [7u64, 40, 300] {
            let (b, bl) = pieces(&mu, lim, u);
            for l in 1..=lim {
                let mut s = 0.0f64;
                let mut sl = 0.0f64;
                for (d, &sign) in mu.iter().enumerate().take(l + 1).skip(1) {
                    if l % d == 0 && (d as u64) <= u {
                        s += sign as f64;
                        sl += sign as f64 * ((l / d) as f64).ln();
                    }
                }
                assert!((b[l] - s).abs() < 1e-9, "b at l={l} U={u}");
                assert!(
                    (bl[l] - sl).abs() < 1e-9 * (1.0 + sl.abs()),
                    "blog at l={l} U={u}"
                );
            }
        }
    }

    #[test]
    fn the_vaughan_form_matches_the_pair_definition() {
        let (q, digits, l, mm, nn) = (3u64, vec![0u64, 1], 8usize, 8u64, 16u64);
        let x = qpow(q, l);
        let vals = column(q, &digits, l);
        let g = vals.len() as u128;
        let bits = bits_of(&vals, x);
        let lim = 2 * nn as usize;
        let spf = spf_table(lim);
        let mu = mobius(lim, &spf);
        let u = iroot(x, 3);
        let (b, blog) = pieces(&mu, lim, u);
        let sg = signs(&b);
        let coefs = vec![b, blog, sg];
        let row = vaughan_row(&bits, mm, nn, g, x, &coefs);
        let delta = g as f64 / x as f64;
        let psi = |m: u64, l: u64| {
            if bits.get((m * l) as usize) {
                1.0 - delta
            } else {
                -delta
            }
        };
        for (j, c) in coefs.iter().enumerate() {
            let mut full = 0.0f64;
            let mut dg = 0.0f64;
            for l1 in nn..2 * nn {
                for l2 in nn..2 * nn {
                    let mut w = 0.0f64;
                    for m in mm..2 * mm {
                        w += psi(m, l1) * psi(m, l2);
                    }
                    full += c[l1 as usize] * c[l2 as usize] * w;
                    if l1 == l2 {
                        dg += c[l1 as usize] * c[l1 as usize] * w;
                    }
                }
            }
            assert!(
                (full - row.form[j]).abs() < 1e-9 * (1.0 + full.abs()),
                "form {j}: {full} against {}",
                row.form[j]
            );
            assert!(
                (dg - row.diag[j]).abs() < 1e-9 * (1.0 + dg.abs()),
                "diag {j}: {dg} against {}",
                row.diag[j]
            );
        }
    }

    #[test]
    fn the_mu_piece_form_is_the_exact_integer() {
        let (q, digits, l, mm, nn) = (5u64, vec![0u64, 1, 2, 3], 5usize, 16u64, 32u64);
        let x = qpow(q, l);
        let vals = column(q, &digits, l);
        let g = vals.len() as u128;
        let bits = bits_of(&vals, x);
        let lim = 2 * nn as usize;
        let spf = spf_table(lim);
        let mu = mobius(lim, &spf);
        let u = iroot(x * x, 5);
        let (b, _) = pieces(&mu, lim, u);
        let row = vaughan_row(&bits, mm, nn, g, x, std::slice::from_ref(&b));
        let bi: Vec<i128> = (0..=lim).map(|i| b[i].round() as i128).collect();
        for i in 1..=lim {
            assert_eq!(bi[i] as f64, b[i], "piece is not integral at {i}");
        }
        let xi = x as i128;
        let gi = g as i128;
        let bs: i128 = (nn..2 * nn).map(|l| bi[l as usize]).sum();
        let mut exact = 0i128;
        for m in mm..2 * mm {
            let mut um = 0i128;
            for l in nn..2 * nn {
                if bits.get((m * l) as usize) {
                    um += bi[l as usize];
                }
            }
            let v = xi * um - gi * bs;
            exact += v * v;
        }
        let want = exact as f64 / (xi as f64 * xi as f64);
        assert!(
            (want - row.form[0]).abs() < 1e-9 * (1.0 + want.abs()),
            "{want} against {}",
            row.form[0]
        );
    }
}
