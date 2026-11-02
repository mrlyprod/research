use mrly_math::counts::ladder::{
    characteristic, determinant, digit_polynomial, even_block, ladder,
};

type Out<T> = Result<T, String>;

fn ok<T, E: std::fmt::Display>(r: Result<T, E>) -> Out<T> {
    r.map_err(|e| e.to_string())
}

// ROUTES

fn automaton(dim: usize, levels: usize) -> Out<Vec<i128>> {
    ok(ladder(3, dim, levels))
}

fn product(dim: usize, levels: usize) -> Out<Vec<i128>> {
    let poly = ok(digit_polynomial(3, dim))?;
    let mut held = vec![1i128];
    let mut out = vec![1i128];
    let mut step = 1usize;
    for _ in 0..levels {
        let mut next = vec![0i128; held.len() + (poly.len() - 1) * step];
        for (i, &h) in held.iter().enumerate() {
            if h == 0 {
                continue;
            }
            for (k, &p) in poly.iter().enumerate() {
                let slot = &mut next[i + k * step];
                *slot = h
                    .checked_mul(p)
                    .and_then(|v| slot.checked_add(v))
                    .ok_or("overflow")?;
            }
        }
        held = next;
        step *= 3;
        let target = dim * (step - 1) / 2;
        out.push(held[target]);
    }
    Ok(out)
}

fn kept(dim: usize) -> Vec<usize> {
    (0..3usize.pow(dim as u32))
        .filter_map(|code| {
            let (mut c, mut sum, mut middles) = (code, 0, 0);
            for _ in 0..dim {
                let d = c % 3;
                middles += usize::from(d == 1);
                sum += d;
                c /= 3;
            }
            (middles <= 1).then_some(sum)
        })
        .collect()
}

fn cells(sums: &[usize], levels: usize, weight: usize, height: usize, target: usize) -> u64 {
    if levels == 0 {
        return u64::from(height == target);
    }
    sums.iter()
        .filter(|&&s| height + s * weight <= target)
        .map(|&s| cells(sums, levels - 1, weight * 3, height + s * weight, target))
        .sum()
}

fn enumeration(dim: usize, budget: f64) -> Vec<i128> {
    let sums = kept(dim);
    let mut out = vec![1i128];
    let mut level = 1;
    while (sums.len() as f64).powi(level as i32) <= budget {
        let target = dim * (3usize.pow(level as u32) - 1) / 2;
        out.push(cells(&sums, level, 1, 0, target) as i128);
        level += 1;
    }
    out
}

// RECURRENCE

fn signature(dim: usize) -> Out<Vec<i128>> {
    let block = ok(even_block(3, dim))?;
    Ok(ok(characteristic(&block))?[1..]
        .iter()
        .map(|c| -c)
        .collect())
}

fn obeys(terms: &[i128], sig: &[i128]) -> Out<usize> {
    let r = sig.len();
    let mut checked = 0;
    for n in r..terms.len() {
        let mut sum = 0i128;
        for (k, &s) in sig.iter().enumerate() {
            sum = s
                .checked_mul(terms[n - 1 - k])
                .and_then(|v| sum.checked_add(v))
                .ok_or("overflow")?;
        }
        if sum != terms[n] {
            return Err(format!("recurrence breaks at level {n}"));
        }
        checked += 1;
    }
    Ok(checked)
}

fn hankel(terms: &[i128], start: usize, size: usize) -> Out<i128> {
    let mut m: Vec<Vec<i128>> = (0..size)
        .map(|i| (0..size).map(|j| terms[start + i + j]).collect())
        .collect();
    let mut sign = 1i128;
    let mut prev = 1i128;
    for k in 0..size {
        if m[k][k] == 0 {
            match (k + 1..size).find(|&r| m[r][k] != 0) {
                Some(r) => {
                    m.swap(k, r);
                    sign = -sign;
                }
                None => return Ok(0),
            }
        }
        for i in k + 1..size {
            for j in k + 1..size {
                let a = m[i][j].checked_mul(m[k][k]).ok_or("overflow")?;
                let b = m[i][k].checked_mul(m[k][j]).ok_or("overflow")?;
                m[i][j] = (a - b) / prev;
            }
        }
        prev = m[k][k];
    }
    Ok(sign * m[size - 1][size - 1])
}

fn join(v: &[i128]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

// VERBS

fn ladders(dims: &[usize]) -> Out<()> {
    for &dim in dims {
        let block = ok(even_block(3, dim))?;
        let terms = automaton(dim, 200)?;
        let sig = signature(dim)?;
        let size = sig.len();
        println!("dim {dim}");
        println!("  even block {block:?}");
        println!("  det {}", ok(determinant(&block))?);
        println!("  signature ({})", join(&sig));
        println!("  hankel from level 0 {}", hankel(&terms, 0, size)?);
        println!("  hankel from level 1 {}", hankel(&terms, 1, size)?);
        println!("  terms 0..{} {}", terms.len() - 1, join(&terms));
    }
    Ok(())
}

fn check(dims: &[usize]) -> Out<()> {
    let mut bad = 0;
    for &dim in dims {
        let terms = automaton(dim, 200)?;
        let poly = product(dim, 9)?;
        let brute = enumeration(dim, 3.0e7);
        let sig = signature(dim)?;
        let agree_poly = poly.iter().zip(&terms).all(|(a, b)| a == b);
        let agree_brute = brute.iter().zip(&terms).all(|(a, b)| a == b);
        let checked = obeys(&terms, &sig)?;
        bad += usize::from(!agree_poly || !agree_brute);
        println!(
            "dim {dim}: automaton levels 0..{}, polynomial product levels 0..{} {}, cell enumeration levels 0..{} {}, recurrence checked on {checked} levels",
            terms.len() - 1,
            poly.len() - 1,
            if agree_poly { "agrees" } else { "DISAGREES" },
            brute.len() - 1,
            if agree_brute { "agrees" } else { "DISAGREES" },
        );
    }
    if bad > 0 {
        return Err(format!("{bad} dimensions disagree"));
    }
    println!("check green");
    Ok(())
}

fn main() -> Out<()> {
    let verb = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    let dims = [2, 3, 4, 5, 6];
    match verb.as_str() {
        "ladders" => ladders(&dims),
        "check" => check(&dims),
        "all" => {
            check(&dims)?;
            ladders(&dims)
        }
        other => Err(format!("unknown verb {other}: ladders, check or all")),
    }
}
