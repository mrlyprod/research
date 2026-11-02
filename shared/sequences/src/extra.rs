use crate::row::{timed, Row, TERMS};
use mrly_core::error::Result;
use mrly_life::source::{counts as source_counts, sequence as source_sequence};
use mrly_life::Source;
use mrly_math::bang::{baseq, factory, Code};
use mrly_math::counts::{self, ladder, six};
use mrly_math::{arcs, press};

const LIMIT: usize = 4000;
const LONG: usize = 24;

fn wide(value: u128) -> Option<i128> {
    i128::try_from(value).ok()
}

pub fn ladders() -> Vec<Row> {
    let mut rows = Vec::new();
    for base in [3, 5, 7] {
        for dimension in 2..=6 {
            let mut terms = Vec::new();
            for levels in (1..TERMS).rev() {
                if let Ok(found) = ladder::ladder(base, dimension, levels) {
                    terms = found;
                    break;
                }
            }
            let stated = (base == 3 && dimension == 3).then(|| "A299916".to_string());
            rows.push(
                Row::new(
                    format!("ladder_base={base}_dim={dimension}"),
                    "mrly_math::counts::ladder::ladder",
                    "ladder",
                    terms,
                )
                .stated(stated),
            );
        }
        rows.push(Row::new(
            format!("ladder_fill_base={base}"),
            "mrly_math::counts::ladder::fill",
            "ladder-fill",
            timed(TERMS, |i| ladder::fill(base, i + 1).ok()),
        ));
        rows.push(Row::new(
            format!("ladder_sign_base={base}"),
            "mrly_math::counts::ladder::sign",
            "ladder-sign",
            timed(TERMS, |i| ladder::sign(base, i + 1).ok().map(i128::from)),
        ));
    }
    rows
}

pub fn hexes() -> Vec<Row> {
    type Slice = fn(usize) -> Result<u128>;
    let slices: [(&str, Slice); 7] = [
        ("core_nodes", six::solid_slice_core_nodes),
        ("core_edges", six::solid_slice_core_edges),
        ("triangles", six::solid_slice_triangles),
        ("boundary", six::solid_slice_boundary),
        ("vertices", six::solid_slice_vertices),
        ("interior", six::solid_slice_interior),
        ("edges", six::solid_slice_edges),
    ];
    let mut rows: Vec<Row> = slices
        .iter()
        .map(|(slug, slice)| {
            Row::new(
                format!("solid_slice_{slug}"),
                &format!("mrly_math::counts::six::solid_slice_{slug}"),
                &format!("solid-slice-{slug}"),
                timed(TERMS, |i| slice(2 * i + 1).ok().and_then(wide)),
            )
        })
        .collect();
    rows.push(
        Row::new(
            "centered_hexagonal",
            "mrly_math::counts::centered_hexagonal",
            "centered-hexagonal",
            timed(TERMS, |i| wide(counts::centered_hexagonal(i + 1))),
        )
        .stated(Some("A003215".to_string())),
    );
    rows.push(Row::new(
        "grid_triangles_side=3_axis=level",
        "mrly_math::counts::six::grid_triangles",
        "grid-triangles/level",
        timed(TERMS, |i| wide(six::grid_triangles(3, i as u32 + 1))),
    ));
    rows.push(Row::new(
        "grid_triangles_axis=side",
        "mrly_math::counts::six::grid_triangles",
        "grid-triangles/side",
        timed(TERMS, |i| wide(six::grid_triangles(2 * i + 3, 1))),
    ));
    rows
}

pub fn profiles() -> Vec<Row> {
    let mut rows = Vec::new();
    for (code, dimension, level) in [(7u128, 2usize, 3u32), (23, 3, 2), (126, 3, 2), (9, 2, 3)] {
        let Ok(tile) = factory::create(Code::from(code), 3, dimension, 2, 1) else {
            continue;
        };
        let Ok(profile) = counts::profile_of_tile(&tile, level) else {
            continue;
        };
        rows.push(Row::new(
            format!("profile_dim={dimension}_code={code}_level={level}"),
            "mrly_math::counts::profile_of_tile",
            "profile",
            profile.into_iter().filter_map(wide).collect(),
        ));
    }
    rows
}

fn burnside_fits(base: usize, dimension: usize) -> bool {
    let cells = (base as f64).powi(dimension as i32);
    let maps = if base == 2 { 2.0 } else { 2.0 * base as f64 };
    let factorial: f64 = (1..=dimension).map(|i| i as f64).product();
    cells + (maps.powi(dimension as i32) * factorial).log2().ceil() <= 127.0
}

fn classes_fit(dimension: usize) -> bool {
    let mut row = vec![1f64];
    for _ in 0..dimension {
        let mut next = vec![1f64; row.len() + 1];
        for k in 1..row.len() {
            next[k] = row[k - 1] + row[k];
        }
        row = next;
    }
    row.iter().map(|c| (c + 1.0).log2()).sum::<f64>() < 127.0
}

pub fn designs() -> Vec<Row> {
    let mut rows = Vec::new();
    for base in 2..=5 {
        let stated = (base == 2).then(|| "A000616".to_string());
        rows.push(
            Row::new(
                format!("design_counts_base={base}"),
                "mrly_math::bang::baseq::sequence",
                "design-counts",
                timed(TERMS, |i| {
                    burnside_fits(base, i + 1)
                        .then(|| baseq::distinct_designs(base, i + 1).ok())
                        .flatten()
                        .and_then(wide)
                }),
            )
            .stated(stated),
        );
    }
    rows.push(
        Row::new(
            "bracelets",
            "mrly_math::bang::baseq::bracelets",
            "bracelets",
            timed(TERMS, |i| {
                baseq::distinct_designs(i + 1, 1).ok().and_then(wide)
            }),
        )
        .stated(Some("A000029".to_string())),
    );
    rows.push(
        Row::new(
            "fill_classes",
            "mrly_math::bang::baseq::class_sequence",
            "fill-classes",
            timed(TERMS, |i| {
                classes_fit(i + 1)
                    .then(|| wide(baseq::classes(i + 1)))
                    .flatten()
            }),
        )
        .stated(Some("A129824".to_string())),
    );
    rows
}

pub fn presses() -> Vec<Row> {
    let mut rows = Vec::new();
    let mut spaces: Vec<(usize, usize, u128)> = Vec::new();
    for base in 2..=6usize {
        for code in 1..(1u128 << base) {
            spaces.push((1, base, code));
        }
    }
    for code in 1..16u128 {
        spaces.push((2, 2, code));
    }
    if let Ok(reps) = baseq::representatives(3, 2) {
        spaces.extend(reps.into_iter().map(|(c, _)| (2, 3, c.get())));
    }
    for (dimension, base, code) in spaces {
        let Ok(members) = press::members(Code::from(code), dimension, base, LONG) else {
            continue;
        };
        rows.push(Row::new(
            format!("press_dim={dimension}_base={base}_code={code}"),
            "mrly_math::press::members",
            "press",
            members.into_iter().filter_map(wide).collect(),
        ));
    }
    rows
}

pub fn loops() -> Vec<Row> {
    let mut rows = Vec::new();
    let mut spaces: Vec<(usize, u128, usize)> = (1..16u128).map(|c| (2usize, c, 10usize)).collect();
    if let Ok(reps) = baseq::representatives(3, 2) {
        spaces.extend(reps.into_iter().map(|(c, _)| (3, c.get(), 6)));
    }
    for (base, code, top) in spaces {
        rows.push(Row::new(
            format!("arcs_base={base}_code={code}"),
            "mrly_math::arcs::draw",
            "arcs-loops",
            timed(top, |i| {
                arcs::draw(Code::from(code), base, i + 1)
                    .ok()
                    .map(|a| a.loops as i128)
            }),
        ));
        let Ok(Some(law)) = arcs::law(Code::from(code), base, 1) else {
            continue;
        };
        rows.push(
            Row::new(
                format!("arcs_law_base={base}_code={code}"),
                "mrly_math::arcs::law",
                "arcs-law",
                timed(TERMS, |i| {
                    arcs::law(Code::from(code), base, i + 1)
                        .ok()
                        .flatten()
                        .map(|l| l.loops as i128)
                }),
            )
            .closed(Some(law.formula)),
        );
    }
    rows
}

pub fn sources() -> Vec<Row> {
    let mut rows = Vec::new();
    for source in Source::all() {
        let take = |raw: Vec<usize>| -> Vec<i128> {
            raw.into_iter().take(LONG).map(|v| v as i128).collect()
        };
        if let Ok(raw) = source_sequence(source, LIMIT) {
            rows.push(
                Row::new(
                    format!("source_{}", source.name()),
                    "mrly_life::source::sequence",
                    "source",
                    take(raw),
                )
                .stated(source.oeis()),
            );
        }
        if let Ok(raw) = source_counts(source, LIMIT, false, false) {
            rows.push(Row::new(
                format!("source_counts_{}", source.name()),
                "mrly_life::source::counts",
                "source-counts",
                take(raw),
            ));
        }
    }
    rows
}

pub fn all() -> Vec<Row> {
    let mut rows = ladders();
    rows.extend(hexes());
    rows.extend(profiles());
    rows.extend(designs());
    rows.extend(presses());
    rows.extend(loops());
    rows.extend(sources());
    rows
}
