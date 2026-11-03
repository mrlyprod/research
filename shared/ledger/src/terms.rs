use super::{Axis, Key, Measure};
use mrly_core::error::{value_error, Result};
use mrly_core::tensor::Tensor;
use mrly_math::bang::{code_to_corners, factory, Code};
use mrly_math::counts::{fill_polynomial, odd_power, Closed, Placements};
use mrly_math::{arcs, counts};
use mrly_math::{six, three, two};

fn grid(number: usize, dimension: usize, level: u32) -> Option<u128> {
    (number as u128).checked_pow((dimension as u32).checked_mul(level)?)
}

fn tile(key: &Key, number: usize) -> Result<Tensor> {
    factory::create(Code::from(key.code), number, key.dimension, key.base, 1)
}

fn complex(key: &Key, number: usize, level: u32) -> Result<i128> {
    let Key {
        code,
        dimension,
        base,
        measure,
        ..
    } = *key;
    let value = match dimension {
        2 => {
            let tally = two::census(&two::create(
                Code::from(code),
                number,
                level as usize,
                0,
                base,
            )?)?;
            match measure {
                Measure::Vertices => tally.vertices as i128,
                Measure::Edges => tally.edges as i128,
                Measure::Euler => tally.euler as i128,
                _ => return value_error("the plane carries no faces."),
            }
        }
        _ => {
            let tally = three::census(&three::create(
                Code::from(code),
                number,
                level as usize,
                base,
            )?)?;
            match measure {
                Measure::Vertices => tally.vertices as i128,
                Measure::Edges => tally.edges as i128,
                Measure::Faces => tally.faces as i128,
                _ => tally.euler as i128,
            }
        }
    };
    Ok(value)
}

fn face(measure: Measure) -> Option<Option<usize>> {
    match measure {
        Measure::Vertices => Some(Some(0)),
        Measure::Edges => Some(Some(1)),
        Measure::Faces => Some(Some(2)),
        Measure::Euler => Some(None),
        _ => None,
    }
}

fn skeletal(key: &Key) -> Result<Option<(Placements, Vec<i128>)>> {
    let Some(k) = face(key.measure).filter(|_| key.axis == Axis::Level) else {
        return Ok(None);
    };
    let placements = Placements::of_tile(&tile(key, key.number())?)?;
    Ok(placements
        .recurrence(k)
        .map(|coefficients| (placements, coefficients)))
}

fn skeleton_terms(key: &Key, placements: &Placements, count: usize) -> (Vec<i128>, bool) {
    let mut out = Vec::with_capacity(count);
    for level in 1..=count as u32 {
        let value = placements.at(level).and_then(|skeleton| match key.measure {
            Measure::Euler => Some(skeleton.euler),
            Measure::Vertices => i128::try_from(skeleton.vertices()).ok(),
            Measure::Edges => i128::try_from(skeleton.edges()).ok(),
            _ => i128::try_from(skeleton.faces()).ok(),
        });
        match value {
            Some(value) => out.push(value),
            None => return (out, true),
        }
    }
    (out, false)
}

fn term(key: &Key, number: usize, level: u32, cells: u128) -> Result<Option<i128>> {
    let Key {
        code,
        dimension,
        base,
        measure,
        ..
    } = *key;
    let fill = || -> Result<Option<u128>> {
        Ok(counts::fill(Code::from(code), number, dimension, 1, base)?.checked_pow(level))
    };
    let within = |cost: Option<u128>| cost.is_some_and(|cost| cost <= cells);
    let value = match measure {
        Measure::Fills => fill()?.and_then(|fill| i128::try_from(fill).ok()),
        Measure::Voids => grid(number, dimension, level)
            .zip(fill()?)
            .map(|(grid, fill)| (grid - fill) as i128),
        Measure::Surface => {
            let corners = code_to_corners(Code::from(code), dimension, base)?;
            counts::Exposure::from_corners(&corners, number, dimension, base)
                .at(level)
                .map(|faces| faces as i128)
        }
        Measure::Peak | Measure::Heights => {
            let span = (number as u128)
                .checked_pow(level)
                .map(|side| dimension as u128 * (side - 1) + 1);
            if !within(span) {
                return Ok(None);
            }
            let counts = counts::profile_of_tile(&tile(key, number)?, level)?;
            Some(if measure == Measure::Peak {
                counts.iter().max().copied().unwrap_or(0) as i128
            } else {
                counts.iter().filter(|&&count| count > 0).count() as i128
            })
        }
        Measure::Vertices | Measure::Edges | Measure::Faces | Measure::Euler => {
            if !within(grid(number, dimension, level)) {
                return Ok(None);
            }
            Some(complex(key, number, level)?)
        }
        Measure::Triangles | Measure::Gaps => {
            let cost = (number as u128)
                .checked_pow(2 * level)
                .and_then(|side| side.checked_mul(16));
            if !within(cost) {
                return Ok(None);
            }
            Some(if measure == Measure::Triangles {
                counts::cut_fills(Code::from(code), number, level)?
            } else {
                counts::cut_voids(Code::from(code), number, level)?
            } as i128)
        }
        Measure::Loops => {
            if !within(grid(number, dimension, level)) {
                return Ok(None);
            }
            let side = number.pow(level);
            let render = factory::create(Code::from(code), number, 2, base, level as usize)?;
            let on: Vec<bool> = (0..side * side).map(|cell| render.at(cell) == 1).collect();
            Some(arcs::trace(side, &on)?.loops as i128)
        }
        Measure::Windows => return value_error("the windows run along the window axis."),
        Measure::Holes | Measure::Pieces => {
            if !within(grid(number, dimension, level)) {
                return Ok(None);
            }
            let slice = six::cut(&three::create(
                Code::from(code),
                number,
                level as usize,
                base,
            )?)?;
            Some(if measure == Measure::Holes {
                six::holes(&slice)? as i128
            } else {
                six::components(&slice)? as i128
            })
        }
    };
    Ok(value)
}

pub fn terms(key: &Key, count: usize, cells: u128) -> Result<(Vec<i128>, bool)> {
    code_to_corners(Code::from(key.code), key.dimension, key.base)?;
    if !key.measure.applies(key.dimension, key.base) {
        return value_error(format!(
            "{} does not read dimension {} base {}.",
            key.measure.slug(),
            key.dimension,
            key.base
        ));
    }
    if !key.measure.axes().contains(&key.axis) {
        return value_error(format!(
            "{} does not run along the {} axis.",
            key.measure.slug(),
            key.axis.slug()
        ));
    }
    if key.axis == Axis::Window {
        let number = key.number();
        let tile = tile(key, number)?;
        return Ok(crate::windows::census(tile.bytes()?, number, count, cells));
    }
    if let Some((placements, _)) = skeletal(key)? {
        return Ok(skeleton_terms(key, &placements, count));
    }
    rendered(key, count, cells)
}

fn rendered(key: &Key, count: usize, cells: u128) -> Result<(Vec<i128>, bool)> {
    let mut out = Vec::with_capacity(count);
    for index in 0..count {
        let (number, level) = key.axis.place(index, key.number());
        match term(key, number, level, cells)? {
            Some(value) => out.push(value),
            None => return Ok((out, true)),
        }
    }
    Ok((out, false))
}

pub fn closed(key: &Key) -> Result<Option<Closed>> {
    let Key {
        code,
        dimension,
        base,
        measure,
        axis,
    } = *key;
    let corners = code_to_corners(Code::from(code), dimension, base)?;
    let fill = || counts::fill(Code::from(code), key.number(), dimension, 1, base);
    let form = match (measure, axis) {
        (Measure::Fills, Axis::Level) => Closed::Power(fill()?),
        (Measure::Voids, Axis::Level) => {
            Closed::Difference((key.number() as u128).pow(dimension as u32), fill()?)
        }
        (Measure::Fills, Axis::Side) if base == 2 => {
            Closed::Polynomial(fill_polynomial(&corners, dimension))
        }
        (Measure::Voids, Axis::Side) if base == 2 => Closed::Polynomial(
            odd_power(dimension)
                .iter()
                .zip(fill_polynomial(&corners, dimension))
                .map(|(all, filled)| all - filled)
                .collect(),
        ),
        (Measure::Surface, Axis::Level) => Closed::Recurrence(
            counts::Exposure::from_corners(&corners, key.number(), dimension, base).recurrence(),
        ),
        (Measure::Vertices | Measure::Edges | Measure::Faces | Measure::Euler, Axis::Level) => {
            match skeletal(key)? {
                Some((_, coefficients)) => Closed::Recurrence(coefficients),
                None => return Ok(None),
            }
        }
        _ => return Ok(None),
    };
    Ok(Some(form))
}

#[cfg(test)]
mod tests {
    use super::super::BUDGET;
    use super::*;

    fn read(code: u128, dimension: usize, measure: Measure, axis: Axis, count: usize) -> Vec<i128> {
        terms(&Key::new(code, dimension, 2, measure, axis), count, BUDGET)
            .unwrap()
            .0
    }

    #[test]
    fn the_closed_measures_read_the_classics() {
        assert_eq!(read(7, 2, Measure::Voids, Axis::Level, 3), [1, 17, 217]);
        assert_eq!(
            read(7, 2, Measure::Surface, Axis::Level, 4),
            [16, 80, 496, 3536]
        );
        assert_eq!(
            read(23, 3, Measure::Fills, Axis::Side, 4),
            [20, 81, 208, 425]
        );
        assert_eq!(
            read(23, 3, Measure::Voids, Axis::Side, 4),
            [7, 44, 135, 304]
        );
        assert_eq!(read(1, 2, Measure::Fills, Axis::Side, 3), [4, 9, 16]);
        assert_eq!(read(9, 2, Measure::Fills, Axis::Side, 3), [5, 13, 25]);
        assert_eq!(read(11, 2, Measure::Fills, Axis::Side, 3), [7, 19, 37]);
        assert_eq!(read(15, 2, Measure::Fills, Axis::Side, 3), [9, 25, 49]);
        assert_eq!(read(129, 3, Measure::Fills, Axis::Side, 3), [9, 35, 91]);
        assert_eq!(read(255, 3, Measure::Fills, Axis::Side, 3), [27, 125, 343]);
    }

    #[test]
    fn the_grid_measures_agree_with_the_censuses() {
        assert_eq!(read(7, 2, Measure::Euler, Axis::Level, 3), [0, -8, -72]);
        assert_eq!(read(23, 3, Measure::Euler, Axis::Level, 2), [-4, -80]);
        assert_eq!(read(23, 3, Measure::Faces, Axis::Level, 1), [96]);
        assert_eq!(read(23, 3, Measure::Vertices, Axis::Level, 1), [64]);
        assert_eq!(read(23, 3, Measure::Edges, Axis::Level, 1), [144]);
        assert_eq!(read(7, 2, Measure::Vertices, Axis::Side, 2), [16, 36]);
        assert_eq!(read(23, 3, Measure::Triangles, Axis::Level, 2), [42, 306]);
        assert_eq!(read(23, 3, Measure::Pieces, Axis::Side, 2), [1, 7]);
        assert_eq!(read(255, 3, Measure::Holes, Axis::Level, 2), [0, 0]);
        assert_eq!(read(23, 3, Measure::Peak, Axis::Level, 2), [6, 42]);
        assert_eq!(read(23, 3, Measure::Heights, Axis::Level, 2), [7, 25]);
    }

    #[test]
    fn the_budget_caps_honestly() {
        let key = Key::new(23, 3, 2, Measure::Euler, Axis::Side);
        assert_eq!(terms(&key, 8, 124).unwrap(), (vec![-4], true));
        assert_eq!(terms(&key, 8, 26).unwrap(), (vec![], true));
        let closed = Key::new(23, 3, 2, Measure::Euler, Axis::Level);
        assert_eq!(
            terms(&closed, 4, 26).unwrap(),
            (vec![-4, -80, -1408, -26432], false)
        );
        let deep = Key::new(23, 3, 2, Measure::Fills, Axis::Level);
        let (fills, capped) = terms(&deep, 40, BUDGET).unwrap();
        assert!(capped && fills.len() == 29);
        let base3 = Key::new(100, 2, 3, Measure::Surface, Axis::Side);
        assert_eq!(terms(&base3, 2, BUDGET).unwrap().0.len(), 2);
        assert!(terms(&Key::new(7, 2, 2, Measure::Faces, Axis::Level), 1, BUDGET).is_err());
        assert!(terms(&Key::new(16, 2, 2, Measure::Fills, Axis::Level), 1, BUDGET).is_err());
    }

    #[test]
    fn the_arcs_and_the_windows_read_the_carpet() {
        let loops = Key::new(495, 2, 3, Measure::Loops, Axis::Level);
        assert_eq!(
            terms(&loops, 8, BUDGET).unwrap(),
            (vec![0, 3, 50, 509, 4444], true)
        );
        let carpet = [2, 10, 40, 74, 152, 242, 344, 442];
        for (code, base) in [(495, 3), (7, 2)] {
            let windows = Key::new(code, 2, base, Measure::Windows, Axis::Window);
            assert_eq!(
                terms(&windows, 8, BUDGET).unwrap(),
                (carpet.to_vec(), false)
            );
        }
        assert!(terms(&Key::new(7, 2, 2, Measure::Windows, Axis::Level), 1, BUDGET).is_err());
        assert!(terms(&Key::new(7, 2, 2, Measure::Fills, Axis::Window), 1, BUDGET).is_err());
    }

    #[test]
    fn the_closed_forms_match_the_terms() {
        for code in 0..16u128 {
            for measure in [Measure::Fills, Measure::Voids] {
                let key = Key::new(code, 2, 2, measure, Axis::Side);
                let Some(Closed::Polynomial(poly)) = closed(&key).unwrap() else {
                    panic!("no polynomial");
                };
                let (terms, _) = terms(&key, 5, BUDGET).unwrap();
                for (index, &term) in terms.iter().enumerate() {
                    let k = index as i128 + 2;
                    let value: i128 = poly
                        .iter()
                        .enumerate()
                        .map(|(p, &c)| c * k.pow(p as u32))
                        .sum();
                    assert_eq!(value, term, "code={code} k={k}");
                }
            }
        }
        assert_eq!(
            closed(&Key::new(23, 3, 2, Measure::Fills, Axis::Level)).unwrap(),
            Some(Closed::Power(20))
        );
        assert_eq!(
            closed(&Key::new(7, 2, 2, Measure::Voids, Axis::Level)).unwrap(),
            Some(Closed::Difference(9, 8))
        );
        assert_eq!(
            closed(&Key::new(7, 2, 3, Measure::Fills, Axis::Side)).unwrap(),
            None
        );
        assert_eq!(
            closed(&Key::new(23, 3, 2, Measure::Holes, Axis::Level)).unwrap(),
            None
        );
        assert_eq!(
            closed(&Key::new(7, 2, 2, Measure::Euler, Axis::Side)).unwrap(),
            None
        );
    }

    #[test]
    fn the_skeleton_closes_the_rendered_face_counts() {
        let faces = [
            Measure::Vertices,
            Measure::Edges,
            Measure::Faces,
            Measure::Euler,
        ];
        let mut rows = 0;
        let mut agreed = 0;
        let mut designs = 0;
        let mut repeated = 0;
        let mut powered = 0;
        for (dimension, base) in crate::SPACES {
            if dimension != 2 && dimension != 3 {
                continue;
            }
            for &code in crate::designs(dimension, base).unwrap() {
                designs += 1;
                let number = base.max(3);
                let measures: Vec<Measure> = faces
                    .into_iter()
                    .filter(|measure| measure.applies(dimension, base))
                    .collect();
                let read: Vec<Vec<i128>> = measures
                    .iter()
                    .map(|&measure| {
                        let key = Key::new(code, dimension, base, measure, Axis::Level);
                        let Some(Closed::Recurrence(coefficients)) = closed(&key).unwrap() else {
                            panic!("no recurrence for {}", key.name());
                        };
                        let (terms, capped) = terms(&key, crate::TERMS, BUDGET).unwrap();
                        assert!(!capped && terms.len() == crate::TERMS, "{}", key.name());
                        for index in coefficients.len()..terms.len() {
                            let replay: i128 = coefficients
                                .iter()
                                .enumerate()
                                .map(|(back, &c)| c * terms[index - back - 1])
                                .sum();
                            assert_eq!(replay, terms[index], "{}", key.name());
                        }
                        let (placements, _) = skeletal(&key).unwrap().unwrap();
                        let k = face(measure).unwrap();
                        let roots = placements.roots(k);
                        if roots.iter().any(|&(_, multiplicity)| multiplicity > 1) {
                            repeated += 1;
                        }
                        if roots
                            .iter()
                            .any(|&(root, multiplicity)| root > 0 && multiplicity > 1)
                        {
                            powered += 1;
                        }
                        rows += 1;
                        terms
                    })
                    .collect();
                for level in 1..=crate::TERMS as u32 {
                    if !grid(number, dimension, level).is_some_and(|cells| cells <= BUDGET) {
                        break;
                    }
                    let census: Vec<i128> = if dimension == 2 {
                        let tally = two::census(
                            &two::create(Code::from(code), number, level as usize, 0, base)
                                .unwrap(),
                        )
                        .unwrap();
                        vec![
                            tally.vertices as i128,
                            tally.edges as i128,
                            tally.euler as i128,
                        ]
                    } else {
                        let tally = three::census(
                            &three::create(Code::from(code), number, level as usize, base).unwrap(),
                        )
                        .unwrap();
                        vec![
                            tally.vertices as i128,
                            tally.edges as i128,
                            tally.faces as i128,
                            tally.euler as i128,
                        ]
                    };
                    for (slot, measure) in measures.iter().enumerate() {
                        assert_eq!(
                            read[slot][level as usize - 1],
                            census[slot],
                            "code={code} dimension={dimension} base={base} {} level={level}",
                            measure.slug()
                        );
                        agreed += 1;
                    }
                }
            }
        }
        assert_eq!(
            (designs, rows, agreed, repeated, powered),
            (859, 2599, 10404, 644, 220)
        );
        let moved = crate::keys(crate::Tier::Closed)
            .into_iter()
            .filter(|key| faces.contains(&key.measure))
            .count();
        assert_eq!(moved, rows);
        assert!(crate::keys(crate::Tier::LevelGrid)
            .iter()
            .all(|key| !faces.contains(&key.measure)));
    }
}
