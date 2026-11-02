use crate::design::Shape;
use crate::mask::{parse_design, MaskSpec};
use mrly_core::error::value_error;
use mrly_life::{counts, Source};
use mrly_core::Result;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Mode {
    Island,
    Wallpaper,
}

impl Mode {
    pub fn token(self) -> &'static str {
        match self {
            Mode::Island => "island",
            Mode::Wallpaper => "wallpaper",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Side {
    Drawn {
        seq: Source,
        zeros: bool,
        ones: bool,
    },
    List(Vec<usize>),
}

impl Side {
    pub fn drawn(seq: Source, zeros: bool, ones: bool) -> Side {
        Side::Drawn { seq, zeros, ones }
    }

    pub fn list(values: &[usize]) -> Side {
        let mut out = values.to_vec();
        out.sort_unstable();
        out.dedup();
        Side::List(out)
    }

    pub fn token(&self) -> String {
        match self {
            Side::Drawn { seq, zeros, ones } => {
                let mut out = seq.name();
                if *zeros {
                    out.push_str("_zeros");
                }
                if *ones {
                    out.push_str("_ones");
                }
                out
            }
            Side::List(values) => values
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        }
    }

    pub fn parse(token: &str) -> Result<Side> {
        if token.starts_with(|c: char| c.is_ascii_digit()) {
            let mut values = Vec::new();
            for part in token.split(',') {
                match part.parse::<usize>() {
                    Ok(v) => values.push(v),
                    Err(_) => return value_error(format!("count {part:?} is not a number.")),
                }
            }
            return Ok(Side::list(&values));
        }
        let Some((seq, tail)) = Source::read(token) else {
            return value_error(format!("no source starts {token:?}."));
        };
        let (zeros, tail) = match tail.strip_prefix("_zeros") {
            Some(rest) => (true, rest),
            None => (false, tail),
        };
        let (ones, tail) = match tail.strip_prefix("_ones") {
            Some(rest) => (true, rest),
            None => (false, tail),
        };
        if !tail.is_empty() {
            return value_error(format!("side {token:?} has a tail {tail:?}."));
        }
        Ok(Side::drawn(seq, zeros, ones))
    }

    pub fn values(&self, budget: usize) -> Result<Vec<usize>> {
        match self {
            Side::Drawn { seq, zeros, ones } => counts(*seq, budget, *zeros, *ones),
            Side::List(values) => Ok(values.iter().copied().filter(|&v| v <= budget).collect()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Recipe {
    pub seed: Shape,
    pub side: usize,
    pub level: usize,
    pub tiling: usize,
    pub mode: Mode,
    pub mask: MaskSpec,
    pub birth: Side,
    pub survive: Side,
    pub x: usize,
    pub cap: usize,
}

impl Recipe {
    pub fn name(&self) -> String {
        let mut out = format!(
            "{}.{}.{} x{} {} mask:{} birth:{} survive:{} x:{}",
            self.seed.token(),
            self.side,
            self.level,
            self.tiling,
            self.mode.token(),
            self.mask.token(),
            self.birth.token(),
            self.survive.token(),
            self.x
        );
        if self.mode == Mode::Island {
            out.push_str(&format!(" cap:{}", self.cap));
        }
        out
    }

    pub fn parse(name: &str) -> Result<Recipe> {
        let words: Vec<&str> = name.split_whitespace().collect();
        if words.len() < 7 {
            return value_error(format!("recipe {name:?} is too short."));
        }
        let (seed, side, level) = parse_design(words[0])?;
        let Some(tiling) = words[1]
            .strip_prefix('x')
            .and_then(|t| t.parse::<usize>().ok())
        else {
            return value_error(format!("tiling {:?} wants x<n>.", words[1]));
        };
        let mode = match words[2] {
            "island" => Mode::Island,
            "wallpaper" => Mode::Wallpaper,
            other => return value_error(format!("mode {other:?} is island or wallpaper.")),
        };
        let mut mask = None;
        let mut birth = None;
        let mut survive = None;
        let mut x = None;
        let mut cap = None;
        for word in &words[3..] {
            let Some((key, value)) = word.split_once(':') else {
                return value_error(format!("word {word:?} wants key:value."));
            };
            match key {
                "mask" => mask = Some(MaskSpec::parse(value)?),
                "birth" => birth = Some(Side::parse(value)?),
                "survive" => survive = Some(Side::parse(value)?),
                "x" => x = value.parse::<usize>().ok(),
                "cap" => cap = value.parse::<usize>().ok(),
                other => return value_error(format!("unknown key {other:?}.")),
            }
        }
        let (Some(mask), Some(birth), Some(survive), Some(x)) = (mask, birth, survive, x) else {
            return value_error(format!("recipe {name:?} misses a part."));
        };
        let cap = match (mode, cap) {
            (Mode::Island, Some(cap)) => cap,
            (Mode::Island, None) => return value_error("an island names its cap."),
            (Mode::Wallpaper, None) => 0,
            (Mode::Wallpaper, Some(_)) => return value_error("a wallpaper has no cap."),
        };
        if mode == Mode::Wallpaper && tiling != 1 {
            return value_error("a wallpaper tiles once.");
        }
        Ok(Recipe {
            seed,
            side,
            level,
            tiling,
            mode,
            mask,
            birth,
            survive,
            x,
            cap,
        })
    }

    pub fn board_side(&self) -> usize {
        self.side.pow(self.level as u32) * self.tiling
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_name_parses_back_to_its_recipe_and_a_recipe_back_to_its_name() {
        let names = [
            "carpet.5.1 x3 island mask:carpet.5.1 birth:carpet_fills survive:carpet_voids x:1024 cap:301",
            "code15.3.2 x1 wallpaper mask:copyinv birth:3 survive:2,3 x:64",
            "star.7.2 x5 island mask:copy birth:random_0_zeros_ones survive:fibonacci_ones x:256 cap:201",
        ];
        for name in names {
            let recipe = Recipe::parse(name).unwrap();
            assert_eq!(recipe.name(), name);
            assert_eq!(Recipe::parse(&recipe.name()).unwrap(), recipe);
        }
        assert!(Recipe::parse("carpet.5.1 x3 wallpaper mask:copy birth:3 survive:3 x:64").is_err());
        assert!(Recipe::parse("carpet.5.1 x1 island mask:copy birth:3 survive:3 x:64").is_err());
        assert!(Side::parse("fibonacci_ones_zeros").is_err());
        assert_eq!(
            Side::parse("binary").unwrap().values(8).unwrap(),
            vec![2, 4, 8]
        );
        assert_eq!(
            Side::parse("evens_zeros").unwrap().values(4).unwrap(),
            vec![0, 2, 4]
        );
    }
}
