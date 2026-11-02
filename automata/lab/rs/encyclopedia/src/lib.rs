pub mod cache;
#[cfg(test)]
mod checks;
pub mod control;
pub mod design;
pub mod engine;
pub mod grammar;
pub mod grid;
pub mod laws;
pub mod library;
pub mod mask;
pub mod readings;
pub mod recipe;
pub mod rows;

pub use design::Shape;
pub use engine::{run, Run, Setup, Stop};
pub use grammar::Grammar;
pub use grid::Grid;
pub use mask::{Kernel, MaskSpec};
pub use readings::{readings, short};
pub use recipe::{Mode, Recipe, Side};

pub fn run_named(name: &str) -> mrly_core::Result<(Recipe, Run)> {
    let recipe = Recipe::parse(name)?;
    let setup = Setup::of(&recipe)?;
    Ok((recipe, run(&setup)))
}
