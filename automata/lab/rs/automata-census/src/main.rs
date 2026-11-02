mod card;
mod dynamics;
mod gasket;
mod groups;
mod identity;
mod rules;
mod seed;
mod table;

fn main() {
    identity::report();
    println!();
    groups::report();
    println!();
    table::report();
    println!();
    dynamics::report();
    println!();
    let census = seed::report();
    println!();
    gasket::report();
    println!();
    card::report(&census);
}
