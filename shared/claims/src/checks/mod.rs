mod moire;
mod roulette;
mod spin;

use crate::registry::{Cost, Entry};

pub fn table() -> &'static [Entry] {
    &[
        (
            "flat-carpet-stack/the-moire-correlation-law-for-odd",
            Cost::Cheap,
            moire::the_moire_correlation_law_for_odd,
        ),
        (
            "flat-carpet-stack/the-stack-is-an-exact-prime",
            Cost::Cheap,
            moire::the_stack_is_an_exact_prime,
        ),
        (
            "the-spirograph-nodes/one-curve-of-a-circle-roulette",
            Cost::Cheap,
            roulette::one_curve_of_a_circle_roulette,
        ),
        (
            "the-spirograph-nodes/two-distinct-curves-of-one-wheel",
            Cost::Cheap,
            roulette::two_distinct_curves_of_one_wheel,
        ),
        (
            "the-spirograph-nodes/a-whole-design-carries-one-node",
            Cost::Cheap,
            roulette::a_whole_design_carries_one_node,
        ),
        (
            "the-spirograph-nodes/a-roulette-cuts-the-plane-into",
            Cost::Cheap,
            roulette::a_roulette_cuts_the_plane_into,
        ),
        (
            "the-spirograph-nodes/the-loop-threshold-is-not-where",
            Cost::Cheap,
            roulette::the_loop_threshold_is_not_where,
        ),
        (
            "the-spirograph-nodes/the-self-law-ends-at-the",
            Cost::Cheap,
            roulette::the_self_law_ends_at_the,
        ),
        (
            "spin/the-coprime-law-survives-the-spin",
            Cost::Cheap,
            spin::the_coprime_law_survives_the_spin,
        ),
    ]
}
