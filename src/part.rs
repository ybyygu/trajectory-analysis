// [[file:../trajectory.note::04fdc2f4][04fdc2f4]]
use crate::common::*;
// 04fdc2f4 ends here

// [[file:../trajectory.note::61ace94b][61ace94b]]
#[derive(Debug, Clone)]
pub struct Image {
    /// The image ID
    image: usize,
    /// The central atom ID
    atom1: usize,
    /// The atom ID nearby the central atom `atom1`
    atom2: usize,
    /// The distance between `atom1` and `atom2`
    distance: f64,
    /// The bond valence for `atom1--atom2` pair
    bond_valence: f64,
}
// 61ace94b ends here
