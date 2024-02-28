// [[file:../../trajectory.note::65c83c1c][65c83c1c]]
#[derive(Debug, Clone)]
/// User options for reaction analysis
pub struct ReactionOptions {
    /// Read trajectory stepping by this number at each frame
    pub read_trajectory_step_by: usize,
    /// The noise event life parameter used in noising removing algorithm.
    pub noise_event_life: usize,
    /// Write reaction species in `reaction-species` and `reactive-frames` directories.
    pub write_reaction_species: bool,
    /// Analysis read in frames in chunk with size of this number.
    pub chunk_size: usize,
    /// Read lattice from xyz title in extxyz format (Lattice=*)
    pub read_lattice_extxyz: bool,
}

impl Default for ReactionOptions {
    fn default() -> Self {
        Self {
            read_trajectory_step_by: 1,
            noise_event_life: 50,
            write_reaction_species: false,
            chunk_size: 150,
            read_lattice_extxyz: true,
        }
    }
}
// 65c83c1c ends here
