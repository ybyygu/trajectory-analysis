// [[file:../trajectory.note::*imports][imports:1]]
#![feature(seek_convenience)]

#[cfg(test)]
#[macro_use] extern crate approx;

use std::collections::HashMap;
// imports:1 ends here

// [[file:../trajectory.note::*mods][mods:1]]
mod adhoc;
mod atoms;
mod graph;
mod lammps;
mod lindemann;

pub mod lammps_;
pub mod xyz;

pub mod common {
    pub use gut::prelude::*;
    pub use gut::cli::*;

    // pub use quicli::prelude::*;
}
// mods:1 ends here

// [[file:../trajectory.note::*exports][exports:1]]
pub use lammps::{analyze_frames, extract_frame};
pub use lindemann::cli::enter_main as lindemann_cli;
// exports:1 ends here

// [[file:../trajectory.note::*frame][frame:1]]
pub struct Frame {
    pub timestep: usize,
    pub natoms: usize,
    pub fragments: HashMap<String, usize>,
    pub positions: HashMap<usize, [f64; 3]>,
    pub cell: [[f64; 3]; 3],
}

impl Frame {
    pub fn new() -> Frame{
        let map1: HashMap<String, usize> = HashMap::new();
        let map2: HashMap<usize, [f64; 3]> = HashMap::new();

        Frame {
            timestep: 0,
            natoms: 0,
            fragments: map1,
            positions: map2,
            cell: [[0_f64; 3]; 3],
        }
    }
}
// frame:1 ends here
