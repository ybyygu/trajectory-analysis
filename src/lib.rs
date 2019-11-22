// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
extern crate petgraph;
extern crate cgmath;

#[cfg(test)]
#[macro_use] extern crate approx;

use std::collections::HashMap;
// imports:1 ends here

// mods

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*mods][mods:1]]
mod atoms;
mod graph;
mod lammps;
mod lindermann;

pub mod xyz;
pub mod lammps_;

pub mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = ::std::result::Result<T, Error>;
}
// mods:1 ends here

// exports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*exports][exports:1]]
pub use lammps::{extract_frame, analyze_frames};
// exports:1 ends here

// frame

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*frame][frame:1]]
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
