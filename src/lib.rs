// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
extern crate petgraph;
extern crate cgmath;

mod atoms;
mod graph;
mod lammps;

use std::collections::HashMap;
pub use lammps::{extract_frame, analyze_frames};

#[cfg(test)]
#[macro_use] extern crate approx;
// imports:1 ends here

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
