// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::68b8f3aa-b3f8-43c0-8b4d-c3165b146535][68b8f3aa-b3f8-43c0-8b4d-c3165b146535]]
extern crate petgraph;
extern crate cgmath;

#[macro_use]
extern crate approx;

mod atoms;
mod graph;
mod lammps;

use std::collections::HashMap;
pub use lammps::{extract_frame, analyze_frames};
// 68b8f3aa-b3f8-43c0-8b4d-c3165b146535 ends here

// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::593685c3-b67d-43c2-95f8-c9fc86f6c5f8][593685c3-b67d-43c2-95f8-c9fc86f6c5f8]]
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
// 593685c3-b67d-43c2-95f8-c9fc86f6c5f8 ends here
