// [[file:../trajectory.note::d3e7a8f3][d3e7a8f3]]
#![feature(test)]
#![feature(result_option_inspect)]
// for lindemann file.stream_len API call
#![feature(seek_convenience)]
#![feature(seek_stream_len)]

#[cfg(test)]
#[macro_use]
extern crate approx;
// d3e7a8f3 ends here

// [[file:../trajectory.note::16fef675][16fef675]]
/// Command line tools
pub mod cli;

mod cp2k;
mod lammps;
mod lindemann;
// mod part;
mod reaction;
mod rings;

// mod atoms;
// mod graph;
// mod topo;

// pub mod adhoc;
// pub mod xyz;

mod common {
    pub use std::collections::HashMap;
    pub use std::path::{Path, PathBuf};

    pub use gut::cli::*;
    pub use gut::prelude::*;
}
use common::*;
// 16fef675 ends here

// [[file:../trajectory.note::61448511][61448511]]
#[cfg(feature = "adhoc")]
/// Docs for local development
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    // export_doc!(part);
    export_doc!(reaction);
    export_doc!(lammps);
    export_doc!(lindemann);
}
// 61448511 ends here

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
