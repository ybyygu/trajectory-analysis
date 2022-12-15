// [[file:../trajectory.note::*imports][imports:1]]
// for linderman
#![feature(seek_convenience)]
#![feature(seek_stream_len)]

#[cfg(test)]
#[macro_use] extern crate approx;
// imports:1 ends here

// [[file:../trajectory.note::16fef675][16fef675]]
pub mod cli;

mod cp2k;
mod part;
mod ring;
mod xyztraj;

// mod atoms;
// mod graph;
// mod lammps;
// mod lindemann;
// mod topo;

// pub mod adhoc;
// pub mod lammps_;
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

    export_doc!(part);
    export_doc!(xyztraj);
}
// 61448511 ends here

// [[file:../trajectory.note::fe8af640][fe8af640]]
#[cfg(feature = "adhoc")]
#[inline]
pub fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

// pub use lammps::{analyze_frames, extract_frame};
// pub use lindemann::cli::enter_main as lindemann_cli;
// pub use topo::topo_cli;
// fe8af640 ends here

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
