// [[file:../trajectory.note::1e216fa1][1e216fa1]]
use crate::common::*;

use gchemol::Lattice;

/// Construct a `Lattice` struct from CP2K input file in `path`.
pub fn read_lattice_from_cp2k_input(path: &Path) -> Result<Lattice> {
    // quick and dirty
    let s = gut::fs::read_file(path)?;

    let mut abc = None;
    let mut alpha_beta_gamma = None;
    for line in s.lines() {
        let line = line.trim();
        if line.starts_with("ABC") {
            let parts: std::result::Result<Vec<f64>, _> = line[3..]
                .split_whitespace()
                .map(|x| x.parse().context("parse cell lengths"))
                .collect();
            abc = parts?.into();
        } else if line.starts_with("ALPHA_BETA_GAMMA") {
            let parts: std::result::Result<Vec<f64>, _> = line[16..]
                .split_whitespace()
                .map(|x| x.parse().context("parse cell angles"))
                .collect();
            alpha_beta_gamma = parts?.into();
        }
    }
    if let Some(abc) = abc {
        if let Some(alpha_beta_gamma) = alpha_beta_gamma {
            if let [a, b, c, ..] = abc[..] {
                if let [alpha, beta, gamma, ..] = alpha_beta_gamma[..] {
                    let lat = Lattice::from_params(a, b, c, alpha, beta, gamma);
                    return Ok(lat);
                }
            }
        }
    }
    bail!("parse cell from cp2k input failure!")
}
// 1e216fa1 ends here
