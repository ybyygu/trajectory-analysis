// [[file:../trajectory.note::*imports][imports:1]]
use super::*;

use gchemol::Molecule;
use gut::fs::*;
// imports:1 ends here

// [[file:../trajectory.note::2698f88c][2698f88c]]
use gchemol::Lattice;

fn get_lattice_from_cp2k_input(file: &Path) -> Result<Lattice> {
    // quick and dirty
    let s = gut::fs::read_file(file)?;

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
// 2698f88c ends here

// [[file:../trajectory.note::*bond][bond:1]]
use gchemol::Atom;

fn is_bonded(atom1: &Atom, atom2: &Atom, distance: f64, bonding_ratio: f64) -> bool {
    let r = bonding_ratio;
    match (atom1.get_vdw_radius(), atom2.get_vdw_radius()) {
        (Some(cr1), Some(cr2)) => {
            let rcut = (cr1 + cr2) * r;
            if distance > rcut {
                false
            } else {
                true
            }
        }
        _ => false,
    }
}
// bond:1 ends here

// [[file:../trajectory.note::*neighbors][neighbors:1]]
use gchemol::neighbors::{Neighbor, Neighborhood};

/// Return a `Neighborhood` struct for probing nearest neighbors in `mol`
///
/// N.B. The neighbor node index is defined using atom serial number
fn create_neighborhood_probe(mol: &Molecule) -> Neighborhood {
    let particles: Vec<_> = mol.atoms().map(|(i, a)| (i, a.position())).collect();
    let mut nh = gchemol::neighbors::Neighborhood::new();
    nh.update(particles);
    if let Some(lat) = mol.lattice {
        nh.set_lattice(lat.matrix().into());
    }

    nh
}

fn calculate_coordination_number(mol: &Molecule, host: usize, neighbors: &[Neighbor], bonding_ratio: f64) -> usize {
    let a1 = mol.get_atom(host).unwrap();
    neighbors
        .into_iter()
        .filter(|n| {
            let a2 = mol.get_atom(n.node).unwrap();
            is_bonded(a1, a2, n.distance, bonding_ratio)
        })
        .count()
}

fn peek_number_of_atoms_from_xyz_file(f: &Path) -> Result<usize> {
    use std::io::prelude::*;
    use std::io::BufRead;

    let fp = File::open(f)?;
    let reader = BufReader::new(fp);
    let line = reader.lines().next().ok_or(anyhow!("empty file?"))?;
    let n = line?.trim().parse().context("parse natoms")?;

    Ok(n)
}
// neighbors:1 ends here

// [[file:../trajectory.note::*core][core:1]]
fn analyze_averaged_coordination_numbers(
    cp2k_input_file: &Path,
    xyz_file: &Path,
    r_cutoff: f64,
    bonding_ratio: f64,
) -> Result<()> {
    use vecfx::*;

    let natoms = peek_number_of_atoms_from_xyz_file(xyz_file.as_ref())?;
    let mols = gchemol::io::read(xyz_file)?;
    let lat = get_lattice_from_cp2k_input(cp2k_input_file.as_ref())?;
    let numbers = (1..(natoms + 1)).collect_vec();
    let mut time_frames = vec![];

    for (i, mut mol) in mols.enumerate() {
        println!("processing frame {}", i + 1);
        mol.set_lattice(lat.clone());
        let nh = create_neighborhood_probe(&mol);
        let cn_map: std::collections::HashMap<_, _> = numbers
            .iter()
            .map(|&x| {
                let neighbors = nh.neighbors(x, r_cutoff).collect_vec();
                (x, calculate_coordination_number(&mol, x, &neighbors, bonding_ratio))
            })
            .collect();
        time_frames.push(cn_map);
    }

    println!("SN\tMean(CN)");
    for i in numbers.iter() {
        let coord_numbers_for_atom_i = time_frames.iter().map(|cn_map| cn_map[i] as f64).collect_vec();
        let averaged_coord_number_for_atom_i = coord_numbers_for_atom_i.mean();
        println!("{}\t{}", i, averaged_coord_number_for_atom_i);
    }

    Ok(())
}
// core:1 ends here

// [[file:../trajectory.note::*cli][cli:1]]
use gut::cli::*;

use structopt::*;

#[derive(Debug, StructOpt)]
/// Calculate averaged coordination numbers over time frames in trajectory.
///
/// # Limitations
///
/// * cell should have no change for each timestep
///
/// * fixed number of atoms, and consistent atom numbering system for each timestep
///
struct Cli {
    #[structopt(flatten)]
    verbose: Verbosity,

    /// The trajectory file in xyz format.
    trjfile: PathBuf,

    /// The main input file for CP2K calculation
    #[structopt(long = "input")]
    inpfile: PathBuf,

    #[structopt(long, default_value = "0.6")]
    /// The bonding ratio for defining a chemical bond
    bonding_ratio: f64,

    #[structopt(short = "r", default_value = "2.0")]
    /// The radius cutoff for searching neighboring atoms
    r_cutoff: f64,
}

pub fn topo_cli() -> Result<()> {
    let args = Cli::from_args();
    args.verbose.setup_logger();

    analyze_averaged_coordination_numbers(&args.inpfile, &args.trjfile, args.r_cutoff, args.bonding_ratio)?;

    Ok(())
}
// cli:1 ends here
