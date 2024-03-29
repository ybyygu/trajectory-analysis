// [[file:../trajectory.note::9cafa605][9cafa605]]
use std::path::{Path, PathBuf};

use crate::lammps::*;

use gut::prelude::*;
use stats::OnlineStats;
// 9cafa605 ends here

// [[file:../trajectory.note::629c873b][629c873b]]
fn calculate_distance_matrix(frame: &LammpsTrajectoryFrame) -> Vec<f64> {
    let natoms = frame.atoms.len();
    debug!("dm: found {} atoms in frame {}", natoms, frame.timestep);

    // atom id counts from 1
    // let coords: Vec<[f64; 3]> = (1..natoms + 1).into_par_iter().map(|i| frame.atoms[&i].xyz.clone()).collect();

    // atom id may be not counted from 1
    let atom_ids = frame.atoms.keys().sorted().collect_vec();
    let coords: Vec<[f64; 3]> = atom_ids.into_par_iter().map(|i| frame.atoms[&i].xyz.clone()).collect();

    (0..natoms)
        .combinations(2)
        .collect_vec()
        .par_iter()
        .map(|p| {
            let i = p[0];
            let j = p[1];
            let [xi, yi, zi] = coords[i];
            let [xj, yj, zj] = coords[j];
            let dij = ((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2)).sqrt();
            dij
        })
        .collect()
}
// 629c873b ends here

// [[file:../trajectory.note::dbce8505][dbce8505]]
fn calculate_distances_center_of_mass(frame: &LammpsTrajectoryFrame, settings: &config::Settings) -> Vec<f64> {
    let natoms = frame.atoms.len();
    debug!("dm: found {} atoms in frame {}", natoms, frame.timestep);

    // atom id may be not counted from 1
    let atom_ids = frame.atoms.keys().sorted().collect_vec();

    let values: Vec<_> = atom_ids
        .iter()
        .map(|i| {
            // atom id counts from 1
            let atom = &frame.atoms[&i];
            let xyz = atom.xyz;
            // type id counts from 1
            // FIXME: not safe
            let j = atom.type_id - 1;
            let data = &settings.atoms[j];
            (data.mass, xyz)
        })
        .collect();

    let coords: Vec<_> = values.iter().map(|(_, xyz)| *xyz).collect();
    let masses: Vec<_> = values.iter().map(|(w, _)| *w).collect();

    calculate_points_radii(&coords, Some(masses))
}

/// Return the distance of points to their geometry center.
fn calculate_points_radii(coords: &[[f64; 3]], weights: Option<Vec<f64>>) -> Vec<f64> {
    let n = coords.len();
    let weights = weights.unwrap_or_else(|| vec![1.0; n]);
    assert_eq!(n, weights.len());

    // calculate the weighted center of points
    let [x0, y0, z0] = {
        let wsum: f64 = weights.iter().sum();
        let com = coords
            .iter()
            .zip(weights.into_iter())
            .fold([0.0; 3], |[acc_x, acc_y, acc_z], ([x, y, z], w)| {
                [acc_x + w * x, acc_y + w * y, acc_z + w * z]
            });
        [com[0] / wsum, com[1] / wsum, com[2] / wsum]
    };

    coords
        .iter()
        .map(|[x, y, z]| ((x - x0).powi(2) + (y - y0).powi(2) + (z - z0).powi(2)).sqrt())
        .collect()
}

#[test]
fn test_points_radii() {
    let points = vec![[0.0; 3]];

    let com = calculate_points_radii(&points, None);
    assert_eq!(com, vec![0.0]);

    let points = vec![[0.0; 3], [2.0, 0.0, 0.0]];

    let com = calculate_points_radii(&points, None);
    assert_eq!(com, vec![1.0; 2]);

    let com = calculate_points_radii(&points, Some(vec![1.008, 12.011]));
    assert_relative_eq!(com[0], 1.845149, epsilon = 1e-4);
    assert_relative_eq!(com[1], 0.154851, epsilon = 1e-4);
}
// dbce8505 ends here

// [[file:../trajectory.note::fa617f7c][fa617f7c]]
mod config {
    use gut::config::*;
    use gut::prelude::*;

    #[derive(Debug, Serialize, Deserialize, Clone)]
    pub(crate) struct Atom {
        /// Element symbol for this Atom
        pub symbol: String,

        /// Atomic mass
        pub mass: f64,
    }

    #[derive(Deserialize, Serialize, Debug)]
    /// User defined parameters for atoms
    pub(crate) struct Settings {
        /// user defined bond valence paramters
        pub atoms: Vec<Atom>,
        /// selected atoms for analysis
        pub selections: Option<Vec<usize>>,
    }

    impl Default for Settings {
        fn default() -> Self {
            let atoms = vec![
                Atom {
                    symbol: "H".into(),
                    mass: 1.008,
                },
                Atom {
                    symbol: "C".into(),
                    mass: 12.011,
                },
            ];

            Settings { atoms, selections: None }
        }
    }

    pub(crate) fn load_settigns_from_config_file() -> Settings {
        Settings::from_toml("lindemann.conf").unwrap()
    }

    #[test]
    fn test_settings() {
        Settings::default().print_toml();
    }
}
// fa617f7c ends here

// [[file:../trajectory.note::ddd39bca][ddd39bca]]
/// Compuate Lindemann index for a single atom
///
/// # Parameters
/// * traj_distances: distances to its neighbors in every frame.
pub fn compute_local_lindemann_index<Neighbors>(traj_distances: impl Iterator<Item = Neighbors>) -> f64
where
    Neighbors: IntoIterator<Item = f64>,
{
    // stats over frame
    let mut stats_array = None;
    for distances in traj_distances {
        let dists = distances.into_iter().collect_vec();
        let stats_ = stats_array.get_or_insert_with(|| vec![OnlineStats::new(); dists.len()]);
        for (i, dij) in dists.into_iter().enumerate() {
            stats_[i].add(dij);
        }
    }

    // stats over neighbors
    stats::mean(stats_array.unwrap().into_iter().map(|o| o.stddev() / o.mean()))
}
// ddd39bca ends here

// [[file:../trajectory.note::7415f651][7415f651]]
/// Compuate Lindemann indices for all atoms
///
/// # Parameters
/// * natoms: the number of atoms per frame.
/// * frames: complete pairwise distances of all atoms in each frame.
pub fn compute_lindemann_indices<Frame>(natoms: usize, frames: impl Iterator<Item = Frame>) -> impl Iterator<Item = f64>
where
    Frame: IntoIterator<Item = f64>,
{
    let npairs = {
        let n = natoms as f64;
        let m = n * (n - 1.0) / 2.0;
        m as usize
    };

    let mut stats_array = vec![OnlineStats::new(); npairs];
    for distances in frames {
        for (i, dij) in distances.into_iter().enumerate() {
            stats_array[i].add(dij);
        }
    }

    let cv_rij: Vec<_> = (0..npairs)
        .into_par_iter()
        .map(|i| stats_array[i].stddev() / stats_array[i].mean())
        .collect();

    // start calculate mean over atom pairs
    (0..natoms).map(move |i| {
        // find neighbors for atom i in pairs
        let it = (0..natoms).combinations(2).enumerate().filter_map(
            |(k, p)| {
                if p[0] == i || p[1] == i {
                    Some(cv_rij[k])
                } else {
                    None
                }
            },
        );
        stats::mean(it)
    })
}
// 7415f651 ends here

// [[file:../trajectory.note::*quick check][quick check:1]]
fn quick_check_natoms_nframes(trjfile: &Path) -> Result<(usize, usize)> {
    use std::fs::File;
    use std::io::prelude::*;
    use std::io::BufReader;
    use std::io::SeekFrom;

    let f = File::open(trjfile)?;
    let mut fp = BufReader::new(f);

    // Read in 4 lines to get the number of atoms.
    let mut buf = String::new();
    for _ in 0..4 {
        buf.clear();
        let size = fp.read_line(&mut buf)?;
        assert_ne!(size, 0);
    }
    let natoms: usize = buf.trim().parse().unwrap();

    // Read in the following lines for atoms
    for _ in 0..natoms {
        let size = fp.read_line(&mut buf)?;
        assert_ne!(size, 0);
    }

    let file_size = fp.stream_len()?;
    let part_size = fp.stream_position()?;

    let estimated_nframes = (file_size / part_size) as usize;

    Ok((natoms, estimated_nframes))
}

#[test]
fn test_quick_check() {
    let fname = "tests/files/lammps-test.dump";
    let x = quick_check_natoms_nframes(fname.as_ref()).unwrap();
    dbg!(x);
}
// quick check:1 ends here

// [[file:../trajectory.note::07c944a2][07c944a2]]
pub mod cli {
    use super::*;

    use gut::cli::*;
    use gut::config::*;

    use indicatif::ProgressBar;

    fn lindemann_process_frames(
        // path to a LAMMPS trajectory file
        trjfile: &Path,
        // number of atoms per frame
        natoms: usize,
        // estimated number of frames in trajectory file
        estimated_nframes: usize,
    ) -> Result<Vec<f64>> {
        // setup progress bar
        let bar =
            ProgressBar::new(estimated_nframes as u64).with_style(indicatif::ProgressStyle::default_bar().progress_chars("#>-"));
        let frames = parse_lammps_dump_file(trjfile)?.map(|frame| calculate_distance_matrix(&frame));
        let indices = compute_lindemann_indices(natoms, frames).inspect(|_| bar.inc(1)).collect();
        bar.finish();

        Ok(indices)
    }

    /// Calculate Lindemann indices for LAMMPS trajectory file (.dump)
    ///
    /// * Current limitations:
    ///
    /// 1. PBC blind (treat as nano-particles)
    /// 2. Required dump fields: x, y, z (Cartesian coordinates)
    #[derive(Debug, Parser)]
    pub struct LindemannCli {
        /// The trajectory file in xyz format.
        trjfile: Option<PathBuf>,

        #[command(flatten)]
        verbose: Verbosity,
    }

    impl LindemannCli {
        pub fn enter_main() -> Result<()> {
            let args = Self::parse();

            if let Some(trjfile) = args.trjfile {
                let (natoms, nframes) = quick_check_natoms_nframes(&trjfile)?;
                let indices = lindemann_process_frames(&trjfile, natoms, nframes)?;

                // FIXME: print with real atom id
                println!("{:^8}\t{:^18}", "atom index", "lindemann index");
                for (i, li) in indices.into_iter().enumerate() {
                    println!("{:^8}\t{:^-18.8}", i + 1, li);
                }
            } else {
                Self::command().print_help();
            }

            Ok(())
        }
    }
}
// 07c944a2 ends here

// [[file:../trajectory.note::f32cd037][f32cd037]]
#[test]
#[ignore]
fn test_lindemann() -> Result<()> {
    use approx::*;

    let fname = "tests/files/lammps-test.dump";
    let natoms = 537;
    let settings = config::Settings::default();

    let nneighbors = 536;
    let distances_traj =
        parse_lammps_dump_file(fname.as_ref())?.map(|frame| calculate_distance_matrix(&frame).into_iter().take(nneighbors));
    let q0 = compute_local_lindemann_index(distances_traj);
    assert_relative_eq!(q0, 0.01002696, epsilon = 1e-4);

    let frames = parse_lammps_dump_file(fname.as_ref())?.map(|frame| calculate_distance_matrix(&frame));
    let indices_ = compute_lindemann_indices(natoms, frames).collect_vec();
    let q0 = indices_[0];
    assert_relative_eq!(q0, 0.01002696, epsilon = 1e-4);

    // test against old version
    //
    // let indices = lindemann_process_frames_(fname.as_ref(), natoms, 0, &settings).unwrap();
    // for i in 0..indices.len() {
    //     assert_relative_eq!(indices[i].1, indices_[i], epsilon = 1e-4);
    // }

    Ok(())
}
// f32cd037 ends here
