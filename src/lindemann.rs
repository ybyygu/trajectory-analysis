// [[file:../trajectory.note::*imports][imports:1]]
use std::path::{Path, PathBuf};

use crate::lammps_::*;

use gut::prelude::*;
// imports:1 ends here

// [[file:../trajectory.note::a90d113a][a90d113a]]
fn calculate_distance_matrix(frame: &LammpsTrajectoryFrame, settings: &config::Settings) -> Vec<f64> {
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
// a90d113a ends here

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

    impl Configure for Settings {}

    pub(crate) fn load_settigns_from_config_file() -> Settings {
        Settings::load_from_file("lindemann.conf")
    }

    #[test]
    fn test_settings() {
        Settings::default().print_toml();
    }
}
// fa617f7c ends here

// [[file:../trajectory.note::*core][core:1]]
use indicatif::ProgressBar;

fn lindemann_process_frames(
    // path to a LAMMPS trajectory file
    trjfile: &Path,
    // number of atoms per frame
    natoms: usize,
    // estimated number of frames in trajectory file
    estimated_nframes: usize,
    // user settings for each atom such as element symbol, mass, etc
    settings: &config::Settings,
) -> Vec<(f64, f64)> {
    let frames = parse_lammps_dump_file(trjfile);
    let npairs = {
        let n = natoms as f64;
        let m = n * (n - 1.0) / 2.0;
        m as usize
    };

    let mut average_com = [0.0; 3];
    // array of mean distances to average center of mass.
    let mut array_mean_dist_com = vec![0.0; natoms];
    let mut array_mean = vec![0.0; npairs];
    let mut array_delta = vec![0.0; npairs];
    let mut array_var = vec![0.0; npairs];
    let mut nframes: f64 = 0.0;

    // setup progress bar
    let bar = ProgressBar::new(estimated_nframes as u64)
        .with_style(indicatif::ProgressStyle::default_bar().progress_chars("#>-"));
    for frame in frames {
        nframes += 1.0;
        let distances1 = calculate_distance_matrix(&frame, settings);
        let distances2 = calculate_distances_center_of_mass(&frame, settings);
        for i in 0..npairs {
            let xn = distances1[i];
            let mean = array_mean[i];
            let var = array_var[i];
            let delta = xn - mean;
            // update mean
            let im = mean + delta / nframes;
            array_mean[i] = im;
            // update variance
            array_var[i] = var + delta * (xn - im);
        }

        // calculate mean distance to com
        for i in 0..natoms {
            let xn = distances2[i];
            let mean = array_mean_dist_com[i];
            let delta = xn - mean;
            array_mean_dist_com[i] = mean + delta / nframes;
        }

        // update progress bar
        bar.inc(1);
    }
    bar.finish();

    // dbg!(array_mean_dist_com);

    let cv_rij: Vec<_> = (0..npairs)
        .into_par_iter()
        .map(|i| {
            let mean = array_mean[i];
            let var = array_var[i];
            let cv = (var / nframes).sqrt() / mean;
            cv
        })
        .collect();

    // start calculate mean over atom pairs
    let pairs: Vec<_> = (0..natoms).combinations(2).collect();
    (0..natoms)
        .zip(array_mean_dist_com.into_iter())
        .map(|(i, di)| {
            // find neighbors for atom i in pairs
            let it = pairs
                .iter()
                .enumerate()
                .filter_map(|(k, p)| if p[0] == i || p[1] == i { Some(cv_rij[k]) } else { None });
            let qi = stats::mean(it);
            (di, qi)
        })
        .collect()
}
// core:1 ends here

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

// [[file:../trajectory.note::*test][test:1]]
#[test]
#[ignore]
fn test_linermann() {
    use approx::*;

    let fname = "tests/files/lammps-test.dump";
    let natoms = 537;
    let settings = config::Settings::default();
    let indices = lindemann_process_frames(fname.as_ref(), natoms, 0, &settings);

    let (d0, q0) = indices[0];
    assert_relative_eq!(q0, 0.01002696, epsilon = 1e-4);
}
// test:1 ends here
