// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use std::path::{Path, PathBuf};

use crate::lammps_::*;

use guts::prelude::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*core][core:1]]

// core:1 ends here

// distances
// aperiodic

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*distances][distances:1]]
fn calculate_distance_matrix(frame: &LammpsTrajectoryFrame) -> Vec<f64> {
    let natoms = frame.atoms.len();
    debug!("dm: found {} atoms in frame {}", natoms, frame.timestep);

    // atom id counts from 1
    let coords: Vec<[f64; 3]> = (1..natoms + 1)
        .into_par_iter()
        .map(|i| frame.atoms[&i].xyz.clone())
        .collect();

    let pairs: Vec<_> = (0..natoms).combinations(2).collect();
    pairs
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
// distances:1 ends here

// pub

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*pub][pub:1]]
use indicatif::ProgressBar;

fn lindemann_process_frames(trjfile: &Path, natoms: usize, estimated_nframes: usize) -> Vec<f64> {
    let frames = parse_lammps_dump_file(trjfile);
    let npairs = {
        let n = natoms as f64;
        let m = n * (n - 1.0) / 2.0;
        m as usize
    };

    let mut array_mean = vec![0.0; npairs];
    let mut array_delta = vec![0.0; npairs];
    let mut array_var = vec![0.0; npairs];
    let mut nframes: f64 = 0.0;

    // setup progress bar
    let bar = ProgressBar::new(estimated_nframes as u64);
    for frame in frames {
        nframes += 1.0;
        // dbg!(frame.timestep);
        let dm = calculate_distance_matrix(&frame);
        // dbg!(dm.len());
        for i in 0..npairs {
            let xn = dm[i];
            let mean = array_mean[i];
            let var = array_var[i];
            let delta = xn - mean;
            // update mean
            let im = mean + delta / nframes;
            array_mean[i] = im;
            // update variance
            array_var[i] = var + delta * (xn - im);
        }

        // update progress bar
        bar.inc(1);
    }
    bar.finish();

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
    let lindermann_indices: Vec<_> = (0..natoms)
        .into_par_iter()
        .map(|i| {
            // find neighbors for atom i in pairs
            let it = pairs
                .iter()
                .enumerate()
                .filter_map(|(k, p)| if p[0] == i || p[1] == i { Some(cv_rij[k]) } else { None });
            let qi = stats::mean(it);
            qi
        })
        .collect();

    lindermann_indices
}
// pub:1 ends here

// quick check
// 快速读取trajectory文件, 获取原子数和帧数信息?

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*quick%20check][quick check:1]]
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

// cli

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*cli][cli:1]]
use structopt::StructOpt;

use guts::cli::*;

/// Calculate Lindermann indices for LAMMPS trajectory file (.dump)
///
/// * Current limitations:
///
/// 1. PBC blind (treat as nano-particles)
/// 2. Required dump fields: x, y, z (Cartesian coordinates)
///
#[derive(Debug, StructOpt)]
struct Cli {
    /// The trajectory file in xyz format.
    #[structopt(parse(from_os_str))]
    trjfile: PathBuf,
}

pub fn lindermann_cli() -> CliResult {
    let args = Cli::from_args();

    let (natoms, nframes) = quick_check_natoms_nframes(&args.trjfile)?;
    let indices = lindemann_process_frames(&args.trjfile, natoms, nframes);

    println!("{:^8}\t{:^18}", "index(i)", "q_i");
    for (i, qi) in indices.into_iter().enumerate() {
        println!("{:^8}\t{:^-18.8}", i + 1, qi);
    }

    Ok(())
}
// cli:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*test][test:1]]
#[test]
#[ignore]
fn test_linermann() {
    use approx::*;

    let fname = "tests/files/lammps-test.dump";
    let natoms = 537;
    let indices = lindemann_process_frames(fname.as_ref(), natoms, 0);

    assert_relative_eq!(indices[0], 0.01002696, epsilon = 1e-4);
}
// test:1 ends here
