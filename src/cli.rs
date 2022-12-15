// [[file:../trajectory.note::9fcb47b8][9fcb47b8]]
use crate::common::*;

use gut::cli::*;
// 9fcb47b8 ends here

// [[file:../trajectory.note::b5926f52][b5926f52]]
/// Calculate averaged coordination numbers over time frames in trajectory.
///
/// # Limitations
/// * cell should be fixed for each timestep
/// * fixed number of atoms, and consistent atom numbering system for each timestep
#[derive(Parser, Debug, Clone)]
pub struct PartCli {
    #[command(flatten)]
    verbose: Verbosity,

    /// The trajectory file in xyz format.
    trjfile: PathBuf,

    /// The main input file for CP2K calculation for reading lattice data.
    #[arg(long = "input")]
    inpfile: PathBuf,

    /// The output file for writing bond valence data in Parquet format.
    #[arg(long = "output")]
    outfile: PathBuf,

    /// Only take one time frame every `n` frames in trajectory. The
    /// default is to take all time frames.
    #[arg(long = "step", default_value = "1")]
    time_frame_stepping: usize,

    /// The radius cutoff for searching neighboring atoms
    #[arg(short, default_value = "2.0")]
    r_cutoff: f64,
}

impl PartCli {
    pub fn enter_main() -> Result<()> {
        let args = PartCli::parse();
        let lat = crate::cp2k::read_lattice_from_cp2k_input(&args.inpfile)?;
        // let mols = gchemol::io::read(&args.trjfile)?.map(|mut mol| {
        if args.time_frame_stepping == 0 {
            bail!("invalid time frame step");
        }
        let n = args.time_frame_stepping;
        println!("analysis time frame stepping by the {n}");
        let mols = crate::xyztraj::read_xyz_trajectory(&args.trjfile)?
            .map(|mut mol| {
                mol.set_lattice(lat.clone());
                mol
            })
            .step_by(n);

        crate::part::write_connection_dataframe_parquet(mols, &args.outfile, args.r_cutoff)?;

        Ok(())
    }
}
// b5926f52 ends here

// [[file:../trajectory.note::3f98b8f6][3f98b8f6]]
pub use crate::lindemann::cli::LindemannCli;
// 3f98b8f6 ends here
