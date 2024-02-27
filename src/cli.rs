// [[file:../trajectory.note::9fcb47b8][9fcb47b8]]
use crate::common::*;

use gut::cli::*;
// 9fcb47b8 ends here

// [[file:../trajectory.note::3f98b8f6][3f98b8f6]]
pub use crate::lindemann::cli::LindemannCli;
// 3f98b8f6 ends here

// [[file:../trajectory.note::2a2c5538][2a2c5538]]
/// Count rings in trajectory file (xyz format only)
#[derive(Debug, Parser)]
pub struct CountRingsCli {
    /// The trajectory file in xyz format.
    trjfile: PathBuf,

    /// The output file for final results.
    #[arg(short)]
    outfile: PathBuf,

    /// The max rings size to be detected.
    #[arg(long = "max", short, default_value = "7")]
    maxsize: usize,

    #[command(flatten)]
    verbosity: Verbosity,
}

impl CountRingsCli {
    pub fn enter_main() -> Result<()> {
        use gchemol::prelude::*;

        let args = Self::parse();
        args.verbosity.setup_logger();

        let txt = crate::rings::count_rings_in_trajectory(args.trjfile, args.maxsize)?;
        txt.to_file(&args.outfile)?;
        println!("Done. Results saved to: {:#?}", args.outfile.display());

        Ok(())
    }
}
// 2a2c5538 ends here

// [[file:../trajectory.note::86912f49][86912f49]]
pub use crate::reaction::cli::*;
// 86912f49 ends here
