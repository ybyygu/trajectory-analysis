// [[file:../../trajectory.note::*imports][imports:1]]
use gut::cli::*;
use gut::prelude::*;
use std::path::{Path, PathBuf};

use structopt::StructOpt;
// imports:1 ends here

// [[file:../../trajectory.note::73f15884][73f15884]]
/// Count rings in trajectory file (xyz format only)
#[derive(Debug, StructOpt)]
struct Cli {
    /// The trajectory file in xyz format.
    #[structopt(parse(from_os_str))]
    trjfile: PathBuf,

    /// The output file for final results.
    #[structopt(parse(from_os_str), short = "o")]
    outfile: PathBuf,

    /// The max rings size to be detected.
    #[structopt(long = "max", short = "m", default_value = "7")]
    maxsize: usize,

    #[structopt(flatten)]
    verbosity: Verbosity,
}

fn main() -> Result<()> {
    use gchemol::prelude::*;

    let args = Cli::from_args();
    args.verbosity.setup_logger();

    use trajectory_analysis::xyz::*;

    let txt = count_rings_in_trajectory(args.trjfile, args.maxsize)?;

    txt.to_file(&args.outfile)?;
    println!("Done. Results saved to: {:#?}", args.outfile.display());

    Ok(())
}
// 73f15884 ends here
