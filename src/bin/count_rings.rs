// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use quicli::prelude::*;
use std::path::{Path, PathBuf};

type Result<T> = ::std::result::Result<T, Error>;
use structopt::StructOpt;
// imports:1 ends here

// cmdline

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*cmdline][cmdline:1]]
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

fn main() -> CliResult {
    use gchemol_old::prelude::*;

    let args = Cli::from_args();
    args.verbosity.setup_env_logger(&env!("CARGO_PKG_NAME"))?;

    use trajectory_analysis::xyz::*;

    let txt = count_rings_in_trajectory(args.trjfile, args.maxsize)?;

    txt.to_file(&args.outfile)?;
    println!("Done. Results saved to: {:#?}", args.outfile.display());

    Ok(())
}
// cmdline:1 ends here
