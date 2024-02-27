// [[file:../../trajectory.note::d74a391a][d74a391a]]
use crate::common::*;

use gut::cli::*;
// d74a391a ends here

// [[file:../../trajectory.note::7a9dfc6b][7a9dfc6b]]
use super::*;
use gut::cli::*;
use gut::config::*;

/// Analysis of reactive trajectory in xyz/extxyz format.
#[derive(Debug, Parser)]
pub struct ReactionCli {
    /// The trajectory file in xyz format.
    trjfile: PathBuf,

    /// Write reaction species (if eached, these files can be found in
    /// the same dir as trajectory file).
    #[clap(short = 'w')]
    write_reaction_species: bool,

    #[command(flatten)]
    verbose: Verbosity,
}

impl ReactionCli {
    pub fn enter_main() -> Result<()> {
        let args = Self::parse();
        args.verbose.setup_logger();

        process(&args)?;

        Ok(())
    }
}
// 7a9dfc6b ends here

// [[file:../../trajectory.note::093b2b9c][093b2b9c]]
fn process(cli: &ReactionCli) -> Result<()> {
    use crate::reaction::algo::find_chemical_reactions_in_trajectory;

    find_chemical_reactions_in_trajectory(&cli.trjfile, cli.write_reaction_species)?;

    Ok(())
}
// 093b2b9c ends here
