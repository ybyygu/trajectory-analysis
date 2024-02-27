// [[file:../../trajectory.note::94912fd0][94912fd0]]
use crate::common::*;

use gosh_dataset::SimpleParquetFileWriter;
// 94912fd0 ends here

// [[file:../../trajectory.note::6e775d47][6e775d47]]
#[derive(Debug, Serialize, Clone, Default)]
pub struct Reaction {
    #[serde(rename = "Local frame")]
    pub local_frame: usize,
    #[serde(rename = "Global frame")]
    pub global_frame: String,
    #[serde(rename = "Reactants")]
    pub reactants: Vec<Vec<usize>>,
    #[serde(rename = "Products")]
    pub products: Vec<Vec<usize>>,
    #[serde(rename = "Reactants composition")]
    pub reactants_composition: String,
    #[serde(rename = "Products composition")]
    pub products_composition: String,
    #[serde(rename = "Reactants fingerprints")]
    pub reactants_fingerprints: Vec<String>,
    #[serde(rename = "Products fingeprints")]
    pub products_fingerprints: Vec<String>,
}

pub struct ReactionWriter {
    writer: SimpleParquetFileWriter,
}

impl ReactionWriter {
    pub fn new(f: &Path) -> Result<Self> {
        Ok(Self {
            writer: SimpleParquetFileWriter::new(f),
        })
    }

    pub fn write_reactions(&mut self, reactions: &[Reaction]) -> Result<()> {
        self.writer.write_row_group(reactions)?;

        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.writer.close();
        Ok(())
    }
}
// 6e775d47 ends here

// [[file:../../trajectory.note::7e14952a][7e14952a]]
use gchemol::prelude::*;
use gchemol::Molecule;

/// Writes molecule `mol` to file `f`. Creates leading directories if
/// not exist.
pub fn write_molecule(f: &Path, mol: &Molecule) -> Result<()> {
    use std::fs::{self, File};

    if let Some(dir_path) = f.parent() {
        fs::create_dir_all(dir_path)?; // Create the directory path if it doesn't exist
    }
    mol.to_file(f)?;

    Ok(())
}
// 7e14952a ends here
