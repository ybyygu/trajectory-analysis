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
