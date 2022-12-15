// [[file:../trajectory.note::04fdc2f4][04fdc2f4]]
use crate::common::*;
// 04fdc2f4 ends here

// [[file:../trajectory.note::61ace94b][61ace94b]]
use parquet_derive::ParquetRecordWriter;

/// Row data for writing in Parquet format.
#[derive(Debug, Clone, ParquetRecordWriter)]
pub struct RowData {
    /// The image ID
    image: usize,
    /// The central atom ID
    atom1: usize,
    /// The atom ID nearby the central atom `atom1`
    atom2: usize,
    /// The distance between `atom1` and `atom2`
    distance: f64,
    /// The bond valence between `atom1` and `atom2`
    bond_valence: Option<f64>,
}
// 61ace94b ends here

// [[file:../trajectory.note::c35cdbe7][c35cdbe7]]
use gchemol::neighbors::{Neighbor, Neighborhood};
use gchemol::{Atom, Molecule};

fn eval_bond_valence(atom1: &Atom, atom2: &Atom, distance: f64) -> f64 {
    let r1 = atom1.get_cov_radius().expect("no cov radius");
    let r2 = atom2.get_cov_radius().expect("no cov radius");
    let r0 = r1 + r2;
    let b = 0.37;
    ((r0 - distance) / b).exp()
}

/// Return a `Neighborhood` struct for probing nearest neighbors in `mol`
///
/// N.B. The neighbor node index is defined using atom serial number
fn create_neighborhood_probe(mol: &Molecule) -> Neighborhood {
    let particles: Vec<_> = mol.atoms().map(|(i, a)| (i, a.position())).collect();
    let mut nh = gchemol::neighbors::Neighborhood::new();
    nh.update(particles);
    if let Some(lat) = mol.lattice {
        nh.set_lattice(lat.matrix().into());
    }

    nh
}

/// 从一帧结构数据中提取需要输出为 Parquet 格式的 row 数据
fn get_neighbor_data(mol: &Molecule, image_id: usize, r_cutoff: f64) -> Vec<RowData> {
    let nh = create_neighborhood_probe(mol);
    let n = mol.natoms();
    (1..=n)
        .flat_map(|i| {
            nh.neighbors(i, r_cutoff).map(move |n| RowData {
                image: image_id,
                atom1: i,
                atom2: n.node,
                distance: n.distance,
                bond_valence: eval_bond_valence(mol.get_atom_unchecked(i), mol.get_atom_unchecked(n.node), n.distance).into(),
            })
        })
        .collect()
}

/// 将轨迹`mols` 中的键连信息写入 parquet 文件 `path` 中.
pub fn write_connection_dataframe_parquet(mols: impl Iterator<Item = Molecule>, path: &Path, r_cutoff: f64) -> Result<()> {
    use parquet_tools::SimpleParquetFileWriter;

    let mut writer = SimpleParquetFileWriter::new(path);
    for (i, mol) in mols.enumerate() {
        let row_group = get_neighbor_data(&mol, i, r_cutoff);
        writer.write_row_group(row_group.as_slice())?;
    }
    writer.close()?;

    Ok(())
}
// c35cdbe7 ends here
