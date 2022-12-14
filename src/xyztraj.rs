// [[file:../trajectory.note::e875348d][e875348d]]
use crate::common::*;

use gchemol::prelude::*;
use gchemol::Molecule;
// e875348d ends here

// [[file:../trajectory.note::2155de6b][2155de6b]]
/// Read an iterator over `Molecule` from trajectory in `path` in xyz
/// format. NOTE: Each image in trajectory is assumed has the same
/// number of atoms.
pub fn read_xyz_trajectory(path: &Path) -> Result<impl Iterator<Item = Molecule>> {
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    let mut file = BufReader::new(File::open(path)?);
    let mut buf = String::new();

    // read the number of atoms from the frist line
    let _ = file.read_line(&mut buf)?;
    let natoms: usize = buf
        .trim()
        .parse()
        .map_err(|e| format_err!("The first line should be an integer: {:?}!", buf))?;

    let mols = (1..).filter_map(move |i| {
        let _ = file.read_line(&mut buf).ok()?;
        if buf.is_empty() {
            return None;
        }
        // read the remaining natoms+1 lines
        for _ in 0..natoms + 1 {
            let _ = file.read_line(&mut buf).ok()?;
        }

        // construct molecule from text stream
        let mut mol = Molecule::from_str(&buf, "text/xyz").ok()?;
        Some(mol)
    });

    Ok(mols)
}
// 2155de6b ends here

// [[file:../trajectory.note::6b46ac30][6b46ac30]]
#[test]
fn test_read_xyz() {
    let file = "data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz";
    let mols = gchemol::io::read(file).unwrap();

    // let mols = read_xyz_trajectory(file.as_ref()).unwrap();

    assert_eq!(mols.count(), 8);
}
// 6b46ac30 ends here
