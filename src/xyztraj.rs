// [[file:../trajectory.note::e875348d][e875348d]]
use crate::common::*;

use gchemol::prelude::*;
use gchemol::Molecule;
// e875348d ends here

// [[file:../trajectory.note::2155de6b][2155de6b]]
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

fn read_line(file: &mut BufReader<File>, buf: &mut String) -> Option<()> {
    match file.read_line(buf) {
        Ok(0) => None,
        Ok(_) => Some(()),
        Err(_) => None,
    }
}

/// Read an iterator over `Molecule` from trajectory in `path` in xyz
/// format. NOTE: Each image in trajectory is assumed has the same
/// number of atoms, and lattice vectors (TV atom) are not parsed.
pub fn read_xyz_trajectory(path: &Path) -> Result<impl Iterator<Item = Molecule>> {
    use gchemol::Atom;

    let mut file = BufReader::new(File::open(path)?);

    // read the number of atoms from the frist line
    let mut buf = String::new();
    let _ = file.read_line(&mut buf)?;
    let natoms: usize = buf
        .trim()
        .parse()
        .map_err(|_| format_err!("The first line should be an integer: {:?}!", buf))?;

    let mols = (0..).map_while(move |i| {
        // skip natom line
        if i > 0 {
            read_line(&mut file, &mut buf)?;
        }
        // skip title line
        read_line(&mut file, &mut buf)?;
        buf.clear();
        // read the remaining natoms lines
        for _ in 0..natoms {
            read_line(&mut file, &mut buf)?;
        }

        // construct molecule from text stream
        // let mol = Molecule::from_str(&buf, "text/pxyz").ok()?;
        let atoms: Vec<Atom> = buf.lines().filter_map(|line| line.parse().ok()).collect();
        let mol = Molecule::from_atoms(atoms);
        buf.clear();
        Some(mol)
    });

    Ok(mols)
}
// 2155de6b ends here

// [[file:../trajectory.note::6b46ac30][6b46ac30]]
#[test]
fn test_read_xyz() {
    let file = "data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz";
    // let mols = gchemol::io::read(file).unwrap();

    let mols = read_xyz_trajectory(file.as_ref()).unwrap();

    assert_eq!(mols.count(), 8);
}
// 6b46ac30 ends here
