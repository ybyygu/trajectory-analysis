// [[file:../trajectory.note::*imports][imports:1]]
use super::*;

use gut::fs::*;
// imports:1 ends here

// [[file:../trajectory.note::*cp2k][cp2k:1]]
fn get_lattice_from_cp2k_input(file: &Path) -> Result<gchemol::Lattice> {
    // quick and dirty
    let s = read_file(file)?;

    for line in s.lines() {
        let line = line.trim();
        if line.starts_with("ABC") {
            let parts: Result<Vec<_>> = line.split_whitespace().map(|x| x.parse()).collect();
        }
    }

    todo!();
}
// cp2k:1 ends here

// [[file:../trajectory.note::*core][core:1]]
#[test]
fn test_coordination_number() -> Result<()> {
    let input_file = "./data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlOnvtpbc.inp";
    let xyz_file = "./data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz";

    // let mols = gchemol::io::read(xyz_file)?;
    let lat = get_lattice_from_cp2k_input(input_file)??;

    Ok(())
}
// core:1 ends here
