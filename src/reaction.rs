// [[file:../trajectory.note::c0770bc4][c0770bc4]]
use crate::common::*;
use gchemol::Molecule;
// c0770bc4 ends here

// [[file:../trajectory.note::38d6eaa6][38d6eaa6]]
fn get_bonding_changes(mol1: &Molecule, mol2: &Molecule) -> Result<[HashSet<[usize; 2]>; 2]> {
    // NOTE: we ignore bond type difference
    let bonds1: HashSet<[usize; 2]> = mol1.bonds().map(|(u, v, bond)| [u, v]).collect();
    let bonds2: HashSet<[usize; 2]> = mol1.bonds().map(|(u, v, bond)| [u, v]).collect();
    let forming = bonds2.difference(&bonds1).copied().collect();
    let breaking = bonds1.difference(&bonds2).copied().collect();
    Ok([forming, breaking])
}

// /// Detects reaction between `mol1` and `mol2` from bond connectivity
// /// changes. Returns reactants and products in list of `Molecule`
// /// objects
// pub fn get_reaction(
//     mol1: &Molecule,
//     mol2: Molecule,
// ) -> Result<(Vec<Molecule>, Vec<Molecule>, Vec<(usize, usize)>, Vec<(usize, usize)>)> {
//     ensure!(mol1.matching_configuration(mol2), "invalid molecule pair!");

//     let (forming, breaking) = get_bonding_changes(mol1, mol2);
//     let mut reactants = HashSet::new();
//     let mut products = HashSet::new();

//     if !forming.is_empty() {
//         for &(u, v) in &forming {
//             let au = mol1.get_atom(u).symbol;
//             let av = mol1.get_atom(v).symbol;
//             let label = format!("{au}{u}--{av}{v}");
//             let f1: HashSet<usize> = mol1.connected_fragment_atoms(u).collect();
//             let f2: HashSet<usize> = mol1.connected_fragment_atoms(v).collect();
//             let f3: HashSet<usize> = mol2.connected_fragment_atoms(u).collect();
//             reactants.insert(f1);
//             reactants.insert(f2);
//             products.insert(f3);
//         }
//     }

//     if !breaking.is_empty() {
//         for &(u, v) in &breaking {
//             let au = mol1.get_atom(u).symbol;
//             let av = mol1.get_atom(v).symbol;
//             let f1: HashSet<usize> = mol2.connected_fragment_atoms(u).collect();
//             let f2: HashSet<usize> = mol2.connected_fragment_atoms(v).collect();
//             let f3: HashSet<usize> = mol1.connected_fragment_atoms(u).collect();
//             products.insert(f1);
//             products.insert(f2);
//             reactants.insert(f3);
//         }
//     }

//     let reactants: Vec<Molecule> = reactants.iter().map(mol1.get_sub_molecule).collect();
//     let products: Vec<Molecule> = products.iter().map(mol2.get_sub_molecule).collect();
//     let bond_forming: Vec<(usize, usize)> = forming.into_iter().collect();
//     let bond_breaking: Vec<(usize, usize)> = breaking.into_iter().collect();

//     (reactants, products, bond_forming, bond_breaking)
// }
// 38d6eaa6 ends here
