// [[file:../trajectory.note::c0770bc4][c0770bc4]]
use crate::common::*;
use gchemol::Molecule;
// c0770bc4 ends here

// [[file:../trajectory.note::707e344d][707e344d]]
/// Return reactants or products composition in string formulation. For
/// example, for two molecules in `mols`, returns "H2 + CH4"
fn get_composition<'a>(mols: impl IntoIterator<Item = &'a Molecule>) -> String {
    let formula_list: HashSet<_> = mols.into_iter().map(|mol| mol.formula()).collect();
    formula_list.into_iter().join(" ")
}
// 707e344d ends here

// [[file:../trajectory.note::38d6eaa6][38d6eaa6]]
use std::collections::BTreeSet as HashSet;

fn get_bonding_changes(mol1: &Molecule, mol2: &Molecule) -> [HashSet<[usize; 2]>; 2] {
    // NOTE: we ignore bond type difference
    let bonds1: HashSet<[usize; 2]> = mol1.bonds().map(|(u, v, bond)| [u, v]).collect();
    let bonds2: HashSet<[usize; 2]> = mol1.bonds().map(|(u, v, bond)| [u, v]).collect();
    let forming = bonds2.difference(&bonds1).copied().collect();
    let breaking = bonds1.difference(&bonds2).copied().collect();
    [forming, breaking]
}

/// Detects reaction between `mol1` and `mol2` from bond connectivity
/// changes. Returns reactants and products in list of `Molecule`
/// objects
pub(crate) fn get_reaction(mol1: &Molecule, mol2: &Molecule) -> Result<([Vec<Molecule>; 2], [HashSet<[usize; 2]>; 2])> {
    ensure!(mol1.matching_configuration(mol2), "invalid molecule pair for reaction!");

    let [forming, breaking] = get_bonding_changes(mol1, mol2);
    let mut reactants = HashSet::new();
    let mut products = HashSet::new();

    if !forming.is_empty() {
        for &[u, v] in &forming {
            let au = mol1.get_atom(u).map(|a| a.symbol()).expect("symbol u");
            let av = mol1.get_atom(v).map(|a| a.symbol()).expect("symbol v");
            let f1: HashSet<usize> = mol1.connected_fragment_atoms(u).collect();
            let f2: HashSet<usize> = mol1.connected_fragment_atoms(v).collect();
            let f3: HashSet<usize> = mol2.connected_fragment_atoms(u).collect();
            reactants.insert(f1);
            reactants.insert(f2);
            products.insert(f3);
        }
    }

    if !breaking.is_empty() {
        for &[u, v] in &breaking {
            let au = mol1.get_atom(u).map(|a| a.symbol()).expect("symbol u");
            let av = mol1.get_atom(v).map(|a| a.symbol()).expect("symbol v");
            let f1: HashSet<usize> = mol2.connected_fragment_atoms(u).collect();
            let f2: HashSet<usize> = mol2.connected_fragment_atoms(v).collect();
            let f3: HashSet<usize> = mol1.connected_fragment_atoms(u).collect();
            products.insert(f1);
            products.insert(f2);
            reactants.insert(f3);
        }
    }

    let reactants: Option<Vec<Molecule>> = reactants.iter().map(|x| mol1.get_sub_molecule(x)).collect();
    let products: Option<Vec<Molecule>> = products.iter().map(|x| mol2.get_sub_molecule(x)).collect();
    let reactants = reactants.unwrap_or_default();
    let products = products.unwrap_or_default();

    Ok(([reactants, products], [forming, breaking]))
}
// 38d6eaa6 ends here

// [[file:../trajectory.note::2972f0e5][2972f0e5]]
mod noise {
    use crate::common::*;

    /// Find start and end positions for noise event pattern (in pair) in `code`
    fn find_noise_codes(code: &str, n_space: usize) -> Vec<(usize, usize)> {
        use regex::Regex;

        let pattern = format!(r#"B(N{{0,{n_space}}})F|F(N{{0,{n_space}}})B"#);
        let re = Regex::new(&pattern).unwrap();
        let mut positions = vec![];

        for mat in re.find_iter(code) {
            positions.push((mat.start(), mat.end()));
        }

        positions
    }

    /// Find and remove matched noise event pattern
    fn remove_noise_bonding_events(bonding_events: &mut Vec<i32>, n_space: usize) -> (bool, Vec<(usize, i32)>) {
        // Mapping rules as a dictionary
        let event_to_code: HashMap<i32, &str> = [(0, "N"), (1, "F"), (-1, "B")].iter().cloned().collect();

        // Convert events to string based on mapping rules
        let event_codes: String = bonding_events.iter().map(|&event| event_to_code[&event]).collect();

        // Find noise codes
        let matched_positions = find_noise_codes(&event_codes, n_space);

        // Replace matched F/B (+1/-1) with N (0)
        let mut positions_with_noises = vec![];
        for (start, end) in matched_positions {
            for i in [start, end - 1] {
                positions_with_noises.push((i, bonding_events[i]));
                bonding_events[i] = 0;
            }
        }

        // NOTE: buffer region in both sides should not be affected by noise removing
        let length = bonding_events.len();
        let mut has_reaction = false;
        for &x in &bonding_events[n_space..length - n_space] {
            if x != 0 {
                has_reaction = true;
                break;
            }
        }

        // Excluding frames in buffer region
        let frames: Vec<(usize, i32)> = positions_with_noises
            .iter()
            .filter(|&&(i, _)| i > n_space || i < length - n_space)
            .cloned()
            .collect();

        (has_reaction, frames)
    }

    #[test]
    fn test_find_noise_events() {
        let s = "NNNBNFNNNNBNNFNNNN";
        let pairs = find_noise_codes(s, 2);
        assert_eq!(pairs.len(), 2);
        let (u, v) = pairs[0];
        assert_eq!(&s[u..v], "BNF");
        let (u, v) = pairs[1];
        assert_eq!(&s[u..v], "BNNF");
        let pairs = find_noise_codes(s, 1);
        assert_eq!(pairs.len(), 1);
        let pairs = find_noise_codes(s, 0);
        assert_eq!(pairs.len(), 0);
    }
}
// 2972f0e5 ends here

// [[file:../trajectory.note::85db89af][85db89af]]
/// Return chemical reaction composition bewteen `mol1` and `mol2` by
/// analysis of bond changes between `mol1` and `mol2`.
pub fn get_reaction_composition(mol1: &Molecule, mol2: &Molecule) -> Result<String> {
    let ([reactants, products], _) = get_reaction(mol1, mol2)?;
    let mut reaction_composition = String::new();
    if !reactants.is_empty() && !products.is_empty() {
        let reactants_composition = get_composition(&reactants);
        let products_composition = get_composition(&products);
        reaction_composition = format!("{} => {}", reactants_composition, products_composition);
    }
    Ok(reaction_composition)
}
// 85db89af ends here
