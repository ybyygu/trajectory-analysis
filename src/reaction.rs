// [[file:../trajectory.note::c0770bc4][c0770bc4]]
use crate::common::*;
use gchemol::Molecule;
// c0770bc4 ends here

// [[file:../trajectory.note::29d234b7][29d234b7]]
mod algo;
mod base;
mod cli;
mod io;
// 29d234b7 ends here

// [[file:../trajectory.note::707e344d][707e344d]]
/// Return reactants or products composition in string formulation. For
/// example, for two molecules in `mols`, returns "H2 + CH4"
fn get_composition<'a>(mols: impl IntoIterator<Item = &'a Molecule>) -> String {
    let formula_list: HashSet<_> = mols.into_iter().map(|mol| mol.formula()).collect();
    formula_list.into_iter().join(" + ")
}
// 707e344d ends here

// [[file:../trajectory.note::38d6eaa6][38d6eaa6]]
use io::Reaction;
use std::collections::BTreeSet as HashSet;

fn get_bonding_changes(mol1: &Molecule, mol2: &Molecule) -> [HashSet<[usize; 2]>; 2] {
    // NOTE: we ignore bond type difference
    let bonds1: HashSet<[usize; 2]> = mol1.bonds().map(|(u, v, _)| [u, v]).collect();
    let bonds2: HashSet<[usize; 2]> = mol2.bonds().map(|(u, v, _)| [u, v]).collect();
    let forming = bonds2.difference(&bonds1).copied().collect();
    let breaking = bonds1.difference(&bonds2).copied().collect();
    [forming, breaking]
}

/// Detects reaction between `mol1` and `mol2` from bond connectivity
/// changes. Returns reactants and products in list of `Molecule`
/// objects
fn get_reaction_mechanism(mol1: &Molecule, mol2: &Molecule) -> Result<([Vec<Molecule>; 2], [HashSet<[usize; 2]>; 2])> {
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

/// Return chemical reaction composition bewteen `mol1` and `mol2` by
/// analysis of bond changes between `mol1` and `mol2`.
pub fn get_reaction_composition(mol1: &Molecule, mol2: &Molecule) -> Result<String> {
    let ([reactants, products], _) = get_reaction_mechanism(mol1, mol2)?;
    let mut reaction_composition = String::new();
    if !reactants.is_empty() && !products.is_empty() {
        let reactants_composition = get_composition(&reactants);
        let products_composition = get_composition(&products);
        reaction_composition = format!("{} => {}", reactants_composition, products_composition);
    }
    Ok(reaction_composition)
}

pub fn get_reaction(mol1: &Molecule, mol2: &Molecule) -> Result<Reaction> {
    // for Molecule.fingerprint method
    use spdkit::prelude::*;

    let mut reaction = Reaction::default();

    let ([reactants, products], _) = get_reaction_mechanism(mol1, mol2)?;
    let mut reaction_composition = String::new();
    if !reactants.is_empty() && !products.is_empty() {
        reaction.reactants_composition = get_composition(&reactants);
        reaction.products_composition = get_composition(&products);
        reaction.reactants = reactants.iter().map(|mol| mol.numbers().collect_vec()).collect();
        reaction.products = products.iter().map(|mol| mol.numbers().collect_vec()).collect();
        reaction.reactants_fingerprints = reactants.iter().map(|mol| mol.fingerprint()).collect();
        reaction.products_fingerprints = products.iter().map(|mol| mol.fingerprint()).collect();
        // write reactants/products to files
        for mol in reactants.iter().chain(&products) {
            let fp = mol.fingerprint();
            let f = format!("/tmp/reaction-species/{fp}.mol2");
            io::write_molecule(f.as_ref(), &mol)?;
        }
    }

    Ok(reaction)
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
    pub fn remove_noise_bonding_events(bonding_events: &mut Vec<i32>, n_space: usize) -> (bool, Vec<(usize, i32)>) {
        // Mapping rules as a dictionary
        let event_to_code: HashMap<i32, &str> = [(0, "N"), (1, "F"), (-1, "B")].into_iter().collect();

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

// [[file:../trajectory.note::6f57ef8b][6f57ef8b]]
use gchemol::trajectory::Trajectory;

/// Bonding events for each pair of atoms
#[derive(Debug, Default, Clone)]
pub struct BondingEvents {
    inner: HashMap<[usize; 2], Vec<i32>>,
}

impl BondingEvents {
    /// Print chemical events for manual inspection
    pub fn print(&self) {
        print_bonding_events(&self.inner);
    }

    /// Insert bonding event value of atom pair `u` and `v`
    pub fn insert(&mut self, u: usize, v: usize, value: Vec<i32>) {
        if u < v {
            self.inner.insert([u, v], value);
        } else {
            self.inner.insert([v, u], value);
        }
    }

    /// Return the bonding event value for atom pair `u` and `v`
    pub fn get(&self, u: usize, v: usize) -> Option<&Vec<i32>> {
        if u < v {
            self.inner.get(&[u, v])
        } else {
            self.inner.get(&[v, u])
        }
    }

    /// The bonding events size
    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Print chemical events for manual inspection
///
/// # Parameters
/// * events: hash map with a pair of atom nodes as the key, and the
///   bonding event codes along time axis.
fn print_bonding_events(events: &HashMap<[usize; 2], Vec<i32>>) {
    if events.is_empty() {
        println!("No chemical events");
        return;
    }

    let map: HashMap<_, _> = [(0, "-"), (1, "↑"), (-1, "↓")].into_iter().collect();
    let mut keys: Vec<_> = events.keys().copied().collect();
    keys.sort();

    println!("bond pair: reaction codes");
    for [u, v] in keys {
        let signals = &events[&[u, v]];
        let codes: String = signals.iter().map(|&s| map[&s]).collect();
        println!("{:03}-{:03}: {}", u, v, codes);
    }

    let jmol_selection: Vec<_> = events.keys().flatten().map(|x| x.to_string()).collect();
    let jmol_console_command = format!("select atomno=[{}]\nselectionhalo\nlabel %i", jmol_selection.join(","));
    println!("If view reaction atoms in jmol, you can use below commands:");
    println!("{}", jmol_console_command);
}
// 6f57ef8b ends here

// [[file:../trajectory.note::*repair noising bonds][repair noising bonds:1]]

// repair noising bonds:1 ends here

// [[file:../trajectory.note::85db89af][85db89af]]
pub fn get_bonding_events_trajectory(mols: &[Molecule]) -> Result<Trajectory> {
    assert!(!mols.is_empty());

    let mut m = 0;
    let mut frames = vec![];
    for (i, pair) in mols.windows(2).enumerate() {
        let mol1 = &pair[0];
        let mol2 = &pair[1];
        // calculate bonding changes and record them
        let [forming, breaking] = get_bonding_changes(mol1, mol2);
        if !forming.is_empty() || !breaking.is_empty() {
            let mut mol2 = mol2.to_owned();
            mol2.properties.store("Bonds forming", forming);
            mol2.properties.store("Bonds breaking", breaking);
            frames.push(mol2);
        }
        m = i + 1;
    }
    let n = frames.len();
    println!("processed {m} frames, found {n} bonding events");
    let traj = Trajectory::try_from(frames)?;

    Ok(traj)
}
// 85db89af ends here

// [[file:../trajectory.note::b91fb579][b91fb579]]
// use polars::prelude::*;
pub fn get_bonding_events(traj: &Trajectory, noise_event_life: impl Into<Option<usize>>) -> Result<BondingEvents> {
    let mut bonding_event_codes = HashMap::new();
    let mut bonding_event_rows = Vec::new();
    let mut bonding_event_columns = HashSet::new();

    let nrows = traj.nframes();
    for i in 0..nrows {
        bonding_event_rows.push(i);
        let forming: HashSet<[usize; 2]> = traj.frames[i].properties.load("Bonds forming")?;
        let breaking: HashSet<[usize; 2]> = traj.frames[i].properties.load("Bonds breaking")?;
        let entry = bonding_event_codes.entry(i).or_insert(HashMap::new());
        for pair in forming {
            bonding_event_columns.insert(pair);
            entry.insert(pair, 1);
        }
        for pair in breaking {
            bonding_event_columns.insert(pair);
            entry.insert(pair, -1);
        }
    }

    let mut d = BondingEvents::default();
    let noise_event_life = noise_event_life.into();
    for col in bonding_event_columns {
        let mut values: Vec<_> = bonding_event_rows
            .iter()
            .map(|row| bonding_event_codes[row].get(&col).unwrap_or(&0))
            .copied()
            .collect();
        if let Some(noise_event_life) = noise_event_life {
            let (has_reaction, affected_frames) = noise::remove_noise_bonding_events(&mut values, noise_event_life);
            if has_reaction {
                d.insert(col[0], col[1], values);
            }
        } else {
            d.insert(col[0], col[1], values);
        }
    }
    Ok(d)
}
// b91fb579 ends here

// [[file:../trajectory.note::1276c516][1276c516]]
pub use base::BondingStates;
// 1276c516 ends here

// [[file:../trajectory.note::ec1621a3][ec1621a3]]
#[test]
fn test_reaction_xyz() -> Result<()> {
    let f = "tests/files/lty.xyz";

    let mut mols: Vec<_> = gchemol::io::read(f)?.collect();
    // create bonds before create trajectory
    for (i, mol) in mols.iter_mut().enumerate() {
        mol.rebond();
        mol.set_title(format!("frame {i}"));
    }

    let traj = get_bonding_events_trajectory(&mols)?;
    let events = get_bonding_events(&traj, None)?;
    assert_eq!(events.len(), 76);
    // events.print();
    let events = get_bonding_events(&traj, 5)?;
    assert_eq!(events.len(), 3);
    // events.print();

    Ok(())
}
// ec1621a3 ends here
