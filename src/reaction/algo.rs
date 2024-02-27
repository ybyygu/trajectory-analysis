// [[file:../../trajectory.note::0ae3c448][0ae3c448]]
use super::base::BondingStates;
use crate::common::*;
// 0ae3c448 ends here

// [[file:../../trajectory.note::22ecaa51][22ecaa51]]
use gchemol::Molecule;

fn get_active_molecules(mols: &[Molecule]) -> Result<Vec<Molecule>> {
    assert!(!mols.is_empty());

    let mut m = 0;
    let mut frames = vec![];
    for (i, pair) in mols.windows(2).enumerate() {
        let mol1 = &pair[0];
        let mol2 = &pair[1];
        // calculate bonding changes and record them
        let [forming, breaking] = super::get_bonding_changes(mol1, mol2);
        if !forming.is_empty() || !breaking.is_empty() {
            let mut mol2 = mol2.to_owned();
            mol2.properties.store("Bonds forming", forming)?;
            mol2.properties.store("Bonds breaking", breaking)?;
            frames.push(mol2);
        }
        m = i + 1;
    }
    let n = frames.len();
    println!("processed {m} frames, found {n} reactive frames");

    Ok(frames)
}
// 22ecaa51 ends here

// [[file:../../trajectory.note::f367b105][f367b105]]
fn remove_inactive_bonding_pairs(mols: &[Molecule]) -> BondingStates {
    let mut states = BondingStates::from_molecules(mols);
    let n_pairs = states.len();
    let n = states.remove_inactive_bonding_pairs();
    println!("Removed {n} inactive bonding paris from {n_pairs} in total.");
    states
}
// f367b105 ends here

// [[file:../../trajectory.note::eb075c08][eb075c08]]
fn remove_noise_bonding_events(states: &mut BondingStates, noise_event_life: usize) -> Vec<([usize; 2], Vec<usize>)> {
    let keys: Vec<_> = states.bonding_pairs().collect();
    let mut bonds_to_repair = Vec::new();
    for key in keys {
        // the frames with noise bonding states which should be
        // inverted (bonded => unbonded, unbonded => bonded)
        let frames = states.remove_noise_events(key, noise_event_life);
        bonds_to_repair.push((key, frames));
    }
    bonds_to_repair
}
// eb075c08 ends here

// [[file:../../trajectory.note::e92233b2][e92233b2]]
fn repair_bonding_states(mols: &mut [Molecule], bonds_to_repair: &[([usize; 2], Vec<usize>)]) {
    use gchemol::Bond;

    for ([u, v], frames) in bonds_to_repair {
        let u = *u;
        let v = *v;
        for &i in frames {
            if mols[i].has_bond(u, v) {
                mols[i].remove_bond(u, v);
            } else {
                mols[i].add_bond(u, v, Bond::default());
            }
        }
    }
}
// e92233b2 ends here

// [[file:../../trajectory.note::c617a958][c617a958]]
fn find_reactions(mols: &[Molecule], states: &BondingStates, noise_event_life: usize) {
    let mol_indices = states.find_reactive_bonds(noise_event_life);
    for [i, j] in mol_indices {
        let mi = &mols[dbg!(i)];
        let mj = &mols[dbg!(j)];
        let x = super::get_reaction_composition(mi, mj);
        dbg!(x);
    }
}
// c617a958 ends here

// [[file:../../trajectory.note::dd2f60bb][dd2f60bb]]
#[test]
fn test_reaction_algo() -> Result<()> {
    let f = "tests/files/lty.xyz";
    let noise_event_life = 5;

    let mut mols: Vec<_> = gchemol::io::read(f)?.collect();
    // create bonds before create trajectory
    for (i, mol) in mols.iter_mut().enumerate() {
        mol.rebond();
        let frame = i + 1;
        mol.set_title(format!("{frame}"));
    }
    let mut mols = get_active_molecules(&mols)?;
    let mut states = remove_inactive_bonding_pairs(&mols);
    let bonds_to_repair = remove_noise_bonding_events(&mut states, noise_event_life);
    assert_eq!(states.len(), 74);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    repair_bonding_states(&mut mols, &bonds_to_repair);
    let n = mols.len();
    let mols_ = get_active_molecules(&mols[noise_event_life..n - noise_event_life])?;
    assert_eq!(mols_.len(), 3);

    let x = find_reactions(&mols, &states, noise_event_life);
    dbg!(x);
    // let traj = get_bonding_events_trajectory(&mols)?;
    // let events = get_bonding_events(&traj, None)?;
    // assert_eq!(events.len(), 76);
    // let events = get_bonding_events(&traj, 5)?;
    // assert_eq!(events.len(), 3);
    // events.print();

    Ok(())
}
// dd2f60bb ends here
