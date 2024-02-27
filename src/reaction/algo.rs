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
    println!("processed {m} frames, found {n} bonding events");

    Ok(frames)
}
// 22ecaa51 ends here

// [[file:../../trajectory.note::f367b105][f367b105]]
fn remove_inactive_bonding_pairs(mols: &[Molecule]) -> BondingStates {
    let mut states = BondingStates::from_molecules(mols);
    let n = states.remove_inactive_bonding_pairs();
    println!("Removed {n} inactive bonding paris.");
    states
}
// f367b105 ends here

// [[file:../../trajectory.note::*step3][step3:1]]
fn remove_noise_bonding_events(states: &mut BondingStates) -> Result<()> {
    todo!()
}
// step3:1 ends here

// [[file:../../trajectory.note::dd2f60bb][dd2f60bb]]
#[test]
fn test_reaction_algo() -> Result<()> {
    let f = "tests/files/lty.xyz";

    let mut mols: Vec<_> = gchemol::io::read(f)?.collect();
    // create bonds before create trajectory
    for (i, mol) in mols.iter_mut().enumerate() {
        mol.rebond();
        let frame = i + 1;
        mol.set_title(format!("{frame}"));
    }

    let mols = get_active_molecules(&mols)?;
    let mut states = remove_inactive_bonding_pairs(&mols);
    remove_noise_bonding_events(&mut states);
    // let mut states = BondingStates::from_molecules(&mols);
    // dbg!(states.nframes());
    // dbg!(states.len());
    // let keys: Vec<_> = states.bonding_pairs().collect();
    // for &[u, v] in &keys {
    // println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    // }

    // let traj = get_bonding_events_trajectory(&mols)?;
    // let events = get_bonding_events(&traj, None)?;
    // assert_eq!(events.len(), 76);
    // let events = get_bonding_events(&traj, 5)?;
    // assert_eq!(events.len(), 3);
    // events.print();

    Ok(())
}
// dd2f60bb ends here
