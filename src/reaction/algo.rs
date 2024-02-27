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
fn find_reactions(mols: &[Molecule], states: &BondingStates, noise_event_life: usize) -> Result<Vec<Reaction>> {
    let mol_indices = states.find_reactive_bonds(noise_event_life);
    let mut reactions = vec![];
    for [i, j] in mol_indices {
        let mi = &mols[i];
        let mj = &mols[j];
        let mut reaction = super::get_reaction(mi, mj)?;
        reaction.local_frame = j;
        reaction.global_frame = mj.title();
        reactions.push(reaction);
    }
    Ok(reactions)
}
// c617a958 ends here

// [[file:../../trajectory.note::2ebc3172][2ebc3172]]
use super::io::{Reaction, ReactionWriter};

pub fn find_chemical_reactions_in_trajectory(f: &str) -> Result<()> {
    use std::collections::VecDeque;

    let n_chunk: usize = 100;
    let overlap_size: usize = 5;

    let mut mols = gchemol::io::read(f)?;
    let mut window = VecDeque::new();
    let mut ichunk = 0;
    // write reactions in parquet format
    let mut writer = ReactionWriter::new("/tmp/reaction.pq".as_ref())?;
    for (i, mut mol) in mols.enumerate() {
        mol.set_title(format!("{i}"));
        window.push_back(mol);
        if window.len() == n_chunk {
            // Process the current window
            // Create one contiguous slice of `Molecule`
            println!("Processing chunk {ichunk}");
            let reactions = process_mol_chunk(window.make_contiguous())?;
            writer.write_reactions(&reactions);
            // Prepare for the next window: keep the last `overlap_size` elements
            while window.len() > overlap_size {
                window.pop_front();
            }
            ichunk += 1;
        }
    }
    writer.close()?;

    // NOTE: we ignore the last chunk
    // process_mol_chunk(window.make_contiguous());

    Ok(())
}

fn get_chemical_reactions(mols: &[Molecule], noise_event_life: usize) -> Result<Vec<Reaction>> {
    let mut mols = get_active_molecules(&mols)?;
    let mut states = remove_inactive_bonding_pairs(&mols);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    println!("When noise events removed:");
    let bonds_to_repair = remove_noise_bonding_events(&mut states, noise_event_life);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }

    repair_bonding_states(&mut mols, &bonds_to_repair);
    let n = mols.len();
    let mols_ = get_active_molecules(&mols[noise_event_life..n - noise_event_life])?;
    find_reactions(&mols, &states, noise_event_life)
}

/// Create bonds and find chemical reactions in this chunk
fn process_mol_chunk(chunk: &mut [Molecule]) -> Result<Vec<Reaction>> {
    chunk.into_par_iter().for_each(|mut mol| {
        // ignore molecules already `rebond` in overlap region
        if mol.nbonds() == 0 {
            mol.rebond();
        }
    });

    get_chemical_reactions(chunk, 5)
}

#[test]
fn test_xyz_chunks() -> Result<()> {
    let f = "/tmp/a.xyz";
    find_chemical_reactions_in_trajectory(f)?;

    Ok(())
}
// 2ebc3172 ends here

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
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    let bonds_to_repair = remove_noise_bonding_events(&mut states, noise_event_life);
    assert_eq!(states.len(), 74);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    repair_bonding_states(&mut mols, &bonds_to_repair);
    let n = mols.len();
    let mols_ = get_active_molecules(&mols[noise_event_life..n - noise_event_life])?;
    assert_eq!(mols_.len(), 2);

    let x = find_reactions(&mols, &states, noise_event_life);
    dbg!(x);

    Ok(())
}
// dd2f60bb ends here
