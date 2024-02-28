// [[file:../../trajectory.note::0ae3c448][0ae3c448]]
use super::base::BondingStates;
use super::options::ReactionOptions;

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
fn find_reactions(
    mols: &[Molecule],
    states: &BondingStates,
    noise_event_life: usize,
    // the root dir for writing reaction species
    reaction_species_dir: Option<&Path>,
    reactive_frames_dir: Option<&Path>,
) -> Result<Vec<Reaction>> {
    let mol_indices = states.find_reactive_bonds(noise_event_life);
    let mut reactions = vec![];
    for [i, j] in mol_indices {
        let mi = &mols[i];
        let mj = &mols[j];
        let mut reaction = super::get_reaction(mi, mj, reaction_species_dir)?;
        reaction.local_frame = j;
        reaction.global_frame = mj.title();
        reactions.push(reaction);
        // write reactive frames for checking
        if let Some(dir) = reactive_frames_dir {
            let f = dir.join(format!("{j:03}.mol2"));
            super::io::write_molecules(&f, &[mi.clone(), mj.clone()])?;
        }
    }
    Ok(reactions)
}
// c617a958 ends here

// [[file:../../trajectory.note::2ebc3172][2ebc3172]]
use super::io::{Reaction, ReactionWriter};

pub fn find_chemical_reactions_in_trajectory(trjfile: &Path, options: &ReactionOptions) -> Result<()> {
    use std::collections::VecDeque;

    let noise_event_life = options.noise_event_life;
    let chunk_size = options.chunk_size;
    ensure!(
        chunk_size > 2 * noise_event_life + 1,
        "invalid chunk_size/noise_event_life parameters"
    );

    let overlap_size = noise_event_life;

    let step_by = options.read_trajectory_step_by;
    ensure!(step_by != 0, "invalid read_trajectory_step_by option!");
    if step_by > 1 {
        println!("Reading trajectory stepping by {step_by} at each frame.");
    }
    // Set lattice from extxyz title
    let mut mols = gchemol::io::read(trjfile)?.step_by(step_by).map(|mut mol| {
        if options.read_lattice_extxyz {
            if let Some(lat) = gchemol::io::formats::ExtxyzFile::read_lattice(&mol.title()) {
                debug!("Set lattice from extxyz title");
                mol.set_lattice(lat);
            }
        }
        mol
    });
    // write reactions in parquet format
    let pqfile = trjfile.with_file_name("reaction.pq");
    let mut writer = ReactionWriter::new(&pqfile)?;

    let mut reaction_species_dir = None;
    let mut reactive_frames_dir = None;
    let mut window = VecDeque::new();
    let mut ichunk = 0;
    if options.write_reaction_species {
        if let Some(p) = trjfile.parent() {
            reaction_species_dir = Some(p.join("reaction-species"));
            reactive_frames_dir = Some(p.join("reactive-frames"));
        }
    }
    for (i, mut mol) in mols.enumerate() {
        mol.set_title(format!("{i}"));
        window.push_back(mol);
        if window.len() == chunk_size {
            // Process the current window
            // Create one contiguous slice of `Molecule`
            println!("Processing chunk {ichunk}");
            let reactions = process_mol_chunk(
                window.make_contiguous(),
                &options,
                reaction_species_dir.as_deref(),
                reactive_frames_dir.as_deref(),
            )?;
            writer.write_reactions(&reactions);
            // Prepare for the next window: keep the last `overlap_size` elements
            while window.len() > overlap_size {
                window.pop_front();
            }
            ichunk += 1;
        }
    }

    // NOTE: ignore the last chunk or not?
    if window.len() > 2 * noise_event_life + 1 {
        process_mol_chunk(
            window.make_contiguous(),
            &options,
            reaction_species_dir.as_deref(),
            reactive_frames_dir.as_deref(),
        );
    }

    writer.close()?;

    Ok(())
}

fn get_chemical_reactions(
    mols: &[Molecule],
    noise_event_life: usize,
    // root dir for writing reaction species
    reaction_species_dir: Option<&Path>,
    reactive_frames_dir: Option<&Path>,
) -> Result<Vec<Reaction>> {
    // NOTE: this is bugging
    // let mut mols = get_active_molecules(&mols)?;
    let mut mols = mols.to_vec();
    let mut states = remove_inactive_bonding_pairs(&mols);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        info!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    info!("When noise events removed:");
    let bonds_to_repair = remove_noise_bonding_events(&mut states, noise_event_life);
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        info!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    repair_bonding_states(&mut mols, &bonds_to_repair);

    find_reactions(&mols, &states, noise_event_life, reaction_species_dir, reactive_frames_dir)
}

/// Create bonds and find chemical reactions in this chunk
fn process_mol_chunk(
    chunk: &mut [Molecule],
    options: &ReactionOptions,
    reaction_species_dir: Option<&Path>,
    reactive_frames_dir: Option<&Path>,
) -> Result<Vec<Reaction>> {
    let noise_event_life = options.noise_event_life;

    chunk.into_par_iter().for_each(|mut mol| {
        // ignore molecules already `rebond` in overlap region
        if mol.nbonds() == 0 {
            mol.rebond();
        }
    });

    get_chemical_reactions(chunk, noise_event_life, reaction_species_dir, reactive_frames_dir)
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

    Ok(())
}
// dd2f60bb ends here
