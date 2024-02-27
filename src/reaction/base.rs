// [[file:../../trajectory.note::74e5e10a][74e5e10a]]
use crate::common::*;

use std::collections::BTreeSet as HashSet;
use std::collections::BTreeMap;
// 74e5e10a ends here

// [[file:../../trajectory.note::eb6afa69][eb6afa69]]
/// Bonding states (bonded or unbonded) for each pair of atoms in each frame
#[derive(Debug, Default, Clone)]
pub struct BondingStates {
    // the total number of frames with bonding states
    nframes: usize,
    // key: [u, v] as atom pair u--v
    // value: a map with key for frame index, and value for bonding state
    inner: BTreeMap<[usize; 2], BTreeMap<usize, bool>>,
}

// for bond: u-v == v-u
fn ordered(key: [usize; 2]) -> [usize; 2] {
    let mut key = key;
    key.sort();
    key
}

impl BondingStates {
    /// Insert bonding state of frame `frame_index` for bond `key` between atom `u` and `v`
    pub fn set_frame(&mut self, frame_index: usize, key: [usize; 2], state: bool) {
        let mut entry = self.inner.entry(ordered(key)).or_default();
        entry.insert(frame_index, state);

        // NOTE: record the largest frame index for the number of frames
        if frame_index + 1 > self.nframes {
            self.nframes = frame_index + 1;
        }
    }

    /// Get bonding state of frame `frame_index` for bond `key` between atom `u` and `v`
    pub fn get_frame(&self, frame_index: usize, key: [usize; 2]) -> bool {
        // if no such bonding pair, we should return false
        if let Some(states) = self.inner.get(&ordered(key)) {
            states.get(&frame_index).copied().unwrap_or_default()
        } else {
            false
        }
    }

    /// The number of bonding pairs.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// The number of frames for each bonding pair.
    pub fn nframes(&self) -> usize {
        self.nframes
    }
}
// eb6afa69 ends here

// [[file:../../trajectory.note::45ddc634][45ddc634]]
impl BondingStates {
    /// Return bonding states code for pair u--v
    pub fn bonding_states_code(&self, key: [usize; 2]) -> impl Iterator<Item = bool> + '_ {
        (0..self.nframes).map(move |iframe| self.get_frame(iframe, key))
    }

    /// Return human readable bonding events code for pair of atom `key` u--v
    pub fn bonding_events_code(&self, key: [usize; 2]) -> String {
        let frames: Vec<_> = (0..self.nframes).collect();

        frames
            .windows(2)
            .map(|pair| {
                let iframe1 = pair[0];
                let iframe2 = pair[1];
                let state1 = self.get_frame(iframe1, key);
                let state2 = self.get_frame(iframe2, key);
                match (state1, state2) {
                    (false, true) => "↑",
                    (true, false) => "↓",
                    (true, true) => "-",
                    (false, false) => "-",
                }
            })
            .collect()
    }
}
// 45ddc634 ends here

// [[file:../../trajectory.note::a4eba9be][a4eba9be]]
impl BondingStates {
    /// Remove bonding pairs no bond breaking or forming
    /// changes. Returns the number of removed bonding pairs.
    pub fn remove_inactive_bonding_pairs(&mut self) -> usize {
        let keys: Vec<_> = self.bonding_pairs().collect();

        // find bonding pairs that have no changes
        let inactive_keys1 = keys.iter().filter(|&key| self.bonding_states_code(*key).all(|bonded| bonded));
        let inactive_keys2 = keys
            .iter()
            .filter(|&key| self.bonding_states_code(*key).all(|bonded| !bonded));

        // remove these bonding pairs
        let mut n_removed = 0;
        for key in inactive_keys1.chain(inactive_keys2).collect_vec() {
            let _ = self.inner.remove(key);
            n_removed += 1;
        }
        n_removed
    }
}
// a4eba9be ends here

// [[file:../../trajectory.note::1d27bbdc][1d27bbdc]]
impl BondingStates {
    /// Remove noise events for atom pair in `key` with a short
    /// `noise_event_life`. Returns the positions of noise bonding
    /// events occured in the frames of bonding states for pair of
    /// atoms in `key`.
    ///
    /// # NOTE
    /// These states in begin and end sides with `noise_event_life` items are
    /// not affected by the removing.
    pub fn remove_noise_events(&mut self, key: [usize; 2], noise_event_life: usize) -> Vec<usize> {
        let mut states: Vec<_> = (0..self.nframes).map(|iframe| self.get_frame(iframe, key)).collect();
        let positions_inverted = fix_noise_states(&mut states, noise_event_life);
        for &iframe in positions_inverted.iter() {
            let state = !self.get_frame(iframe, key);
            self.set_frame(iframe, key, state);
        }
        positions_inverted
    }

    /// Find reactive bonds in valid frame region in context of
    /// `noise_event_life`. Panics if `noise_event_life` too large to
    /// fit the whole frame region.
    pub fn find_reactive_bonds(&self, noise_event_life: usize) -> Vec<[usize; 2]> {
        let frames: Vec<_> = (0..self.nframes).collect();

        let istart = noise_event_life;
        let iend = self.nframes - noise_event_life;
        assert!(istart < iend, "invalid istart..iend {istart}..{iend}");

        let mut reactive_bonds = HashSet::new();
        for &bond in self.inner.keys() {
            let states: Vec<_> = self.bonding_states_code(bond).collect();
            let changes = find_reactive_changes(&states[istart..iend]);
            if !changes.is_empty() {
                for [u, v] in changes {
                    // fix index shift issue
                    reactive_bonds.insert([u + istart, v + istart]);
                }
            }
        }
        reactive_bonds.into_iter().collect()
    }
}

/// Find start and end positions for noise event pattern (in pair) in `code`
fn find_noise_codes(code: &str, n_space: usize) -> Vec<(usize, usize)> {
    use regex::Regex;

    // NOTE: here we use B for Breaking, F for Forming. These chars
    // are 1 char long, and will not break regex search logical.
    let pattern1 = format!(r#"B(-{{0,{n_space}}}?)F"#);
    let pattern2 = format!(r#"F(-{{0,{n_space}}}?)B"#);
    let re1 = Regex::new(&pattern1).unwrap();
    let re2 = Regex::new(&pattern2).unwrap();
    let mut positions1 = vec![];
    let mut positions2 = vec![];

    for mat in re1.find_iter(code) {
        positions1.push((mat.start(), mat.end()));
    }

    for mat in re2.find_iter(code) {
        positions2.push((mat.start(), mat.end()));
    }

    // we select the one that matches more
    if positions1.len() < positions2.len() {
        positions2
    } else {
        positions1
    }
}

fn states_to_events(states: &[bool]) -> String {
    let nframes = states.len();
    let frames: Vec<_> = (0..nframes).collect();
    frames
        .windows(2)
        .map(|pair| {
            let iframe1 = pair[0];
            let iframe2 = pair[1];
            let state1 = states[iframe1];
            let state2 = states[iframe2];
            match (state1, state2) {
                // (false, true) => "↑",
                // (true, false) => "↓",
                (false, true) => "F",
                (true, false) => "B",
                (true, true) => "-",
                (false, false) => "-",
            }
        })
        .collect()
}

pub(self) fn fix_noise_states(states: &mut [bool], noise_event_life: usize) -> Vec<usize> {
    let ecode = states_to_events(states);
    let matched_positions = find_noise_codes(&ecode, noise_event_life);

    let mut positions_inverted = vec![];
    // the valid region of active frames considering `noise_event_life`
    let istart = noise_event_life;
    let iend = states.len() - noise_event_life;
    assert!(istart < iend, "invalid active frame region: {istart}, {iend}");
    for (start, end) in matched_positions {
        for i in start + 1..end {
            if i >= istart && i < iend {
                states[i] = !states[i];
                positions_inverted.push(i);
            }
        }
    }
    positions_inverted
}

/// Find positions that have chemical reactions
fn find_reactive_changes(states: &[bool]) -> Vec<[usize; 2]> {
    assert!(!states.is_empty());

    let ecode = states_to_events(states);
    ecode.match_indices(&['F', 'B']).map(|(p, _)| [p, p + 1]).collect()
}

#[test]
fn test_find_noise_codes() {
    let xx = find_noise_codes("-B---F-B-F--B--F---B", 5);
    dbg!(xx);

    let mut states: Vec<_> = "+++--+++---+++++"
        .chars()
        .map(|x| if x == '+' { true } else { false })
        .collect();
    let ecode = states_to_events(&states);
    assert_eq!(ecode, "--B-F--B--F----");

    let pairs = find_noise_codes(&ecode, 1);
    assert_eq!(pairs.len(), 1);
    let pairs = find_noise_codes(&ecode, 0);
    assert_eq!(pairs.len(), 0);

    let pairs = find_noise_codes(&ecode, 2);
    assert_eq!(pairs.len(), 2);
    let (u, v) = pairs[0];
    assert_eq!(&ecode[u..v], "B-F");
    let (u, v) = pairs[1];
    assert_eq!(&ecode[u..v], "B--F");

    fix_noise_states(&mut states, 1);
    let ecode = states_to_events(&states);
    assert_eq!(ecode, "-------B--F----");

    let frames = find_reactive_changes(&states);
    assert_eq!(frames.len(), 2);
    for [i, j] in frames {
        assert_ne!(states[i], states[j]);
    }
}
// 1d27bbdc ends here

// [[file:../../trajectory.note::3b9c65f4][3b9c65f4]]
use gchemol::Molecule;

impl BondingStates {
    /// Construct `BondingStates` from bonds of each `Molecule` in `mols`.
    pub fn from_molecules<'a>(mols: impl IntoIterator<Item = &'a Molecule>) -> Self {
        let mut states = Self::default();
        for (iframe, mol) in mols.into_iter().enumerate() {
            for (u, v, _) in mol.bonds() {
                states.set_frame(iframe, [u, v], true);
            }
        }
        states
    }

    /// Returns an iterator over bonding pairs.
    pub fn bonding_pairs(&self) -> impl Iterator<Item = [usize; 2]> + '_ {
        self.inner.keys().copied()
    }
}
// 3b9c65f4 ends here

// [[file:../../trajectory.note::ace163f1][ace163f1]]
#[test]
fn test_bonding_states() -> Result<()> {
    let f = "tests/files/lty.xyz";
    let mols = gchemol::io::read(f)?;

    let mut mols: Vec<_> = mols.take(100).collect();
    for mol in mols.iter_mut() {
        mol.rebond();
    }

    let states = BondingStates::from_molecules(&mols[..10]);
    // bonded atom 208-212
    assert_eq!(states.get_frame(2, [208, 212]), true);
    assert_eq!(states.get_frame(20, [208, 212]), false);
    // unbonded atom 208-209
    assert_eq!(states.get_frame(2, [208, 209]), false);
    assert_eq!(states.get_frame(20, [208, 209]), false);
    assert_eq!(states.get_frame(3, [208, 207]), false);
    assert_eq!(states.get_frame(4, [208, 207]), true);

    // test noise removing with more frames
    let mut states = BondingStates::from_molecules(&mols);
    // bonding events code
    let keys: Vec<_> = states.bonding_pairs().collect();
    for &[u, v] in &keys {
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }

    let noise_event_life = 5;
    let reactive_bonds = states.find_reactive_bonds(noise_event_life);
    assert_eq!(reactive_bonds.len(), 72);

    for &[u, v] in &keys {
        states.remove_noise_events([u, v], noise_event_life);
        println!("{u:03}-{v:03}: {}", states.bonding_events_code([u, v]));
    }
    let reactive_bonds = states.find_reactive_bonds(noise_event_life);
    assert_eq!(reactive_bonds.len(), 3);

    Ok(())
}
// ace163f1 ends here
