// [[file:../trajectory.note::77290756][77290756]]
//! Find rings in a `Molecule`
//!
//! Credit:
//!
//! Heavily inspired by the codes developed by vitroid: https://github.com/vitroid/CountRings
// 77290756 ends here

// [[file:../trajectory.note::0221ebf7][0221ebf7]]
use crate::common::*;

use gchemol::Molecule;

use std::collections::HashSet;
// 0221ebf7 ends here

// [[file:../trajectory.note::ca262281][ca262281]]
fn find_ring(
    mol: &Molecule,    // parent molecule
    members: &[usize], // current node list
    max: usize,        // max ring size?
) -> (usize, Rings) {
    let n = members.len();
    if n > max {
        return (max, vec![]);
    }

    let mut results = vec![];
    let last = members[n - 1];
    let mut max = max;
    for adj in mol.connected(last) {
        if members.contains(&adj) {
            if adj == members[0] {
                // Ring is closed.
                // It is the best and unique answer.
                let s: HashSet<_> = members.to_vec().into_iter().collect();
                if !shortcuts(mol, members) {
                    return (members.len(), vec![s]);
                }
            } else {
                // Shortcut ring
            }
        } else {
            let mut ms = members.to_vec();
            ms.push(adj);
            let (newmax, newres) = find_ring(mol, &ms, max);
            if newmax < max {
                max = newmax;
                results = newres;
            } else if newmax == max {
                results.extend(newres);
            }
        }
    }

    (max, results)
}

fn shortcuts(mol: &Molecule, members: &[usize]) -> bool {
    let n = members.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (j - i).min(n - (j - i));
            if d > shortest_pathlen(mol, members[i], members[j]) {
                return true;
            }
        }
    }
    return false;
}

// FIXME: caching
fn shortest_pathlen(mol: &Molecule, i: usize, j: usize) -> usize {
    if let Some(n) = mol.nbonds_between(i, j) {
        n
    } else {
        0
    }
}

fn find_rings(mol: &Molecule, max_ring_size: usize) -> Rings {
    let mut rings = vec![];
    for x in mol.numbers() {
        let mut neis = mol.connected(x).collect_vec();
        neis.sort();
        for p in neis.iter().combinations(2) {
            let y = p[0];
            let z = p[1];
            let triplet = [*y, x, *z];
            let (max, mut results) = find_ring(mol, &triplet, max_ring_size);
            for i in results {
                if !rings.contains(&i) {
                    rings.push(i);
                }
            }
        }
    }

    rings
}
// ca262281 ends here

// [[file:../trajectory.note::92cea8ed][92cea8ed]]
pub type Rings = Vec<HashSet<usize>>;

pub trait FindRings {
    fn find_rings(&self, max: usize) -> Rings;
}

impl FindRings for Molecule {
    /// Find rings up to `max` atoms in a `Molecule`.
    ///
    /// # Arguments
    /// * nmax: max ring size to be searched.
    fn find_rings(&self, nmax: usize) -> Rings {
        let rings = find_rings(self, nmax);

        rings
    }
}
// 92cea8ed ends here

// [[file:../trajectory.note::c9d04499][c9d04499]]
#[test]
fn test_find_rings_xyz() -> Result<()> {
    use gchemol::prelude::*;

    let path = "data/8e/0f1aaf-d18d-4893-aa5b-4b3cd736ca99/a.mol2";
    let mol = Molecule::from_file(path)?;

    // output in NGPH for test against vitroid/CountRings
    // println!("@NGPH");
    // println!("{}", mol.natoms());
    // for (i, j, _) in mol.view_bonds() {
    //     println!("{} {}", i - 1, j - 1);
    // }
    // println!("-1 -1");

    let rings = mol.find_rings(6);

    // one 3-member-ring, one 4-member-ring, and one 6-member-ring
    let ns = rings.iter().map(|ri| ri.len()).sorted().collect_vec();
    assert_eq!(ns, [3, 4, 6]);

    Ok(())
}
// c9d04499 ends here
