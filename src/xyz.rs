// header

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*header][header:1]]
/// Credit:
///
/// Heavily inspired by the codes developed by vitroid: https://github.com/vitroid/CountRings
// header:1 ends here

// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use std::path::Path;

use gchemol::prelude::*;
use gchemol::Molecule;

use crate::common::*;
// imports:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*base][base:1]]
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

pub fn count_rings_in_trajectory<P: AsRef<Path>>(path: P, max: usize) -> Result<()> {
    let f = File::open(path)?;
    let mut f = BufReader::new(f);
    let mut buf = String::new();

    for iframe in 0.. {
        let _ = f.read_line(&mut buf)?;
        if buf.is_empty() {
            println!("done.");
            break;
        }

        let natoms: usize = buf.trim().parse().map_err(|e| {
            eprintln!("The first line should be an integer: {}!", buf);
            e
        })?;

        // read the remaining natoms+1 lines
        for _ in 0..natoms + 1 {
            let _ = f.read_line(&mut buf)?;
        }

        // construct molecule from text stream
        let mut mol = Molecule::parse_from(&buf, "text/xyz")?;

        // build bonding connectivity
        mol.rebond();
        let rings = mol.find_rings(max);

        // for count the numbers rings in same size
        let mut map: HashMap<usize, usize> = HashMap::new();

        for ri in rings {
            let n = ri.len();
            *map.entry(n).or_insert(0) += 1;
        }
        let mut keys: Vec<_> = map.keys().collect();
        keys.sort();
        if ! keys.is_empty() {
            println!("frame {:?}", iframe);
            for k in keys {
                println!("{}, {:}", k, map[k]);
            }
        }

        // reset frame buf
        buf.clear();
    }
    Ok(())
}
// base:1 ends here

// entry/find rings

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*entry/find%20rings][entry/find rings:1]]
use std::collections::HashSet;

pub type Rings = Vec<HashSet<AtomIndex>>;

pub trait FindRings {
    fn find_rings(&self, max: usize) -> Rings;
}

impl FindRings for Molecule {
    /// # Arguments
    ///
    /// * max: max ring size.
    ///
    fn find_rings(&self, max: usize) -> Rings {
        let rings = find_rings(self, max);

        rings
    }
}
// entry/find rings:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*core][core:1]]
use gchemol::AtomIndex;
use itertools::Itertools;

pub(crate) fn find_rings(mol: &Molecule, max_ring_size: usize) -> Rings {
    let mut rings = vec![];
    for x in mol.sites() {
        let mut neis = mol.neighbors(x);
        neis.sort();
        for p in neis.iter().combinations(2) {
            let y = p[0];
            let z = p[1];
            let triplet = [*y, x, *z];
            let (max, mut results) = find_ring(mol, &triplet, max_ring_size);
            for i in results {
                // let j: HashSet<_> = results.into_iter().collect();
                if !rings.contains(&i) {
                    rings.push(i);
                }
            }
        }
    }

    rings
}

fn find_ring(
    mol: &Molecule,        // parent molecule
    members: &[AtomIndex], // current node list
    max: usize,            // max ring size?
) -> (usize, Rings) {
    let n = members.len();
    if n > max {
        return (max, vec![]);
    }

    let mut results = vec![];
    // let s: HashSet<_> = members.to_vec().into_iter().collect();
    let last = members[n - 1];
    let mut max = max;
    for adj in mol.neighbors(last) {
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

fn shortcuts(mol: &Molecule, members: &[AtomIndex]) -> bool {
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
fn shortest_pathlen(mol: &Molecule, i: AtomIndex, j: AtomIndex) -> usize {
    if let Some(n) = mol.nbonds_between(i, j) {
        n
    } else {
        0
    }
}
// core:1 ends here

// tests

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*tests][tests:1]]
#[test]
fn test_read_xyz() -> Result<()> {
    // let path = "./test-files/nanoreactor/90ps~101ps.xyz";
    // let mol = read_frame(path)?;

    let path = "data/8e/0f1aaf-d18d-4893-aa5b-4b3cd736ca99/a.mol2";
    let path = "/home/ybyygu/a.mol2";
    let mols = gchemol::io::read(path)?;
    let mol = &mols[0];

    // output in NGPH for test against vitroid/CountRings
    // println!("@NGPH");
    // println!("{}", mol.natoms());
    // for (i, j, _) in mol.view_bonds() {
    //     println!("{} {}", i - 1, j - 1);
    // }
    // println!("-1 -1");

    let rings = mol.find_rings(6);
    for ri in rings {
        println!("{:}, {:?}", ri.len(), ri);
    }

    Ok(())
}
// tests:1 ends here
