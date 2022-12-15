// [[file:../trajectory.note::2d0268ed][2d0268ed]]
use crate::common::*;

use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

use gchemol::prelude::*;
use gchemol::Molecule;
use indicatif::ProgressBar;

pub fn count_rings_in_trajectory<P: AsRef<Path>>(path: P, max: usize) -> Result<String> {
    let f = File::open(&path)?;
    let mut f = BufReader::new(f);
    let mut buf = String::new();

    let mut out = String::new();

    println!("Working on {:?} ...", path.as_ref().display());

    // set up progress bar
    let metadata = std::fs::metadata(path)?;
    let bar = ProgressBar::new(metadata.len());
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
        let mut mol = Molecule::from_str(&buf, "text/xyz")?;

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
        out.push_str(&format!("frame {}\n", iframe));
        for k in 3..(max + 1) {
            if let Some(n) = map.get(&k) {
                out.push_str(&format!("{}, {:}\n", k, n));
            } else {
                out.push_str(&format!("{}, {:}\n", k, 0));
            }
        }

        // update progress bar
        bar.inc(buf.len() as u64);

        // reset frame buf
        buf.clear();
    }
    bar.finish();

    Ok(out)
}
// 2d0268ed ends here
