// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use std::path::Path;
use std::collections::HashMap;

use guts::prelude::*;
use text_parser::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*core][core:1]]
/// Minimal representation for LAMMPS trajectory frame.
#[derive(Debug, Default)]
pub struct LammpsTrajectoryFrame {
    /// ITEM: TIMESTEP
    pub timestep: usize,
    /// ITEM: ATOMS
    pub atoms: HashMap<usize, LammpsAtom>,
}

/// Minimal Atom representation for LAMMPS.
#[derive(Debug)]
pub struct LammpsAtom {
    // Element number
    pub n: usize,
    // Cartesian coordinates
    pub xyz: [f64; 3],
}

impl LammpsAtom {
    pub fn new(n: usize, xyz: [f64; 3]) -> Self {
        Self { n, xyz }
    }
}

struct LammpsTrajectoryParser {
    //
}

impl LammpsTrajectoryParser {
    /// Process trajectory file that could be very large in size.
    fn process(&self, trjfile: &Path) {
        unimplemented!()
    }
}

fn read_lammps_dump_file() {
    unimplemented!()
}
// core:1 ends here

// meta
// 读入Frame元数据.

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*meta][meta:1]]
#[derive(Debug)]
struct FrameData {
    timestep: usize,
    natoms: usize,
}

named!(read_meta_data<&str, FrameData>,
    do_parse!(
                  tag!("ITEM: TIMESTEP")        >> eol >>
        timestep: read_usize                    >> // current timestep in this frame
                  tag!("ITEM: NUMBER OF ATOMS") >> eol >>
        natoms  : read_usize                    >> // number of atoms
        (
            FrameData {
                timestep, natoms
            }
        )
    )
);

#[test]
fn test_read_meta_data() {
    let txt = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
537
ITEM: BOX BOUNDS pp pp pp
-200.487 200.487
-200.487 200.487
-200.487 200.487
";
    let (_, x) = read_meta_data(txt).unwrap();
    assert_eq!(x.timestep, 0);
    assert_eq!(x.natoms, 537);
}
// meta:1 ends here

// src

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*src][src:1]]
#[derive(Debug)]
struct BoxData {
    t: String,
    a: String,
    b: String,
    c: String,
}

named!(read_box_data<&str, BoxData>,
   do_parse!(
          tag!("ITEM: BOX BOUNDS")   >>
       t: read_until_eol >> //
       a: read_until_eol >> //
       b: read_until_eol >> //
       c: read_until_eol >> //
       (
           BoxData {
               t: t.into(), a: a.into(), b: b.into(), c: c.into()
           }
       )
   )
);

#[test]
fn test_read_box_data() {
    let txt = "ITEM: BOX BOUNDS pp pp pp
-200.487 200.487
-200.487 200.487
-200.487 200.487
";
    let (_, x) = read_box_data(txt).unwrap();
}
// src:1 ends here

// src

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*src][src:1]]
// ITEM: ATOMS id type x y z c_eng c_cn c_cnt c_cna
named!(read_atom_header<&str, &str>, do_parse!(
            tag!("ITEM: ATOMS")     >>
    header: read_until_eol          >> // the line including column headers
    (
        header
    )
));

#[test]
fn test_read_atom_header() {
    let txt = "ITEM: ATOMS id type x y z c_eng c_cn c_cnt c_cna
1 1 0.542919 3.02575 -9.81464 -3.53666 9 9.42606 5
";
    let (_, x) = read_atom_header(txt).unwrap();
    assert_eq!(x, " id type x y z c_eng c_cn c_cnt c_cna");
}

fn read_atoms(input: &str, natoms: usize) -> IResult<&str, HashMap<usize, LammpsAtom>> {
    let (rest, header_line) = read_atom_header(input)?;
    let (rest, atom_lines) = many_m_n(natoms, natoms, read_until_eol)(rest)?;

    // collect column headers
    let headers: Vec<_> = header_line.trim().split_whitespace().collect();
    let nheaders = headers.len();

    // parse atom properties
    let atoms: HashMap<_, _> = atom_lines
        .into_iter()
        .map(|line| {
            let items: Vec<_> = line.trim().split_whitespace().collect();
            assert_eq!(items.len(), nheaders);

            let mut x: Option<f64> = None;
            let mut y: Option<f64> = None;
            let mut z: Option<f64> = None;
            let mut n: Option<usize> = None;
            let mut id: Option<usize> = None;
            for (k, v) in headers.iter().zip(items.into_iter()) {
                match k {
                    &"x" | &"xu" => x = v.parse().ok(),
                    &"y" | &"yu" => y = v.parse().ok(),
                    &"z" | &"zu" => z = v.parse().ok(),
                    &"id" => id = v.parse().ok(),
                    &"type" => n = v.parse().ok(),
                    _ => (),
                }
            }
            let xyz = [
                x.expect("invalid x data"),
                y.expect("invalid y data"),
                z.expect("invalid z data"),
            ];
            let id = id.expect("invalid atom id data");
            let n = n.expect("invalid atom number data");
            let atom = LammpsAtom::new(n, xyz);
            (id, atom)
        })
        .collect();

    Ok((rest, atoms))
}
// src:1 ends here

// tests

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*tests][tests:1]]
#[test]
fn test_read_atoms() {
    let txt = "ITEM: ATOMS id type x y z c_eng c_cn c_cnt c_cna
1 1 0.542919 3.02575 -9.81464 -3.53666 9 9.42606 5
3 2 0.274667 2.58508 -6.92683 -4.05727 12 4.9903 5
2 1 -0.566175 4.44199 -8.2352 -3.78452 11 5.16628 5
4 1 -1.0285 4.21504 -5.29927 -4.12135 10 4.87744 5
5 3 -0.209709 2.20696 -4.20742 -3.99243 10 0.410664 3
6 1 -1.26547 3.85937 -2.59029 -4.01771 10 0.241619 5
7 1 -0.638752 1.99003 -1.2366 -4.03166 11 0.141509 3
8 1 -1.75827 3.6228 0.296721 -4.14419 13 0.28179 5
9 1 -0.927052 1.51729 1.59026 -3.98979 11 0.889157 5
10 1 -2.18531 3.375 2.91709 -4.08537 12 4.52081 5
11 1 -1.33645 1.21535 4.42074 -4.04849 12 0.307786 3
";

    let (_, m) = read_atoms(txt, 5).unwrap();
    assert_eq!(m.len(), 5);
    assert_eq!(m[&5].xyz[0], -0.209709);
    assert_eq!(m[&5].n, 3);

    let txt = "ITEM: ATOMS x y z type id
0.1832399964 1.684999943 3.850500107 1 1
4.906760216 5.054999828 0.6794999838 1 2
2.361759901 5.054999828 1.585500002 1 3
2.728240013 1.684999943 2.944499969 1 4
0.9467399716 0.4246200025 1.485839963 1 5 ";
    let (_, m) = read_atoms(txt, 3).unwrap();
    assert_eq!(m.len(), 3);
    assert_eq!(m[&3].xyz[0], 2.361759901);
    assert_eq!(m[&3].n, 1);
}
// tests:1 ends here

// frame

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*frame][frame:1]]
fn read_lammps_dump_frame(input: &str) -> IResult<&str, LammpsTrajectoryFrame> {
    let (rest, frame_data) = read_meta_data(input)?;
    let (rest, box_data) = read_box_data(rest)?;
    let (rest, atoms) = read_atoms(rest, frame_data.natoms)?;

    // assign frame data
    let frame = {
        let mut frame = LammpsTrajectoryFrame::default();
        frame.timestep = frame_data.timestep;
        frame.atoms = atoms;
        frame
    };

    Ok((rest, frame))
}

#[test]
fn test_parser() -> Result<()> {
    let fname = "tests/files/lammps-test.dump";
    let parser = TextParser::default();
    let frames: Vec<_> = parser.parse_file(fname, read_lammps_dump_frame).collect();
    assert_eq!(frames.len(), 3);

    Ok(())
}
// frame:1 ends here

// pub

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*pub][pub:1]]
/// Parse LAMMPS trajectory file (.dump), returning a list of frames.
pub fn parse_lammps_dump_file(trjfile: &Path) -> impl Iterator<Item = LammpsTrajectoryFrame> + '_ {
    let parser = TextParser::default();
    parser.parse_file(trjfile, read_lammps_dump_frame)
}
// pub:1 ends here
