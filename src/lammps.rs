// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use std::fs::File;
use std::error::Error;
use std::io::{self, BufReader, BufWriter};
use std::io::prelude::*;
use std::collections::HashMap;
use std::path::Path;

use petgraph::prelude::*;
use petgraph as pg;

use crate::atoms::{AtomData, TrajectoryFrame, write_as_cif};
use crate::Frame;
use crate::graph::fragments_from_atoms;
// imports:1 ends here

// extract frame structure
// 提出特定frame对应的结构, 生成cif文件.
// - 元素信息: symbols
// - 键连信息: bonds
// - 结构信息: dump, cell


// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*extract%20frame%20structure][extract frame structure:1]]
pub fn extract_frame(filename: &str, target_timestep: usize, ciffile: &str) -> Result<(), Box<Error>>
{
    // 1. guess required lammps files from input filename
    let path = Path::new(filename);
    let path_lammps_data = path.with_extension("data");
    let path_lammps_dump = path.with_extension("dump");
    let path_lammps_bonds = path.with_extension("bonds");
    let path_lammps_bonds_terse = path.with_extension("bonds-terse");

    if ! path_lammps_data.is_file() {
        let msg = format!("data file not found: {:}", path_lammps_data.display());
        Err(msg)?;
    }
    if ! path_lammps_bonds_terse.is_file() {
        eprintln!("bonds-terse file not found, creating now ...");
        create_terse_copy_of_lammps_bonds_file(&path_lammps_bonds, &path_lammps_bonds_terse)?;
    }

    if ! path_lammps_dump.is_file() {
        let msg = format!("dump file not found: {:}", path_lammps_dump.display());
        Err(msg)?;
    }

    // get positions from dump file
    let mut frame = get_frame_from_lammps_dump_file(&path_lammps_dump, target_timestep)?;

    // get symbols from data file
    frame.symbols = parse_lammps_data_file(&path_lammps_data)?;

    // assign connectivity
    frame.neighbors = get_connectivity_from_terse_bonds_file(&path_lammps_bonds_terse, target_timestep)?;

    let path = Path::new(ciffile);
    write_as_cif(frame, &path);

    Ok(())
}
// extract frame structure:1 ends here

// fragment analysis

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*fragment%20analysis][fragment analysis:1]]
pub fn analyze_frames(filename: &str, outfile: &str, maxcols: usize) -> Result<(), Box<Error>>{
    let frames = parse_lammps_files(filename)?;
    write_formated_text(&frames, outfile, maxcols)?;

    Ok(())
}

fn parse_lammps_files(filename: &str) -> Result<Vec<Frame>, Box<Error>> {
    // 1. guess required lammps files from input filename
    let path = Path::new(filename);
    let path_lammps_data = path.with_extension("data");
    let path_lammps_dump = path.with_extension("dump");
    let path_lammps_bonds = path.with_extension("bonds");
    let path_lammps_bonds_terse = path.with_extension("bonds-terse");

    if ! path_lammps_data.is_file() {
        let msg = format!("data file not found: {:}", path_lammps_data.display());
        Err(msg)?;
    }
    if ! path_lammps_bonds_terse.is_file() {
        eprintln!("bonds-terse file not found, creating now ...");
        create_terse_copy_of_lammps_bonds_file(&path_lammps_bonds, &path_lammps_bonds_terse)?;
    }

    // read atom indices and symbols
    let symbols = parse_lammps_data_file(&path_lammps_data)?;

    // assign connectivity
    let frames = parse_terse_bonds_file(&path_lammps_bonds_terse, &symbols);

    frames
}

fn write_formated_text(frames: &Vec<Frame>, outfile: &str, max_columns: usize) -> Result<(), Box<Error>>{
    // create output file
    let f = File::create(outfile)?;
    let mut writer = BufWriter::new(f);

    let mut species:HashMap<String, usize> = HashMap::new();
    for frame in frames {
        for (k, v) in &frame.fragments {
            let x = species.entry(k.to_string()).or_insert(0_usize);
            *x += v;
        }
    }

    let mut count_vec: Vec<_> = species.iter().collect();
    count_vec.sort_by_key(|k| k.1);
    count_vec.reverse();

    let vs:Vec<String> = count_vec.iter().map(|x| x.0.to_string()).collect();
    writer.write("Timestep ".as_bytes());

    let mut mc = vs.len();
    if max_columns < mc {
        mc = max_columns;
    }
    let vs = &vs[..mc];
    writer.write(format!("{:}\n", vs.join(" ")).as_bytes());

    for frame in frames {
        let s = format!("{:^width$}", frame.timestep, width="Timestep ".len());
        writer.write(s.as_bytes());
        let mut lst = Vec::new();
        for k in vs.iter() {
            let count = frame.fragments.get(k).unwrap_or(&0_usize);
            lst.push(format!("{:^width$}", count, width=k.len()));
        }
        let s = format!("{}\n", lst.join(" "));
        writer.write(s.as_bytes());
    }

    Ok(())
}
// fragment analysis:1 ends here

// src
// - 关键信息: 所有原子对应的元素类型.
// - 数据类型选择HashMap, key为index, value为元素符号


// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*src][src:1]]
/// read data from lammps .data file
fn parse_lammps_data_file(path: &Path) -> Result<HashMap<usize, String>, Box<Error>>
{
    eprintln!("reading data file: {}", path.display());

    let fp = File::open(path)?;
    let mut reader = BufReader::new(fp);
    let mut lines_iter = reader.lines().peekable();

    // sanity check
    if let Some(&Ok(ref firstline)) = lines_iter.peek() {
        if ! firstline.starts_with("LAMMPS data file") {
            let msg = format!("read in a wrong file: {}", firstline);
            Err(msg)?;
        }
    } else {
        let msg = format!("Expect more lines: {}", path.display());
        Err(msg)?;
    }

    // skip the first two lines
    for _ in 0..2 {
        lines_iter.next();
    }

    // 1. read number of total atoms
    // 684  atoms
    let mut natoms = 0;
    if let Some(line) = lines_iter.next() {
        let line = line?;
        assert!(line.contains(" atoms"), format!("cannot find number of atoms: {}", line));
        let mut attrs = line.split_whitespace();
        if let Some(s) = attrs.nth(0) {
            natoms = s.parse().unwrap();
        } else {
            let msg = format!("failed to get natoms: {}", line);
            Err(msg)?;
        }
    } else {
        Err("data file is incomplete: failed to get natoms!")?;
    }

    // 2. read in number of atom types
    let mut natom_types = 0_usize;
    loop {
        if let Some(line) = lines_iter.next() {
            let line = line?;
            if line.ends_with("atom types") {
                if let Some(n) = line.split_whitespace().nth(0) {
                    natom_types = n.parse().unwrap();
                }
                break;
            }
        } else {
            Err("cannot find atom types lines in lammps data file")?;
        }
    }

    // 3. parse atom types
    // NOTE: element symbol is supposed to be after `#`
    //     1  50.941500   # V
    assert!(natom_types > 0);
    let mut mapping_symbols = HashMap::new();
    loop {
        if let Some(line) = lines_iter.next() {
            let line = line?;
            if line.starts_with("Masses") {
                // skip one blank line
                lines_iter.next();
                // mapping: atom_index ==> atom_symbol
                for _ in 0..natom_types {
                    if let Some(line) = lines_iter.next() {
                        let line = line?;
                        let mut attrs = line.split_whitespace();
                        let k = attrs.nth(0).unwrap();
                        let v = attrs.last().unwrap();
                        mapping_symbols.insert(k.to_string(), v.to_string());
                    }
                }
                break;
            }
        } else {
            Err("failed to read Masses section")?;
        }
    }

    // 4. read in atom index and atom type
    assert!(natoms > 0);

    let mut symbols = HashMap::new();
    loop {
        if let Some(line) = lines_iter.next() {
            let line = line?;
            if line.starts_with("Atom") {
                // skip one blank line
                lines_iter.next();
                for _ in 0..natoms {
                    if let Some(line) = lines_iter.next() {
                        let line = line?;
                        let mut attrs = line.split_whitespace();
                        let index = attrs.next().unwrap();
                        let t = attrs.next().unwrap();
                        let index = index.parse().unwrap();
                        let symbol = mapping_symbols.get(t).unwrap().to_string();
                        symbols.insert(index, symbol);
                    } else {
                        Err("Atom records are incomplete.")?;
                    }
                }
                break;
            }
        } else {
            Err("cannot find Atom lines in lammps data file")?;
        }
    }

    Ok(symbols)
}
// src:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*test][test:1]]
#[test]
#[ignore]
fn test_parse_data_file() {
    let filename = "/home/ybyygu/Incoming/FeC reaxff tests/FeCO/terse-tests/test.data";
    let path = Path::new(&filename);
    let symbols = parse_lammps_data_file(&path);
    println!("{:?}", symbols);
}
// test:1 ends here

// src
// #+name: daedfe6b-34ed-4dd1-94a2-4e698a00a42c

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::daedfe6b-34ed-4dd1-94a2-4e698a00a42c][daedfe6b-34ed-4dd1-94a2-4e698a00a42c]]
fn get_frame_from_lammps_dump_file
    (path: &Path, target_timestep: usize) -> Result<TrajectoryFrame, Box<Error>>
{
    eprintln!("reading dump file: {}", path.display());

    let fp = File::open(path)?;
    let mut reader = BufReader::new(fp);

    let mut frame = TrajectoryFrame::new();
    let mut natoms = 0_usize;
    let mut timestep = 0_usize;
    let mut buf = String::new();
    loop {
        // 0. sanity check
        buf.clear();
        let nb = reader.read_line(&mut buf)?;
        if nb <= 0 {
            eprintln!("reached the end of the file: {}", path.display());
            break;
        }
        assert!(buf.starts_with("ITEM: TIMESTEP"), format!("Expect the frame header, but: {}", buf));

        // 1. get current timestep
        buf.clear();
        let nb = reader.read_line(&mut buf)?;
        if nb <= 0 {
            Err("Expect more lines: failed to read timestep!")?;
        }
        timestep = buf.trim_right().parse()?;

        // 2. get number of atoms
        for i in 0..2 {
            buf.clear();
            let nb = reader.read_line(&mut buf)?;
            if nb > 0 {
                if i == 1 {
                    natoms = buf.trim_right().parse()?;
                }
            } else {
                Err("Expect more lines: failed to read number of atoms!")?;
            }
        }

        // 3. get lammps box and atoms
        assert!(natoms > 0, buf);
        println!("current timestep = {:?}", timestep);
        if timestep < target_timestep {
            for _ in 0..(natoms + 5) {
                let nb = reader.read_line(&mut buf)?;
                if nb <= 0 {
                    Err("Expect more lines: failed to read lammps box and atoms")?;
                }
            }
        } else if timestep == target_timestep {
            frame.timestep = timestep;
            frame.natoms = natoms;

            buf.clear();
            // 3.1 the lammps box
            for _ in 0..4 {
                let nb = reader.read_line(&mut buf)?;
                if nb <= 0 {
                    Err("Expect more lines: failed to read lammps box!")?;
                }
            }
            let (cell, origin) = get_lammps_dump_box(&buf)?;
            frame.cell = cell;
            frame.cell_origin = origin;

            // 3.2 the atom records
            buf.clear();
            for _ in 0..(natoms+1) {
                let nb = reader.read_line(&mut buf)?;
                if nb <= 0 {
                    Err("Expect more lines: failed to read all atom records!")?;
                }
            }
            let positions = get_lammps_dump_positions(&buf, natoms)?;
            frame.positions = positions;
            // ignore remaining lines
            break;
        } else {
            let msg = format!("Requested timestep {} not found in {}", target_timestep, path.display());
            Err(msg)?;
        }
    }

    Ok(frame)
}
// daedfe6b-34ed-4dd1-94a2-4e698a00a42c ends here

// tests
// #+name: 69313c0d-d969-40b4-a338-ff274407d54d

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::69313c0d-d969-40b4-a338-ff274407d54d][69313c0d-d969-40b4-a338-ff274407d54d]]
#[test]
#[ignore]
fn test_parse_lammps_dump_file() {
    let filename = "/home/ybyygu/Workspace/Programming/reaction-analysis/tests/FeCO/Fe100_8816_50CO_500_500K.dump";
    let path = Path::new(&filename);

    let frame = get_frame_from_lammps_dump_file(&path, 200_usize).unwrap();
    println!("{:?}", frame.positions);
}
// 69313c0d-d969-40b4-a338-ff274407d54d ends here

// src
// #+name: 32374d0e-1d81-4ec0-a8b3-0fb7950a625a

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::32374d0e-1d81-4ec0-a8b3-0fb7950a625a][32374d0e-1d81-4ec0-a8b3-0fb7950a625a]]
use std::f64;

enum BoxStyle {
    Orthogonal,
    Triclinic,
}

fn get_lammps_dump_box(txt: &str) -> Result<([[f64; 3]; 3], [f64; 3]), Box<Error>>{
    use self::BoxStyle::{Orthogonal, Triclinic};

    let mut lines_iter = txt.lines();

    let mut hi    = [0.0; 3];
    let mut lo    = [0.0; 3];
    let mut tilt  = [0.0; 3];
    let mut style = Orthogonal;

    if let Some(line) = lines_iter.next() {
        if line.starts_with("ITEM: BOX BOUNDS") {
            let attrs = line.split_whitespace();
            style = match attrs.count() {
                6 => Orthogonal,
                9 => Triclinic,
                _ => Err(format!("unexpected box style: {}", &line))?,
            };

            for i in 0..3 {
                if let Some(line) = lines_iter.next() {
                    let mut attrs:Vec<f64> = line.split_whitespace().map(|x| x.parse::<f64>().unwrap()).collect();
                    match style {
                        Orthogonal => {
                            lo[i] = attrs[0];
                            hi[i] = attrs[1];
                        },
                        Triclinic => {
                            lo[i] = attrs[0];
                            hi[i] = attrs[1];
                            tilt[i] = attrs[2];
                        },
                    }
                } else {
                    Err("lammps box is incomplete!")?;
                }
            }
        } else {
            let msg = format!("expect LAMMPS BOX header, but found: {}", &line);
            Err(msg)?;
        }
    } else {
        Err("why")?;
    }

    let mut va = [0.0; 3];
    let mut vb = [0.0; 3];
    let mut vc = [0.0; 3];
    let mut origin = lo;
    match style {
        Orthogonal => {
            va[0] = hi[0] - lo[0];
            vb[1] = hi[1] - lo[1];
            vc[2] = hi[2] - lo[2];
            origin = lo;
        },
        Triclinic  => {
            let xy = tilt[0];
            let xz = tilt[1];
            let yz = tilt[2];

            // x vector
            let xlo = lo[0] - [0.0, xy, xz, xy+xz].iter().fold(f64::MAX, |a, &b| a.min(b));
            let xhi = hi[0] - [0.0, xy, xz, xy+xz].iter().fold(f64::MIN, |a, &b| a.max(b));
            va[0] = xhi - xlo;
            // y vector
            let ylo = lo[1] - [0.0, yz].iter().fold(f64::MAX, |a, &b| a.min(b));
            let yhi = hi[1] - [0.0, yz].iter().fold(f64::MIN, |a, &b| a.max(b));
            vb[0] = xy;
            vb[1] = yhi - ylo;
            // z vector
            let zlo = lo[2];
            let zhi = hi[2];
            vc[0] = xz;
            vc[1] = yz;
            vc[2] = zhi - zlo;
            origin = [xlo, ylo, zlo];
        },
    }

    Ok(([va, vb, vc], origin))
}
// 32374d0e-1d81-4ec0-a8b3-0fb7950a625a ends here

// tests
// #+name: 398e563b-0ad5-4845-a5c3-97c115748e74

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::398e563b-0ad5-4845-a5c3-97c115748e74][398e563b-0ad5-4845-a5c3-97c115748e74]]
#[test]
fn test_lammps_box() {
    let box1 = "ITEM: BOX BOUNDS pp pp pp
-0.195983 11.329
-0.195983 11.329
-0.195983 11.329";

    // results from ovito
    let cell_vector1 = [11.525, 0.0, 0.0];
    let cell_vector2 = [0.0, 11.525, 0.0];
    let cell_vector3 = [0.0, 0.0, 11.525];
    let cell_origin  = [-0.195983, -0.195983, -0.195983];

    let (vts, origin) = get_lammps_dump_box(&box1).unwrap();

    assert_relative_eq!(vts[0][0], cell_vector1[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[1][1], cell_vector2[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[2][2], cell_vector3[2] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[0], cell_origin[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[1], cell_origin[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[2], cell_origin[2] as f64, epsilon=1.0e-4);

    let box2 = "ITEM: BOX BOUNDS xy xz yz pp pp pp
-0.08189 15.3282 -0.045807
0.072939 15.5755 0
0.001924 17.4877 0";

    // results from ovito
    let cell_vector1 = [15.3643, 0.0, 0.0];
    let cell_vector2 = [-0.045807, 15.5026, 0.0];
    let cell_vector3 = [0.0, 0.0, 17.4858];
    let cell_origin  = [-0.036083, 0.072939, 0.001924];

    let (vts, origin) = get_lammps_dump_box(&box2).unwrap();
    assert_relative_eq!(vts[0][0], cell_vector1[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[1][0], cell_vector2[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[1][1], cell_vector2[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[2][2], cell_vector3[2] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[0], cell_origin[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[1], cell_origin[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[2], cell_origin[2] as f64, epsilon=1.0e-4);

    let box3 = "ITEM: BOX BOUNDS pp pp ff
0.0000000000000000e+00 2.2931000000000001e+01
0.0000000000000000e+00 2.2931000000000001e+01
-1.0000000000000000e+00 5.0497999999999998e+01
";
    let cell_vector1 = [22.931, 0.0, 0.0];
    let cell_vector2 = [0.0, 22.931, 0.0];
    let cell_vector3 = [0.0, 0.0, 51.498];
    let cell_origin = [0.0, 0.0, -1.0];
    let (vts, origin) = get_lammps_dump_box(&box3).unwrap();
    assert_relative_eq!(vts[0][0], cell_vector1[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[1][1], cell_vector2[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(vts[2][2], cell_vector3[2] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[0], cell_origin[0] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[1], cell_origin[1] as f64, epsilon=1.0e-4);
    assert_relative_eq!(origin[2], cell_origin[2] as f64, epsilon=1.0e-4);
}
// 398e563b-0ad5-4845-a5c3-97c115748e74 ends here

// src
// #+name: dd0f9789-eb7f-492a-aa0d-24a76d346f76

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::dd0f9789-eb7f-492a-aa0d-24a76d346f76][dd0f9789-eb7f-492a-aa0d-24a76d346f76]]
fn get_lammps_dump_positions(txt: &str, natoms: usize) -> Result<HashMap<usize, [f64; 3]>, Box<Error>>{
    let mut lines_iter = txt.lines();
    let prefix = "ITEM: ATOMS";

    // 1. sanity check and get header labels
    let mut labels = Vec::new();
    if let Some(line) = lines_iter.next() {
        if line.starts_with(&prefix) {
            labels = line[prefix.len()..].split_whitespace().collect();
        } else {
            Err(format!("Expected header not found: {}", &line))?;
        }
    } else {
        Err("failed to read atoms header.")?;
    }

    // 2. assign values according to header labels
    let mut positions: HashMap<usize, [f64; 3]> = HashMap::new();
    for _ in 0..natoms {
        if let Some(line) = lines_iter.next() {
            let mut attrs:Vec<_> = line.split_whitespace().collect();
            assert!(attrs.len() == labels.len(), line.to_string());
            let mut index = 0_usize;
            let mut x = 0_f64;
            let mut y = 0_f64;
            let mut z = 0_f64;
            for (k, v) in labels.iter().zip(attrs.iter()) {
                match k {
                    &"x"  => x = v.parse().unwrap(),
                    &"y"  => y = v.parse().unwrap(),
                    &"z"  => z = v.parse().unwrap(),
                    &"id" => index = v.parse().unwrap(),
                    _     => (),
                }
            }
            positions.insert(index, [x, y, z]);
        } else {
            Err("atom records are incomplete")?;
        }
    }

    Ok(positions)
}
// dd0f9789-eb7f-492a-aa0d-24a76d346f76 ends here

// tests

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*tests][tests:1]]
#[test]
fn test_parse_dump_positions() {
    let txt = "ITEM: ATOMS id type x y z
1 1 3.77622 3.9054 0.009267
2 1 3.77622 3.9054 2.90503
3 1 1.89072 1.66683 4.37145
4 1 4.97252 1.29462 4.37145
5 2 4.73984 2.09641 0.131493";

    let natoms = 5_usize;
    let positions = get_lammps_dump_positions(&txt, natoms).unwrap();
    assert_relative_eq!(3.77622, positions[&1_usize][0], epsilon=1e-4);
    assert_relative_eq!(0.131493, positions[&5_usize][2], epsilon=1e-4);

    let txt = "ITEM: ATOMS x y z type id
0.1832399964 1.684999943 3.850500107 1 1
4.906760216 5.054999828 0.6794999838 1 2
2.361759901 5.054999828 1.585500002 1 3
2.728240013 1.684999943 2.944499969 1 4
0.9467399716 0.4246200025 1.485839963 1 5 ";

    let positions = get_lammps_dump_positions(&txt, 5_usize).unwrap();
    assert_relative_eq!(0.183239996, positions[&1_usize][0], epsilon=1e-4);
    assert_relative_eq!(4.906760216, positions[&2_usize][0], epsilon=1e-4);
}
// tests:1 ends here

// on file basis
// #+name: 8497df51-3eb3-41fc-a87d-c1688c94c29f

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::8497df51-3eb3-41fc-a87d-c1688c94c29f][8497df51-3eb3-41fc-a87d-c1688c94c29f]]
fn get_connectivity_from_terse_bonds_file
    (path: &Path, target_timestep: usize) -> Result<HashMap<usize, Vec<usize>>, Box<Error>>
{
    eprintln!("reading bonds file: {}", path.display());

    let fp = File::open(path)?;
    let mut reader = BufReader::new(fp);

    // parse bonds
    let mut neighbors: HashMap<usize, Vec<usize>> = HashMap::new();

    let mut timestep = 0;
    let mut natoms = 0;
    let mut buf = String::new();
    loop {
        // 0. sanity check and get current timestep
        buf.clear();
        let nb = reader.read_line(&mut buf)?;
        if nb <= 0 {
            eprintln!("reached the end of the file: {}", path.display());
            break;
        }

        let label = "# Timestep";
        timestep = get_int_data_from_comment_line(&buf, &label)?;

        // 1. get natoms
        buf.clear();
        let nb = reader.read_line(&mut buf)?;
        if nb <= 0 {
            Err("Expect more lines: failed to read number of atoms!")?;
        }
        let label = "# Number of particles";
        natoms = get_int_data_from_comment_line(&buf, &label)?;

        // 2. read atom records
        buf.clear();
        if timestep < target_timestep {
            for _ in 0..natoms {
                let nb = reader.read_line(&mut buf)?;
                if nb <= 0 {
                    Err("Expect more lines: failed to read all atom records!")?;
                }
            }
        } else if timestep == target_timestep {
            for _ in 0..natoms {
                let nb = reader.read_line(&mut buf)?;
                if nb <= 0 {
                    Err("Expect more lines: failed to read all atom records!")?;
                }
            }
            neighbors = get_connectivity_from_terse_bonds_file_frame(&buf, natoms)?;
            // ignore other parts
            break;
        } else {
            let msg = format!("Requested timestep {} not found in {}", target_timestep, path.display());
            Err(msg)?;
        }
    }

    Ok(neighbors)
}
// 8497df51-3eb3-41fc-a87d-c1688c94c29f ends here



// #+name: 1cb4fcf1-d093-41ce-a011-f88a95c9bf7b

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::1cb4fcf1-d093-41ce-a011-f88a95c9bf7b][1cb4fcf1-d093-41ce-a011-f88a95c9bf7b]]
fn parse_terse_bonds_file (path: &Path, symbols: &HashMap<usize, String>)
                           -> Result<Vec<Frame>, Box<Error>>
{
    // create file handler and a buffle reader
    let fp = File::open(path)?;
    let mut reader = BufReader::new(fp);
    let mut lines_iter = reader.lines().peekable();

    // parse data
    let mut timestep = 0;
    let mut natoms = 0;
    let mut frames = Vec::new();
    loop {
        if lines_iter.peek().is_none() {
            eprintln!("reached the end of the file: {}", path.display());
            break;
        }

        // process a single frame
        let frame = parse_terse_bonds_file_single_frame(&mut lines_iter, &symbols)?;
        eprintln!("timestep {:}, done.", frame.timestep);
        eprintln!("fragments {:?}", frame.fragments);
        // // for test
        // if frame.timestep > 5000_usize {
        //     break;
        // }
        frames.push(frame);
    }

    Ok(frames)
}
// 1cb4fcf1-d093-41ce-a011-f88a95c9bf7b ends here

// on frame basis
// 读入不同原子对应的电荷及connectivity.
// 简化版中:
// - 以行数来定原子index.
// - 第一列为浮点数, 对应partial charge
// - 其后几列为与当前原子直接成键的所有原子对应的编号数据. 该数据为键连原子编号与当
//   前原子编号之差, 即index - current

// #+name: dd3a4020-2ed5-4c62-a6f7-1d80e0fc6198

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::dd3a4020-2ed5-4c62-a6f7-1d80e0fc6198][dd3a4020-2ed5-4c62-a6f7-1d80e0fc6198]]
fn get_connectivity_from_terse_bonds_file_frame(
    txt: &str,
    natoms: usize) -> Result<HashMap<usize, Vec<usize>>, Box<Error>>
{
    let mut neighbors = HashMap::new();

    let mut lines_iter = txt.lines();

    for n in 1..natoms+1 {
        if let Some(line) = lines_iter.next() {
            let (charge, nns) = parse_terse_bonds_file_single_line(&line);
            let mut connected = vec![];
            for x in nns {
                connected.push(x+n);
            }
            neighbors.insert(n, connected);
        } else {
            Err("Atom data is incomplete.")?;
        }
    }

    Ok(neighbors)
}
// dd3a4020-2ed5-4c62-a6f7-1d80e0fc6198 ends here



// #+name: b88e850d-3754-49fe-a0a6-cdf0ba8e2169

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::b88e850d-3754-49fe-a0a6-cdf0ba8e2169][b88e850d-3754-49fe-a0a6-cdf0ba8e2169]]
#[test]
fn test_get_connectivity_from_terse_bonds_frame() {
    let txt = "0.007 8 9 16 896 1017 1016 905 904 120 121 1 112 128
0.008 128 112 120 121 1 1016 904 896 16 8 1017 905 9
0.009 112 120 128 896 1016 1017 904 905 121 1 8 16 9
0.008 120 121 112 1017 1016 905 904 896 8 9 128 16 1
0.008 896 1017 1016 1 112 904 120 8 905 16 128 121 9";

    let neighbors = get_connectivity_from_terse_bonds_file_frame(&txt, 5_usize).unwrap();
    assert_eq!(neighbors.len(), 5);

    let connected = neighbors.get(&1).unwrap();
    assert_eq!(connected.len(), 13);
    assert_eq!(connected[0], 9);
    assert_eq!(connected[1], 10);
}
// b88e850d-3754-49fe-a0a6-cdf0ba8e2169 ends here



// #+name: 5f005858-636b-4d77-aca5-e6be1baca10a

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::5f005858-636b-4d77-aca5-e6be1baca10a][5f005858-636b-4d77-aca5-e6be1baca10a]]
fn parse_terse_bonds_file_single_frame<I> (
    lines_iter: &mut I,
    symbols: &HashMap<usize, String>) -> Result<Frame, Box<Error>>
    where I: Iterator<Item=io::Result<String>>,
{
    // 1. read current timestep
    let mut timestep = 0;
    let label = "# Timestep";
    if let Some(line) = lines_iter.next() {
        let line = line?;
        timestep = get_int_data_from_comment_line(&line, &label)?;
    } else {
        Err("Failed to read timestep!")?;
    }

    // 2. read number of atoms
    let mut natoms = 0;
        let label = "# Number of particles";
    if let Some(line) = lines_iter.next() {
        let line = line?;
        natoms = get_int_data_from_comment_line(&line, &label)?;
    } else {
        Err("Failed to read number of atoms!")?;
    }

    // 3. read connectivity for each atom
    let mut atoms = Vec::new();
    assert!(natoms > 0);
    for n in 1..natoms+1 {
        let mut data = AtomData::new();
        if let Some(line) = lines_iter.next() {
            let line = line?;
            let (charge, neighbors) = parse_terse_bonds_file_single_line(&line);
            data.index = n;
            data.charge = charge;
            data.symbol = symbols.get(&data.index).unwrap().to_string();
            for x in neighbors {
                data.neighbors.push(x+n);
            }
            atoms.push(data);
        } else {
            Err("Atom data is incomplete.")?;
        }
    }

    assert!(atoms.len() == natoms);

    // 4. create frame
    let mut frame = Frame::new();
    frame.timestep = timestep;
    let fragments = fragments_from_atoms(&atoms);
    frame.fragments = fragments;

    Ok(frame)
}
// 5f005858-636b-4d77-aca5-e6be1baca10a ends here

// on line basis
// #+name: 39637608-12b2-4724-ac38-cfc5d4f9c990

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::39637608-12b2-4724-ac38-cfc5d4f9c990][39637608-12b2-4724-ac38-cfc5d4f9c990]]
fn parse_terse_bonds_file_single_line(line: &str) -> (f64, Vec<usize>) {
    let mut attrs = line.split_whitespace();
    let first = attrs.nth(0).unwrap();
    let charge:f64 = first.parse().unwrap();
    let neighbors:Vec<usize> = attrs.map(|x| x.parse::<usize>().unwrap()).collect();

    (charge, neighbors)
}
// 39637608-12b2-4724-ac38-cfc5d4f9c990 ends here



// #+name: b0d25b51-fbb3-460a-882f-3cdb7c6f6619

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::b0d25b51-fbb3-460a-882f-3cdb7c6f6619][b0d25b51-fbb3-460a-882f-3cdb7c6f6619]]
#[test]
fn test_parse_terse_bonds_line() {
    let s = "0.007 8 9 16 896 1017 1016 905 904 120 121 1 112 128";
    let (charge, neighbors) = parse_terse_bonds_file_single_line(&s);
    assert_eq!(0.007, charge);
    assert_eq!(13, neighbors.len());
    assert_eq!(8, neighbors[0]);
}
// b0d25b51-fbb3-460a-882f-3cdb7c6f6619 ends here

// on file basis
// 处理文件级别的逻辑
// #+name: 7ea4d968-ef79-438d-8fed-ce109e391554

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::7ea4d968-ef79-438d-8fed-ce109e391554][7ea4d968-ef79-438d-8fed-ce109e391554]]
/// Parameters
/// ----------
/// inputfile: the path of a .bonds file
/// outputfile: the path of the output file
fn create_terse_copy_of_lammps_bonds_file
    (inputfile: &Path, outputfile: &Path) -> Result<String, Box<Error>>
{
    // 1.0 prepare input file stream
    let mut fp = File::open(inputfile)?;
    let reader = BufReader::new(fp);
    let mut lines_iter = reader.lines().peekable();

    // 1.1 prepare output file stream
    if outputfile.exists() {
        panic!("output file already exists: {}", outputfile.display());
    }

    let mut fp = File::create(outputfile)?;
    let mut writer = BufWriter::new(fp);

    // 2 create a terse copy for each frame
    // loop over frames in trajectory file
    loop {
        if lines_iter.peek().is_none() {
            println!("reached the end of the file.");
            break;
        }

        let output = parse_lammps_bonds_single_snapshot(&mut lines_iter).unwrap();
        writer.write_all(&output.as_bytes());
        if let Some(line) = lines_iter.next() {
            ;
        } else {
            println!("Warning: missing final blank comment.");
            break;
        }
    }

    Ok("done".to_string())
}
// 7ea4d968-ef79-438d-8fed-ce109e391554 ends here

// on frame basis
// 处理轨迹中单一帧
// #+name: 2bf6281c-b33e-4181-b7fb-8a3eb6ad88b2

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::2bf6281c-b33e-4181-b7fb-8a3eb6ad88b2][2bf6281c-b33e-4181-b7fb-8a3eb6ad88b2]]
use std::iter::Peekable;

fn parse_lammps_bonds_single_snapshot<I>(lines_iter: &mut I) -> Result<String, String>
    where I: Iterator<Item=io::Result<String>>,
{
    let mut timestep = 0 as usize;
    let mut natoms = 0 as usize;

    let mut output:Vec<String> = Vec::new();

    // 1. read in meta data from comments
    // expected => Some(Ok("# Timestep 0 "))
    for n in 0..7 {
        let line = lines_iter.next().unwrap().unwrap();
        assert!(line.starts_with("# "), line);
        match n {
            0 => {
                timestep = get_int_data_from_comment_line(&line, "# Timestep").unwrap();
                output.push(line);
            },
            2 => {
                natoms = get_int_data_from_comment_line(&line, "# Number of particles").unwrap();
                output.push(line);
            },
            _ => (),
        }
    }

    println!("processing timestep {:} ...", timestep);
    let mut data = Vec::new();
    // 2. read atom data
    for _ in 0..natoms {
        if let Some(line) = lines_iter.next() {
            let line = line.unwrap();
            let x = get_terse_line_from_bonds_data_line(&line).unwrap();
            data.push(x);
        } else {
            let msg = format!("atom data is incomplete\ntimestep = {}", timestep);
            return Err(msg);
        }
    }

    // 3. append to output
    // atom indices are sorted in ascending order
    data.sort();
    for (_, s) in data {
        output.push(s);
    }

    Ok(output.join("\n") + "\n")
}

// #[test]
// fn test_snapshot_parse() {
//     let filename = "/home/ybyygu/Incoming/FeC reaxff tests/FeCO/Fe100_8816_50CO_500_500K.bonds";

//     let f = File::open(filename).unwrap();
//     let mut reader = BufReader::new(f);
//     let mut lines_iter = reader.lines().peekable();

//     let output = parse_lammps_bonds_single_snapshot(&mut lines_iter).unwrap();
//     let path = Path::new("tt");
//     let mut file = File::create(&path).unwrap();
//     file.write_all(output.as_bytes());
// }
// 2bf6281c-b33e-4181-b7fb-8a3eb6ad88b2 ends here

// 有陨版: 保留连接表, 不保留bond order

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*%E6%9C%89%E9%99%A8%E7%89%88:%20%E4%BF%9D%E7%95%99%E8%BF%9E%E6%8E%A5%E8%A1%A8,%20%E4%B8%8D%E4%BF%9D%E7%95%99bond%20order][有陨版: 保留连接表, 不保留bond order:1]]
fn get_terse_line_from_bonds_data_line(line: &str) -> Result<(usize, String), String>{
    if line.starts_with("# ") {
        let msg = format!("incorrect line: {}", line);
        return Err(msg);
    }

    let mut attrs:Vec<&str> = line.split_whitespace().collect();
    let cur = attrs[0].parse::<usize>().unwrap();
    let nbonds = attrs[2].parse::<usize>().unwrap();
    let neighbors = &attrs[3..3+nbonds];

    // partial charge
    let charge = &attrs.last().unwrap();

    // using a string to store the result
    let mut result = format!("{}", charge);

    // 1. keep the neighbor which is larger than current
    // 2. store the shift relative to current
    for n in neighbors {
        let n:usize = n.parse().unwrap();
        if n > cur {
            result.push_str(format!(" {}", n - cur).as_str());
        }
    }

    Ok((cur, result))
}

#[test]
fn test_terse_line() {
    // test1
    let s = " 1108 3 1 1107 0 1.479 1.479 2.000 -0.313";
    let new = get_terse_line_from_bonds_data_line(&s);
    assert!(new.is_ok());
    let (cur, result) = new.unwrap();
    assert!(cur == 1108);
    assert!(result.starts_with("-0.313"));

    // test2
    let s = " 137 1 9 9 249 257 138 145 153 265 273 129 0 0.450 0.450 0.710 0.410 0.709 0.450 0.450 0.709 0.709 5.047 -1.000 0.025";
    let new = get_terse_line_from_bonds_data_line(&s);
    assert!(new.is_ok());
    let (cur, result) = new.unwrap();
    assert!(cur == 137);
    assert_eq!(result, "0.025 112 120 1 8 16 128 136");
}
// 有陨版: 保留连接表, 不保留bond order:1 ends here

// functions parsing line

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*functions%20parsing%20line][functions parsing line:1]]
fn get_int_data_from_comment_line(line: &str, prefix: &str) -> Result<usize, String> {
    if line.starts_with(prefix) {
        let s = line[prefix.len()..].trim().parse::<usize>();
        match s {
            Ok(v) => return Ok(v),
            Err(why) => return Err(format!("{:?}", why)),
        }
    } else {
        let msg = format!(
            "Failed to get value {:?} for current frame: {:?}",
            prefix, line
        );
        Err(msg)
    }
}

#[test]
fn test_get_int_data_from_comment_line() {
    // get timestep
    let r = get_int_data_from_comment_line("# Timestep 10", "# Timestep");
    assert_eq!(r, Ok(10));
    // get number of atoms
    let r = get_int_data_from_comment_line("# Number of particles 684", "# Number of particles");
    assert_eq!(r, Ok(684));

    let r = get_int_data_from_comment_line("# Timestep 0.0", "# Timestep");
    assert!(r.is_err());
    let r = get_int_data_from_comment_line("12 22\n", "# Timestep");
    assert!(r.is_err());
}

// fn get_atom_data_from_line(line: &str) -> Result<(AtomData, &[usize]), String> {
fn get_atom_data_from_line(line: &str) -> Result<AtomData, String> {
    let mut data = AtomData::new();

    let error = format!("Failed to parse atom data from: {}", line);

    // 1. get index
    let mut attrs = line.split_whitespace();
    if let Some(v) = attrs.next() {
        match v.parse::<usize>() {
            Ok(v) => {
                data.index = v;
            }
            Err(why) => {
                return Err(format!("{:?}", why));
            }
        }
    } else {
        return Err(error);
    }

    // 2. get particle type
    if let Some(v) = attrs.next() {
        data.symbol = v.to_string();
    } else {
        return Err("failed to read particle type.".to_string());
    }

    // 3. get number of neighbors
    let mut nneighbors = 0;
    if let Some(v) = attrs.next() {
        match v.parse::<usize>() {
            Ok(v) => {
                nneighbors = v;
            }
            Err(why) => {
                return Err(format!("{:?}", why));
            }
        }
    } else {
        return Err("failed to read number of neighbors.".to_string());
    }

    // 4. get neighbors
    // let mut neighbors = vec![];
    for _ in 0..nneighbors {
        if let Some(v) = attrs.next() {
            match v.parse::<usize>() {
                Ok(v) => {
                    // neighbors.push(v);
                    data.neighbors.push(v);
                }
                Err(why) => {
                    return Err(format!("{:?}", why));
                }
            }
        } else {
            return Err(error);
        }
    }

    Ok(data)
}

#[test]
fn test_get_atom_data_from_line() {
    let line =
        " 121 3 2 301 28 0         0.978         0.978         1.956         2.000        -0.736 ";
    let r = get_atom_data_from_line(&line);
    assert!(r.is_ok());
    // let (data, _) = r.unwrap();
    let data = r.unwrap();
    assert!(data.index == 121);
    assert!(data.symbol == "3");
}
// functions parsing line:1 ends here
