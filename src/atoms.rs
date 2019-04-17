// globals

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*globals][globals:1]]
use std::collections::HashMap;
use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::io::Write;
// globals:1 ends here

// atom
// 定义与Atom相关的元信息

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*atom][atom:1]]
#[derive(Debug, Default, Clone)]
pub struct AtomData{
    pub index: usize,
    pub symbol: String,
    pub neighbors: Vec<usize>,
    pub charge: f64,
    pub position: [f64; 3],
}

impl AtomData {
    pub fn new() -> Self {
        let ns:Vec<usize> = Vec::new();

        AtomData {
            index: 0,
            charge: 0.0,
            symbol: "H".to_string(),
            neighbors: ns,
            position: [0.0; 3],
        }
    }
}
// atom:1 ends here

// trajectory frame
// #+name: ec36c72e-4d04-447e-bd3d-bd4c6c3c49bb

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::ec36c72e-4d04-447e-bd3d-bd4c6c3c49bb][ec36c72e-4d04-447e-bd3d-bd4c6c3c49bb]]
// data structure for a single frame
pub struct TrajectoryFrame {
    pub timestep    : usize,            // current timestep
    pub natoms      : usize,
    pub cell        : [[f64; 3]; 3],
    pub cell_origin : [f64; 3],
    pub symbols     : HashMap<usize, String>,
    pub positions   : HashMap<usize, [f64; 3]>,
    pub neighbors   : HashMap<usize, Vec<usize>>,
}
// ec36c72e-4d04-447e-bd3d-bd4c6c3c49bb ends here



// #+name: 8a0cf3ff-227f-493b-aa5e-20ac3cac8f5c

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::8a0cf3ff-227f-493b-aa5e-20ac3cac8f5c][8a0cf3ff-227f-493b-aa5e-20ac3cac8f5c]]
impl TrajectoryFrame {
    pub fn new() -> Self {
        let positions: HashMap<usize, [f64; 3]> = HashMap::new();
        let symbols  : HashMap<usize, String> = HashMap::new();
        let cell = [[0_f64; 3]; 3];
        let neighbors: HashMap<usize, Vec<usize>> = HashMap::new();

        TrajectoryFrame {
            timestep    : 0,
            natoms      : 0,
            positions   : positions,
            symbols     : symbols,
            cell        : cell,
            cell_origin : [0_f64; 3],
            neighbors   : neighbors,
        }
    }
}
// 8a0cf3ff-227f-493b-aa5e-20ac3cac8f5c ends here

// src
// #+name: 6746df7d-0871-43bf-98cd-2e00e15020a5

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::6746df7d-0871-43bf-98cd-2e00e15020a5][6746df7d-0871-43bf-98cd-2e00e15020a5]]
use cgmath::{Vector3, Matrix3, Point3, Deg};
use cgmath::prelude::*;

fn cart_to_frac(matrix: Matrix3<f64>,
                coordinates: Vec<Vector3<f64>>) -> Vec<Vector3<f64>>
{
    let mut fractional = Vec::new();
    let inv = matrix.transpose().invert().unwrap();
    for v in coordinates {
        fractional.push(inv*v);
    }

    fractional
}
// 6746df7d-0871-43bf-98cd-2e00e15020a5 ends here



// 最近邻镜像原子
// #+name: c57f4ca0-4e68-4e30-a91b-7cbd47b7071c

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::c57f4ca0-4e68-4e30-a91b-7cbd47b7071c][c57f4ca0-4e68-4e30-a91b-7cbd47b7071c]]
use std::f64;

fn get_nearest_image(
    cell: Matrix3<f64>,
    position1: Point3<f64>,
    position2: Point3<f64>) -> (Vector3<f64>, f64)
{
    let d = position1.distance(position2);

    // loop 27 possible point images
    let relevant_images = [-1, 0, 1];
    let mut distance = f64::MAX;
    let mut image = Vector3::from_value(0_f64);
    for x in relevant_images.iter() {
        for y in relevant_images.iter() {
            for z in relevant_images.iter() {
                let p = position2 + (*x as f64)*cell.x + (*y as f64)*cell.y + (*z as f64)*cell.z;
                let d = position1.distance(p);
                if d < distance {
                    distance = d;
                    image.x = *x as f64;
                    image.y = *y as f64;
                    image.z = *z as f64;
                }
            }
        }
    }

    (image, distance)
}

#[test]
fn test_get_nearest_image() {
    let mat1 = Matrix3::new(5.09, 0.00, 0.00,
                            0.00, 6.74, 0.00,
                            0.00, 0.00, 4.53);

    let p1  = Point3::new(0.18324000,   1.68500000,   3.85050000);
    let p13 = Point3::new(4.53010000,   1.68500000,   2.03850000);
    let p10 = Point3::new(0.94674000,   2.94538000,   1.48584000);
    let dp1_13 = 1.95847;
    let dp1_10 = 2.61920;

    let (image, d) = get_nearest_image(mat1, p1, p13);
    assert_relative_eq!(d, dp1_13, epsilon=1e-4);
    assert_relative_eq!(image.x, -1.0, epsilon=1e-4);
    assert_relative_eq!(image.y, 0.0, epsilon=1e-4);
    assert_relative_eq!(image.z, 0.0, epsilon=1e-4);

    let (image, d) = get_nearest_image(mat1, p1, p10);
    assert_relative_eq!(d, dp1_10, epsilon=1e-4);
}
// c57f4ca0-4e68-4e30-a91b-7cbd47b7071c ends here



// #+name: c90b9d01-096f-47c1-bbfd-2649862e61dc

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::c90b9d01-096f-47c1-bbfd-2649862e61dc][c90b9d01-096f-47c1-bbfd-2649862e61dc]]
fn cell_vectors_to_parameters(matrix: Matrix3<f64>) -> (f64, f64, f64, f64, f64, f64) {
    let a = matrix.x.magnitude();
    let b = matrix.y.magnitude();
    let c = matrix.z.magnitude();

    let alpha: Deg<_> = matrix.y.angle(matrix.z).into();
    let beta: Deg<_> = matrix.x.angle(matrix.z).into();
    let gamma: Deg<_> = matrix.x.angle(matrix.y).into();

    (a, b, c, alpha.0, beta.0, gamma.0)
}
// c90b9d01-096f-47c1-bbfd-2649862e61dc ends here



// #+name: 9cab3b07-9781-48cd-a6fb-e6ee248b93dd

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::9cab3b07-9781-48cd-a6fb-e6ee248b93dd][9cab3b07-9781-48cd-a6fb-e6ee248b93dd]]
#[test]
fn test_cell() {
    // ovito/tests/files/LAMMPS/multi_sequence_1.dump
    let mat1 = Matrix3::new(5.09, 0.00, 0.00,
                            0.00, 6.74, 0.00,
                            0.00, 0.00, 4.53);
    let inv = mat1.transpose().invert().unwrap();

    let v1 = Vector3::new(2.1832, 1.6850, 3.8505);
    let v2 = Vector3::new(6.9068, 5.0550, 0.6795);
    let v3 = Vector3::new(4.3618, 5.0550, 1.5855);

    let fracs = cart_to_frac(mat1, vec![v1, v2, v3]);
    assert_relative_eq!(fracs[0].x, 0.4289, epsilon=1e-3);
    assert_relative_eq!(fracs[0].y, 0.2500, epsilon=1e-3);
    assert_relative_eq!(fracs[0].z, 0.8500, epsilon=1e-3);
    assert_relative_eq!(fracs[1].x, 1.3569, epsilon=1e-3);
    assert_relative_eq!(fracs[2].z, 0.3500, epsilon=1e-3);

    let mat2 = Matrix3::new(15.3643, 0.0, 0.0,
                            4.5807, 15.5026, 0.0,
                            0.0, 0.0, 17.4858);

    let (a, b, c, alpha, beta, gamma) = cell_vectors_to_parameters(mat2);
    assert_relative_eq!(a, 15.3643, epsilon=1e-4);
    assert_relative_eq!(b, 16.1652, epsilon=1e-4);
    assert_relative_eq!(c, 17.4858, epsilon=1e-4);

    assert_relative_eq!(alpha, 90.0, epsilon=1e-4);
    assert_relative_eq!(beta, 90.0, epsilon=1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon=1e-4);
}
// 9cab3b07-9781-48cd-a6fb-e6ee248b93dd ends here

// src
// #+name: fb3fb6d8-9ae3-4b29-a7c6-92c61776c867

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::fb3fb6d8-9ae3-4b29-a7c6-92c61776c867][fb3fb6d8-9ae3-4b29-a7c6-92c61776c867]]
pub fn write_as_cif(frame: TrajectoryFrame, path: &Path) -> Result<(), Box<Error>>{
    let mut lines = String::new();

    // meta inforation
    lines.push_str("data_test\n");
    lines.push_str("_audit_creation_method            'gosh'\n");
    lines.push_str("_symmetry_space_group_name_H-M    'P1'\n");
    lines.push_str("_symmetry_Int_Tables_number       1\n");
    lines.push_str("_symmetry_cell_setting            triclinic\n");
    lines.push_str("\n");

    // cell parameters
    lines.push_str("loop_\n");
    lines.push_str("_symmetry_equiv_pos_as_xyz\n");
    lines.push_str(" x,y,z\n");
    let cell = Matrix3::new(frame.cell[0][0],
                            frame.cell[0][1],
                            frame.cell[0][2],
                            frame.cell[1][0],
                            frame.cell[1][1],
                            frame.cell[1][2],
                            frame.cell[2][0],
                            frame.cell[2][1],
                            frame.cell[2][2]);
    let (a, b, c, alpha, beta, gamma) = cell_vectors_to_parameters(cell);

    lines.push_str(&format!("_cell_length_a     {:10.4}\n", a));
    lines.push_str(&format!("_cell_length_b     {:10.4}\n", b));
    lines.push_str(&format!("_cell_length_c     {:10.4}\n", c));
    lines.push_str(&format!("_cell_angle_alpha  {:10.4}\n", alpha));
    lines.push_str(&format!("_cell_angle_beta   {:10.4}\n", beta));
    lines.push_str(&format!("_cell_angle_gamma  {:10.4}\n", gamma));
    lines.push_str("\n");

    // atom fractional coordinates
    lines.push_str("loop_\n");
    lines.push_str("_atom_site_type_symbol\n");
    lines.push_str("_atom_site_label\n");
    lines.push_str("_atom_site_fract_x\n");
    lines.push_str("_atom_site_fract_y\n");
    lines.push_str("_atom_site_fract_z\n");

    let cell_origin = Vector3::new(frame.cell_origin[0], frame.cell_origin[1], frame.cell_origin[2]);
    for index in 1..(frame.natoms+1) {
        let position = frame.positions.get(&index).unwrap();
        let symbol = frame.symbols.get(&index).unwrap();
        let name = format!("{}{}", symbol, index);
        let coords = Vector3::new(position[0], position[1], position[2]) - cell_origin;
        let v = cell.transpose().invert().unwrap()*coords;
        let s = format!("{:4}{:6}{:12.5}{:12.5}{:12.5}\n", symbol, name, v.x, v.y, v.z);
        lines.push_str(&s);
    }

    // bonds
    if frame.neighbors.len() > 0 {
        lines.push_str("loop_\n");
        lines.push_str("_geom_bond_atom_site_label_1\n");
        lines.push_str("_geom_bond_atom_site_label_2\n");
        lines.push_str("_geom_bond_distance\n");
        lines.push_str("_geom_bond_site_symmetry_2\n");
        lines.push_str("_ccdc_geom_bond_type\n");
        for current in 1..(frame.natoms+1) {
            let symbol1 = frame.symbols.get(&current).unwrap();
            let name1 = format!("{}{}", symbol1, current);
            let p1 = frame.positions.get(&current).unwrap();
            let p1 = Point3::new(p1[0], p1[1], p1[2]) - cell_origin;

            let connected = frame.neighbors.get(&current).unwrap();
            for other in connected {
                if *other > current {
                    let symbol2 = frame.symbols.get(&other).unwrap();
                    let name2 = format!("{}{}", symbol2, other);
                    let p2 = frame.positions.get(&other).unwrap();
                    let p2 = Point3::new(p2[0], p2[1], p2[2]) - cell_origin;
                    let (image, distance) = get_nearest_image(cell, p1, p2);
                    if image.x == 0. && image.y == 0. && image.z == 0. {
                        lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name1, name2, distance, "."));
                    } else {
                        let symcode = get_image_symcode(image);
                        lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name1, name2, distance, symcode));
                        let (image, distance) = get_nearest_image(cell, p2, p1);
                        let symcode = get_image_symcode(image);
                        lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name2, name1, distance, symcode));
                    }
                }
            }
        }
    }

    // save as a file
    let f = File::create(path)?;
    let mut writer = BufWriter::new(&f);
    writer.write_all(&lines.as_bytes())?;

    Ok(())
}
// fb3fb6d8-9ae3-4b29-a7c6-92c61776c867 ends here



// #+name: bacd2858-ded3-4982-a986-dbee08aa6a51

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::bacd2858-ded3-4982-a986-dbee08aa6a51][bacd2858-ded3-4982-a986-dbee08aa6a51]]
fn get_image_symcode(image: Vector3<f64>) -> String {
    let mut symcode = String::new();

    symcode.push_str("1_");
    symcode.push_str(&format!("{:1}", image.x + 5_f64));
    symcode.push_str(&format!("{:1}", image.y + 5_f64));
    symcode.push_str(&format!("{:1}", image.z + 5_f64));

    symcode
}


#[test]
fn test_get_image_symcode() {
    let v = Vector3::new(1., 0., 0.);
    let x = get_image_symcode(v);
    assert_eq!(x, "1_655");
    let v = Vector3::new(0., 0., 0.);
    let x = get_image_symcode(v);
    assert_eq!(x, "1_555");
}
// bacd2858-ded3-4982-a986-dbee08aa6a51 ends here

// tests

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*tests][tests:1]]
#[test]
fn test_wirte_cif() {
    let mut frame = TrajectoryFrame::new();
    let path = Path::new("/tmp/test.cif");

    frame.timestep = 0;

    let cell = [[15.3643, 0.0, 0.0], [4.5807, 15.5026, 0.0], [0.0, 0.0, 17.4858]];
    frame.cell = cell;

    let mut positions = HashMap::new();
    positions.insert(1, [2.1832, 1.6850, 3.8505]);
    positions.insert(2, [6.9068, 5.0550, 0.6795]);
    frame.positions = positions;

    let mut symbols = HashMap::new();
    symbols.insert(1, "C".to_string());
    symbols.insert(2, "H".to_string());
    frame.symbols = symbols;

    write_as_cif(frame, &path);
}
// tests:1 ends here

// 定义最简单的原子结构信息.

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*%E5%AE%9A%E4%B9%89%E6%9C%80%E7%AE%80%E5%8D%95%E7%9A%84%E5%8E%9F%E5%AD%90%E7%BB%93%E6%9E%84%E4%BF%A1%E6%81%AF.][定义最简单的原子结构信息.:1]]
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

#[derive (Default, Debug, Clone, Copy)]
/// simple atom data structure
pub struct Atom {
    pub index: u64,
    pub symbol: &'static str,
}

impl PartialEq for Atom {
    fn eq(&self, other: &Atom) -> bool {
        self.index == other.index
    }
}

impl Eq for Atom {}

impl Hash for Atom {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}

impl Ord for Atom {
    fn cmp(&self, other: &Atom) -> Ordering {
        self.index.cmp(&other.index)
    }
}

impl PartialOrd for Atom {
    fn partial_cmp(&self, other: &Atom) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[test]
fn test_atom() {
    let a = Atom{
        index: 1,
        symbol: "H",
    };

    let b = Atom {
        index: 2,
        symbol: "H",

    };
    let mut c = Atom {
        index: 1,
        symbol: "H",
    };

    assert!(a != b);
    assert!(a == c);

    assert!(a.index == 1);
    assert!(a.symbol == "H");

    c.symbol = "C";
    assert!(c.symbol == "C");
}
// 定义最简单的原子结构信息.:1 ends here

// from symbols to formula
// 实现从元素符号列表到分子式.

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*from%20symbols%20to%20formula][from symbols to formula:1]]
// it is better to use generics function,
// but it is really difficult for me now
pub fn get_reduced_formula(symbols: &[&str]) -> String {
    // 1. count symbols: CCCC ==> C 4
    let mut counts = HashMap::new();
    for x in symbols {
        let c = counts.entry(x).or_insert(0);
        *c += 1;
    }

    let mut syms: Vec<String> = Vec::new();
    let mut to_append = String::new();
    // 2. format the formula
    for (k, v) in counts {
        // 2.1 omit number if the count is 1: C1H4 ==> CH4
        let mut s = String::new();
        if v > 1 {
            s = v.to_string();
        }
        // 2.2 special treatments for C and H
        let reduced = format!("{}{}", k, s);
        if *k == "C" {
            syms.insert(0, reduced);
        } else if *k == "H" {
            to_append = reduced;
        } else {
            syms.push(reduced);
        }
    }
    // 3. final output
    syms.push(to_append);
    let formula = syms.join("");
    formula
}

#[test]
fn test_formula() {
    let symbols   = vec!["C", "H", "C", "H", "H", "H"];
    let formula = get_reduced_formula(&symbols);
    assert!(formula == "C2H4".to_string());
    let symbols   = vec!["C", "H", "C", "H", "H", "O", "H", "O"];
    let formula = get_reduced_formula(&symbols);
    println!("{:?}", formula);
    assert!(formula == "C2O2H4".to_string());
}
// from symbols to formula:1 ends here
