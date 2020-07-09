// [[file:../trajectory.note::*imports][imports:1]]
use gchemol::prelude::*;
use structopt::*;
// imports:1 ends here

// [[file:../trajectory.note::*gaussian function][gaussian function:1]]
#[inline]
fn gaussian_factor(w: f64, ca: f64, ca_ij: f64) -> f64 {
    (-w * (ca - ca_ij).powi(2)).exp()
}
// gaussian function:1 ends here

// [[file:../trajectory.note::*variables][variables:1]]
use std::iter::FromIterator;

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Variables {
    /// The lower and upper limit of d variable.
    d_range: (f64, f64),
    /// The lower and upper limit of lambda variable.
    l_range: (f64, f64),
    /// Grid size for d and l specified in n_grid x n_grid format
    n_grid: usize,
    /// The W parameter in Gaussian function.
    w: f64,
    /// Gas phase water molecules defined in oxygen and the associated hydrogon
    /// atoms in H2O.
    h2o_gas_map: std::collections::HashMap<usize, (usize, usize)>,
    /// Surface oxygen atoms
    o_atoms_surface: Vec<usize>,
}

impl Default for Variables {
    fn default() -> Self {
        let waters = vec![
            (37, (91, 92)),
            (38, (93, 94)),
            (39, (95, 96)),
            (40, (97, 98)),
            (41, (99, 100)),
            (42, (101, 102)),
            (43, (103, 104)),
            (44, (105, 106)),
            (45, (107, 108)),
            (46, (109, 110)),
        ];

        let o_atoms_surface = vec![25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36];
        let h2o_gas_map = std::collections::HashMap::from_iter(waters);

        Self {
            d_range: (-5.0, 5.0),
            l_range: (-5.0, 5.0),
            n_grid: 10,
            w: 1.0,
            o_atoms_surface,
            h2o_gas_map,
        }
    }
}

fn write_default_vars() -> Result<()> {
    let var = Variables::default();
    let s = serde_json::to_string_pretty(&var)?;
    println!("{}", s);

    Ok(())
}
// variables:1 ends here

// [[file:../trajectory.note::*cli][cli:1]]
use gchemol::Molecule;
use gut::prelude::*;

use gut::cli::*;
use gut::config::*;

#[derive(Debug, StructOpt)]
struct Cli {
    #[structopt(flatten)]
    verbose: gut::cli::Verbosity,

    /// The trajectory file in xyz format.
    trjfile: std::path::PathBuf,
    /// Variables for calculation in json format.
    json_file: std::path::PathBuf,

    /// plot using plotly
    #[structopt(long, short)]
    plot: bool,
}

pub fn cli() -> Result<()> {
    if std::env::args().count() <= 1 {
        write_default_vars()?;

        eprintln!("default json variables outputed.");
        eprintln!("Input 'adhoc -h' for help.");

        return Ok(());
    }

    let args = Cli::from_args();
    args.verbose.setup_logger();

    let s = gut::fs::read_file(args.json_file)?;
    let vars: Variables = serde_json::from_str(&s)?;

    let mols = gchemol::io::read_all(&args.trjfile)?;
    let m = vars.h2o_gas_map.len() as f64;
    let n = mols.len() as f64;
    let (d_lower, d_upper) = vars.d_range;
    let (l_lower, l_upper) = vars.l_range;
    let w = vars.w;
    let mut zz = vec![];

    let ng = vars.n_grid;
    let xx: Vec<_> = (0..=ng)
        .map(|i| d_lower + (d_upper - d_lower) * (i as f64) / (ng as f64))
        .collect();
    let yy: Vec<_> = (0..=ng)
        .map(|i| l_lower + (l_upper - l_lower) * (i as f64) / (ng as f64))
        .collect();

    for i in 0..=ng {
        let mut zi = vec![];
        for j in 0..ng {
            // let d = d_lower + (d_upper - d_lower) * (i as f64) / (vars.n_grid as f64);
            // let l = l_lower + (l_upper - l_lower) * (j as f64) / (vars.n_grid as f64);
            let d = xx[i];
            let l = yy[j];

            // calculate free energy etc
            let mut sum = 0.0;
            for mol in mols.iter() {
                let f = calculate_proton_diffusion_distribution(mol, w, d, l, &vars);
                sum += f;
            }

            let a = w / (n * m * std::f64::consts::PI);
            let pdl = a * sum;
            println!("{:-10.4} {:-10.4} {:-10.4}", d, l, pdl);
            zi.push(pdl);
        }
        zz.push(zi);
    }

    plot_3d(zz, xx, yy);

    Ok(())
}

fn calculate_proton_diffusion_distribution(mol: &Molecule, w: f64, d: f64, lambda: f64, vars: &Variables) -> f64 {
    // 气相水中的氧原子及该气相水中的氧对应的 2个H原子
    let h2o_gas_map = &vars.h2o_gas_map;

    // 表面的氧原子
    let o_atoms_surface = &vars.o_atoms_surface;

    let mut inner_sum = 0.0;
    info!("{:^4} {:^4} {:^10} {:^10}", "o_g", "o_s", "d_ij", "lambda_ij");
    for (&o_g, &(h1, h2)) in h2o_gas_map.iter() {
        let a_o_g = mol.get_atom(o_g).unwrap();
        let a_h1 = mol.get_atom(h1).unwrap();
        let a_h2 = mol.get_atom(h2).unwrap();

        // 只选择距离气相O最近的那个表面O原子
        let mut d_ij = std::f64::MAX;
        let mut o_s_min = 0;
        for &o_s in o_atoms_surface.iter() {
            let a_o_s = mol.get_atom(o_s).unwrap();
            let d_ij_ = a_o_g.distance(a_o_s);
            if d_ij_ < d_ij {
                d_ij = d_ij_;
                o_s_min = o_s;
            }
        }

        // 找到距离表面O最近的H
        let (d_o_s_h, a_h) = {
            let a_o_s = mol.get_atom(o_s_min).unwrap();
            let d1 = a_o_s.distance(a_h1);
            let d2 = a_o_s.distance(a_h2);
            if d1 < d2 {
                (d1, a_h1)
            } else {
                (d2, a_h2)
            }
        };
        let d_o_g_h = a_o_g.distance(a_h);
        let lambda_ij = d_o_g_h - d_o_s_h;

        info!("{:4} {:4} {:-10.5} {:-10.5}", o_g, o_s_min, d_ij, lambda_ij);
        inner_sum += gaussian_factor(w, d, d_ij) * gaussian_factor(w, lambda, lambda_ij);
    }

    inner_sum
}
// cli:1 ends here

// [[file:../trajectory.note::*plot][plot:1]]
fn plot_3d(z: Vec<Vec<f64>>, x: Vec<f64>, y: Vec<f64>) {
    use plotly::*;
    use plotly::layout::Axis;

    let mut plot = Plot::new();

    // layout tuning
    // let xaxis = Axis::new().tick_format("-0.2f");
    // let yaxis = Axis::new().tick_format("-0.2f");
    // let layout = Layout::new().xaxis(xaxis).yaxis(yaxis);
    // plot.set_layout(layout);

    let trace1 = Surface::new(z).x(x).y(y).cauto(false);
    plot.add_trace(trace1);
    plot.show();
}
// plot:1 ends here

// [[file:../trajectory.note::*test][test:1]]
#[test]
fn test_xyz_traj() -> Result<()> {
    let f = "data/ba/cd22e1-0409-4aae-b545-956390eba4b5/10h2o-al2o3-ni-400k.xyz";

    for mol in gchemol::io::read(f)? {
        // calculate_proton_diffusion_distribution(&mol, 1.0, 5.0, -5.0);
    }

    Ok(())
}
// test:1 ends here
