// [[file:../trajectory.note::2db75e1a][2db75e1a]]
use gut::prelude::*;
use trajectory_analysis::docs::xyztraj::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

const input_file: &'static str = "data/55/798fcc-6c7a-424f-8c87-7e8b11300345/输入和输出轨迹/15Ca300k-pos-1.xyz";
// const input_file: &'static str = "data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz";

fn read_xyz_traj_raw() {
    let mols = gchemol::io::read(input_file).unwrap();
    let _ = mols.count();
}

fn read_xyz_traj_opt() {
    let mols = read_xyz_trajectory(input_file.as_ref()).unwrap();
    let _ = mols.count();
}


pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("read-xyz-traj-raw", |b| b.iter(|| read_xyz_traj_raw()));
    c.bench_function("read-xyz-traj-opt", |b| b.iter(|| read_xyz_traj_opt()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
// 2db75e1a ends here
