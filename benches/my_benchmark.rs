// [[file:../trajectory.note::2db75e1a][2db75e1a]]
use gut::prelude::*;
use trajectory_analysis::docs::xyztraj::*;
use trajectory_analysis::fibonacci;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn read_xyz_traj_raw() {
    let mols = gchemol::io::read("data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz")
        .unwrap()
        .collect_vec();
    assert_eq!(mols.len(), 8);
}

fn read_xyz_traj_opt() {
    let mols = read_xyz_trajectory("data/55/798fcc-6c7a-424f-8c87-7e8b11300345/SiAlCaO1800k.xyz".as_ref())
        .unwrap()
        .collect_vec();
    assert_eq!(mols.len(), 8);
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("read-xyz-traj", |b| b.iter(|| read_xyz_traj_raw));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
// 2db75e1a ends here
