// [[file:../trajectory.note::*imports][imports:1]]
use gchemol::prelude::*;
// imports:1 ends here

// [[file:../trajectory.note::*test][test:1]]
use gut::prelude::*;

#[test]
#[ignore]
fn test_xyz_traj() -> Result<()> {
    let f = "data/ba/cd22e1-0409-4aae-b545-956390eba4b5/10h2o-al2o3-ni-400k.xyz";
    for mol in gchemol::io::read(f)? {
        let positions_o = mol
            .atoms()
            .filter_map(|(_, a)| if a.symbol() == "O" { Some(a.position()) } else { None })
            .collect_vec();
        dbg!(positions_o.len());
    }

    Ok(())
}
// test:1 ends here
