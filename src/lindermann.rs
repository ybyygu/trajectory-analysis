// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use stats::mean;
// imports:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*test][test:1]]
#[test]
fn test_mean() {
    let data = vec![1.0, 3.0, 9.2];
    let m = mean(data.into_iter());
    dbg!(m);
}
// test:1 ends here
