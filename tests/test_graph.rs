// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::3ac4ebc7-6b99-4984-8571-f0cc7ab27990][3ac4ebc7-6b99-4984-8571-f0cc7ab27990]]
extern crate petgraph;

use std::collections::HashSet;

use petgraph::prelude::*;
use petgraph as pg;
use petgraph::graph::node_index as n;

#[test]
fn test_get_connected_components() {
    let mut graph = Graph::new_undirected();

    // weight
    let w = 0;

    let mut v = Vec::new();

    // Add 8 vertices to G
    for i in 0..8 {
        v.push(graph.add_node(i));
    }

    // create bonds for CH4
    for i in 1..5 {
        graph.add_edge(v[0], v[i], w);
    }

    graph.add_edge(v[5], v[6], w);
    graph.add_edge(v[5], v[7], w);

    let ncc = pg::algo::connected_components(&graph);
    assert_eq!(2, ncc);

    for x in pg::algo::kosaraju_scc(&graph) {
        println!("{:?}", x);
    }

    // store connected components
    let mut ccs1 = HashSet::new();
    ccs1.insert(v[0]);
    ccs1.insert(v[1]);
    ccs1.insert(v[2]);
    ccs1.insert(v[3]);
    ccs1.insert(v[4]);

    let mut ccs2 = HashSet::new();
    ccs2.insert(v[5]);
    ccs2.insert(v[6]);
    ccs2.insert(v[7]);
}
// 3ac4ebc7-6b99-4984-8571-f0cc7ab27990 ends here
