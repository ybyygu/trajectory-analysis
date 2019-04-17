// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
use std::collections::HashMap;

use petgraph::prelude::*;
use petgraph as pg;

use crate::atoms::get_reduced_formula;
use crate::atoms::AtomData;
// imports:1 ends here

// fragments

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*fragments][fragments:1]]
pub fn fragments_from_atoms(atoms: &Vec<AtomData>) -> HashMap<String, usize>
{
    let mut graph = Graph::new_undirected();

    // add nodes
    for d in atoms {
        graph.add_node(&d.symbol);
    }

    // update edges
    for d in atoms {
        let icurrent = NodeIndex::new(&d.index-1);
        for x in &d.neighbors {
            let iother = NodeIndex::new(x - 1);
            graph.update_edge(icurrent, iother, 1);
        }
    }

    // get fragments from connected components
    let sccs = pg::algo::kosaraju_scc(&graph);
    let mut counts = HashMap::new();
    let mut symbols = vec![];
    for x in sccs {
        symbols.clear();
        for index in x {
            let sym = &graph[index];
            symbols.push(sym.as_str());
        }
        // count formulas
        let formula = get_reduced_formula(&symbols);
        let c = counts.entry(formula).or_insert(0);
        *c += 1;
    }

    counts
}
// fragments:1 ends here
