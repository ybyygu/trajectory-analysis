// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::8da0ca06-6b0c-4d1b-8eec-827c7459cf2b][8da0ca06-6b0c-4d1b-8eec-827c7459cf2b]]
use std::collections::HashMap;

use petgraph::prelude::*;
use petgraph as pg;

use atoms::get_reduced_formula;
use atoms::AtomData;
// 8da0ca06-6b0c-4d1b-8eec-827c7459cf2b ends here

// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::080d17e9-5876-4808-8e5a-a181129da4fd][080d17e9-5876-4808-8e5a-a181129da4fd]]
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
// 080d17e9-5876-4808-8e5a-a181129da4fd ends here
