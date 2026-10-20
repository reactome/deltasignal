"""
Recycle edges (specs/081): a modified protein's return to the complex it starts in.

Reactome curates p130 (RBL2) leaving RBL2:E2F4,E2F5:TFDP1,TFDP2 when cyclin D:CDK4
phosphorylates it, and PP2A making free p130 again, but not free p130 rebinding E2F4.
So the complex is a root, the generator builds it from its own starting copy of p130,
and the free p130 PP2A makes leads nowhere. The cycle is never closed.

The generator (`LNG_RECYCLE_POOLS`) writes the cycle as a conserved pool: the resting
complex copies and the modified form are its states, the kinase reaction the forward
step, and the dephosphorylation the return step. It also writes a `recycle` edge from
the free copy to the starting copy. That edge is structural only:

- it joins the cycle's nodes into one cyclic component, where pools are solved
  (`index_reactions`, keyword `structural`);
- it is never an input. `convert_to_reaction_network` drops it, so the starting copy
  keeps its own value, and the solve removes it before any other rule reads the network;
- the benchmark does not count it when deciding roots, so a p130 knockout still pins
  the starting copy.

The solve reports `recycle_edges`.
"""

const RECYCLE_EDGE = "recycle"

"""`(network without recycle edges, the (source, target) pairs they joined)`."""
function split_recycle_edges(network::ReactionNetwork)
    any(e -> e.edge_type == RECYCLE_EDGE, network.edges) || return network, Tuple{String, String}[]
    pairs = sort!(unique([(e.parent_uuid, e.child_uuid) for e in network.edges if e.edge_type == RECYCLE_EDGE]))
    kept = [e for e in network.edges if e.edge_type != RECYCLE_EDGE]
    out = ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, network.pools,
                          network.pathogen_stids, network.carrier_loops, network.node_leaves)
    return out, pairs
end
