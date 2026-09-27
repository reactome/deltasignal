"""
Drug-derived entities (specs/032).

Reactome curates drug actions inside signalling pathways: MAP2K and MAPK
inhibitors bind the RAF:scaffold:MAP2K:MAPK complex in the RAF/MAP kinase
cascade, and PARP inhibitors bind PARP in DNA repair. The generator keeps them,
because they are curated, and ships `drugs.csv` naming which nodes are
drug-derived (a Drug; a complex with any drug component; a set whose every
member is one).

A benchmark case, or any question about an untreated cell, describes a cell
WITHOUT the drug. Left to propagate, a drug is a root held at baseline, so
raising RAS raises the drug-bound complexes, and their inhibitory edges then
shut the kinase cascade off: every RAS/RAF perturbation in RAF/MAP kinase read
DOWN (experimental 28.6% against MP-BioPath's 93.9%, specs/032).

  DS_DRUG_MODE=propagate  (default) — previous behaviour.
  DS_DRUG_MODE=inert      — drug-derived nodes are pinned at baseline, as
                            cofactors are: participants of fold 1.0 that
                            never carry a perturbation. Not 0, because under
                            `divide` an inhibitor at 0 de-represses tenfold,
                            which would model drug WITHDRAWAL. An explicit
                            observation still wins.
"""

const DS_VALID_DRUG_MODES = Set(["propagate", "inert"])

function drug_mode()::String
    value = get(ENV, "DS_DRUG_MODE", "propagate")
    if !(value in DS_VALID_DRUG_MODES)
        throw(ArgumentError(
            "DS_DRUG_MODE must be one of $(join(sort(collect(DS_VALID_DRUG_MODES)), ", ")); got $(repr(value))"))
    end
    return value
end

"""uuids in this network whose Reactome stable id the bundle declares drug-derived."""
function drug_uuids(network::ReactionNetwork)::Set{String}
    out = Set{String}()
    stids = network.drug_stids
    stids === nothing && return out
    for (uuid, node) in network.nodes
        rid = node.reactome_id
        if rid !== nothing && rid in stids
            push!(out, uuid)
        end
    end
    return out
end
