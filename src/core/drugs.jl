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

"""
Pathogen-derived entities (specs/048), held exactly as drugs are.

Reactome curates host-pathogen interactions inside human pathways: the
SARS-CoV-2 N:M:PDPK1 complex negatively regulates "PDPK1 phosphorylates AKT at
T308". A PDPK1 knockout zeroes that complex, de-represses the reaction, and the
AKT loop rails to 100x where the truth is DOWN. In a cell without the virus the
complex does not exist; held at baseline it is a constant that cancels in
fold-change. The generator lists the nodes in `pathogens.csv` (a leaf with a
species, none human; a complex with any such component; a set whose every
member is one).

  DS_PATHOGEN_MODE=inert      (default since specs/048 amendment 2) —
                                pathogen-derived nodes pinned at baseline, as
                                DS_DRUG_MODE=inert pins drugs. An explicit
                                observation still wins.
  DS_PATHOGEN_MODE=propagate  — the previous behaviour.

What the bundle lists is the generator's call. Since LNG_PATHOGEN_PROTEIN
(default on) only entities carrying a pathogen PROTEIN are listed: viral RNA
is the ligand DDX58/IFIH1 senses, and holding it cut 12 curator cases.
"""
const DS_VALID_PATHOGEN_MODES = Set(["propagate", "inert"])

function pathogen_mode()::String
    value = get(ENV, "DS_PATHOGEN_MODE", "inert")
    if !(value in DS_VALID_PATHOGEN_MODES)
        throw(ArgumentError(
            "DS_PATHOGEN_MODE must be one of $(join(sort(collect(DS_VALID_PATHOGEN_MODES)), ", ")); got $(repr(value))"))
    end
    return value
end

"""uuids in this network whose Reactome stable id the bundle declares pathogen-derived."""
pathogen_uuids(network::ReactionNetwork)::Set{String} = _uuids_with_stid(network, network.pathogen_stids)

function _uuids_with_stid(network::ReactionNetwork, stids)::Set{String}
    out = Set{String}()
    stids === nothing && return out
    for (uuid, node) in network.nodes
        rid = node.reactome_id
        if rid !== nothing && rid in stids
            push!(out, uuid)
        end
    end
    return out
end
