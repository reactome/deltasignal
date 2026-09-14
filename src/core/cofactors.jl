"""Metabolic cofactors — represented in the networks, not propagated through.

The logic networks stay true to the curated data: the generator emits every
participant Reactome records, including ATP, ADP, H2O and Pi, so the networks
remain reusable by anyone. Whether a cofactor should carry a perturbation is a
MODELLING decision, and modelling decisions belong here (constitution I).

As one shared node a cofactor connects everything to everything — ATP reaches
degree 262 in DNA_Double-Strand_Break_Repair against 16 across MP-BioPath's
entire 85-network corpus — so a knockout anywhere can reach a readout by
travelling through the cell's energy currency. That is not a pathway.

But the fix is "participant, not conduit", NOT deletion, and the two measure
oppositely:

  DELETING them from the networks (an LNG-side filter, since abandoned) cost
  84 cases on 23,622 curator cases, macro-F1 0.7835 -> 0.7766. Only 25% of the
  losses became unreachable; 137 kept a path and simply weakened. Removing the
  node removes real structure with it.

  HOLDING them at baseline gains 37 on 21,450 curator cases (71 pathways,
  Release97), accuracy 0.8430 -> 0.8447, macro-F1 0.8095 -> 0.8108. 70
  predictions change; the experiment moves in NONE of the 21,986 shared cases,
  because nothing leaves the network. The mechanism is the one diagnosed:
  53 spurious change calls removed against 16 true ones lost — 27 UP->NO_CHANGE
  and 26 DOWN->NO_CHANGE correct, against 7 and 7 wrong.

  (The two numbers come from different pipelines and case sets and are not
  directly comparable to each other; each is a paired A/B against its own
  baseline on one shared catalog build.)

Per-pathway the gain concentrates in Signaling_by_MET (+18), DAP12 (+12) and
SCF-KIT (+12), against one regression, GPVI-mediated_activation_cascade (-12),
which is PIK3CA and PTPN11 into two readouts.

On the ten-pathway experimental set the mode is neutral: 13 predicted values
move and no case crosses a class boundary. Those networks carry 147 cofactor
nodes in total, so they cannot settle this question either way.

  DS_COFACTOR_MODE=inert      (default) — a cofactor stays an AND input, so
                                reaction completeness is unaffected, but it is
                                pinned at baseline and cannot carry a
                                perturbation.
  DS_COFACTOR_MODE=propagate  — previous behaviour, for reproducing earlier
                                results or re-testing the decision.

The list is chemistry, not curation: energy, redox, phosphate, one-carbon and
bulk ions, every compartment variant in Release97. It deliberately EXCLUDES
Ca2+, PI(3,4,5)P3, PI(4,5)P2, cAMP, cGMP, DAG and I(1,4,5)P3, which are second
messengers and are the signal in a signalling pathway — an earlier attempt to
exclude the SimpleEntity class cost 163 cases by deleting exactly those. It
also excludes ubiquitin and SUMO, whose transfer IS the regulatory event.

260 stable ids across 28 molecules.
"""

const COFACTOR_STIDS = Set([
    # ADP
    "R-ALL-113581", "R-ALL-113582", "R-ALL-114564", "R-ALL-114565",
    "R-ALL-196180", "R-ALL-211606", "R-ALL-29370", "R-ALL-5632457",
    "R-ALL-5696026", "R-ALL-8931884",
    # AMP
    "R-ALL-109275", "R-ALL-159448", "R-ALL-164121", "R-ALL-389620",
    "R-ALL-76577", "R-ALL-9931265",
    # ATP
    "R-ALL-113592", "R-ALL-113593", "R-ALL-114570", "R-ALL-114571",
    "R-ALL-139836", "R-ALL-211579", "R-ALL-29358", "R-ALL-389573",
    "R-ALL-5632460", "R-ALL-5696069", "R-ALL-8931885", "R-ALL-9935785",
    # AdoHcy
    "R-ALL-2162252", "R-ALL-5278409", "R-ALL-71285", "R-ALL-77502",
    # AdoMet
    "R-ALL-2162265", "R-ALL-5279190", "R-ALL-71284", "R-ALL-77087",
    # CDP
    "R-ALL-110094", "R-ALL-110638", "R-ALL-111809", "R-ALL-8851220",
    "R-ALL-8851368", "R-ALL-8851506",
    # CO2
    "R-ALL-1132345", "R-ALL-113528", "R-ALL-1237009", "R-ALL-159751",
    "R-ALL-159942", "R-ALL-189480", "R-ALL-29376", "R-ALL-389536",
    "R-ALL-5668565", "R-ALL-8863760",
    # CTP
    "R-ALL-110577", "R-ALL-110614", "R-ALL-29470", "R-ALL-8851086",
    "R-ALL-8851239", "R-ALL-8851515",
    # Cl-
    "R-ALL-188972", "R-ALL-2730999", "R-ALL-2731008", "R-ALL-2731020",
    "R-ALL-29572", "R-ALL-352022", "R-ALL-6788973", "R-ALL-879867",
    "R-ALL-9861558",
    # CoA-SH
    "R-ALL-162743", "R-ALL-1678675", "R-ALL-193514", "R-ALL-2485002",
    "R-ALL-29374", "R-ALL-76194", "R-ALL-8939024", "R-ALL-9766372",
    # FAD
    "R-ALL-113596", "R-ALL-113597", "R-ALL-141335", "R-ALL-141707",
    "R-ALL-189484", "R-ALL-2160490", "R-ALL-217259", "R-ALL-29386",
    "R-ALL-9861429",
    # FADH2
    "R-ALL-164934", "R-ALL-31649",
    # GDP
    "R-ALL-111349", "R-ALL-113525", "R-ALL-114549", "R-ALL-114622",
    "R-ALL-114623", "R-ALL-1467290", "R-ALL-1996292", "R-ALL-205689",
    "R-ALL-29420", "R-ALL-5617810", "R-ALL-8851528", "R-ALL-9628545",
    "R-ALL-9942379",
    # GMP
    "R-ALL-113578", "R-ALL-113579", "R-ALL-29626", "R-ALL-5696194",
    "R-ALL-744239",
    # GTP
    "R-ALL-113571", "R-ALL-113573", "R-ALL-114625", "R-ALL-114626",
    "R-ALL-1806221", "R-ALL-1996291", "R-ALL-2130170", "R-ALL-2213216",
    "R-ALL-29438", "R-ALL-5617813", "R-ALL-8851242", "R-ALL-8851508",
    "R-ALL-983318",
    # H+
    "R-ALL-1132304", "R-ALL-113529", "R-ALL-1470067", "R-ALL-156540",
    "R-ALL-1614597", "R-ALL-163953", "R-ALL-193465", "R-ALL-194688",
    "R-ALL-2000349", "R-ALL-2429673", "R-ALL-2872447", "R-ALL-351626",
    "R-ALL-372511", "R-ALL-374900", "R-ALL-425969", "R-ALL-425978",
    "R-ALL-425999", "R-ALL-427899", "R-ALL-428040", "R-ALL-428548",
    "R-ALL-5228597", "R-ALL-5244410", "R-ALL-5339575", "R-ALL-70106",
    "R-ALL-74722", "R-ALL-8858135", "R-ALL-9631150", "R-ALL-9668967",
    "R-ALL-9683057",
    # H2O
    "R-ALL-109276", "R-ALL-113518", "R-ALL-113519", "R-ALL-113521",
    "R-ALL-1222475", "R-ALL-141343", "R-ALL-1605715", "R-ALL-189422",
    "R-ALL-2022884", "R-ALL-2429665", "R-ALL-29356", "R-ALL-351603",
    "R-ALL-5278291", "R-ALL-5668574", "R-ALL-5685882", "R-ALL-5693747",
    "R-ALL-6781870", "R-ALL-8851517", "R-ALL-9923826", "R-ALL-9926997",
    # K+
    "R-ALL-29804", "R-ALL-5626313", "R-ALL-74126", "R-ALL-9859159",
    # NAD+
    "R-ALL-1132064", "R-ALL-113526", "R-ALL-192307", "R-ALL-194653",
    "R-ALL-29360", "R-ALL-352330", "R-ALL-427523", "R-ALL-5688282",
    "R-ALL-9912885",
    # NADH
    "R-ALL-1130844", "R-ALL-192305", "R-ALL-194697", "R-ALL-29362",
    "R-ALL-73473",
    # NADP+
    "R-ALL-1130860", "R-ALL-113563", "R-ALL-113564", "R-ALL-194668",
    "R-ALL-2000348", "R-ALL-29366", "R-ALL-351628", "R-ALL-389556",
    "R-ALL-5623650", "R-ALL-9749714", "R-ALL-9861421",
    # NADPH
    "R-ALL-1132417", "R-ALL-113600", "R-ALL-113601", "R-ALL-113602",
    "R-ALL-194725", "R-ALL-2000347", "R-ALL-29364", "R-ALL-351627",
    "R-ALL-5623644", "R-ALL-9749712",
    # Na+
    "R-ALL-2872443", "R-ALL-2889072", "R-ALL-2889083", "R-ALL-2892454",
    "R-ALL-425958", "R-ALL-425971", "R-ALL-425977", "R-ALL-74113",
    "R-ALL-83910",
    # O2
    "R-ALL-1131511", "R-ALL-113533", "R-ALL-113534", "R-ALL-113535",
    "R-ALL-113685", "R-ALL-1222561", "R-ALL-1236709", "R-ALL-189461",
    "R-ALL-2230955", "R-ALL-29368", "R-ALL-351593", "R-ALL-352327",
    "R-ALL-5668566",
    # PPi
    "R-ALL-111294", "R-ALL-113541", "R-ALL-113542", "R-ALL-114654",
    "R-ALL-159450", "R-ALL-2046049", "R-ALL-389593", "R-ALL-6806656",
    # Pi
    "R-ALL-109277", "R-ALL-113548", "R-ALL-113550", "R-ALL-113551",
    "R-ALL-114640", "R-ALL-2255331", "R-ALL-29372", "R-ALL-5228339",
    "R-ALL-8851226", "R-ALL-8851513", "R-ALL-947590", "R-ALL-9839058",
    # UDP
    "R-ALL-110096", "R-ALL-110732", "R-ALL-158602", "R-ALL-205687",
    "R-ALL-417888", "R-ALL-8851514", "R-ALL-9683078",
    # UTP
    "R-ALL-110731", "R-ALL-113562", "R-ALL-29494", "R-ALL-417877",
    "R-ALL-8851238", "R-ALL-8851510",
])

const DS_VALID_COFACTOR_MODES = Set(["inert", "propagate"])

"""Cofactor handling for this solve, validated loudly."""
function cofactor_mode()::String
    value = get(ENV, "DS_COFACTOR_MODE", "inert")
    if !(value in DS_VALID_COFACTOR_MODES)
        throw(ArgumentError(
            "DS_COFACTOR_MODE must be one of $(join(sort(collect(DS_VALID_COFACTOR_MODES)), ", ")); got $(repr(value))"))
    end
    return value
end

"""uuids in this network that are metabolic cofactors."""
function cofactor_uuids(network)::Set{String}
    out = Set{String}()
    for (uuid, node) in network.nodes
        rid = node.reactome_id
        if rid !== nothing && rid in COFACTOR_STIDS
            push!(out, uuid)
        end
    end
    return out
end
