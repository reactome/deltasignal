# specs/061 census: every depletion edge D -| T whose depleter's entity contains
# T's entity (the specs/060 amendment-1 pairs), and the SOURCE copies the rule
# would read: nodes with T's entity reached walking backward from D along
# non-depletion edges, stopping at each such copy. Run per catalog build:
#   CAT=<build dir> julia --project=/app own_source_census.jl > out.tsv
using DeltaSignal; const DS = DeltaSignal
ENV["DS_DEPLETION_OWN_PRODUCT"] = "suppress_only_entity"
cat_dir = ENV["CAT"]
println("pathway\tdepleter\ttarget\ttarget_stid\tn_sources\ttarget_is_source\tsources")
for pw in sort(readdir(cat_dir))
    dir = joinpath(cat_dir, pw)
    isfile(joinpath(dir, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"),
                                    joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
    pairs = DS.own_entity_depletion_pairs(net)
    isempty(pairs) && continue
    back = Dict{String, Vector{String}}()
    for e in net.edges
        e.edge_type == "depletion" && continue
        push!(get!(back, e.child_uuid, String[]), e.parent_uuid)
    end
    for (d, t) in sort(collect(pairs))
        tid = net.nodes[t].reactome_id
        seen = Set{String}([d]); stack = copy(get(back, d, String[])); src = String[]
        while !isempty(stack)
            n = pop!(stack)
            n in seen && continue
            push!(seen, n)
            if net.nodes[n].reactome_id == tid
                push!(src, n); continue
            end
            append!(stack, get(back, n, String[]))
        end
        sort!(src)
        println(pw, "\t", d, "\t", t, "\t", tid, "\t", length(src), "\t", t in src, "\t", join(src, ","))
    end
end
