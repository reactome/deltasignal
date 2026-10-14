# specs/054 T5: baseline elasticity (d log F / d log x) of every input, the cyclic
# components it closes, and which edges carry the leading eigenvalue.
using DeltaSignal; const DS = DeltaSignal
using LinearAlgebra, SparseArrays
dir = "/catalog/R-HSA-5673001"
net0 = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
net, = DS.apply_carrier_loops(net0, DS.carrier_loop_mode())
rx = DS.convert_to_reaction_network(net)
nodes = sort(collect(keys(net.nodes))); idx = Dict(u => i for (i, u) in enumerate(nodes))
bl = Dict(u => n.baseline for (u, n) in net.nodes)
ss = first(DS.self_inhibitor_setup(net))
rxns, comp, ncomp = DS.index_reactions(rx, idx, bl, Set{String}(); self_shared = ss)
cfg = DS.resolve_reaction_eval_config()
held = Set(idx[u] for u in union(DS.cofactor_uuids(net), DS.pathogen_uuids(net)) if haskey(idx, u))
x = [bl[u] for u in nodes]
nm(i) = string(something(net.nodes[nodes[i]].reactome_id, "?"), " ", first(net.nodes[nodes[i]].display_name, 50))
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
outs = String.(split(first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")[9], "|"))
pins = String.(split(first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")[8], "|"))
# Elasticity d log F_t / d log x_s, both sides, over every input of every reaction.
h = 1e-4
I = Int[]; J = Int[]; V = Float64[]; kink = 0
for r in rxns
    t = r.target_idx; t in held && continue
    base = DS.compute_reaction_output_vec(x, r; config = cfg)
    for s in unique(vcat(r.activator_indices, r.inhibitor_indices, r.depletion_indices))
        s in held && continue
        xs = x[s]
        x[s] = xs * exp(h);  up = DS.compute_reaction_output_vec(x, r; config = cfg)
        x[s] = xs * exp(-h); dn = DS.compute_reaction_output_vec(x, r; config = cfg)
        x[s] = xs
        gu = (log(max(up, 1e-300)) - log(max(base, 1e-300))) / h
        gd = (log(max(base, 1e-300)) - log(max(dn, 1e-300))) / h
        abs(gu - gd) > 1e-3 && (global kink += 1)
        g = (gu + gd) / 2
        abs(g) > 1e-9 && (push!(I, t); push!(J, s); push!(V, g))
    end
end
n = length(nodes)
println("reactions ", length(rxns), "; nodes ", n, "; elastic edges ", length(V), "; kinked at baseline ", kink,
        "; baseline consistent: ", all(abs(DS.compute_reaction_output_vec(x, r; config = cfg) - x[r.target_idx]) < 1e-12 for r in rxns))
adj = [Set{Int}() for _ in 1:n]; for k in eachindex(V); push!(adj[J[k]], I[k]); end
cid, nc = DS.tarjan_scc_components(adj)
sz = zeros(Int, nc); for c in cid; sz[c] += 1; end
E = sparse(I, J, V, n, n)
outc = unique(cid[idx[u]] for u in outs if haskey(idx, u))
upc = Set{Int}()  # components upstream of the readouts
let seen = Set(idx[u] for u in outs), st = collect(seen)
    radj = [Int[] for _ in 1:n]; for k in eachindex(V); push!(radj[I[k]], J[k]); end
    while !isempty(st); v = pop!(st); push!(upc, cid[v]); for w in radj[v]; w in seen || (push!(seen, w); push!(st, w)); end; end
end
pinc = Set(cid[idx[u]] for u in pins if haskey(idx, u))
println("readout components ", [(c, sz[c]) for c in outc], "; pin components ", [(c, sz[c]) for c in pinc])

sid(i) = something(net.nodes[nodes[i]].reactome_id, "?")
sid0(i) = first(split(sid(i), "::"))
for c in sort([c for c in 1:nc if sz[c] > 1 && c in upc], by = c -> -sz[c])
    m = findall(==(c), cid); Ec = Matrix(E[m, m])
    F = eigen(Ec); k = argmax(real.(F.values)); lam = real(F.values[k])
    lam > 1.0001 || continue
    v = real.(F.vectors[:, k]); w = real.(eigen(Matrix(Ec')).vectors[:, argmax(real.(eigen(Matrix(Ec')).values))])
    sc = dot(w, v)
    println("\n== component size ", sz[c], ", leading eigenvalue ", round(lam; digits = 4))
    # row sums of in-component elasticity: reactions where the loop enters more than once
    agg = Dict{Tuple{String,String},Float64}()
    rs = Dict{String,Tuple{Float64,Int}}()
    for a in 1:length(m), b in 1:length(m)
        g = Ec[a, b]; g == 0 && continue
        s = w[a] * v[b] * g / sc          # share of lambda carried by edge b -> a
        key = (sid0(m[b]), sid0(m[a])); agg[key] = get(agg, key, 0.0) + s
    end
    for a in 1:length(m)
        r = sum(Ec[a, :]); rs0 = get(rs, sid0(m[a]), (0.0, 0)); rs[sid0(m[a])] = (max(rs0[1], r), rs0[2] + 1)
    end
    tot = sum(values(agg))
    println("   edge shares of lambda (sum ", round(tot; digits = 4), "), top 20 by entity pair:")
    for (key, s) in sort(collect(agg), by = x -> -abs(x[2]))[1:min(20, length(agg))]
        println("   ", round(s / tot; digits = 3), "  ", key[1], " -> ", key[2])
    end
    println("   targets whose in-loop elasticities sum > 1.0001 (max over copies, copies):")
    for (t, (r, n)) in sort(collect(rs), by = x -> -x[2][1])
        r > 1.0001 || break
        println("   ", round(r; digits = 3), "  ", t, "  (", n, " copies)")
    end
end
