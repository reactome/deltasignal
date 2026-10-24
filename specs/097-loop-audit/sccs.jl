# Iterative Tarjan: the non-trivial strongly connected components of a directed graph.
function sccs(nodes::Vector{String}, adj::Dict{String, Vector{String}})
    # iterative Tarjan
    index = Dict{String, Int}(); low = Dict{String, Int}(); onstack = Set{String}()
    stack = String[]; out = Vector{Vector{String}}(); i = 0
    for root in nodes
        haskey(index, root) && continue
        work = [(root, 1)]
        while !isempty(work)
            v, k = pop!(work)
            if k == 1
                i += 1; index[v] = i; low[v] = i; push!(stack, v); push!(onstack, v)
            end
            nbrs = get(adj, v, String[])
            recurse = false
            for j in k:length(nbrs)
                w = nbrs[j]
                if !haskey(index, w)
                    push!(work, (v, j + 1)); push!(work, (w, 1)); recurse = true; break
                elseif w in onstack
                    low[v] = min(low[v], index[w])
                end
            end
            recurse && continue
            if low[v] == index[v]
                comp = String[]
                while true
                    w = pop!(stack); delete!(onstack, w); push!(comp, w)
                    w == v && break
                end
                length(comp) > 1 && push!(out, comp)
            end
            if !isempty(work)
                u = work[end][1]
                low[u] = min(low[u], low[v])
            end
        end
    end
    out
end

