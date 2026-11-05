# specs/136: slots the leaf clause adds (DS_SELF_ACTIVATOR_LEAVES=1), per pathway.
using DeltaSignal; const DS = DeltaSignal
ENV["DS_SELF_ACTIVATOR_LEAVES"] = "1"
root = ENV["BUILD"]; tot = 0; strict = 0
for d in sort(readdir(root))
    isfile(joinpath(root, d, "logic_network.csv")) || continue
    net = redirect_stdout(devnull) do
        DS.parse_complete_network(joinpath(root, d, "logic_network.csv"), joinpath(root, d, "stid_to_uuid_mapping.csv"), nothing)
    end
    m, rule, n = DS.self_activator_setup(net)
    global tot += n; global strict += length(m) - n
    n > 0 && println(d, "\t", n)
end
println("LEAF_SLOTS\t", tot, "\tALL_FLAGGED\t", tot + strict)
