Diagnostics for specs/057 (ubiquitination substrate).

- `ub_reactions.cypher`: every ubiquitination (Ub an input) in a catalog
  pathway, with input leaves (stId, genes) and output leaves. Run with
  cypher-shell `--format plain`; the result is the first argument to
- `ub_census.py RESULT BUILD_DIR`: the depletion edges aimed at a partner
  that passes through its ubiquitination unchanged.
- `met_trace.jl`: Signaling by MET solved on the canonical bundle and on a
  copy without the p-Y-CBL pool -| GRB2 depletion edge (`/s/metx`), pinning
  p-Y1003,4Y-MET, it plus plain MET (93ce1da2), or HGF. Reads `/s/met_cases.tsv`,
  the MET rows of the canonical curator case log in the specs/054 layout:

      awk -F'\t' 'NR>1 && $1=="Signaling_by_MET" && $7==1 {print "cur\t"$2"\t"$3"\t"$4"\t"$6"\t"$5"\t"$11"\t"$12"\t"$13}'

  Run in the dev container with the build at /app/pathway_catalog.
- `ub_complex_substrates.cypher`: per ubiquitination, the output complexes
  that hold Ub as a direct component and their other direct components (the
  real substrate in Reactome's "Ub-X" complex representation; research R1).
- `dump.jl`: one solve of a bundle (env PW = bundle dir, PINS = uuids joined
  by |, V = pinned value, OUT = file of uuid, activity on the 0-100 scale).
