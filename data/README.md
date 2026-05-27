# Data Directory

This directory is for storing input files, output results, and cached data during development.

## Structure

```
data/
├── input/     # Input network files, observations, etc.
├── output/    # Generated results, analysis outputs
└── cache/     # Temporary cached computations
```

## Usage

**Examples:**
```bash
# Save parsed network
julia cli/deltasignal.jl parse \
  --logic examples/sample_logic_network.tsv \
  --uuid-map examples/sample_uuid_mapping.tsv \
  --output data/output/my_network.json

# Save analysis results
julia cli/deltasignal.jl solve \
  --network data/output/my_network.json \
  --observations examples/sample_observations.csv \
  --output data/output/results_$(date +%Y%m%d).json
```

## Note

This directory is **gitignored** - files here won't be committed to version control.
Use `examples/` for sample data that should be versioned.
