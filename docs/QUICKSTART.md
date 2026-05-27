# 🚀 DeltaSignal Quick Start Guide

This guide will get you running DeltaSignal in under 5 minutes.

## Prerequisites

- Julia 1.10 or later (LTS recommended)
- Basic familiarity with command line

## 1. Install Julia

```bash
# Linux/macOS
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.7-linux-x86_64.tar.gz
tar -xzf julia-1.10.7-linux-x86_64.tar.gz
export PATH="$PWD/julia-1.10.7/bin:$PATH"

# Verify installation
julia --version
```

## 2. Setup DeltaSignal

```bash
# Clone repository (or download files)
cd deltasignal

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 3. Test Installation

```bash
# Test basic functionality
julia test/test_basic.jl

# Test steady-state solver
julia test/test_steady_state.jl
```

Expected output:
```
🧬 Testing DeltaSignal TSV parsing...
✅ Parsed 6 edges
✅ Parsed 9 nodes  
✅ Parsed 2 set mappings
✨ Basic parsing test completed successfully!
```

## 4. Run Example Analysis

### Parse a pathway network:
```bash
julia cli/deltasignal.jl parse \
  --logic examples/sample_logic_network.tsv \
  --uuid-map examples/sample_uuid_mapping.tsv \
  --set-map examples/sample_set_mappings.tsv \
  --output my_network.json \
  --validate
```

### Solve steady-state:
```bash
julia cli/deltasignal.jl solve \
  --network my_network.json \
  --observations examples/sample_observations.csv \
  --output my_results.json
```

## 5. Understanding the Results

The solver will output node activities and influence scores:

```
📊 Final node activities:
  parent-001: 67.5% [OBS] ← Your observation was 75%
  child-002: 34.0% [OBS]  ← Your observation was 25%
  
🎯 Top influential nodes:
  1. parent-002: 269.93 ← Most important driver
  2. parent-004: 49.07  ← Secondary driver
```

## 6. Next Steps

- **Modify observations**: Edit `examples/sample_observations.csv` with your data
- **Try different networks**: Create your own TSV files following the format
- **Explore parameters**: Use `--help` to see solver options
- **Check validation**: Use `--validate` flag to verify network consistency

## Troubleshooting

**Julia packages fail to install:**
```bash
# Try updating package registry
julia -e 'using Pkg; Pkg.Registry.update()'
```

**CLI command not found:**
```bash
# Make sure you're in the project directory
cd deltasignal
julia cli/deltasignal.jl --help
```

**Test failures:**
```bash
# Check Julia version
julia --version  # Should be 1.10+

# Check project activation
julia --project=. -e 'using Pkg; Pkg.status()'
```

## What's Next?

- Read the [full README](README.md) for detailed documentation
- Explore the mathematical model in `src/core/`
- Check out planned features in the roadmap
- Join development - contributions welcome!

---

**Need help?** Open an issue or check the documentation.