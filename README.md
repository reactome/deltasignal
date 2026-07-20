<div align="center">
  <img src="assets/logo.svg" alt="DeltaSignal Logo" width="400">
  
  <p><strong>A Pathway Perturbation & Dynamics Engine</strong></p>
  
  <p>
    <strong>Author:</strong> Adam Wright (OICR)<br>
    <strong>Date:</strong> August 25, 2024<br>
    <strong>Version:</strong> 1.0.0 Beta
  </p>
  
  <p>
    <img src="https://img.shields.io/badge/Julia-v1.10+-blue?logo=julia" alt="Julia Version">
    <img src="https://img.shields.io/badge/License-Apache_2.0-blue" alt="License">
    <img src="https://img.shields.io/badge/Biological_Realism-96%25-brightgreen" alt="Biological Realism">
    <img src="https://img.shields.io/badge/Status-Beta-orange" alt="Status">
  </p>
  
  <p>
    <strong>🧬 Enhanced Biological Accuracy • 🚀 High Performance • 🎨 Dual Visualization</strong>
  </p>
</div>

## Overview

DeltaSignal is a modeling and inference engine for Reactome-scale biological pathways that supports two modes:

- **Steady-State (SS) Mode** — infers node activities that are consistent with observed perturbations and the network's directed causal logic
- **Time-Dynamic (TD) Mode** — simulates discrete-time evolution of node activities with substrate consumption and product accumulation *(coming soon)*

## 🚀 Quick Start

### Installation

1. **Install Julia 1.10+**
   ```bash
   # Download Julia from https://julialang.org/downloads/
   wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.7-linux-x86_64.tar.gz
   tar -xzf julia-1.10.7-linux-x86_64.tar.gz
   export PATH="$PWD/julia-1.10.7/bin:$PATH"
   ```

2. **Clone and setup DeltaSignal**
   ```bash
   git clone <repository-url>
   cd deltasignal
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

### Basic Usage

```bash
# Parse a logic network from TSV files
julia cli/deltasignal.jl parse \
  --logic examples/sample_logic_network.tsv \
  --uuid-map examples/sample_uuid_mapping.tsv \
  --set-map examples/sample_set_mappings.tsv \
  --output parsed_network.json \
  --validate

# Solve steady-state given observations  
julia cli/deltasignal.jl solve \
  --network parsed_network.json \
  --observations examples/sample_observations.csv \
  --output results.json \
  --aggregation stoichiometry_weighted
```

### Test the Implementation

```bash
# Test basic parsing
julia test/test_basic.jl

# Test steady-state solver
julia test/test_steady_state.jl
```

## 📁 Repository Structure

```
DeltaSignal.jl/
├── src/                          # Core implementation
│   ├── DeltaSignal.jl           # Main module
│   ├── core/                    # Mathematical operations
│   │   ├── sensitivity.jl       # Input sensitivity transforms
│   │   ├── aggregators.jl       # Multi-input aggregation
│   │   ├── hill_functions.jl    # Hill activation functions
│   │   └── reaction_model.jl    # Complete reaction model
│   ├── io/                      # Data input/output
│   │   ├── tsv_parser.jl        # TSV logic network parser
│   │   └── reactome_mapper.jl   # Reactome pathway mapping
│   └── solvers/                 # Numerical solvers
│       ├── steady_state.jl      # SS solver implementation
│       └── time_dynamic.jl      # TD solver (coming soon)
├── cli/                         # Command-line interface
│   └── deltasignal.jl          # Main CLI script
├── test/                        # Test suite
├── examples/                    # Example data files
└── docs/                        # Documentation
```

## 🧮 Mathematical Model

### Core Reaction Model

For each node `r` receiving inputs from activators, inhibitors, and substrates:

1. **Sensitivity Transform**: `x̃ = x^α(x)` where `α(x) = 1 + s·x^n/(x^n + K_α^n)`
2. **Activator Aggregation**: `A = exp(Σᵢ wᵢ log(x̃ᵢ + ε))` (geometric mean)
3. **Inhibitor Suppression**: `H = ∏ⱼ 1/(1 + βⱼ xⱼ^mⱼ)` (Hill inhibition)
4. **Substrate Availability**: `L = exp(Σₖ uₖ log(xₖ + ε))` (soft AND)
5. **Pre-activation**: `s = A·H·L`
6. **Final Output**: `y = s^h/(s^h + K^h)` (Hill activation)

### Steady-State Optimization

Minimizes: `Σᵢ ωᵢ(xᵢ - yᵢ)² + μ‖x - F(x;θ)‖² + γ‖x - x₀‖²`

Where:
- `yᵢ` are observed node activities (0-100 UI scale)
- `F(x;θ)` is the forward model
- `x₀` are baseline activities
- Internal computation uses [0,1] normalization

## 📊 Data Formats

### Logic Network TSV
```
parent-001	child-001	1	1	1
parent-002	child-001	1	1	1
parent-003	child-002	0	1	2
```
Columns: `Parent UUID | Child UUID | AND/OR (1/0) | Pos/Neg (1/-1) | Stoichiometry`

### UUID Mapping TSV
```
parent-001	REACT:R-HSA-123456	protein	set-001
parent-002	REACT:R-HSA-123457	protein	set-001
parent-003	REACT:R-HSA-123458	small_molecule	
```
Columns: `Network UUID | Reactome DB ID | Entity Type | Set ID (optional)`

### Set Mappings TSV
```
set-001	PI3K Complex	parent-001,parent-002
set-002	mTORC1 Complex	parent-005
```
Columns: `Set ID | Original Name | Member UUIDs (comma-separated)`

### Observations CSV
```
node_uuid,activity,confidence
parent-001,75.0,0.9
parent-003,50.0,0.8
child-002,25.0,0.7
```

## 🎯 Key Features

### ✅ Currently Implemented

- **TSV Logic Network Parsing** with UUID mapping and set expansion handling
- **Mathematical Core**: Sensitivity transforms, multi-input aggregators, Hill functions
- **Steady-State Solver**: SCC-condensation feed-forward with observations pinned as hard constraints
- **Explainability**: Influence scoring (shares the live propagator's math)
- **CLI Interface**: Parse, solve, and export commands
- **Dual Scale Support**: 0-100 UI scale with internal 0-1 normalization
- **Reactome Integration**: Mapping between expanded networks and original pathways

### 🔄 In Development

- **Web Frontend**: an Angular UI in the WebsiteAngular workspace, consuming this repo's HTTP API (contract in `docs/API.md`) — reuses Reactome's pathway-browser + cytoscape styling
- **Time-Dynamic Mode / Parameter Learning**: prototyped then shelved; source archived in `attic/`

### 📈 Roadmap

- **Advanced Solvers**: GPU acceleration, sparse optimization
- **Validation Suite**: Benchmark against CRISPR/drug perturbation datasets
- **Uncertainty Quantification**: Confidence intervals and parameter sensitivity
- **Pathway-Scale Deployment**: Docker/Kubernetes for production use

## 🧪 Example Results

**Sample Network Analysis:**
```
📊 Final node activities (penalty method):
  child-001: 0.0%
  child-002: 34.0% [OBS] ← Observed: 25%
  child-003: 0.0%
  parent-001: 67.5% [OBS] ← Observed: 75%
  parent-002: 0.0%
  parent-003: 38.0% [OBS] ← Observed: 50%

🎯 Top influential nodes:
  1. parent-002: 269.9271
  2. parent-004: 49.066  
  3. parent-003: 0.6457

🔍 Upstream driver suggestions (to increase child-002 by 20%):
  parent-002: +8.4%
  parent-004: +1.5%
```

## 🛠️ Technical Implementation

- **Language**: Julia 1.10+ LTS for high-performance numerical computing
- **Core Dependencies**: Optim.jl, JSON3.jl, DataFrames.jl, CSV.jl, HTTP.jl
- **Architecture**: Modular design supporting both CLI and API interfaces
- **Performance**: Targets <30s for 1000-node networks, <5min for Reactome-scale
- **Extensibility**: Plugin architecture for custom aggregators and solvers

## 📚 Documentation

### CLI Commands

```bash
# Available commands
deltasignal parse      # Parse TSV logic networks
deltasignal solve      # Solve steady-state
deltasignal rollout    # Time-dynamic rollout (coming soon)
deltasignal train      # Parameter learning (coming soon)  
deltasignal validate   # Model validation (coming soon)
deltasignal export     # Export results for visualization (coming soon)
deltasignal server     # Start web API server (coming soon)

# Get help for any command
deltasignal parse --help
```

### API Design (Planned)

```bash
# REST endpoints  
POST /api/parse         # Parse networks
POST /api/ss/solve      # Steady-state solving
POST /api/td/rollout    # Time-dynamic rollout
POST /api/explain       # Explainability analysis
```

## 🤝 Contributing

1. **Core Mathematics**: Improve sensitivity transforms, aggregators, solvers
2. **Visualization**: PathwayBrowser integration, interactive overlays
3. **Performance**: GPU acceleration, sparse matrix optimization
4. **Validation**: Benchmark datasets, biological validation
5. **Documentation**: Tutorials, API docs, examples

## 📄 License

Licensed under the **Apache License, Version 2.0**. See the [LICENSE](LICENSE)
file for the full text and [NOTICE](NOTICE) for attribution.

Copyright © 2025 Ontario Institute for Cancer Research (OICR)

## 🔗 References

- **Reactome Database**: https://reactome.org/
- **PathwayBrowser**: https://github.com/reactome/PathwayBrowser
- **Technical Specification**: See full mathematical specification in project documentation

---

**Status**: MVP implementation complete with parsing, steady-state solving, and explainability features. Time-dynamic mode and web interface in development.