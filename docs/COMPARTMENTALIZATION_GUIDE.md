# 🏠 Compartmentalization in DeltaSignal

## Understanding Reactome Compartmentalization

### How Reactome Handles Compartments

In Reactome, **the same biological entity gets different identifiers for each compartment**:

```
Example: p53 protein
- Cytoplasmic p53:     REACT:R-HSA-391251  
- Nuclear p53:         REACT:R-HSA-69620
- Mitochondrial p53:   REACT:R-HSA-109606
```

### Logic Network Structure

The **logic network files already reflect this compartmentalization**:

```tsv
# Different UUIDs for the same protein in different compartments
uuid-p53-cyto    uuid-p53-nucleus    1    1    1    # Transport to nucleus
uuid-p53-nucleus uuid-gene-p21       1    1    1    # Nuclear transcription
```

This design **prevents non-biological loops** by properly representing spatial separation.

## DeltaSignal Implementation

### 1. **Reactome ID-Based Compartment Detection**

```julia
# Check Reactome ID patterns for compartment-specific identifiers
const REACTOME_COMPARTMENT_PATTERNS = Dict{Compartment, Vector{String}}(
    NUCLEUS => ["R-HSA-.*[Nn]ucleus", "R-HSA-.*[Nn]uclear", "R-HSA-.*[Cc]hromatin"],
    MITOCHONDRIA => ["R-HSA-.*[Mm]itochondria", "R-HSA-.*[Mm]itochondrial"],
    MEMBRANE => ["R-HSA-.*[Mm]embrane", "R-HSA-.*[Pp]lasma"],
    EXTRACELLULAR => ["R-HSA-.*[Ee]xtracellular", "R-HSA-.*[Ss]ecret"]
    # ... etc
)

function infer_node_compartment(node::NetworkNode)::Compartment
    reactome_id = get(node.reactome_id, "", "")
    
    # Check Reactome ID patterns first
    for (compartment, patterns) in REACTOME_COMPARTMENT_PATTERNS
        for pattern in patterns
            if occursin(Regex(pattern), reactome_id)
                return compartment
            end
        end
    end
    
    # Fallback to name-based inference
    return infer_compartment_from_name(node)
end
```

### 2. **Cross-Compartment Transport Detection**

```julia
# Detect when reaction inputs come from different compartments
for input_uuid in all_input_uuids
    input_compartment = infer_node_compartment(input_node)
    
    if input_compartment != target_compartment
        # Create transport parameters for cross-compartment reaction
        transport_params = create_transport_params(input_compartment, target_compartment)
        transport_map[input_uuid] = transport_params
    end
end
```

### 3. **Compartment-Specific Environments**

Each compartment has different biochemical conditions that affect kinetics:

```julia
# Nuclear environment (transcription)
CompartmentEnvironment(
    7.3,    # pH - slightly basic
    0.12,   # ionic strength 
    0.4,    # high protein crowding (chromatin)
    3e-3,   # lower ATP than cytoplasm
    -0.2    # slightly oxidizing
)

# Mitochondrial environment (energy metabolism)  
CompartmentEnvironment(
    7.8,    # more basic (matrix)
    0.2,    # high ionic strength
    0.5,    # very high protein density
    8e-3,   # high ATP synthesis
    -0.3    # very reducing (NADH/FADH2)
)
```

### 4. **Transport Kinetics**

Different compartment pairs have different transport characteristics:

```julia
# Nuclear import - active, energy-dependent
CYTOPLASM → NUCLEUS: TransportParams(0.1, 3.0, 0.3, 0.5, 0.9)

# Membrane receptor binding - high selectivity  
EXTRACELLULAR → MEMBRANE: TransportParams(0.8, 0.0, 0.1, 0.0, 0.9)

# Mitochondrial import - highly regulated
CYTOPLASM → MITOCHONDRIA: TransportParams(0.05, 1.5, 0.2, 0.8, 0.95)
```

## Example: Growth Factor Signaling

```
Extracellular Growth Factor (R-HSA-190236-extracellular)
    ↓ [receptor binding]
Membrane Receptor (R-HSA-109582-membrane)  
    ↓ [signal transduction]
Cytoplasmic PI3K (R-HSA-109704-cytoplasm)
    ↓ [nuclear import]
Nuclear FOXO (R-HSA-112399-nucleus)
    ↓ [transcription]
Nuclear p21 Gene (R-HSA-69620-nucleus)
```

Each arrow represents a properly compartmentalized reaction with:
- **Appropriate environmental parameters** for the target compartment
- **Transport kinetics** for cross-compartment steps  
- **Different Reactome IDs** for the same entity in different locations

## Benefits of Correct Compartmentalization

### 1. **Prevents Non-Biological Loops**
- Cytoplasmic protein A → Nuclear protein A (different UUIDs)
- No artificial feedback from nuclear back to cytoplasmic version

### 2. **Realistic Kinetics**
- Nuclear reactions: slower, higher cooperativity (chromatin effects)
- Mitochondrial reactions: different pH, high energy availability
- Membrane reactions: lipid environment effects

### 3. **Transport Modeling**  
- Cross-compartment steps have realistic delays
- Energy costs for active transport
- Selectivity factors for different molecules

### 4. **Proper Spatial Organization**
- Transcription factors must be in nucleus to regulate genes
- Metabolic enzymes have compartment-specific kinetics
- Signaling cascades follow proper spatial organization

## Usage in DeltaSignal

```julia
# Create compartmentalized network
comp_reactions = create_compartmentalized_reactions(network)

# Analysis shows proper compartment distribution
🏠 Compartment distribution:
  CYTOPLASM: 45 reactions (60%)
  NUCLEUS: 12 reactions (16%)  
  MITOCHONDRIA: 8 reactions (11%)
  MEMBRANE: 6 reactions (8%)
  EXTRACELLULAR: 4 reactions (5%)

🚚 Cross-compartment transport: 15 reactions
```

This approach respects the **Reactome compartment structure** while providing **realistic biochemical modeling** of spatial organization in cells.