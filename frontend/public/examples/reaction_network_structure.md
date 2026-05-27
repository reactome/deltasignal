# Reaction Network Data Structure

## Overview
This document describes the data structure for reaction-based networks (like Reactome), which differs from the current regulatory/causal networks.

## Data Structure

### Nodes
Two types of nodes:
1. **Entity Nodes** (circular): Proteins, metabolites, complexes, etc.
2. **Reaction Nodes** (diamond): Chemical reactions, transport, binding events

### Network Format
```json
{
  "network_type": "reaction",
  "entities": [
    {
      "uuid": "entity_glucose",
      "name": "Glucose",
      "type": "metabolite",
      "compartment": "cytosol"
    },
    {
      "uuid": "entity_g6p",
      "name": "Glucose-6-phosphate",
      "type": "metabolite",
      "compartment": "cytosol"
    },
    {
      "uuid": "entity_hexokinase",
      "name": "Hexokinase",
      "type": "protein",
      "compartment": "cytosol"
    }
  ],
  "reactions": [
    {
      "uuid": "reaction_hexokinase",
      "name": "Glucose phosphorylation",
      "type": "biochemical_reaction",
      "ec_number": "2.7.1.1",
      "compartment": "cytosol"
    }
  ],
  "connections": [
    {
      "source": "entity_glucose",
      "target": "reaction_hexokinase",
      "role": "substrate"
    },
    {
      "source": "entity_atp",
      "target": "reaction_hexokinase",
      "role": "substrate"
    },
    {
      "source": "entity_hexokinase",
      "target": "reaction_hexokinase",
      "role": "catalyst"
    },
    {
      "source": "reaction_hexokinase",
      "target": "entity_g6p",
      "role": "product"
    },
    {
      "source": "reaction_hexokinase",
      "target": "entity_adp",
      "role": "product"
    }
  ]
}
```

## Roles
- **substrate**: Input to reaction
- **product**: Output from reaction
- **catalyst**: Enzyme/protein that catalyzes reaction
- **modifier**: Regulatory molecule (activator/inhibitor)
- **cofactor**: Required cofactor

## Reaction Types
- **biochemical_reaction**: Standard enzymatic reaction
- **transport**: Movement between compartments
- **binding**: Protein-protein or protein-DNA binding
- **complex_formation**: Assembly of protein complexes
- **degradation**: Protein/metabolite breakdown

## Visualization Mapping
- **Entity nodes**: Circle shape, colored by type
- **Reaction nodes**: Diamond shape, colored by reaction type
- **Edges**: Different styles based on role (substrate→reaction, reaction→product, etc.)