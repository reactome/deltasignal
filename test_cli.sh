#!/bin/bash

# Test CLI commands

set -e

echo "Testing DeltaSignal CLI..."
echo "=========================="

# Create output directory
mkdir -p data/output

# Test 1: Parse command
echo ""
echo "Test 1: Parse network"
echo "---------------------"
docker compose -f docker-compose.dev.yml run --rm julia-api \
  julia --project=/app /app/cli/deltasignal.jl \
  parse \
  --logic /app/examples/sample_logic_network.tsv \
  --uuid-map /app/examples/sample_uuid_mapping.tsv \
  --set-map /app/examples/sample_set_mappings.tsv \
  --output /app/data/output/test_network.json \
  --validate

# Test 2: Solve command
echo ""
echo "Test 2: Solve steady-state"
echo "--------------------------"
docker compose -f docker-compose.dev.yml run --rm julia-api \
  julia --project=/app /app/cli/deltasignal.jl \
  solve \
  --network /app/data/output/test_network.json \
  --observations /app/examples/sample_observations.csv \
  --output /app/data/output/test_results.json

# Test 3: Export command
echo ""
echo "Test 3: Export results"
echo "---------------------"
docker compose -f docker-compose.dev.yml run --rm julia-api \
  julia --project=/app /app/cli/deltasignal.jl \
  export \
  --results /app/data/output/test_results.json \
  --format csv \
  --output /app/data/output/test_export.csv

echo ""
echo "✅ All CLI tests passed!"
