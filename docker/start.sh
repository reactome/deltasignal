#!/bin/bash
set -e

echo "🧬 Starting DeltaSignal container..."

# Create necessary directories
mkdir -p /var/log/supervisor /var/log/nginx /var/run

# Initialize Julia if needed (warm up compilation)
echo "Warming up Julia compilation..."
julia --project=/app -e "using DeltaSignal; println(\"DeltaSignal module loaded successfully\")" || {
    echo "Warning: DeltaSignal module failed to load, continuing anyway..."
}

# Check if API server can start
echo "Testing API server startup..."
timeout 30s julia --project=/app /app/src/api/server.jl --port=8081 --test || {
    echo "Warning: API server test failed, will try in supervisor..."
}

# Start supervisor to manage all services
echo "Starting supervisor..."
exec /usr/bin/supervisord -c /etc/supervisor/conf.d/supervisord.conf