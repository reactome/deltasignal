# DeltaSignal — engine + HTTP API image (single stage).
#
# The UI lives in the WebsiteAngular workspace and talks to this API over HTTP
# (see docs/API.md), so there is no frontend build here. Both compose files
# build the `julia-builder` target and run src/api/server.jl.
FROM julia:1.10-bullseye as julia-builder

WORKDIR /app

# System deps: curl for the container healthcheck.
RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

# Julia dependencies (Project.toml first for layer caching).
COPY Project.toml .
COPY src/ ./src/
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Application code.
COPY cli/ ./cli/
COPY examples/ ./examples/
COPY test/ ./test/

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost:8080/api/health || exit 1

# Default command runs the API bound to all interfaces; compose overrides as needed.
CMD ["julia", "--project=/app", "/app/src/api/server.jl", "--port=8080", "--host=0.0.0.0"]
