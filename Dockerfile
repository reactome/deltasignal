# DeltaSignal - Multi-stage Docker build
# Stage 1: Julia backend build
FROM julia:1.10-bullseye as julia-builder

WORKDIR /app

# Copy Julia project files
COPY Project.toml .
COPY src/ ./src/

# Install Julia dependencies (will generate Manifest.toml)
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy remaining code
COPY cli/ ./cli/
COPY examples/ ./examples/
COPY test/ ./test/

# Stage 2: Node.js frontend build
FROM node:18-alpine as frontend-builder

WORKDIR /app/frontend

# Copy package files
COPY frontend/package.json frontend/package-lock.json ./

# Install dependencies
RUN npm ci --only=production

# Copy frontend source
COPY frontend/ ./

# Build frontend
RUN npm run build

# Stage 3: Production runtime
FROM julia:1.10-bullseye

# Install system dependencies
RUN apt-get update && apt-get install -y \
    nginx \
    supervisor \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Copy Julia application from builder
COPY --from=julia-builder /app ./

# Copy frontend build from builder
COPY --from=frontend-builder /app/frontend/dist ./public

# Copy configuration files
COPY docker/nginx.conf /etc/nginx/nginx.conf
COPY docker/supervisord.conf /etc/supervisor/conf.d/supervisord.conf
COPY docker/start.sh /start.sh

# Make start script executable
RUN chmod +x /start.sh

# Create necessary directories
RUN mkdir -p /var/log/supervisor /var/log/nginx /var/run

# Expose ports
EXPOSE 80 8080

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost/api/health || exit 1

# Start supervisor to manage services
CMD ["/start.sh"]