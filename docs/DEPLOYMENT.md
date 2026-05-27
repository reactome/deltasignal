# 🚀 DeltaSignal Deployment Guide

Complete deployment guide for DeltaSignal pathway analysis system.

## Prerequisites

### System Requirements
- **Docker** 20.10+ and Docker Compose 2.0+
- **8GB+ RAM** (recommended for pathway analysis)
- **4+ CPU cores** (for parallel Julia computation)

### Optional (for local development)
- **Julia 1.10+** LTS
- **Node.js 18+** and npm

## Quick Start

### 1. Development Environment

```bash
# Clone repository
git clone <repository-url>
cd deltasignal

# Start development environment
./scripts/deploy.sh development
```

Access at:
- **Frontend**: http://localhost:3000
- **API**: http://localhost:8081

### 2. Production Deployment

```bash
# Build and deploy production
./scripts/deploy.sh production
```

Access at:
- **Application**: http://localhost
- **API**: http://localhost:8080

## Deployment Options

### Docker Compose (Recommended)

#### Development
```bash
# Start with hot reload
docker-compose --profile dev up --build

# View logs
docker-compose logs -f deltasignal-dev

# Stop
docker-compose down
```

#### Production
```bash
# Build and start
docker-compose up --build -d

# Scale API instances
docker-compose up --scale deltasignal=3 -d

# View logs
docker-compose logs -f deltasignal
```


### Manual Docker Build

```bash
# Build image
./scripts/build.sh

# Run container
docker run -d \
  --name deltasignal \
  -p 80:80 \
  -p 8080:8080 \
  -e JULIA_NUM_THREADS=4 \
  deltasignal:latest

# View logs
docker logs -f deltasignal
```

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `JULIA_NUM_THREADS` | 4 | Julia computation threads |
| `NODE_ENV` | production | Node.js environment |
| `API_BASE_URL` | http://localhost:8080 | Backend API URL |

### Resource Requirements

#### Minimum
- **Memory**: 2GB
- **CPU**: 2 cores
- **Disk**: 5GB

#### Recommended
- **Memory**: 4GB
- **CPU**: 4 cores
- **Disk**: 10GB

#### High-performance
- **Memory**: 8GB+
- **CPU**: 8+ cores
- **Disk**: 20GB+ SSD

## Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   Web Browser   │────│      Nginx       │────│  Julia API      │
│                 │    │   (Port 80)      │    │  (Port 8080)    │
│  - Cytoscape.js │    │                  │    │                 │
│  - D3.js        │    │  - Static files  │    │  - DeltaSignal  │
│  - React UI     │    │  - API proxy     │    │  - HTTP.jl      │
└─────────────────┘    └──────────────────┘    └─────────────────┘
```

## Health Checks

### Application Health
```bash
# Frontend health
curl http://localhost/health

# API health
curl http://localhost:8080/api/health

# Docker health check
docker inspect --format='{{.State.Health.Status}}' deltasignal
```


## Monitoring

### Logs

#### Docker Compose
```bash
# All services
docker-compose logs -f

# Specific service
docker-compose logs -f deltasignal
```


### Metrics

The application exposes metrics at:
- `/api/health` - Health status and basic info
- `/metrics` - Prometheus metrics (if enabled)

## Security

### HTTPS Setup

#### Docker Compose with Let's Encrypt
```yaml
# Add to docker-compose.yml
services:
  reverse-proxy:
    image: traefik:v2.10
    command:
      - --certificatesresolvers.letsencrypt.acme.email=your-email@domain.com
      - --certificatesresolvers.letsencrypt.acme.storage=/acme.json
      - --certificatesresolvers.letsencrypt.acme.httpchallenge.entrypoint=web
    volumes:
      - /var/run/docker.sock:/var/run/docker.sock
      - ./acme.json:/acme.json
```


### Security Headers

The application includes:
- Content Security Policy (CSP)
- X-Frame-Options: DENY
- X-Content-Type-Options: nosniff
- CORS configuration

## Performance Tuning

### Julia Performance
```bash
# Increase Julia threads
export JULIA_NUM_THREADS=8

# Use faster BLAS
export OPENBLAS_NUM_THREADS=4
```

### Frontend Performance
- Static assets are cached for 1 year
- Gzip compression enabled
- CSS/JS minification in production

### Database (if added)
```bash
# For large pathway databases, consider:
# - PostgreSQL with proper indexing
# - Redis for caching
# - Neo4j for graph queries
```

## Backup and Recovery

### Data Backup
```bash
# Backup examples and results
docker run --rm -v deltasignal_data:/data -v $(pwd):/backup alpine \
  tar czf /backup/deltasignal-data.tar.gz /data

# Restore
docker run --rm -v deltasignal_data:/data -v $(pwd):/backup alpine \
  tar xzf /backup/deltasignal-data.tar.gz -C /
```

### Configuration Backup
```bash
# Version control all configuration files:
# - docker-compose.yml
# - nginx.conf
# - package.json
```

## Troubleshooting

### Common Issues

#### Julia compilation errors
```bash
# Clear Julia cache
docker exec deltasignal rm -rf /root/.julia/compiled

# Rebuild image
docker-compose build --no-cache
```

#### Frontend build failures
```bash
# Check Node.js version
node --version  # Should be 18+

# Clear npm cache
npm cache clean --force

# Rebuild
cd frontend && npm ci && npm run build
```

#### Memory issues
```bash
# Increase Docker memory limit
# Docker Desktop -> Settings -> Resources -> Memory: 8GB

# Monitor usage
docker stats deltasignal
```

### Log Analysis

#### High CPU usage
```bash
# Check Julia threads
docker exec deltasignal julia -e "println(Threads.nthreads())"

# Profile computation
# Add profiling to Julia code
```

#### Network connectivity
```bash
# Test API from container
docker exec deltasignal curl http://localhost:8080/api/health

# Test external access
curl http://localhost/api/health
```

## Scaling

### Horizontal Scaling

#### Docker Compose
```bash
# Scale API instances
docker-compose up --scale deltasignal=3 -d

# Add load balancer
# Configure nginx upstream
```


### Vertical Scaling
```yaml
# Increase resources in docker-compose.yml
services:
  deltasignal:
    deploy:
      resources:
        limits:
          cpus: '4'
          memory: 8G
```

## Updates

### Rolling Updates

#### Docker Compose
```bash
# Build new version
VERSION=1.1.0 ./scripts/build.sh

# Update with zero downtime
docker-compose up -d --no-deps deltasignal
```


## Support

- **Documentation**: See README.md and inline code comments
- **Issues**: Check logs first, then create GitHub issue
- **Performance**: Monitor resource usage and scale appropriately

---

**Next Steps**: After deployment, see [QUICKSTART.md](QUICKSTART.md) for usage instructions.