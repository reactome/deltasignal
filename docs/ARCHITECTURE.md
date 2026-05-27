# DeltaSignal Architecture

## Overview

DeltaSignal follows a **microservices architecture** with Docker-first deployment, implementing best practices for service separation and scalability.

## Service Architecture

```
                    ┌─────────────────┐
                    │   Nginx Proxy   │ :80
                    │ (Load Balancer) │
                    └─────┬───────────┘
                          │
              ┌───────────┼───────────┐
              │           │           │
              ▼           ▼           ▼
    ┌─────────────┐ ┌─────────────┐ ┌────────────┐
    │  Frontend   │ │  Julia API  │ │   Health   │
    │  (Nginx)    │ │  Backend    │ │  Monitor   │
    │    :3000    │ │    :8080    │ │    ---     │
    └─────────────┘ └─────────────┘ └────────────┘
```

## Services

### 1. Julia API Backend (`src/`)
- **Purpose**: Core computation engine and REST API
- **Technology**: Julia 1.10, HTTP.jl
- **Port**: 8080 (internal)
- **Health Check**: `/api/health`
- **Features**:
  - Pathway analysis algorithms (`src/core/`)
  - Network solving and optimization (`src/solvers/`)
  - Data processing and validation (`src/io/`)
  - CLI interface integration (`cli/`)
  - REST API server (`src/api/server.jl`)

### 2. Frontend Service (`frontend/`)
- **Purpose**: User interface and static asset serving
- **Technology**: Vite + Node.js, served by Nginx
- **Port**: 3000 (internal), 80 (via proxy)
- **Health Check**: `/health`
- **Features**:
  - Interactive network visualization
  - File upload and processing
  - Results visualization
  - Responsive design

### 3. Reverse Proxy (`proxy/`)
- **Purpose**: Load balancing, SSL termination, routing
- **Technology**: Nginx Alpine
- **Port**: 80 (public), 443 (HTTPS)
- **Features**:
  - Request routing to appropriate services
  - Static asset caching
  - Security headers
  - Rate limiting
  - WebSocket support for development

### 4. Health Monitor (Optional)
- **Purpose**: Service health aggregation
- **Technology**: Alpine Linux + shell scripts
- **Features**:
  - Periodic health checks
  - Service status reporting
  - Alerting capabilities

## Deployment Modes

### Development Environment
**File**: `docker-compose.dev.yml`

```yaml
services:
  julia-api:     # Direct port 8080
  frontend-dev:  # Vite dev server, port 3000-3001
  test-runner:   # On-demand testing
  cli:          # Interactive Julia shell
```

**Features**:
- Hot reload for all components
- Volume mounts for live code changes
- Development API server
- Direct service access

### Production Environment  
**File**: `docker-compose.prod.yml`

```yaml
services:
  api:          # Julia backend (internal)
  frontend:     # Static Nginx (internal)
  proxy:        # Reverse proxy (public :80)
  health:       # Monitoring (optional)
```

**Features**:
- Service isolation with internal networking
- Production-optimized builds
- Health checks and restart policies
- Security hardening

## Network Architecture

### Development Networks
- **Single network**: All services communicate directly
- **External access**: Direct port exposure (3000, 8080, 3001)

### Production Networks
- **Frontend Network**: Proxy ↔ Frontend communication
- **Backend Network**: Proxy ↔ API ↔ Health monitoring
- **Public Access**: Only through reverse proxy (:80, :443)

## Data Flow

### Request Flow (Production)
1. **Client** → `Nginx Proxy` (:80)
2. **Nginx Proxy** → Route decision:
   - `/api/*` → `Julia API Backend` (:8080)
   - `/*` → `Frontend Service` (:3000)
3. **Response** flows back through proxy with caching/compression

### Development Flow
1. **Client** → Direct service access
   - Frontend: `localhost:3000` (Vite dev server)
   - API: `localhost:8080` (Julia HTTP server)
   - Dev API: `localhost:3001` (Express server)

## Security Features

### Production Security
- **Non-root containers**: All services run as unprivileged users
- **Network isolation**: Services only communicate through defined networks
- **Security headers**: CSP, HSTS, X-Frame-Options, etc.
- **Rate limiting**: API and general request throttling
- **Health checks**: Automated service monitoring

### Development Security
- **Volume mounts**: Read-only where possible
- **Environment isolation**: Separate from production

## Scalability

### Horizontal Scaling
```bash
# Scale API backend
docker compose -f docker-compose.prod.yml up --scale api=3

# Scale with load balancer
docker compose -f docker-compose.prod.yml up --scale api=3 --scale proxy=2
```

### Vertical Scaling
```yaml
services:
  api:
    deploy:
      resources:
        limits:
          cpus: '4.0'
          memory: 8G
```

## Monitoring & Health

### Health Endpoints
- **API**: `GET /api/health` → Julia service status
- **Frontend**: `GET /health` → Nginx service status  
- **Proxy**: `GET /health` → Load balancer status

### Logging
- **Development**: `docker compose logs -f`
- **Production**: Centralized logging to volumes
- **Monitoring**: Optional health-monitor service

## File Structure

```
/
├── src/                    # Julia source code (core)
│   ├── DeltaSignal.jl     # Main module
│   ├── core/              # Mathematical operations
│   ├── io/                # Data I/O
│   ├── solvers/           # Numerical solvers
│   ├── learning/          # Parameter learning
│   └── api/               # HTTP API server
├── cli/                   # Command-line tools
│   └── deltasignal.jl    # Main CLI script
├── test/                  # Julia test suite
├── examples/              # Example data files
├── frontend/              # Frontend service
│   ├── Dockerfile.prod    # Production frontend container
│   ├── Dockerfile.dev     # Development container
│   ├── nginx.conf         # Frontend nginx config
│   └── src/               # Frontend source
├── docker/                # Docker configuration
│   ├── nginx.conf         # Reverse proxy config
│   ├── supervisord.conf   # Process manager
│   └── start.sh           # Container startup
├── proxy/                 # Nginx reverse proxy
│   ├── nginx.conf         # Main proxy config
│   └── conf.d/            # Virtual host configs
├── scripts/               # Deployment automation
│   ├── build.sh           # Multi-service builds
│   └── deploy.sh          # Environment deployment
├── docs/                  # Documentation
├── docker-compose.dev.yml # Development environment
├── docker-compose.prod.yml# Production environment
├── Dockerfile             # Multi-stage production build
└── Project.toml           # Julia dependencies
```

## Best Practices Implemented

1. **Service Separation**: One process per container
2. **12-Factor App**: Environment config, stateless services
3. **Security**: Non-root users, network isolation
4. **Observability**: Health checks, structured logging
5. **Scalability**: Horizontal and vertical scaling support
6. **Development Experience**: Hot reload, volume mounts
7. **Production Ready**: Reverse proxy, caching, compression
8. **DevOps**: Automated builds, deployment scripts

## Migration from Legacy

The legacy monolithic `Dockerfile` has been replaced with:
- ✅ **Separate service Dockerfiles**
- ✅ **Multi-environment compose files**
- ✅ **Reverse proxy architecture**  
- ✅ **Updated deployment scripts**

This provides better separation of concerns, improved scalability, and easier maintenance.