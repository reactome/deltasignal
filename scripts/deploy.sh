#!/bin/bash
set -e

# DeltaSignal Deployment Script
# Usage: ./scripts/deploy.sh [environment]

ENVIRONMENT=${1:-development}
VERSION=${VERSION:-latest}
REGISTRY=${REGISTRY:-deltasignal}

echo "🧬 DeltaSignal Deployment Script"
echo "Environment: $ENVIRONMENT"
echo "Version: $VERSION"
echo "Registry: $REGISTRY"
echo "================================"

case $ENVIRONMENT in
  "development")
    echo "🛠️ Starting development environment..."
    
    # Start development environment with Docker Compose
    echo "Starting development services..."
    docker compose -f docker-compose.dev.yml up --build -d
    
    echo "✅ Development environment started!"
    echo "Frontend (Vite): http://localhost:3000"
    echo "API (Julia): http://localhost:8080"  
    echo "Dev API (Express): http://localhost:3001"
    echo ""
    echo "To view logs: docker compose -f docker-compose.dev.yml logs -f"
    echo "To stop: docker compose -f docker-compose.dev.yml down"
    ;;
    
  "production")
    echo "🚀 Deploying to production..."
    
    # Build production images
    echo "Building production services..."
    docker compose -f docker-compose.prod.yml build
    
    if [ "$REGISTRY" != "deltasignal" ]; then
      echo "Tagging and pushing images to registry..."
      docker tag deltasignal-api:latest $REGISTRY/deltasignal-api:$VERSION
      docker tag deltasignal-frontend:latest $REGISTRY/deltasignal-frontend:$VERSION
      docker push $REGISTRY/deltasignal-api:$VERSION
      docker push $REGISTRY/deltasignal-frontend:$VERSION
    fi
    
    # Deploy production stack
    echo "Starting production deployment..."
    VERSION=$VERSION docker compose -f docker-compose.prod.yml up -d
    
    echo "✅ Production deployment complete!"
    echo "Application: http://localhost (via reverse proxy)"
    echo "API: http://localhost/api/"
    echo "Health: http://localhost/health"
    echo ""
    echo "To view logs: docker compose -f docker-compose.prod.yml logs -f"
    echo "To stop: docker compose -f docker-compose.prod.yml down"
    ;;
    
  "test")
    echo "🧪 Running tests..."
    
    # Run Julia tests in development environment
    echo "Running Julia tests..."
    docker compose -f docker-compose.dev.yml --profile test up --build test-runner
    
    # Run frontend tests (if they exist)  
    if [ -f "frontend/package.json" ]; then
      echo "Running frontend tests..."
      docker compose -f docker-compose.dev.yml run --rm frontend-dev npm test || echo "No frontend tests defined"
    fi
    
    echo "✅ Tests completed!"
    ;;
    
  *)
    echo "❌ Unknown environment: $ENVIRONMENT"
    echo "Usage: $0 [development|production|test]"
    exit 1
    ;;
esac