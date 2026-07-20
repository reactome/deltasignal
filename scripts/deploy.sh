#!/bin/bash
set -e

# DeltaSignal Deployment Script — engine + HTTP API only.
# The UI lives in the WebsiteAngular workspace (see docs/API.md).
# Usage: ./scripts/deploy.sh [development|production|test]

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
    echo "🛠️ Starting development API..."
    docker compose -f docker-compose.dev.yml up --build -d julia-api

    echo "✅ API started!"
    echo "API (Julia): http://localhost:8080  (health: /api/health)"
    echo ""
    echo "To view logs: docker compose -f docker-compose.dev.yml logs -f julia-api"
    echo "To stop: docker compose -f docker-compose.dev.yml down"
    ;;

  "production")
    echo "🚀 Deploying API to production..."
    docker compose -f docker-compose.prod.yml build

    if [ "$REGISTRY" != "deltasignal" ]; then
      echo "Tagging and pushing image to registry..."
      docker tag deltasignal-api:latest $REGISTRY/deltasignal-api:$VERSION
      docker push $REGISTRY/deltasignal-api:$VERSION
    fi

    echo "Starting production deployment..."
    VERSION=$VERSION docker compose -f docker-compose.prod.yml up -d

    echo "✅ Production deployment complete!"
    echo "API: http://localhost:8080/api/  (health: /api/health)"
    echo ""
    echo "To view logs: docker compose -f docker-compose.prod.yml logs -f"
    echo "To stop: docker compose -f docker-compose.prod.yml down"
    ;;

  "test")
    echo "🧪 Running Julia tests..."
    docker compose -f docker-compose.dev.yml --profile test up --build test-runner
    echo "✅ Tests completed!"
    ;;

  *)
    echo "❌ Unknown environment: $ENVIRONMENT"
    echo "Usage: $0 [development|production|test]"
    exit 1
    ;;
esac
