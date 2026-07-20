#!/bin/bash
set -e

# DeltaSignal Build Script — engine + HTTP API only.
# The UI lives in the WebsiteAngular workspace (see docs/API.md).
# Usage: ./scripts/build.sh [api|dev|prod] [--push]

SERVICE=${1:-dev}
PUSH_IMAGES=${2:-false}
VERSION=${VERSION:-latest}
REGISTRY=${REGISTRY:-deltasignal}

echo "🏗️ DeltaSignal Build Script"
echo "Service: $SERVICE   Version: $VERSION   Registry: $REGISTRY   Push: $PUSH_IMAGES"
echo "================================"

case $SERVICE in
  "api")
    image_name="$REGISTRY/deltasignal-api:$VERSION"
    echo "Building API image ($image_name)..."
    docker build -t "$image_name" --target julia-builder .
    if [ "$PUSH_IMAGES" = "--push" ]; then
      echo "Pushing $image_name..."
      docker push "$image_name"
    fi
    echo "✅ API image built"
    ;;

  "dev")
    echo "Building development environment..."
    docker compose -f docker-compose.dev.yml build
    echo "✅ Development environment built"
    ;;

  "prod"|"production")
    echo "Building production stack..."
    docker compose -f docker-compose.prod.yml build
    echo "✅ Production stack built"
    ;;

  *)
    echo "❌ Unknown service: $SERVICE"
    echo "Usage: $0 [api|dev|prod] [--push]"
    exit 1
    ;;
esac

echo ""
echo "🎉 Build completed successfully!"
