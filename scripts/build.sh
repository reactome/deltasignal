#!/bin/bash
set -e

# DeltaSignal Multi-Service Build Script
# Usage: ./scripts/build.sh [service] [--push]

SERVICE=${1:-all}
PUSH_IMAGES=${2:-false}
VERSION=${VERSION:-latest}
REGISTRY=${REGISTRY:-deltasignal}

echo "🏗️ DeltaSignal Build Script"
echo "Service: $SERVICE"
echo "Version: $VERSION"
echo "Registry: $REGISTRY"
echo "Push: $PUSH_IMAGES"
echo "================================"

build_service() {
    local service=$1
    local context=$2
    local dockerfile=$3
    local image_name="$REGISTRY/deltasignal-$service:$VERSION"
    
    echo "Building $service service..."
    docker build -t "$image_name" -f "$dockerfile" "$context"
    
    if [ "$PUSH_IMAGES" = "--push" ]; then
        echo "Pushing $image_name..."
        docker push "$image_name"
    fi
    
    echo "✅ $service service built successfully"
}

case $SERVICE in
    "backend"|"api")
        build_service "api" "./backend" "./backend/Dockerfile"
        ;;
        
    "frontend"|"web")
        build_service "frontend" "./frontend" "./frontend/Dockerfile.prod"
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
        
    "all")
        echo "Building all services..."
        build_service "api" "./backend" "./backend/Dockerfile"
        build_service "frontend" "./frontend" "./frontend/Dockerfile.prod"
        echo "✅ All services built successfully"
        ;;
        
    *)
        echo "❌ Unknown service: $SERVICE"
        echo "Usage: $0 [backend|frontend|dev|prod|all] [--push]"
        echo ""
        echo "Examples:"
        echo "  ./scripts/build.sh all                    # Build all services"
        echo "  ./scripts/build.sh backend --push         # Build and push backend"
        echo "  ./scripts/build.sh dev                    # Build dev environment"
        echo "  VERSION=1.2.0 ./scripts/build.sh prod     # Build prod with version"
        exit 1
        ;;
esac

echo ""
echo "🎉 Build completed successfully!"

if [ "$PUSH_IMAGES" = "--push" ]; then
    echo "📦 Images pushed to registry"
fi