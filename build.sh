#!/bin/bash
# Modern Python package build script using uv and pyproject.toml

set -e  # Exit on error

echo "🧹 Cleaning previous builds..."
rm -rf build dist *.egg-info

echo "📦 Building package with uv..."
uv build

echo "🚀 Uploading to PyPI..."
uv publish

echo "✅ Build and publish complete!"