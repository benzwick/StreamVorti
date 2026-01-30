#!/bin/bash
set -e

echo "Building and installing Google Test..."

cd /usr/src/gtest
sudo cmake .
sudo make
sudo cp lib/*.a /usr/lib/ 2>/dev/null || sudo cp *.a /usr/lib/ 2>/dev/null || true

echo "Google Test installed successfully"
