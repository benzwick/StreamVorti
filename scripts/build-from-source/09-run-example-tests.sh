#!/bin/bash
set -e

MFEMRUN_PATH=${1:-"build/MfemRun"}
STREAMVORTI_PATH=${2:-"build/StreamVorti"}

echo "Running test examples..."

# Test MfemRun in its own directory
mkdir -p test_mfemrun
cd test_mfemrun
echo "========================================="
echo "Test 1: MfemRun (5x5 mesh, 10 neighbors)"
echo "========================================="
${MFEMRUN_PATH} -dim 2 -nx 5 -ny 5 -nn 10 -sd -sn
echo ""
echo "MfemRun outputs:"
ls -lh *.dat 2>/dev/null || echo "No .dat files found"
echo "Files generated: $(ls -1 *.dat 2>/dev/null | wc -l)"

# Test StreamVorti in its own directory
cd ..
mkdir -p test_streamvorti
cd test_streamvorti
echo ""
echo "========================================="
echo "Test 2: StreamVorti (5x5 mesh, 10 neighbors, 10 timesteps)"
echo "========================================="
${STREAMVORTI_PATH} -dim 2 -nx 5 -ny 5 -nn 10 -tf 0.01 -dt 1e-3 -sd -sn -pv
echo ""
echo "StreamVorti outputs:"
ls -lh output_dat/*.dat 2>/dev/null || echo "No .dat files found"
ls -lh ParaView/ 2>/dev/null || echo "No ParaView directory found"
echo "DAT files generated: $(ls -1 output_dat/*.dat 2>/dev/null | wc -l)"
echo "ParaView files generated: $(find ParaView -type f 2>/dev/null | wc -l)"

echo ""
echo "All test examples completed successfully"
