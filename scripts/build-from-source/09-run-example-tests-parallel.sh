#!/bin/bash
set -e

STREAMVORTI_PAR_PATH=${1:-"build/StreamVorti_par"}
NP=${2:-2}

echo "Running parallel test with np=${NP}..."

mkdir -p test_streamvorti_par_np${NP}
cd test_streamvorti_par_np${NP}
echo "========================================="
echo "StreamVorti_par np=${NP} (5x5 mesh, 10 neighbors, 10 timesteps)"
echo "========================================="
mpirun --oversubscribe -np ${NP} ${STREAMVORTI_PAR_PATH} -dim 2 -nx 5 -ny 5 -nn 10 -ft 0.01 -dt 1e-3 -sd -sn -pv
echo ""
echo "StreamVorti_par np=${NP} outputs:"
ls -lh output_dat/*.dat 2>/dev/null || echo "No .dat files found"
ls -lh ParaView/ 2>/dev/null || echo "No ParaView directory found"
echo "DAT files generated: $(ls -1 output_dat/*.dat 2>/dev/null | wc -l)"
echo ""
echo "Parallel test np=${NP} completed successfully"
