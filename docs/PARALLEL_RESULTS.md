# Parallel DCPSE — Validation and Performance Results

Results for `StreamVorti_par` on the 2D lid-driven cavity benchmark
(Re=100, dt=0.001, t=30 s to steady state) using the stream-function /
vorticity formulation with DCPSE derivatives.

---

## Validation: Ghia et al. (1982) benchmark

Comparison of the u-velocity centerline profile (x=0.5) against
Ghia et al. tabulated data at Re=100.

| Ranks | L2 error | L∞ error | L2 < 5% | L∞ < 10% | Result |
|-------|----------|----------|---------|----------|--------|
| np=1  | 0.30%    | 0.59%    | PASS    | PASS     | **PASS** |
| np=4  | 0.39%    | 0.80%    | PASS    | PASS     | **PASS** |

Grid: 40×40 nodes, 30 neighbours per node (DCPSE stencil).

The small increase from np=1 to np=4 is within expected floating-point
variation from the different partition-boundary stencils; both are well
within the 5% / 10% acceptance thresholds.

---

## Strong scaling: 40×40 grid

| Ranks | Time/iter | Speedup | Notes |
|-------|-----------|---------|-------|
| np=1  | 5.2 ms    | 1.00×   | Single rank, no MPI overhead |
| np=2  | 3.5 ms    | 1.49×   | 2 face neighbours per rank |
| np=4  | 3.3 ms    | 1.58×   | 2 face neighbours for corner partitions |

The np=4 gain over np=2 is modest because the 40×40 grid is small: each
partition at np=4 has ~400 nodes, so GMRES communication (one global
all-reduce per iteration, ~80 iterations/step) dominates over the local
compute. Larger grids show better scaling.

---

## Strong scaling: 80×80 grid

| Ranks | Time/iter | Speedup |
|-------|-----------|---------|
| np=1  | ~20 ms    | 1.00×   |
| np=2  | ~19 ms    | 1.04×   |
| np=4  | ~6.7 ms   | 2.98×   |

The np=2 anomaly (1.04× speedup) is caused by MFEM partitioning the 80×80
mesh into two long horizontal strips, maximising the shared-boundary length
relative to interior. The GMRES communication cost is proportionally higher
for strip partitions than for squarish (np=4) partitions.

---

## Environment

- Platform: macOS 15 (Darwin 25.3, x86_64 under Rosetta 2)
- CPU: Apple M-series (x86 emulation via Rosetta)
- MPI: Open MPI 5.0.8 (`/usr/local/opt/open-mpi`)
- MFEM: 4.8 with HYPRE and METIS
- Solver: HypreGMRES (~75–100 iterations, tolerance 1e-8)
- OpenMP: libomp 20.1.8 (Homebrew), used for vorticity update loop only

---

## Reproducing the results

```bash
# Build
mkdir build && cd build
cmake -DMFEM_DIR=/opt/mfem/mfem-4.8 -DCMAKE_BUILD_TYPE=Release ..
make -j6 StreamVorti_par

# Run 30-second simulation (np=4)
mpirun -np 4 ./StreamVorti_par -nx 40 -ny 40 -Re 100 -dt 0.001 -ft 30.0

# Validate against Ghia et al.
sbcl --load ../lisp/compare-ghia.lisp \
     --eval '(streamvorti.validation:run-validation
               :results-file "output_dat/mfem_square10x10_u_centerline_x0.5.dat"
               :reynolds 100)'
```
