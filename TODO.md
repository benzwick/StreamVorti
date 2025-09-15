# TODO

- restructure following similar structure as MFEM
  examples/poisson
  miniapps/dcpse
  src/efm
  src/general

- replace CGAL with nanoflann or similar
  - https://stackoverflow.com/questions/15124900/why-are-kd-trees-so-damn-slow-for-nearest-neighbor-search-in-point-sets
  - https://pointclouds.org/

- replace config manager with MFEM ArgParser (no more config files)

- rename library and replace headers



2025-06-12
- Use attributes for boundary conditons (Cartesian2D creates the automatically)
- save solution to paraview format
  - https://docs.mfem.org/html/classmfem_1_1ParaViewDataCollection.html
- Debug linear solver for streamfunction
