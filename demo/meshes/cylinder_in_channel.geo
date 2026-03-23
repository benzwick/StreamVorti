// Gmsh geometry file for cylinder-in-channel benchmark
//
// DFG benchmark (Schafer & Turek, 1996):
//   Channel: 2.2 x 0.41
//   Cylinder center: (0.2, 0.2), radius: 0.05
//
// Boundary attributes:
//   1 = inlet   (x = 0)
//   2 = outlet  (x = 2.2)
//   3 = top     (y = 0.41)
//   4 = bottom  (y = 0)
//   5 = cylinder
//
// Usage:
//   gmsh -2 cylinder_in_channel.geo -o cylinder_in_channel.msh
//   gmsh -2 cylinder_in_channel.geo -format msh2 -o cylinder_in_channel.msh
//
// For MFEM mesh format conversion:
//   gmsh -2 cylinder_in_channel.geo -format msh2 -o cylinder_in_channel.msh
//   Then use MFEM's mesh converter or load .msh directly

// Mesh size parameters
lc_far = 0.05;      // Far-field mesh size
lc_cylinder = 0.01; // Mesh size near cylinder
lc_wake = 0.02;     // Mesh size in the wake region

// Channel dimensions
L = 2.2;    // Length
H = 0.41;   // Height

// Cylinder parameters
cx = 0.2;   // Center x
cy = 0.2;   // Center y
r = 0.05;   // Radius

// Channel corners
Point(1) = {0, 0, 0, lc_far};
Point(2) = {L, 0, 0, lc_far};
Point(3) = {L, H, 0, lc_far};
Point(4) = {0, H, 0, lc_far};

// Cylinder center and surface points
Point(5) = {cx, cy, 0, lc_cylinder};
Point(6) = {cx + r, cy, 0, lc_cylinder};
Point(7) = {cx, cy + r, 0, lc_cylinder};
Point(8) = {cx - r, cy, 0, lc_cylinder};
Point(9) = {cx, cy - r, 0, lc_cylinder};

// Channel edges
Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Outlet
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Inlet

// Cylinder arcs
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// Channel boundary (counterclockwise)
Curve Loop(1) = {1, 2, 3, 4};

// Cylinder boundary (clockwise = hole)
Curve Loop(2) = {5, 6, 7, 8};

// Surface with hole
Plane Surface(1) = {1, 2};

// Physical groups (boundary attributes for MFEM)
Physical Curve("inlet", 1) = {4};
Physical Curve("outlet", 2) = {2};
Physical Curve("top", 3) = {3};
Physical Curve("bottom", 4) = {1};
Physical Curve("cylinder", 5) = {5, 6, 7, 8};

Physical Surface("fluid", 1) = {1};

// Mesh refinement near cylinder and in wake
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_cylinder;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.5;

// Wake refinement box
Field[3] = Box;
Field[3].VIn = lc_wake;
Field[3].VOut = lc_far;
Field[3].XMin = cx;
Field[3].XMax = L;
Field[3].YMin = cy - 3*r;
Field[3].YMax = cy + 3*r;

// Combine fields
Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

// Meshing options
Mesh.Algorithm = 6;       // Frontal-Delaunay
Mesh.RecombineAll = 0;    // 0 = triangles, 1 = quads
Mesh.Smoothing = 10;
