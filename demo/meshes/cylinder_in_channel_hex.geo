// Gmsh geometry file for cylinder-in-channel with quadrilateral elements
//
// Same geometry as cylinder_in_channel.geo but with structured/unstructured
// quad meshing for tensor-product space-time elements (hexahedra).
//
// Usage:
//   gmsh -2 cylinder_in_channel_hex.geo -o cylinder_in_channel_hex.msh

// Mesh size parameters
lc_far = 0.05;
lc_cylinder = 0.01;

// Channel dimensions
L = 2.2;
H = 0.41;

// Cylinder parameters
cx = 0.2;
cy = 0.2;
r = 0.05;

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
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder arcs
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};

// Physical groups
Physical Curve("inlet", 1) = {4};
Physical Curve("outlet", 2) = {2};
Physical Curve("top", 3) = {3};
Physical Curve("bottom", 4) = {1};
Physical Curve("cylinder", 5) = {5, 6, 7, 8};
Physical Surface("fluid", 1) = {1};

// Refinement fields
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_cylinder;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.5;

Background Field = 2;

// Recombine into quads
Mesh.Algorithm = 8;       // Frontal-Delaunay for Quads
Mesh.RecombineAll = 1;    // Force all quads
Mesh.Smoothing = 10;
Recombine Surface {1};
