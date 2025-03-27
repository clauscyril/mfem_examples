// Gmsh project created on Thu Mar 27 14:52:28 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.05, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
