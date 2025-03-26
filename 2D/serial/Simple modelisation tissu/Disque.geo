// Gmsh project created on Wed Mar 26 09:28:49 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.01, 0, 2*Pi};
//+
Curve Loop(1) = {1};
