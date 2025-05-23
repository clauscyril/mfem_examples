// Mesh.MshFileVersion = 2.2;
// Gmsh project created on Wed Mar 26 15:02:47 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {2};
Physical Surface(1100) = {1};
// Physical Line(1000) = {1};
Mesh.CharacteristicLengthMax = 0.06;