// Gmsh project created on Thu Mar 27 15:34:05 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 0.01, 0.0045, 0};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(2) = {2};
Mesh.CharacteristicLengthMax = 0.0005;
