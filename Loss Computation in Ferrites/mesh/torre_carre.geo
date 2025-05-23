// Gmsh project created on Wed May 21 18:27:37 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.0048, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 0.0101, 0, 2*Pi};
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Surface(1) = {1, 2};
//+
Curve Loop(3) = {2};
//+
Curve Loop(4) = {1};
//+
Plane Surface(1) = {3, 4};
//+
Extrude {0, 0, 0.0076} {
  Surface{1}; Curve{2}; Curve{1}; 
}
//+
Physical Surface("bdrc", 11) = {6, 4, 1, 5};
//+
Physical Volume("volume", 12) = {1};
