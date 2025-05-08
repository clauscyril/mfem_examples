// Gmsh project created on Wed May 07 17:05:47 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.016, -Pi/2, Pi/2, 2*Pi};
Mesh.CharacteristicLengthMax = 0.0009; // Plus petit = raffinement
Mesh.CharacteristicLengthMin = 0.0009;