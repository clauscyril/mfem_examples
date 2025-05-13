// Gmsh project: Parallelepiped with Physical Surfaces
SetFactory("OpenCASCADE");

// Dimensions
Lx = 0.01;
Ly = 0.0045;
Lz = 0.01;

// Mesh size
lc = 0.0003;
Mesh.CharacteristicLengthMax = lc;

// Base rectangle in XY plane
Rectangle(1) = {0, 0, 0, Lx, Ly, 0};

// Extrude along Z to form the 3D volume
// Save all created entities
out[] = Extrude {0, 0, Lz} {
  Surface{1}; Layers{15}; Recombine;
};

// Physical Volume (entire 3D region)
Physical Volume("Volume") = {out[1]};

// The extrusion creates surfaces in the following order:
// out[0] : top surface (extruded)
// out[2] to out[5] : lateral sides
// Surface{1} : bottom surface

// Assign Physical Surfaces (boundary markers) for BCs
Physical Surface("Bottom") = {1};         // z = 0
Physical Surface("Top") = {out[0]};       // z = Lz
Physical Surface("Front") = {out[2]};     // y = 0
Physical Surface("Back") = {out[3]};      // y = Ly
Physical Surface("Left") = {out[4]};      // x = 0
Physical Surface("Right") = {out[5]};     // x = Lx
