// Dimensions du pavé
Lx = 10;
Ly = 10;
Lz = 10;

// Discrétisation (nombre de divisions)
nx = 10;
ny = 10;
nz = 10;

// Points du cube
Point(1) = {-Lx/2, -Ly/2, -Lz/2};
Point(2) = {Lx/2, -Ly/2, -Lz/2};
Point(3) = {Lx/2, Ly/2, -Lz/2};
Point(4) = {-Lx/2, Ly/2, -Lz/2};

Point(5) = {-Lx/2, -Ly/2, Lz/2};
Point(6) = {Lx/2, -Ly/2, Lz/2};
Point(7) = {Lx/2, Ly/2, Lz/2};
Point(8) = {-Lx/2, Ly/2, Lz/2};

// Lignes de la base
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Lignes du haut
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Lignes verticales
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Surfaces latérales
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};  // base z=0

Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};  // haut z=Lz

Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};

Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};

Line Loop(21) = {3, 12, -7, -11};
Plane Surface(22) = {21};

Line Loop(23) = {4, 9, -8, -12};
Plane Surface(24) = {23};

// Surface loop et volume
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

// Maillage structuré
Transfinite Line {1,2,3,4} = ny + 1 Using Progression 1;
Transfinite Line {5,6,7,8} = ny + 1 Using Progression 1;
Transfinite Line {9,10,11,12} = nz + 1 Using Progression 1;

Transfinite Surface {14,16,18,20,22,24};


Recombine Volume {26}; // Maillage en hexaèdres
//+
Physical Surface("Dirichlet", 27) = {24, 22, 20, 18};
//+
Show "*";
//+
Volume(27) = {25};
