// cube.geo - Cube 1x1x1 démarrant à (0,0,0)

lc = 0.1; // taille caractéristique (maillage)

// Définition des points
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {1, 1, 1, lc};
Point(8) = {0, 1, 1, lc};

// Arêtes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Surfaces
Line Loop(13) = {1, 2, 3, 4};       // Face inférieure
Line Loop(14) = {5, 6, 7, 8};       // Face supérieure
Line Loop(15) = {1, 10, -5, -9};    // Face avant
Line Loop(16) = {2, 11, -6, -10};   // Face droite
Line Loop(17) = {3, 12, -7, -11};   // Face arrière
Line Loop(18) = {4, 9, -8, -12};    // Face gauche

Plane Surface(19) = {13};
Plane Surface(20) = {14};
Plane Surface(21) = {15};
Plane Surface(22) = {16};
Plane Surface(23) = {17};
Plane Surface(24) = {18};

// Volume
Surface Loop(25) = {19, 20, 21, 22, 23, 24};
Volume(26) = {25};

// Maillage 3D
Mesh 3;
//+
Physical Surface("dir", 27) = {24, 21, 22, 23};
//+
Physical Surface("libre", 28) = {20, 19};
//+
Physical Volume("vol", 29) = {26};
