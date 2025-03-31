// Script Gmsh pour générer un maillage quadrilatéral homogène d'un disque

// Paramètres du disque
radius = 1.0;
num_radial = 10;
num_circum = 20;

// Création du disque
Point(1) = {0, 0, 0, 1.0};
Point(2) = {radius, 0, 0, 1.0};
Point(3) = {0, radius, 0, 1.0};
Point(4) = {-radius, 0, 0, 1.0};
Point(5) = {0, -radius, 0, 1.0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Maillage
Transfinite Surface {1};
Transfinite Line {1, 2, 3, 4} = num_circum Using Progression 1;
Transfinite Line {1} = num_radial Using Progression 1;

Recombine Surface {1};

// Ajustement de la taille des éléments pour plus d'homogénéité
Mesh.CharacteristicLengthFactor = 1.0;
Mesh.CharacteristicLengthMin = radius / num_radial;
Mesh.CharacteristicLengthMax = radius / num_radial;

// Génération du maillage
Mesh 2;
