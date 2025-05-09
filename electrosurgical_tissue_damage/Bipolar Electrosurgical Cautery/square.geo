SetFactory("OpenCASCADE");

// Crée le rectangle : (coin inférieur gauche) + (largeur, hauteur)
Rectangle(1) = {0, 0, 0, 0.01, 0.0045, 0}; // Tag 1

// Attributs physiques des lignes (bords)
Physical Line("Bas", 1) = {1};      // Bas → y = 0
Physical Line("Droite", 2) = {2};   // Droite → x = 10 mm
Physical Line("Haut", 3) = {3};     // Haut → y = 4.5 mm
Physical Line("Gauche", 4) = {4};   // Gauche → x = 0 mm

// Attribut physique de la surface (utile pour marquer le domaine dans MFEM)
Physical Surface("Domaine", 10) = {1};

// Paramètre de maillage
Mesh.CharacteristicLengthMax = 0.0003;