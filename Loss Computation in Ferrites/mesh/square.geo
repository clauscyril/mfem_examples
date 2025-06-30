SetFactory("OpenCASCADE");

// Crée le rectangle : (coin inférieur gauche) + (largeur, hauteur)
Rectangle(1) = {0.0048, 0, 0, 0.0053 , 0.0076, 0}; // Tag 1

// Attributs physiques des lignes (bords)
Physical Line("Bas", 1) = {1};    
Physical Line("Droite", 2) = {2};   
Physical Line("Haut", 3) = {3};     
Physical Line("Gauche", 4) = {4};   

Physical Surface("Domaine", 10) = {1};

 Mesh.ElementOrder = 2;
// Paramètre de maillage
Mesh.CharacteristicLengthMax = 1;