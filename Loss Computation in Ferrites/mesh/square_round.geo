// Paramètres
lc = 0.0003; // taille de maille
W = 0.0053;  // largeur
H = 0.0076;  // hauteur
r = 0.001;   // rayon des coins

// Points sur les bords droits entre les arrondis
Point(1) = {r, 0, 0, lc};
Point(2) = {W - r, 0, 0, lc};
Point(3) = {W, r, 0, lc};
Point(4) = {W, H - r, 0, lc};
Point(5) = {W - r, H, 0, lc};
Point(6) = {r, H, 0, lc};
Point(7) = {0, H - r, 0, lc};
Point(8) = {0, r, 0, lc};

// Centres des coins arrondis
Point(9)  = {r, r, 0, lc};
Point(10) = {W - r, r, 0, lc};
Point(11) = {W - r, H - r, 0, lc};
Point(12) = {r, H - r, 0, lc};

// Arcs de cercle (coins arrondis)
Circle(1) = {2, 10, 3}; // coin bas droit
Circle(2) = {4, 11, 5}; // coin haut droit
Circle(3) = {6, 12, 7}; // coin haut gauche
Circle(4) = {8, 9, 1};  // coin bas gauche

// Lignes droites entre les arcs
Line(5) = {3, 4}; // droite verticale droite
Line(6) = {5, 6}; // haut horizontal
Line(7) = {7, 8}; // gauche verticale
Line(8) = {1, 2}; // bas horizontal

// Contour fermé (ordre trigonométrique)
Line Loop(9) = {4, 8, 1, 5, 2, 6, 3, 7};

// Surface plane
Plane Surface(10) = {9};

// Maillage
Mesh 2;
