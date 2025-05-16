// Tore plein avec section carrée
SetFactory("OpenCASCADE");

// Paramètres
R = 0.5;         // Rayon du cercle directeur (rayon du grand cercle)
a = 0.2;         // Demi-côté du carré (donc section = 2a x 2a)
angle = 2*Pi;    // Balayage sur 360° (tore complet)

// Création du carré dans le plan YZ (centré en (R, 0, 0))
Point(1) = {R, -a, -a, 1.0};
Point(2) = {R,  a, -a, 1.0};
Point(3) = {R,  a,  a, 1.0};
Point(4) = {R, -a,  a, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Rotation du carré autour de l'axe Z pour générer le tore
// Résultat : tore plein avec section carrée
Surface Loop(100) = {1};
Extrude {0, 0, 0}, {0, 0, 0}, {0, 0, 1}, angle {
    Surface{1}; Layers{40}; Recombine;
}

// Marquage physique
Physical Volume("ToreCarre", 100) = {1};
