// Tore plein (volume)
// Gmsh project modified on 2025-05-15
SetFactory("OpenCASCADE");

// Paramètres du tore
r1 = 0.5;     // Grand rayon (du centre du trou au centre du tube)
r2 = 0.2;     // Petit rayon (rayon du tube)
angle = 2*Pi; // Angle du tore (2*Pi pour un tore complet)

// Création du tore plein (volume)
Torus(1) = {0, 0, 0, r1, r2, angle};  // Crée un tore volumique OpenCASCADE

// Marquage physique des surfaces et du volume
Physical Surface("Dirichlet", 101) = {1}; // Optionnel : marquer surface pour conditions limites
Physical Volume("Volume", 100) = {1};     // Marque le volume pour MFEM ou autres solveurs
