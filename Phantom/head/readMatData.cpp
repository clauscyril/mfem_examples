#include <iostream>
#include "mat.h"   // Permet de lire les fichiers .mat (Données MATLAB)

using namespace std;

int main(){

    // Ouvrir le fichier .mat
    const char *file = "D:/Documents/projets/Alvar/Alvar_v16.mat";
    MATFile *pmat = matOpen(file, "r");
    
    if (pmat == nullptr) {
        cout << "Erreur lors de l'ouverture du fichier .mat" << endl;
        return 1;
    }

    // Obtenir la liste des variables dans le fichier
    int n = 0;
    char **dir = matGetDir(pmat, &n);
    if (dir == nullptr) {
        std::cerr << "Erreur lors de la récupération de la liste des variables" << std::endl;
        matClose(pmat);
        return 1;
    }
    // Afficher les noms des variables
    std::cout << "Liste des variables dans le fichier .mat :" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << dir[i] << std::endl;
        // mxFree(dir[i]); // Libérer la mémoire allouée pour chaque nom de variable
    }
    mxFree(dir); // Libérer la mémoire allouée pour le tableau de noms de variables

    // Fermer le fichier .mat
    matClose(pmat);

    return 0;
}
