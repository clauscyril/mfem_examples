#include <mfem.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;

static real_t e = 0.001f;


int main(int argc, char* argv[])
{   
    // Il faut passer deux arguments lors de l'execution du programme
    int n  = 4;       // Nombre d'éléments du maillage (4 par défaut)
    int order = 1;    // Ordre des fonctions de test (1 par défaut)  

    cout << argc;

    if (argv[1]) {
        n = stoi(argv[1]);
    }
    cout << argv[2];

    if(argv[2]) {
        order = stoi(argv[2]);
    }
    // Création du maillage 
    Mesh mesh =  Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1, 1);
    return 0;
}