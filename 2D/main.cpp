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

    // cout << argc;

    if (argv[1]) {
        n = stoi(argv[1]);
    }
    cout << argv[2];

    if(argv[2]) {
        order = stoi(argv[2]);
    }
    // Création du maillage 
    Mesh mesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1, 1);

    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);




    // Affichage du maillage
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock.precision(8);
    // sol_sock_i.precision(8);
    sol_sock << "solution\n" << mesh << "window_title 'Mesh'";  // Génère une erreur dans glvis à corriger    

    return 0;
}