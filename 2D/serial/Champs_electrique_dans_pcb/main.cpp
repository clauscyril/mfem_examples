#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

double eps_function(const Vector &x) { 
    if (x(1) > 0.001) {
        return 8.85e-12;
    } else {
        return 10 * 8.85e-12;
    }  
}


int main(int argc,char* argv[]) {
    int order;
    if (argv[1]) {
        order = stoi(argv[1]);
    }
    const char *path = "D:/Documents/projets/MFEM/mfem_examples/2D/serial/Champs_electrique_dans_pcb/floating.msh";
    Mesh mesh(path, 1, 1);

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;


    FiniteElementCollection *fec = new H1_FECollection(order, dim);     // On utilise les éléments finis de 
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;

    // Tableau des vrais dof (pour les conditions de dirichlet)
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh.GetNV(); i++) {
        const double *v = mesh.GetVertex(i);
        if (abs(v[1]) < 1e-6 or (v[0] <= 0.00236 && (v[1] >= 0.001 && v[1] <= 0.00135))) {  // Vérifie si y = 0
            ess_tdof_list.Append(i);
        } 
        
    }
    cout << "Nombre de nœuds sur y = 0 : " << ess_tdof_list.Size() << endl;
    cout << fespace->GetTrueVSize() << endl;
    // mesh.GetBdrElementVertices(0, );
    GridFunction v(fespace);
    v = 1.f;
    for (int i = 0; i < ess_tdof_list.Size(); i++){
        const double *u = mesh.GetVertex(ess_tdof_list[i]);
        if (u[1] > 0){
            v(ess_tdof_list[i]) = 0.001f;
        } else {
            v(ess_tdof_list[i]) = 0.f;
        }
    }

    cout << "Nombre de conditions de Dirichlet : " << ess_tdof_list.Size() << endl;

    LinearForm *b = new LinearForm(fespace);
    *b = 0.f;
    // b->AddDomainIntegrator(new VectorFEDomainLFIntegrator());
    b->Assemble();


    // double eps = 1.0;
    FunctionCoefficient eps(eps_function);
    BilinearForm *a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(eps));
    a->Assemble();

    OperatorPtr A;
    Vector B, X;

    a->FormLinearSystem(ess_tdof_list, v, *b, A, X, B);
    cout << "Taille du système linéaire : " << A->Height() << endl;

    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 1000, 1e-12, 0.f);

    // On récupère les solutions 
    a->RecoverFEMSolution(X, *b, v);

    // Affichage de la solution avec glvis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock.precision(8);
    // sol_sock_i.precision(8);
    sol_sock << "solution\n" << mesh << v << "window_title 'Solution'"
            << flush;  // Génère une erreur dans glvis à corriger 

    return 0;
}