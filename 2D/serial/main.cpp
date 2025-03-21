#include <mfem.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;

static real_t e = 0.001f;
real_t freq = 1.0, kappa;
void E_exact(const Vector &x, Vector &E);
void f_exact(const Vector &x, Vector &f);


int main(int argc, char* argv[])
{   
    // Il est nécessaire de passer deux arguments lors de l'execution du programme
    int n;       // Nombre d'éléments du maillage
    int order;    // Ordre des fonctions de test  

    if (argv[1]) {
        n = stoi(argv[1]);
    }

    if(argv[2]) {
        order = stoi(argv[2]);
    }


    kappa = freq * M_PI;

    // Création du maillage (Maillage carré de nxn éléments, [0,1]x[0,1])
    Mesh mesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1, 1);
    // Mesh mesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1, 1);

    mesh.UniformRefinement();
    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;

    cout << "Dimension : " << dim << "\nSpaceDimension : " << spaceDim << endl;


    // Définition des espace d'élements finis
    FiniteElementCollection *fec = new ND_FECollection(order, dim);     // On utilise les éléments finis de Nedelec
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;


    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;

    Array<int> ess_tdof_list;
    if (mesh.bdr_attributes.Size())
    {
        Array<int> ess_bdr(mesh.bdr_attributes.Max());
        ess_bdr = 1; // On impose les conditions comme des conditions de dirichlet
        cout << "ess_bdr : " << ess_bdr[0] << endl;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    // Définition de la forme linéaire B
    VectorFunctionCoefficient f(spaceDim, f_exact);
    LinearForm *b = new LinearForm(fespace);
    b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
    b->Assemble();


    // Définition du vecteur solution comme un grid function liée à l'espace d'éléments finis
    GridFunction x(fespace);
    VectorFunctionCoefficient E(spaceDim, E_exact);
    x.ProjectCoefficient(E);

    // Forme bilinéaire
    Coefficient *muinv = new ConstantCoefficient(1.f);
    Coefficient *sigma = new ConstantCoefficient(1.f);
    BilinearForm *a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
    a->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));
    a->Assemble();


    // Définition du système linéaire à résoudre
    OperatorPtr A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

    cout << "Taille du système linéaire : " << A->Height() << endl;

    // Solver
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 1000, 1e-12, 0.f);

    // On récupère les solutions 
    a->RecoverFEMSolution(X, *b, x);


    // On sauvegarde le mesh 
    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);

    cout << endl;

    // Affichage de la solution avec glvis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock.precision(8);
    // sol_sock_i.precision(8);
    sol_sock << "solution\n" << mesh << x << "window_title 'Solution'"
             << flush;  // Génère une erreur dans glvis à corriger    

    return 0;
}



void E_exact(const Vector &x, Vector &E)
{
    E(0) = sin(kappa * x(1));
    E(1) = sin(kappa * x(0));
    if (x.Size() == 3) { E(2) = 0.0; }
}

void f_exact(const Vector &x, Vector &f)
{
    f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
    f(1) = (1. + kappa * kappa) * sin(kappa * x(0));
    if (x.Size() == 3) { f(2) = 0.0; }  
}