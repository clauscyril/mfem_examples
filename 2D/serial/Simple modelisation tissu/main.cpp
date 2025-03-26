#include <mfem.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

static real_t e = 0.001f;
real_t freq = 1.0, kappa;
void A_func(const Vector &x, Vector &E);



int main(int argc, char* argv[])
{   
    int n;       // Nombre d'éléments du maillage
    int order = 1;    // Ordre des fonctions de test  

    const char *path = "../disque.msh";
    Mesh mesh(path, 1, 1);

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;

    mesh.UniformRefinement();
    
    FiniteElementCollection *fec = new H1_FECollection(order, dim);     // On utilise les éléments finis noeudaux
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;

    mesh.PrintInfo(cout);


    // Conditions de type Dirichlet
    // Array<int> ess_tdof_list;

    cout << mesh.GetNBE();

    ComplexGridFunction v(fespace);
    v = 0.f;


    VectorFunctionCoefficient sigmajwA(dim + 1, A_func);

    ComplexLinearForm *b = new ComplexLinearForm(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    b->AddDomainIntegrator(NULL, new DomainLFGradIntegrator(sigmajwA)); // Forme linéaire imaginaire pur 
    b->Assemble();


    ConstantCoefficient sigma(1.f);
    SesquilinearForm*a = new SesquilinearForm(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    a->AddDomainIntegrator(new DiffusionIntegrator(sigma), NULL);       // Forme Bilinéaire réel
    a->Assemble();

    
    

    // Gestion des conditions aux limites
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0; // Imposer u = 1 sur les bords
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    ConstantCoefficient ZeroCoef(0.f);
    

    v.ProjectBdrCoefficient(ZeroCoef, ZeroCoef, ess_bdr);  


    OperatorHandle A;
    Vector B, U;

    a->FormLinearSystem(ess_tdof_list, v, *b, A, U, B);
    cout << A->Width() << endl << endl;

    cout << "test1" << endl;
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, U, 1, 1000, 1e-12, 0.f);
    cout << "test2" << endl;

    
    // On récupère les solutions 
    a->RecoverFEMSolution(U, *b, v);
    cout << "test3" << endl;
        // Visualisation GLvis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_r(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock_r.precision(8);
    // sol_sock_i.precision(8);
    sol_sock_r << "solution\n" << mesh << v.real()
                << "window_title 'Solution: Real Part'" 
                << "pause\n" << "keys c\n" << flush;
    
    cout << "test" << endl;

    return 0;
}


void A_func(const Vector &x, Vector &A){
    double B = 1;
    double sigma_ = 1;
    double omega = M_PI*2* 50;
    A(0) = 0;
    A(1) = sigma_*omega*B*x(3);
    A(2) = 0;  
}
