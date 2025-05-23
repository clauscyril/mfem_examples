#include <mfem.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

void A_func(const Vector &x, Vector &A);
double A_func_bdr(const Vector &x);

int main(int argc, char* argv[])
{   
    int n;       
    int order = 1;  

    const char *path = "../disque.msh";
    Mesh mesh(path, 1, 1);

    mesh.UniformRefinement();
    // mesh.UniformRefinement();
    

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;

    FiniteElementCollection *fec = new H1_FECollection(order, dim);    
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;
    mesh.PrintInfo(cout);

    // Gestion des conditions aux limites
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh.GetNV(); i++) {
        const double *v = mesh.GetVertex(i);
        if (abs(v[0]) > 0.4 ) {  // 
            ess_tdof_list.Append(i);
        }
        
    }

    GridFunction v(fespace);
    v = 0.f;
    for (int i = 0; i < ess_tdof_list.Size(); i++){
        const double *u = mesh.GetVertex(ess_tdof_list[i]);
        if (u[0] > 0) {
            v(ess_tdof_list[i]) = 1.f;
        } else {
            v(ess_tdof_list[i]) = -1.f;
        }
    }


    cout << mesh.bdr_attributes.Size() << endl << endl;

    VectorFunctionCoefficient sigmajwA(dim + 1, A_func);
    FunctionCoefficient omegaA(A_func_bdr);


    Array<int> surf_attrib;
    surf_attrib.Append(1);  // Numéro de la surface où la charge doit être nulle

    LinearForm *b = new LinearForm(fespace);
    b->Assemble();

    ConstantCoefficient sigma(0.2f);

    // Forme bilinéaire réelle principale
    BilinearForm *a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(sigma)); 
    // a->AddDomainIntegrator(new MassIntegrator(m),NULL);
    a->Assemble();

    OperatorPtr A;
    Vector B, X;

    a->FormLinearSystem(ess_tdof_list, v, *b, A, X, B);
    cout << "Taille du système linéaire : " << A->Height() << endl;

    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 1000, 1e-12, 0.f);

    // On récupère les solutions 
    a->RecoverFEMSolution(X, *b, v);



    ofstream sol_r_ofs("sol_real.gf");

    sol_r_ofs.precision(8);

    v.Save(sol_r_ofs);

    GradientGridFunctionCoefficient grad_v_coeff(&v);

    FiniteElementCollection *fec_grad = new ND_FECollection(order, dim );  // Raviart-Thomas (RT) pour le gradient
    FiniteElementSpace *fespace_grad = new FiniteElementSpace(&mesh, fec_grad);
    
    GridFunction grad_v(fespace_grad);

    grad_v.ProjectCoefficient(grad_v_coeff);

    GridFunction J(fespace_grad);
    VectorFunctionCoefficient A_coeff(dim, A_func);
    J.ProjectCoefficient(A_coeff);
    
    J += grad_v;


    cout << "test" << endl;
   // Visualisation GLVis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_i(vishost, visport);
    socketstream sol_sock_J(vishost, visport);
    sol_sock_i.precision(8);
    sol_sock_J.precision(8);
    sol_sock_i << "solution\n" << mesh << v
                << "window_title 'Solution: Potentiel'" 
                << "pause\n" << "keys c\n" << flush;
    sol_sock_J << "solution\n" << mesh << grad_v
                << "window_title 'Solution: Gradient'" 
                << "pause\n" << "keys c\n" << flush;


    
    cout << "test" << endl;
    
    delete a;
    // delete pcOp;
    delete b;
    delete fespace;
    delete fec;
    delete fespace_grad;

    return 0;
}
void A_func(const Vector &x, Vector &A){
    double B = 1;
    double sigma_ = 0.2;
    double omega = M_PI*2*200;
    // double r = sqrt(x(1)*x(1) + x(0)*x(0));
    A(0) = sigma_*omega*B/2*x(1);
    A(1) = -sigma_*omega*B/2*x(0);
    A(2) = 0.f;  
}
double A_func_bdr(const Vector &x)
{
    double B = 1;
    double omega = M_PI*2* 200;
    real_t norm = sqrt(x(1)*x(1) + x(0)*x(0));
    return B*omega*(-x(0)/norm * omega*B/2*x(1) + x(1)/norm * omega*B/2*x(0));  
}