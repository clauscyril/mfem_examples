#include <mfem.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

void A_func1(const Vector &x, Vector &A);
double A_func_bdr1(const Vector &x);
void A_func2(const Vector &x, Vector &A);
double A_func_bdr2(const Vector &x);

double B = 0.1;
double sigma_ = 0.2;
double omega = M_PI*2*50;
double alpha_ = 1.f;


int main(int argc, char* argv[])
{   
    int n;       
    int order = 1;  

    const char *path = "../disque.msh";
    Mesh mesh(path, 1, 1);

    mesh.UniformRefinement();
    mesh.UniformRefinement();
    // mesh.UniformRefinement();
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
    // mesh.PrintInfo(cout);

    // Gestion des conditions aux limites
    Array<int> ess_tdof_list;

    // for (int i = 0; i < mesh.GetNV(); i++) {
    //     const double *v = mesh.GetVertex(i);
    //     if (abs(v[0]) > 0.45) {  // 
    //         ess_tdof_list.Append(i);
    //     } 
    // }

    GridFunction v(fespace);
    // v = 0.f;
    // for (int i = 0; i < ess_tdof_list.Size(); i++){
    //     const double *u = mesh.GetVertex(ess_tdof_list[i]);
    //     if (u[0] > 0) {
    //         v(ess_tdof_list[i]) = 1.f;
    //     } else {
    //         v(ess_tdof_list[i]) = -1.f;
    //     }
    // }
    // IMPOSE TOUS LES BORDS à 0.
    // // Marquer les DOFs associés aux conditions de Dirichlet
    // Array<int> ess_bdr(mesh.bdr_attributes.Max());
    // ess_bdr = 1; // Appliquer la condition de Dirichlet sur tout le bord
    // fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    cout << mesh.bdr_attributes.Size() << endl << endl;

    VectorFunctionCoefficient sigmajwA(dim, A_func1);
    // VectorFunctionCoefficient sigmajwA(dim, A_func2);

    // FunctionCoefficient omegaA(A_func_bdr1);
    // FunctionCoefficient omegaA(A_func_bdr2);
    // Array<int> surf_attrib;
    // surf_attrib.Append(1);  

    LinearForm *b = new LinearForm(fespace);
    b->AddDomainIntegrator(new DomainLFGradIntegrator(sigmajwA)); 
    // b->AddBoundaryIntegrator(new BoundaryLFIntegrator(omegaA), surf_attrib);
    b->Assemble();

    ConstantCoefficient sigma(sigma_);

    // Forme bilinéaire réelle principale
    BilinearForm *a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a->Assemble();

    OperatorPtr A;
    Vector B, X;

    a->FormLinearSystem(ess_tdof_list, v, *b, A, X, B);
    cout << "Taille du système linéaire : " << A->Height() << endl;

    GSSmoother M((SparseMatrix&)(*A));
    // PCG(*A, M, B, X, 1, 1000, 1e-12, 0.f);

    CG(*A, B, X, 1, 1000, 1e-12, 0.f);




    // On récupère les solutions 
    a->RecoverFEMSolution(X, *b, v);

    ofstream sol_r_ofs("sol_real.gf");

    sol_r_ofs.precision(8);

    v.Save(sol_r_ofs);

    GradientGridFunctionCoefficient grad_v_coeff(&v);

    FiniteElementCollection *fec_grad = new ND_FECollection(order, dim);  // Raviart-Thomas (RT) pour le gradient
    FiniteElementSpace *fespace_grad = new FiniteElementSpace(&mesh, fec_grad);

    GridFunction grad_v(fespace_grad);
    

    grad_v.ProjectCoefficient(grad_v_coeff);
    cout << "Dimension : " << dim << endl;


    GridFunction J(fespace_grad);
    
    VectorFunctionCoefficient A_coeff(dim , A_func1);
    // VectorFunctionCoefficient A_coeff(dim , A_func2);
    J.ProjectCoefficient(A_coeff);
    
    grad_v *= sigma_;
    J -= grad_v;

    ofstream sol_j("sol_j.gf");

    sol_j.precision(8);

    J.Save(sol_j);


    ParaViewDataCollection *pd = NULL;
  
    pd = new ParaViewDataCollection("Solutions", &mesh);
    pd->SetPrefixPath("ParaView");
    pd->RegisterField("solution", &v);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();




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
    sol_sock_J << "solution\n" << mesh << J
                << "window_title 'Solution: Gradient'" 
                << "pause\n" << "keys c\n" << flush;
    cout << "test" << endl;
    
    delete a;
    // delete pcOp;
    delete b;
    delete fespace;
    delete fec;
    // delete fespace_grad;

    return 0;
}
void A_func1(const Vector &x, Vector &A_vect){
    // double r = sqrt(x(1)*x(1) + x(0)*x(0));
    // cout << A_vect.Size();
    // A_vect(0) = -sigma_*omega*B/2*x(1) + 2*sigma_*omega*alpha_ * x(0);
    // A_vect(1) = +sigma_*omega*B/2*x(0) + 2*sigma_*omega*alpha_ * x(1);
    A_vect(0) = -sigma_*omega*B/2*x(1);
    A_vect(1) = +sigma_*omega*B/2*x(0);
    // A_vect(2) = 0.f;  
    // A_vect(0) = 0;
    // A_vect(1) = 0;
    // A_vect(2) = 0.f;  
}


void A_func2(const Vector &x, Vector &A){
    // double r = sqrt(x(1)*x(1) + x(0)*x(0));
    A(0) = 0.f;
    A(1) = sigma_*omega*B*x(0);
    A(2) = 0.f;  
}
