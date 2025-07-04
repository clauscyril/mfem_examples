#include<iostream>
#include "mfem.hpp"

using namespace mfem;

int main() {
    
    const char *mesh_file = "../../../../mesh/square.msh";
    Mesh *mesh = new Mesh(mesh_file, 1 , 1);
    int order = 1;
    int dim = mesh->Dimension();

    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 0;  dir_bdr[2] = 1, dir_bdr[0] = 0;
    Array<int> dir_bdr_fake(mesh->bdr_attributes.Max());
    dir_bdr_fake = 0;
    Array<int> ess_tdof_list;
    Array<int> ess_tdof_list_fake;
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);
    fespace->GetEssentialTrueDofs(dir_bdr_fake, ess_tdof_list_fake);

    auto bdr_func = [](const Vector &x) -> real_t
    {
        return x(0)*x(0);
    };


    FunctionCoefficient bdr_coeff(bdr_func);
    ConstantCoefficient one(1.0);

    GridFunction x(fespace);
    x = 0;
    x.ProjectBdrCoefficient(bdr_coeff, dir_bdr);


    BilinearForm a(fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.EnableStaticCondensation();
    a.Assemble();

    LinearForm b(fespace);
    b = 0;
    b.Assemble();

    SparseMatrix A = a.SpMat();
    Vector X, B;
    B = b;

    std::cout << a.StaticCondensationIsEnabled() << std::endl;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // StaticCondensation static_cond(fespace);
    // static_cond.SetEssentialTrueDofs(ess_tdof_list);
    // std::cout << static_cond.GetNPrDofs() <<std::endl;
    // std::cout << static_cond.GetNExDofs() <<std::endl;  
    
    UMFPackSolver solver;
    solver.SetOperator(A);
    std::cout << "Size Linear System : " << B.Size() << std::endl;
    // B.Print(std::cout);


    solver.Mult(B,X);
    a.RecoverFEMSolution(X,b,x);

    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock.precision(8);
    // sol_sock_i.precision(8);
    sol_sock << "solution\n" << *mesh << x << "window_title 'Solution'" << std::flush;  // Génère une erreur dans glvis à corriger    

    
    return 0;
}