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
    int order = 2;  

    const char *path = "../disque.msh";
    Mesh mesh(path, 1, 1);

    mesh.UniformRefinement();
    // mesh.UniformRefinement();
    // mesh.UniformRefinement();
    // mesh.UniformRefinement();

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;

    // mesh.UniformRefinement();
    
    FiniteElementCollection *fec = new H1_FECollection(order, dim);    
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;
    mesh.PrintInfo(cout);

    // Gestion des conditions aux limites
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(mesh.bdr_attributes.Max());  
    ess_bdr = 0;  
    
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    cout << mesh.bdr_attributes.Size() << endl << endl;

    
    ComplexGridFunction v(fespace);
    v = 0.f;

    // ConstantCoefficient DirichletCoeff(100.f);
    // ConstantCoefficient ZeroCoef(0.f);
    // v.ProjectBdrCoefficient(DirichletCoeff, DirichletCoeff, ess_tdof_list);

    VectorFunctionCoefficient sigmajwA(dim + 1, A_func);
    FunctionCoefficient omegaA(A_func_bdr);
    // Forme linéaire complexe
    // ConstantCoefficient omega(2*M_PI*50.f);


    Array<int> surf_attrib;
    surf_attrib.Append(1);  // Numéro de la surface où la charge doit être nulle

    ComplexLinearForm *b = new ComplexLinearForm(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    b->AddDomainIntegrator(NULL, new DomainLFGradIntegrator(sigmajwA)); 
    b->AddBoundaryIntegrator(NULL, new BoundaryLFIntegrator(omegaA), surf_attrib); 
    // b->AddDomainIntegrator(new DomainLFGradIntegrator(sigmajwA),NULL); 
    // b->AddBoundaryIntegrator(new BoundaryLFIntegrator(omegaA),NULL, surf_attrib); 
    b->Assemble();

    ConstantCoefficient sigma(0.2f);
    ConstantCoefficient m(1);

    // Forme bilinéaire réelle principale
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    a->AddDomainIntegrator(new DiffusionIntegrator(sigma), NULL); 
    // a->AddDomainIntegrator(new MassIntegrator(m),NULL);
    a->Assemble();

    // 🔹 Ajout du préconditionneur `pc0p` (une version simplifiée de `a`)
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(sigma)); 
    // pcOp->AddDomainIntegrator(new MassIntegrator(m)); 
    pcOp->Assemble();  

    OperatorHandle A, Pc;
    Vector B, U;

    // A->PrintMatlab(cout);

    a->FormLinearSystem(ess_tdof_list, v, *b, A, U, B);
    cout << "Taille du système : " << A->Width() << endl << endl;

    SparseMatrix *A_sparse = A.As<SparseMatrix>(); // Convertit A en SparseMatrix
    if (A_sparse)
    {
        ofstream out("matrice_A.txt");
        A_sparse->PrintMatlab(out);
        out.close();
        cout << "Matrice A enregistrée dans 'matrice_A.txt'" << endl;
    }


    // B.Print(cout);



    
    // // 🔹 Utilisation du préconditionneur `pc0p` dans le solveur
    // pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    // GSSmoother M((SparseMatrix&)(*Pc));  

    // cout << "test1" << endl;
    // PCG(*A, M, B, U, 1, 1000, 1e-12, 0.f);
    // // cout << "test2" << endl;

    // // Récupération des solutions
    // a->RecoverFEMSolution(U, *b, v);
    // // cout << "test3" << endl;
    



    // Solver avec Preconditionner et MGRESSolver

    Array<int> blockOffsets;
    blockOffsets.SetSize(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = A->Height() / 2;
    blockOffsets[2] = A->Height() / 2;
    blockOffsets.PartialSum();

    BlockDiagonalPreconditioner BDP(blockOffsets);

    Operator * pc_r = NULL;
    Operator * pc_i = NULL;

    OperatorHandle PCOp;
    
    pcOp->SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
    pcOp->FormSystemMatrix(ess_tdof_list, PCOp);

    pc_r = new DSmoother(*PCOp.As<SparseMatrix>());

    real_t s = -1.0;
    
    pc_i = new ScaledOperator(pc_r, s);


    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);

    BDP.owns_blocks = 1;


    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);
    gmres.SetOperator(*A.Ptr());
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(1000);
    gmres.SetPrintLevel(1);
    gmres.Mult(B, U);

    a->RecoverFEMSolution(U, *b, v);

    ofstream sol_r_ofs("sol_r.gf");
    ofstream sol_i_ofs("sol_i.gf");
    sol_r_ofs.precision(8);
    sol_i_ofs.precision(8);
    v.real().Save(sol_r_ofs);
    v.imag().Save(sol_i_ofs);


    // AFFICHAGE DU GRADIENT Pour afficher j

    // Définition d’un nouvel espace pour le gradient (DG0 ou RT pour H(div))
    FiniteElementCollection *grad_fec = new L2_FECollection(0, mesh.Dimension());
    FiniteElementSpace *grad_fes = new FiniteElementSpace(&mesh, grad_fec, mesh.Dimension());

    // GridFunction pour stocker le gradient
    GridFunction grad_u(grad_fes);

    GridFunction x_real(fespace);
    GridFunction x_imag(fespace);

    std::ifstream real_file("sol_r.gf");
    std::ifstream imag_file("sol_i.gf");
    
    x_real.Load(real_file);
    x_imag.Load(imag_file);
    
    // GridFunction pour stocker le gradient réel et imaginaire
    GridFunction grad_real(grad_fes);
    GridFunction grad_imag(grad_fes);

        // Projection du gradient
    GradientGridFunctionCoefficient grad_real_coeff(&x_real);
    GradientGridFunctionCoefficient grad_imag_coeff(&x_imag);
    grad_real.ProjectCoefficient(grad_real_coeff);
    grad_imag.ProjectCoefficient(grad_imag_coeff);

    // Visualisation GLVis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_r(vishost, visport);
    socketstream sol_sock_i(vishost, visport);
    sol_sock_r.precision(8);
    sol_sock_i.precision(8);
    sol_sock_r << "solution\n" << mesh << v.real()
                << "window_title 'Solution: Real Part'" 
                << "pause\n" << "keys c\n" << flush;
    sol_sock_i << "solution\n" << mesh << v.imag()
                << "window_title 'Solution: Imaginary Part'" 
                << "pause\n" << "keys c\n" << flush;


    // char vishost[] = "localhost";
    // int  visport   = 19916;
    socketstream grad_sock_r(vishost, visport);
    socketstream grad_sock_i(vishost, visport);
    grad_sock_r.precision(8);
    grad_sock_i.precision(8);
    grad_sock_r << "solution\n" << mesh << grad_real
                << "window_title 'Solution: Real Part'" 
                << "pause\n" << "keys c\n" << flush;
    grad_sock_i << "solution\n" << mesh << grad_imag
                << "window_title 'Solution: Imaginary Part'" 
                << "pause\n" << "keys c\n" << flush;
    
    cout << "test" << endl;
    
    delete a;
    delete pcOp;
    delete b;
    delete fespace;
    delete fec;

    return 0;
}

void A_func(const Vector &x, Vector &A){
    double B = 0.01;
    double sigma_ = 0.2;
    double omega = M_PI*2*50;
    A(0) = +sigma_*omega*B/2*x(1);
    A(1) = -sigma_*omega*B/2*x(0);
    A(2) = 0.f;  
}
double A_func_bdr(const Vector &x)
{
    double B = 0.01;
    double omega = M_PI*2* 50;
    real_t norm = sqrt(x(1)*x(1) + x(0)*x(0));
    return B*omega*(-x(0)/norm * omega*B/2*x(1) + x(1)/norm * omega*B/2*x(1));  
}