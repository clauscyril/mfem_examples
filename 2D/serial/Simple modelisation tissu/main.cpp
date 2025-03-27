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
    int order = 4;  

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



    // PERMET D'AFFICHER LES COORDONNEES DES NOEUDS DE BORD : 
    
    set<int> boundary_vertices; // Utilisation d'un set pour éviter les doublons
    
    // Récupération des sommets des éléments de bord
    for (int i = 0; i < mesh.GetNBE(); i++) // GetNBE() : nombre d'éléments de bord
    {
        Array<int> vertices;
        mesh.GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++)
        {
            boundary_vertices.insert(vertices[j]); // Stocker l'indice du sommet
        }
    }
    // Affichage des coordonnées des sommets sur le bord
    for (int v : boundary_vertices)
    {
        const double *x = mesh.GetVertex(v);
        Vector coord(2);
        coord(0) = x[0]; coord(1) = x[1];
        cout << "(" << x[0] << ", " << x[1] << ") : A.n = " << A_func_bdr(coord) << endl;
    }

    
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

    b->Assemble();

    ConstantCoefficient sigma(-0.2f);
    // ConstantCoefficient m(1);

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



    
    // 🔹 Utilisation du préconditionneur `pc0p` dans le solveur
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));  

    cout << "test1" << endl;
    PCG(*A, M, B, U, 1, 1000, 1e-12, 0.f);
    // cout << "test2" << endl;

    // Récupération des solutions
    a->RecoverFEMSolution(U, *b, v);
    // cout << "test3" << endl;
    



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
    gmres.SetMaxIter(10000);
    gmres.SetPrintLevel(1);
    gmres.Mult(B, U);

    a->RecoverFEMSolution(U, *b, v);

    ofstream sol_r_ofs("sol_r.gf");
    ofstream sol_i_ofs("sol_i.gf");
    sol_r_ofs.precision(8);
    sol_i_ofs.precision(8);
    v.real().Save(sol_r_ofs);
    v.imag().Save(sol_i_ofs);

    GridFunction v_i = v.imag();

    GradientGridFunctionCoefficient grad_v_coeff(&v_i);

    FiniteElementCollection *fec_grad = new ND_FECollection(order, dim );  // Raviart-Thomas (RT) pour le gradient
    FiniteElementSpace *fespace_grad = new FiniteElementSpace(&mesh, fec_grad);
    
    GridFunction grad_v(fespace_grad);

    grad_v.ProjectCoefficient(grad_v_coeff);
    GridFunction J(fespace_grad);
    VectorFunctionCoefficient A_coeff(dim , A_func);

    J.ProjectCoefficient(A_coeff);
    
    J += grad_v;


    cout << "test" << endl;
   // Visualisation GLVis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_i(vishost, visport);
    socketstream sol_sock_grad(vishost, visport);
    sol_sock_i.precision(8);
    sol_sock_grad.precision(8);
    sol_sock_i << "solution\n" << mesh << v.imag()
                << "window_title 'Solution: Imaginary Part'" 
                << "pause\n" << "keys c\n" << flush;
    sol_sock_grad << "solution\n" << mesh << J
                << "window_title 'Solution: Gradient'" 
                << "pause\n" << "keys c\n" << flush;


    
    cout << "test" << endl;
    
    delete a;
    delete pcOp;
    delete b;
    delete fespace;
    delete fec;
    delete fespace_grad;

    return 0;
}
void A_func(const Vector &x, Vector &A){
    double B = 1;
    double sigma_ = 0.2;
    double omega = M_PI*2*50;
    // double r = sqrt(x(1)*x(1) + x(0)*x(0));
    A(0) = sigma_*omega*B/2*x(1);
    A(1) = -sigma_*omega*B/2*x(0);
    A(2) = 0.f;  
}
double A_func_bdr(const Vector &x)
{
    double B = 1;
    double omega = M_PI*2* 50;
    real_t norm = sqrt(x(1)*x(1) + x(0)*x(0));
    return B*omega*(-x(0)/norm * omega*B/2*x(1) + x(1)/norm * omega*B/2*x(0));  
}