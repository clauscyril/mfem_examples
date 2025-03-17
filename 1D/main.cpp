#include <mfem.hpp>
#include <iostream>
#include <cmath>
// #include <complex>

using namespace mfem;
using namespace std;

ConstantCoefficient Hs = ConstantCoefficient(1.f);
static real_t e = 0.001;

static real_t sigma = 1e7;
static real_t mu(100*4e-7*M_PI);
static real_t omega = 2*M_PI*400;

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
    Mesh mesh =  Mesh::MakeCartesian1D(n,real_t(e));


    // Définir un espace d'éléments finis
    FiniteElementCollection *fec = NULL;
    fec = new H1_FECollection(order, mesh.Dimension());  // Ordre fonctions de tests et dimension du mesh
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);

    // Vecteur solution
    ComplexGridFunction u(fespace);
    u = 0.0;


    // Modification du nombre de quadrature (pour le calcul des intégrales)
    int quad_order = 10*order;
    IntegrationRules ir;
    const IntegrationRule *irule = &ir.Get(Geometry::SEGMENT, quad_order);

    // Second membre (forme liniéaire)
    ComplexLinearForm b(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    // ComplexLinearForm *b = new ComplexLinearForm()
    b.Vector::operator=(0.0);
    b.Assemble();
    
    // Forme Bilinéaire 

    ConstantCoefficient facteur(sigma*mu*omega);
    ConstantCoefficient one(-1.0);

    // Dans le cas où on a une equation complexe, il faut utiliser SesquilinearForm et non BilinearForm
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::BLOCK_SYMMETRIC);
    a->AddDomainIntegrator(new DiffusionIntegrator(one, irule),NULL);
    a->AddDomainIntegrator(NULL, new MassIntegrator(facteur, irule));
    a->Assemble();

    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(one, irule));
    pcOp->AddDomainIntegrator(new MassIntegrator(facteur, irule));
    pcOp->Assemble();

    
    // Gestion des conditions aux limites
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1; // Imposer u = 1 sur les bords
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    ConstantCoefficient ZeroCoef(0.f);
    
    u.ProjectBdrCoefficient(Hs, ZeroCoef, ess_bdr);  // On impose Hs + j*ZeroCoef sur les bords


    // Définition des Matrices / Vecteurs du système linéaire
    OperatorHandle A;
    Vector B, U;

    a->FormLinearSystem(ess_tdof_list, u, b, A, U, B);
    cout << A->Width() << endl << endl;

    // COPIE DE L'EXEMPLE 22
    // 11. Define and apply a GMRES solver for AU=B with a block diagonal
    //     preconditioner based on the appropriate sparse smoother.

    
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

    a->RecoverFEMSolution(U, b, u);


    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);

    ofstream sol_r_ofs("sol_r.gf");
    ofstream sol_i_ofs("sol_i.gf");
    sol_r_ofs.precision(8);
    sol_i_ofs.precision(8);
    u.real().Save(sol_r_ofs);
    u.imag().Save(sol_i_ofs);


    // Visualisation GLvis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_r(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    sol_sock_r.precision(8);
    // sol_sock_i.precision(8);
    sol_sock_r << "solution\n" << mesh << u.real()
               << "window_title 'Solution: Real Part'" 
               << "pause\n" << "keys c\n" << flush;
    
    return 0;
}

float u_ex(float z){
    return 0.f;
}

// complex<real_t> u0_exact(const Vector &x){
//     // int dim = x.Size();
//     complex<real_t> i(0.0, 1.0);
//     complex<real_t> A = (sqrt(-i * omega * sigma*mu));
//     // return 1/sinh(A*e) * (sinh(A*e*0.5 + A * x));
//     // return 
// }