#include <iostream>
#include <complex>
#include "mfem.hpp"

#define M_PI 3.14159265358979323846

using namespace mfem;

void bdr_func(const Vector &x, Vector &H);

// paramètres 
real_t rho = 1.;
real_t sigma = 1.;
real_t eps = 1.;
real_t mu = 1;

real_t tau = 1;  // Parmaètre déterminé par un fit 

real_t omega = 1.;

// Dimensions Torre
real_t h = 1;
real_t w = 2*M_PI*1000.;

// Paramètre Conditions aux limites
real_t N = 1;
real_t I = 1;
real_t r = 1;


// Coefficients complexes
// j imaginaire pur 
const std::complex<real_t> j(0., 1.);

const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps);  // On est "obliger" de caster 1 en real_t (Au moins metre 1.)
const std::complex<real_t> mu_eq = j * mu / ((real_t)1. + tau * j * omega);


int main(){
    const char *path = "../../mesh/torre.msh";
    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->PrintInfo(std::cout);

    int order = 1;
    int dim = 3;

    FiniteElementCollection *fec = new ND_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    ComplexGridFunction h(fespace);

    // Conditions aux limites :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    // dir_bdr[26] = 1; 
    dir_bdr = 1; 
    
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // On précise dans le fespace quelles sont les "lignes" de Dirichlets
    VectorFunctionCoefficient bdr_coeff(dim, bdr_func);
    Vector zero(dim);
    zero = 0.;
    VectorConstantCoefficient zero_coeff(zero);
    h.ProjectBdrCoefficientTangent(bdr_coeff, zero_coeff, dir_bdr);


    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    b.Assemble();

    // ConstantCoefficient rho_real(sigma/(sigma*sigma + pow(eps*omega, 2))); 
    // ConstantCoefficient rho_imag(-omega/(sigma*sigma + pow(eps*omega, 2))); 
    ConstantCoefficient rho_real(rho_eq.real()); 
    ConstantCoefficient rho_imag(rho_eq.imag()); 
    

    // ConstantCoefficient mu_real(mu*omega*omega*tau/(1+pow(tau*omega,2)));
    // ConstantCoefficient mu_imag(mu*omega/(1+pow(tau*omega,2)));
    ConstantCoefficient mu_real(mu_eq.real());
    ConstantCoefficient mu_imag(mu_eq.imag());


    std::cout << "mu_real = " << mu_eq.real() << ", mu imag = " << mu_eq.imag() << std::endl;

    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN);
    a->AddDomainIntegrator(new CurlCurlIntegrator(rho_real), new CurlCurlIntegrator(rho_imag));
    a->AddDomainIntegrator(new VectorFEMassIntegrator(mu_real), new VectorFEMassIntegrator(mu_imag));
    a->Assemble();

    // Définition du preconditionneur (Uniquement la partie réel choisie ici)
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new CurlCurlIntegrator(rho_real));
    pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(mu_real));
    pcOp->Assemble();


    OperatorHandle A, Pc;
    Vector B, H; 

    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B);
    
    // Utilisation du préconditionneur `pc0p` dans le solveur
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));  

    PCG(*A, M, B, H, 1, 1000, 1e-12, 0.f);

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
    gmres.Mult(B, H);

    a->RecoverFEMSolution(H,b,h);

    GridFunction h_r = h.real();
    GridFunction h_i = h.imag();

   // Visualisation GLVis
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_i(vishost, visport);
    socketstream sol_sock_grad(vishost, visport);
    sol_sock_i.precision(8);
    sol_sock_grad.precision(8);
    sol_sock_i << "solution\n" << *mesh << h_r
                << "window_title 'Solution: Real Part'" 
                << "pause\n" << "keys c\n" << std::flush;
    sol_sock_grad << "solution\n" << *mesh << h_i
                << "window_title 'Solution: Imaginary Part'" 
                << "pause\n" << "keys c\n" << std::flush;

                
    ParaViewDataCollection paraview_dc("torre_solution", mesh);
    paraview_dc.SetPrefixPath("paraview_output"); // Optionnel : dossier de sortie
    paraview_dc.RegisterField("H_real", &h_r);
    paraview_dc.RegisterField("H_imag", &h_i);
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(0.0);
    paraview_dc.Save();


    return 0;
}

void bdr_func(const Vector &x, Vector &H)
{
    // Cas Torre :
    real_t r = sqrt(x[0]*x[0] + x[1]*x[1]);
    real_t module = N*I/(2 * M_PI * r);
    H(0) = -x[1]/r * module;
    H(1) = x[0]/r * module;
    H(2) = 0.;


    // Cas Cube : 
    // real_t Ri = 1;
    // real_t r = Ri + x[0];
    // H(0) = 0;
    // H(1) = 0;
    // H(2) = N*I/(2 * M_PI * r);
}